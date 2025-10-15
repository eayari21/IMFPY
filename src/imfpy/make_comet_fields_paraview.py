# make_comet_fields_paraview.py
# ParaView (pvpython) comet E & B on z=0 with masking + constant-length glyphs.

from paraview.simple import *
paraview.simple._DisableFirstRenderCameraReset()

# ---------------------------- Tunables -----------------------------
B0 = 1.0    # far-field Bx
v0 = 1.0    # far-field speed (along -z)
z0 = 0.5    # subsolar ionopause position
c  = 1.0    # speed of light

L  = 1.0    # domain extent in x,y: [-L, L]
N  = 101    # samples per side (image-data grid)
rmask_rel = 0.06   # extra center mask radius (fraction of L)
subsample = 4      # ExtractSubset sampling (bigger=fewer glyphs)
arrow_len = 0.18   # constant arrow length in data units

# ------------------------ ProgrammableSource ------------------------
src = ProgrammableSource()
src.OutputDataSetType = 'vtkImageData'
src.Script = f"""
import numpy as np
from vtkmodules.vtkCommonDataModel import vtkImageData
from vtkmodules.util import numpy_support

B0 = {B0}; v0 = {v0}; z0 = {z0}; c = {c}
L = {L}; N = {N}; rmask_rel = {rmask_rel}
EPS = 1e-12

nx = ny = N; nz = 1
xmin, xmax = -L, L; ymin, ymax = -L, L; z0plane = 0.0

img = vtkImageData()
img.SetDimensions(nx, ny, nz)
img.SetOrigin(xmin, ymin, z0plane)
img.SetSpacing((xmax-xmin)/(nx-1), (ymax-ymin)/(ny-1), 1.0)

xs = np.linspace(xmin, xmax, nx)
ys = np.linspace(ymin, ymax, ny)
X, Y = np.meshgrid(xs, ys, indexing='ij')
Z = np.zeros_like(X)

r = np.hypot(X, Y); r = np.where(r==0.0, EPS, r)
R = np.hypot(r, Z);  R = np.where(R==0.0, EPS, R)

# S^2 = r^2 + 2 z0^2 (z/R - 1)
S_sq_raw = r**2 + 2.0*z0**2*(Z/R - 1.0)
validS   = S_sq_raw >= 0.0
S_sq     = np.where(validS, S_sq_raw, 0.0)
S        = np.sqrt(S_sq)
S        = np.where((S==0.0) & validS, EPS, S)

center = (np.hypot(X, Y) < (rmask_rel*L))
valid_mask = np.logical_and(validS, ~center).astype(np.uint8)

# ---- E field ----
coeff  = (B0*v0)/(c*r)
common = (-S/(r**2)) + (1.0/S) - (z0**2*Z)/(S*R**3)
Ex = coeff * common * (X*Y)
Ey = coeff * (S + common*(Y**2))
Ez = coeff * ((z0**2 * r**2)/(S*R**3)) * Y
E_mag = np.hypot(Ex, Ey)

# ---- B field (Bz=0 case) ----
den = -1.0 + (z0**2)*Z/(R**3)   # v_z / v0
den = np.where(np.abs(den) < EPS, np.sign(den)*EPS, den)
Bx = B0 * ( (( S/(r**2) - 1.0/S + (z0**2)*Z/(S*R**3) )*(Y**2) - S ) / (r*den) )
By = B0 * ( ( -S/(r**2) + 1.0/S - (z0**2)*Z/(S*R**3) )*(X*Y) / (r*den) )
Bz = np.zeros_like(Bx)
B_mag = np.hypot(Bx, By)

# normalized directions for constant-length glyphs
Edir_x = Ex/(E_mag+EPS); Edir_y = Ey/(E_mag+EPS); Edir_z = np.zeros_like(E_mag)
Bdir_x = Bx/(B_mag+EPS); Bdir_y = By/(B_mag+EPS); Bdir_z = np.zeros_like(B_mag)

# ONES array for constant scaling
ones = np.ones_like(E_mag, dtype='f')

# shove into VTK
from vtkmodules.util.numpy_support import numpy_to_vtk
def add_vec(name, a, b, d):
    vec = np.stack([a,b,d], axis=-1).astype('f').reshape(-1,3)
    v = numpy_to_vtk(vec, deep=True); v.SetName(name); return v
def add_sca(name, s, dtype='f'):
    v = numpy_to_vtk(np.ascontiguousarray(s.astype(dtype)).ravel(), deep=True); v.SetName(name); return v

output.ShallowCopy(img)
pd = output.GetPointData()
pd.AddArray(add_vec('E', Ex, Ey, Ez))
pd.AddArray(add_vec('B', Bx, By, Bz))
pd.AddArray(add_sca('E_mag', E_mag))
pd.AddArray(add_sca('B_mag', B_mag))
pd.AddArray(add_vec('E_dir', Edir_x, Edir_y, Edir_z))
pd.AddArray(add_vec('B_dir', Bdir_x, Bdir_y, Bdir_z))
pd.AddArray(add_sca('ones', ones))
pd.AddArray(add_sca('ValidMask', valid_mask, dtype='u1'))
"""

# ---------------------------- Mask + Subsample -----------------------
thr = Threshold(Input=src)
# Array picker across versions
try:
    thr.Scalars = ['POINTS', 'ValidMask']
except Exception:
    try: thr.SelectInputScalars = ['POINTS', 'ValidMask']
    except Exception: pass

# 5.10+ threshold API
if hasattr(thr, 'LowerThreshold') and hasattr(thr, 'UpperThreshold'):
    thr.ThresholdMethod = 'Between'
    thr.LowerThreshold = 0.5
    thr.UpperThreshold = 1.0
else:
    thr.ThresholdRange = [0.5, 1.0]

sub = ExtractSubset(Input=thr)
sub.SampleRateI = subsample
sub.SampleRateJ = subsample
sub.SampleRateK = 1
sub.IncludeBoundary = 1

# ------------------------------ View --------------------------------
view = GetActiveViewOrCreate('RenderView')
view.OrientationAxesVisibility = 1

# ------------------------------ E glyphs -----------------------------
gE = Glyph(Input=sub, GlyphType='Arrow')
gE.OrientationArray = ['POINTS','E_dir']   # just set the array; no "UseOrientationArray" needed
gE.ScaleArray = ['POINTS','ones']          # scale by ones -> constant length
gE.ScaleFactor = arrow_len
dE = Show(gE, view)
ColorBy(dE, ('POINT_DATA', 'E_mag'))
dE.SetScalarBarVisibility(view, True)

# ------------------------------ B glyphs -----------------------------
gB = Glyph(Input=sub, GlyphType='Arrow')
gB.OrientationArray = ['POINTS','B_dir']
gB.ScaleArray = ['POINTS','ones']          # constant length
gB.ScaleFactor = arrow_len
dB = Show(gB, view)
ColorBy(dB, ('POINT_DATA', 'B_mag'))
dB.Opacity = 0.55
dB.SetScalarBarVisibility(view, True)

view.ResetCamera()
Render()

# ---------------------------- Save screenshots ----------------------------
SaveScreenshot("E_field.png", view, ImageResolution=[1600, 1200])
# if you want just the E glyph, you can Hide(dB, view) before saving.

SaveScreenshot("B_field.png", view, ImageResolution=[1600, 1200])
# if you want just the B glyph, Hide(dE, view) and Show(dB, view) before saving.

# Or save both visible together:
SaveScreenshot("EB_fields.png", view, ImageResolution=[1600, 1200])
