# visualize_EB_simple.py
# Minimal ParaView Python: show E and B vector fields separately as arrow glyphs
# Usage:
#   pvpython visualize_EB_simple.py path/to/data.vti            # uses arrays "E" and "B"
#   pvpython visualize_EB_simple.py path/to/data.vtk Ex B_field # custom names

from paraview.simple import *
import sys, math

# ---------------- user args ----------------
filename = sys.argv[1] if len(sys.argv) > 1 else "data.vti"
E_NAME   = sys.argv[2] if len(sys.argv) > 2 else "E"
B_NAME   = sys.argv[3] if len(sys.argv) > 3 else "B"
# -------------------------------------------

# Open the data (works for most VTK family formats)
src = OpenDataFile(filename)
if src is None:
    raise RuntimeError(f"Could not open: {filename}")

# Create a view
view = CreateView("RenderView")
view.OrientationAxesVisibility = 0
view.Background = [1,1,1]  # white background

# Helper: pick a reasonable glyph scale from dataset bounds
info   = src.GetDataInformation()
bounds = info.GetBounds() if info else (0,1,0,1,0,1)
dx, dy, dz = bounds[1]-bounds[0], bounds[3]-bounds[2], bounds[5]-bounds[4]
diag = math.sqrt(dx*dx + dy*dy + dz*dz) if (dx or dy or dz) else 1.0
GLYPH_SCALE = 0.03 * diag  # quick heuristic

def make_glyph(vec_name):
    g = Glyph(Input=src, GlyphType="Arrow")
    # Vector selection for glyph orientation
    g.Vectors = ["POINTS", vec_name]              # change to "CELLS" if your vectors are cell data
    # Scale arrows by vector magnitude
    g.ScaleArray = ["POINTS", vec_name]
    g.ScaleFactor = GLYPH_SCALE
    # Keep glyph count sane and simple
    g.GlyphMode = "Uniform Spatial Distribution"
    g.MaximumNumberOfSamplePoints = 2000
    return g

def show_and_save(vec_name, out_png, title=""):
    HideAll(view)
    glyph = make_glyph(vec_name)
    disp  = Show(glyph, view)

    # Color by vector magnitude
    ColorBy(disp, ("POINTS", vec_name, "Magnitude"))
    disp.SetScalarBarVisibility(view, True)
    lut = GetColorTransferFunction(vec_name)
    sb  = GetScalarBar(lut, view)
    sb.Title = f"|{title or vec_name}|"
    sb.ComponentTitle = ""

    ResetCamera(view)
    Render()
    SaveScreenshot(out_png, view, ImageResolution=[1600, 1200])
    print(f"Wrote: {out_png}")

# Render E then B
show_and_save(E_NAME, "E_arrows.png", title="E")
show_and_save(B_NAME, "B_arrows.png", title="B")
