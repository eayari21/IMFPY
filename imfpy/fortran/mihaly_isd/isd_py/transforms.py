# isd_py/transforms.py
import numpy as np

def hae_to_hse(xh, yh, zh, alf, slamb):
    """
    Fortran:
      xhsep=xh*cos(slamb)+yh*sin(slamb)
      yhsep=-xh*sin(slamb)+yh*cos(slamb)
      xhse=xhsep
      yhse=zh*sin(alf)+yhsep*cos(alf)
      zhse=zh*cos(alf)-yhsep*sin(alf)
    """
    xhsep = xh*np.cos(slamb) + yh*np.sin(slamb)
    yhsep = -xh*np.sin(slamb) + yh*np.cos(slamb)
    xhse = xhsep
    yhse = zh*np.sin(alf) + yhsep*np.cos(alf)
    zhse = zh*np.cos(alf) - yhsep*np.sin(alf)
    return xhse, yhse, zhse

def hse_to_hae(xx, yy, zz, alf, slamb):
    """
    Fortran output block:
      yheip=zz*sin(-alf)+yy*cos(-alf)
      zheip=zz*cos(-alf)-yy*sin(-alf)
      xhei=xx*cos(-slamb)+yheip*sin(-slamb)
      yhei=-xx*sin(-slamb)+yheip*cos(-slamb)
      zhei=zheip
    """
    yheip = zz*np.sin(-alf) + yy*np.cos(-alf)
    zheip = zz*np.cos(-alf) - yy*np.sin(-alf)
    xhei = xx*np.cos(-slamb) + yheip*np.sin(-slamb)
    yhei = -xx*np.sin(-slamb) + yheip*np.cos(-slamb)
    zhei = zheip
    return xhei, yhei, zhei
