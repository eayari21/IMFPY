# isd_py/main.py
import numpy as np

from .constants import ConstantsCGS
from .startup import read_istart, startup_grain_radius_micron, integration_timestep_seconds
from .grain import dust_setup
from .transforms import hae_to_hse, hse_to_hae
from .field_model import FieldModel
from .merson import MersonIntegrator
from .rhs import rhs

def load_charge_table(path="charge.dat"):
    # Fortran reads 41 lines: date(i), potenc(i)
    data = np.loadtxt(path, dtype=float)
    date = data[:, 0]
    potenc = data[:, 1]
    return date, potenc

def main():
    const = ConstantsCGS()

    # Observation date input (Fortran reads year)
    tobs = float(input("Enter observational date [year] (e.g. 2005): ").strip())
    print(f"\nT_obs={tobs}\n")

    if tobs >= 2011.0:
        tobs = tobs - 22.0

    # Output
    fcoord = open("Coord.dat", "w")
    print("\nOutput file: Coord.dat (t,x,y,z,vx,vy,vz,p,qVolts)\n")

    # Startup inputs
    tpl_days, ifgr, ifrp, ifl, nochar = read_istart("istart.in")
    deltat = integration_timestep_seconds()

    # Grain
    rg_micron = startup_grain_radius_micron(0.3)
    gr_cm = rg_micron * 1e-4
    grain = dust_setup(gr_cm, const, dust_type=2.0, beta_file="beta-2.dat")
    print(f"rg [micron]= {rg_micron}")
    print(f"beta={grain.beta}\n")

    # Charge table
    date_tab, potenc_tab = load_charge_table("charge.dat")
    tdatemin = 1970.0
    deltatd = 1.0

    # Parameters from Fortran main
    omegasun = 2.86533e-6
    b0 = 3.5e-5
    fi0 = 0.0

    visd = 27.0  # km/s
    rstau = 80.0
    alf = 7.25 * const.pi / 180.0
    xstau = rstau * np.cos(alf)
    xstkm = rstau * const.oneau * 1e-5
    toff = xstkm / visd / 86400.0 / 365.0

    t0_years = tobs - 1964.0 - toff
    print(f"Particles start in {t0_years + 1964.0:.6f}")

    # Start grid (x,z lattice) and initial direction
    nxx = 20
    nzz = 20
    xstmin = -20.0
    xstmax = -5.0
    zstmin = 0.0
    zstmax = 10.0
    nstart = nxx * nzz
    print(f"\nNumber of particles at start= {nstart}\n")
    print(f"Initial distance of ISD particles: {rstau} AU")

    # Sun tilt + frame angle (Fortran)
    slamb = 76.0 * const.pi / 180.0

    # ISD arriving direction in HAE
    dustlong = 259.0 * const.pi / 180.0
    dustlat  = 5.0 * const.pi / 180.0

    rst = rstau * const.oneau
    xh0 = rstau * np.cos(dustlat) * np.cos(dustlong)
    yh0 = rstau * np.cos(dustlat) * np.sin(dustlong)
    zh0 = rstau * np.sin(dustlat)
    print(f"x0,y0,z0= {xh0} {yh0} {zh0}")

    # Initial velocity in HAE [cm/s]
    v0 = visd * 1e5
    vxh = v0 * np.cos(dustlat) * np.cos(dustlong - const.pi)
    vyh = v0 * np.cos(dustlat) * np.sin(dustlong - const.pi)
    vzh = -v0 * np.sin(dustlat)

    # Convert velocity to HSE
    vxhsep = vxh*np.cos(slamb) + vyh*np.sin(slamb)
    vyhsep = -vxh*np.sin(slamb) + vyh*np.cos(slamb)
    vxhse = vxhsep
    vyhse = vzh*np.sin(alf) + vyhsep*np.cos(alf)
    vzhse = vzh*np.cos(alf) - vyhsep*np.sin(alf)

    # Initialize particle arrays in HSE coordinates
    rng = np.random.default_rng(12345)
    xx = np.zeros(nstart, dtype=float)
    yy = np.zeros(nstart, dtype=float)
    zz = np.zeros(nstart, dtype=float)
    vxx = np.full(nstart, vxhse, dtype=float)
    vyy = np.full(nstart, vyhse, dtype=float)
    vzz = np.full(nstart, vzhse, dtype=float)
    qq  = np.zeros(nstart, dtype=float)
    pol = np.zeros(nstart, dtype=float)
    alive = np.ones(nstart, dtype=bool)

    i = 0
    for ix in range(1, nxx+1):
        for iz in range(1, nzz+1):
            xh = (xstmin + ix*(xstmax-xstmin)/nxx) * const.oneau
            yh = (yh0 - rng.random()) * const.oneau
            zh = (zstmin + iz*(zstmax-zstmin)/nzz) * const.oneau

            xhse, yhse, zhse = hae_to_hse(xh, yh, zh, alf, slamb)

            xx[i] = xhse
            yy[i] = yhse
            zz[i] = zhse

            # Fortran initial potential +7 V => q = 7 * gr / 300
            qq[i] = 7.0 * grain.gr / 300.0
            i += 1

    # Field model
    field = FieldModel(const=const, t0_years=t0_years, fi0=fi0)

    # Integrator settings (Fortran merson params)
    integrator = MersonIntegrator(acc=1e-5, h=1.0, hmin=1e-10, jtest=0)

    flags = {"ifgr": ifgr, "ifrp": ifrp, "ifl": ifl, "nochar": nochar}

    time = 0.0
    nstep = 0

    # Fortran prints yearly; outputs every 10 days (modstep2=20 for deltat=0.5 days)
    modstep = 365          # steps per year when dt=0.5 day
    modstep2 = 20          # 10 days

    # Convert tpl_days to years-like stop condition in Fortran:
    # Fortran: if tend > 365*86400*tplan then stop, and tplan=tpl (days) in startup print
    # Interpreting tpl as "years" in the Fortran is inconsistent in that print,
    # but the code uses: 365*86400*tplan. We'll treat tpl as years here if you want,
    # OR set tpl to years in istart.in. If your istart.in was built for this code, keep it.
    tplan_years = float(tpl_days)

    while True:
        tend = time + deltat
        tdate = t0_years + 1964.0 + tend/86400.0/365.0

        if tend > (365.0 * 86400.0 * tplan_years):
            print("\nTime > T_plan\n")
            break

        # integrate each particle from time -> tend
        for k in range(nstart):
            if not alive[k]:
                continue

            y0 = np.array([xx[k], yy[k], zz[k], vxx[k], vyy[k], vzz[k], qq[k]], dtype=float)

            def local_rhs(tt, yyvec):
                return rhs(tt, yyvec, const, grain, field, flags)

            try:
                y1 = integrator.integrate(local_rhs, time, tend, y0)
            except RuntimeError:
                alive[k] = False
                continue

            # Evaluate igone/icrash at end by re-calling rhs once (cheap)
            dy, igone, icrash, p, tilt, thetad, theta = rhs(tend, y1, const, grain, field, flags)

            if igone or icrash:
                alive[k] = False
                continue

            xx[k], yy[k], zz[k], vxx[k], vyy[k], vzz[k], qq[k] = y1
            pol[k] = p

            # Update charge based on date, matching Fortran:
            # itd = (tdate - tdatemin)/deltatd   (1-based indexing in Fortran)
            itd = int(np.floor((tdate - tdatemin) / deltatd))
            itd = max(0, min(itd, len(potenc_tab)-1))
            qdate = potenc_tab[itd] * grain.gr / 300.0
            qq[k] = qdate

        # advance
        time = tend
        nstep += 1

        # yearly print
        if (nstep % modstep) == 0:
            print(f"Time [years]: {time/86400.0/365.0:5.1f}")

        # output every 10 days
        if (nstep % modstep2) == 0:
            for k in range(nstart):
                if not alive[k]:
                    continue

                xhei, yhei, zhei = hse_to_hae(xx[k], yy[k], zz[k], alf, slamb)
                vxhei, vyhei, vzhei = hse_to_hae(vxx[k], vyy[k], vzz[k], alf, slamb)

                xau = xhei / const.oneau
                yau = yhei / const.oneau
                zau = zhei / const.oneau

                # qVolts as in Fortran: qq/gr*300
                qvolts = qq[k] / grain.gr * 300.0

                fcoord.write(
                    f"{time/86400.0:16.6e}"
                    f"{xau:16.6e}{yau:16.6e}{zau:16.6e}"
                    f"{(vxhei*1e-5):16.6e}{(vyhei*1e-5):16.6e}{(vzhei*1e-5):16.6e}"
                    f"{pol[k]:16.6e}{qvolts:16.6e}\n"
                )
            fcoord.flush()

        # stop if all gone/crashed
        if not np.any(alive):
            print("\nAll particles lost.\n")
            break

    fcoord.close()
    print("\nProgram completed successfully, Bye!")

if __name__ == "__main__":
    main()
