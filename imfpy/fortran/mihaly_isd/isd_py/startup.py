# isd_py/startup.py
import re

def read_istart(path="istart.in"):
    """
    Fortran-compatible reader for istart.in

    Accepts:
      - headers
      - inline labels
      - arbitrary whitespace
      - key=value formats

    Extracts the FIRST FIVE numeric values in file order.
    """
    numbers = []

    with open(path, "r") as f:
        for line in f:
            # strip comments
            line = line.split("!")[0]
            line = line.split("#")[0]

            # extract all floats/ints from the line
            tokens = re.findall(
                r"[-+]?\d*\.?\d+(?:[eEdD][-+]?\d+)?", line
            )
            for tok in tokens:
                numbers.append(float(tok.replace("D", "E").replace("d", "e")))

    if len(numbers) < 5:
        raise ValueError(
            f"istart.in must contain at least 5 numeric values, found {len(numbers)}"
        )

    tpl, ifgr, ifrp, ifl, nochar = numbers[:5]

    return float(tpl), int(ifgr), int(ifrp), int(ifl), int(nochar)


def startup_grain_radius_micron(default=0.3):
    # Fortran hard-codes rgmicron=0.3 for ns=1
    return float(default)

def integration_timestep_seconds():
    # Fortran: deltat = 86400/2
    return 86400.0 / 2.0
