#import syotools.environment
import os
import sys

import numpy as np
import astropy.units as u

from syotools.spectra.spec_defaults import syn_spectra_library
from syotools.wrappers.common import generate_test
from syotools.models.telescope import Telescope

telescopes = ["EAC5"]
telescope = Telescope()
telescope.set_from_hwome("EAC5")
instruments = list(telescope.instruments.keys())
seds = syn_spectra_library
redshifts = np.logspace(0,5,20)
extinctions = np.linspace(0,6,10)
snrs = np.logspace(0.01,1e2,20)
exptimes = np.linspace(0.1,1e5,30)
magnitudes = np.linspace(5, 30, 20)
bandpasses = ["johnson,v", "galex,fuv", "2mass,j"]
targets = ["exptime", "snr", "magnitude"]

def create_comparisons(reset):
    os.makedirs(f"{os.path.dirname(__file__)}/baselines/magnitudes", exist_ok=True)

    redshift = 0
    extinction = 0
    exptime = 10000
    snr = 10
    sed = "G2V Star"
    target = "magnitude"

    for telescope in telescopes:
        for instrument in instruments:
            if "Imager" not in instrument:
                continue
            for magnitude in magnitudes:
                test_setup = {"telescope": telescope, "instrument": instrument, "sed": sed, "magnitude": magnitude, "snr": snr, "exptime": exptime, "redshift": redshift, "extinction": extinction, "target": target}
                filename = f"{os.path.dirname(__file__)}/baselines/magnitudes/magnitude_{telescope}_{instrument}_{magnitude}.yml"

                generate_test(test_setup, filename, reset)

if __name__ == "__main__":
    reset=False
    if (len(sys.argv) > 1) and (sys.argv[1] == "reset"):
        reset = True
    create_comparisons(reset)
