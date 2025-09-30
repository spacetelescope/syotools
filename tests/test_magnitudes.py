#import syotools.environment
import sys
import pickle

import pytest
import numpy as np
import astropy.units as u

from syotools.spectra.spec_defaults import syn_spectra_library
from syotools.models import Camera, Spectrograph, Telescope, Source, SourcePhotometricExposure, SourceSpectrographicExposure
from syotools.wrappers.common import compute_observation, check_relative_diff

telescopes = ["EAC1", "EAC2", "EAC3"]
instruments = ["hri", "uvi"]
seds = syn_spectra_library
redshifts = np.logspace(0,5,20)
extinctions = np.linspace(0,6,10)
snrs = np.logspace(0.01,1e2,20)
exptimes = np.linspace(0.1,1e5,30)
magnitudes = np.linspace(5, 30, 20)
bandpasses = ["johnson,v", "galex,fuv", "2mass,j"]
targets = ["exptime", "snr", "magnitude"]

def create_comparisons(reset):
    saved = []

    redshift = 0
    extinction = 0
    exptime = 10000
    snr = 10
    sed = "G2V Star"
    target = "exptime"

    for telescope in telescopes:
        for instrument in instruments:
            for magnitude in magnitudes:
                print(telescope, instrument, sed, magnitude, snr, exptime, redshift, extinction, target, end="")
                try:
                    result = compute_observation(telescope, instrument=instrument, sed=sed, magnitude=magnitude, snr=snr, exptime=exptime, redshift=redshift, extinction=extinction, target=target)
                    #result = np.median(result)
                    print(result)
                    saved.append({"telescope": telescope, "instrument": instrument, "sed": sed, "magnitude": magnitude, "snr": snr, "exptime": exptime, "redshift": redshift, "extinction": extinction, "target": target, "expected": result})
                except Exception as err:
                    print(f" Error in calculation: {err}")
    if reset:
        with open("test_magnitudes.pickle", "wb") as picklefile:
            pickle.dump(saved, picklefile)

'''
LOAD IT
'''
try:    
    with open("test_magnitudes.pickle", "rb") as picklefile:
        test_setups = pickle.load(picklefile)
except FileNotFoundError:
    create_comparisons(True)
    with open("test_magnitudes.pickle", "rb") as picklefile:
        test_setups = pickle.load(picklefile)

@pytest.mark.parametrize("inputs", test_setups)
def test_etc_exptimes(inputs):
    result = compute_observation(inputs["telescope"], instrument=inputs["instrument"], sed=inputs["sed"], 
                        magnitude=inputs["magnitude"], snr=inputs["snr"], exptime=inputs["exptime"], 
                        redshift=inputs["redshift"], extinction=inputs["extinction"], target=inputs["target"])
    result = [res.value for res in result]
    assert check_relative_diff(result, [res.value for res in inputs["expected"]], 0.0005) #1e-3)


if __name__ == "__main__":
    reset=False
    if (len(sys.argv) > 1) and (sys.argv[1] == "reset"):
        reset = True
    create_comparisons(reset)
