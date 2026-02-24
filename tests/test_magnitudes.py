#import syotools.environment
import sys

import yaml
import pytest
import numpy as np
import astropy.units as u

from syotools.spectra.spec_defaults import syn_spectra_library
from syotools.utils.yaml_utils import read_yaml, write_yaml
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
                print(telescope, instrument, sed, magnitude, snr, exptime, redshift, extinction, target)
                try:
                    actual = compute_observation(telescope, instrument=instrument, sed=sed, magnitude=magnitude, snr=snr, exptime=exptime, redshift=redshift, extinction=extinction, target=target)
                    #result = np.median(result)
                    result = []
                    if actual is not None:
                        for band in actual:
                            result.append({"mean": np.nanmean(band).value, "median": np.nanmedian(band).value, "std": np.nanstd(band).value, "len": len(band)})
                    saved.append({"telescope": telescope, "instrument": instrument, "sed": sed, "magnitude": magnitude, "snr": snr, "exptime": exptime, "redshift": redshift, "extinction": extinction, "target": target, "expected": result})
                except Exception as err:
                    print(f" Error in calculation: {err}")
    if reset:
        write_yaml(saved, "tests/baselines/test_magnitudes.yml")

'''
LOAD IT
'''
try:
    test_setups = read_yaml("tests/baselines/test_magnitudes.yml")
except (FileNotFoundError, yaml.io.UnsupportedOperation):
    create_comparisons(True)
    test_setups = read_yaml("tests/baselines/test_magnitudes.yml")

@pytest.mark.parametrize("inputs", test_setups)
def test_etc_exptimes(inputs):
    actual = compute_observation(inputs["telescope"], instrument=inputs["instrument"], sed=inputs["sed"], 
                        magnitude=inputs["magnitude"], snr=inputs["snr"], exptime=inputs["exptime"], 
                        redshift=inputs["redshift"], extinction=inputs["extinction"], target=inputs["target"])
    result = []
    if actual is not None:
        for band in actual:
            result.append({"mean": np.nanmean(band).value, "median": np.nanmedian(band).value, "std": np.nanstd(band).value, "len": len(band)})
    assert check_relative_diff(result, inputs["expected"], 0.0005) #1e-3)


if __name__ == "__main__":
    reset=False
    if (len(sys.argv) > 1) and (sys.argv[1] == "reset"):
        reset = True
    create_comparisons(reset)
