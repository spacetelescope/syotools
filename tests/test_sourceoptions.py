#import syotools.environment
import sys
import pytest
import numpy as np


from syotools.spectra.spec_defaults import syn_spectra_library

import pickle

from syotools.models import Camera, Spectrograph, Telescope, Source, SourcePhotometricExposure, SourceSpectrographicExposure
import astropy.units as u

telescopes = ["EAC1", "EAC2", "EAC3"]
instruments = ["hri", "uvi"]
seds = syn_spectra_library
redshifts = np.logspace(0,5,20)
extinctions = np.linspace(0,6,10)
snrs = np.logspace(0.01,1e2,20)
exptimes = np.logspace(0.1,1e6,20)
magnitudes = np.linspace(5, 30, 20)
bandpasses = ["johnson,v", "galex,fuv", "2mass,j"]
targets = ["exptime", "snr", "magnitude"]

def _do_calculation(tel, inst, exp, mode=None, source=None, snr=10.0, exptime=100, bandpass=None, target="magnitude", verbose=True):

    print(target)

    # a setting for spectrographs
    if mode is not None:
        inst.mode = mode

    if target == "magnitude":
        
        exp.exptime = [[exptime, exptime, exptime, exptime, exptime, exptime, exptime, exptime, exptime, exptime], 'hr']
        exp._snr = [snr] * u.Unit('electron(1/2)')  
        tel.add_camera(inst)
        inst.add_exposure(exp)

        exp.unknown = target
        result = exp.magnitude

    elif target == "exptime":
        
        exp._snr = [snr] * u.Unit('electron(1/2)')  
        exp.unknown = target
        tel.add_camera(inst)
        inst.add_exposure(exp)

        if verbose:
            print('-- Computing Exptime as the Unknown --') 
            for bb, ee in zip(inst.bandnames, exp.exptime): print("{}, SNR = {}".format(bb, ee)) 
        result = exp.exptime

    elif target == "snr":

        exp.exptime = [[exptime, exptime, exptime, exptime, exptime, exptime, exptime, exptime, exptime, exptime], 'hr']

        inst.add_exposure(exp)
        exp.unknown = target
        result = exp.snr

    return result



def compute_observation(telescope, instrument="hri", sed="G2V Star", magnitude=20.0, snr=10.0, exptime=100, redshift=0, extinction=0, bandpass="johnson,v", target="magnitude", verbose=False):
    
    
    print("telescope:", telescope)
    print("instrument:", instrument)
    print("sed: ", sed)
    print("magnitude", magnitude)
    print("snr: ", snr)
    print("exptime", exptime)
    print("redshift", redshift)
    print("extinction", extinction)
    print("bandpass", bandpass)
    print("target", target)
    # create the source
    source = Source()
	
    source.set_sed(sed, magnitude, redshift, extinction, bandpass=bandpass)   


    # create a Telescope, Camera, and Exposure 
    tel = Telescope()
    tel.set_from_sei(telescope)
    if instrument in ["camera", "hri"]:
        inst = Camera()
        inst.set_from_sei('HRI')
        exp = SourcePhotometricExposure()

        exp.source = source
        exp.verbose = verbose

        result = [_do_calculation(tel, inst, exp, source=source, snr=snr, exptime=exptime, bandpass=bandpass, target=target, verbose=verbose)]
        
    elif instrument in ["spectroscopy", "uvi"]:
        inst = Spectrograph()
        inst.set_from_sei('UVI')
        inst.bandnames = inst.modes
        exp = SourceSpectrographicExposure() 
        exp.source = source
        exp.verbose = verbose

        result = []
        for mode in inst.modes:
            result.append(_do_calculation(tel, inst, exp, mode=mode, source=source, snr=snr, exptime=exptime, bandpass=bandpass, target=target, verbose=verbose))

    return result


def create_comparisons():
    saved = []
    for telescope in telescopes:
        for instrument in instruments:
            for sed in seds:
                redshift=0
                extinction=0
                print(telescope, instrument, sed, magnitudes[0], snrs[0], exptimes[0], redshifts[0], extinctions[0], targets[0], end="")
                try:
                    result = compute_observation(telescope, instrument=instrument, sed=sed, magnitude=magnitudes[0], snr=snrs[0], exptime=exptimes[0], redshift=redshift, extinction=extinction, target=targets[0])
                    #result = np.median(result)
                    print(result)
                    saved.append({"telescope": telescope, "instrument": instrument, "sed": sed, "magnitude": magnitudes[0], "snr": snrs[0], "exptime": exptimes[0], "redshift": redshift, "extinction": extinction, "target": targets[0], "expected": result})
                except Exception as err:
                    print(f" Error in calculation: {err}")
    with open("test_seds.pickle", "wb") as picklefile:
        pickle.dump(saved, picklefile)

'''
LOAD IT
'''
try:    
    with open("test_seds.pickle", "rb") as picklefile:
        test_setups = pickle.load(picklefile)
except FileNotFoundError:
    create_comparisons()

def check_relative_diff(actual, expected, rel_tol=0.1):
    """
    Simple function to check if two lists are approximately equal within a relative tolerance.
    Reports percentage differences for values that exceed the tolerance.

    Args:
        actual: List of actual values
        expected: List of expected values
        rel_tol: Relative tolerance (default: 0.1 or 10%)

    Returns:
        True if all values are within tolerance, False otherwise
    """
    if len(actual) != len(expected):
        print(f"Lists have different lengths: actual={len(actual)}, expected={len(expected)}")
        return False

    all_within_tolerance = True
    differences = []

    for i, (a, e) in enumerate(zip(actual, expected)):
        for j, (act, exp) in enumerate(zip(actual[i], expected[i])):
            if exp == 0:
                # Can't calculate relative difference if expected is zero
                if act != 0:
                    all_within_tolerance = False
                    differences.append((i, a, e, "inf"))
            else:
                rel_diff = abs(act - exp) / abs(exp)
                pct_diff = 100 * rel_diff

                if rel_diff > rel_tol:
                    all_within_tolerance = False
                    differences.append((i, j, act, exp, pct_diff))

    if not all_within_tolerance:
        print("\nValues exceeding relative tolerance:")
        for i, j, act, exp, pct in differences:
            if pct == "inf":
                print(f"  Index {i},{j}: actual={act}, expected={exp}, difference=infinite (division by zero)")
            else:
                print(f"  Index {i},{j}: actual={act:.6g}, expected={exp:.6g}, difference={pct:.2f}%")

    return all_within_tolerance

@pytest.mark.parametrize("inputs", test_setups)
def test_etc_seds(inputs):
    result = compute_observation(inputs["telescope"], instrument=inputs["instrument"], sed=inputs["sed"], 
                        magnitude=inputs["magnitude"], snr=inputs["snr"], exptime=inputs["exptime"], 
                        redshift=inputs["redshift"], extinction=inputs["extinction"], target=inputs["target"])
    result = [res.value for res in result]
    print(inputs["expected"][0].value)
    assert check_relative_diff(result, [res.value for res in inputs["expected"]], 0.0005) #1e-3)


if __name__ == "__main__":
    print(sys.argv[1] == "reset")
    if len(sys.argv) >= 1 and sys.argv[1] == "reset":
        create_comparisons()
