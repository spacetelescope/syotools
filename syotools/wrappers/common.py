#import syotools.environment
import sys
import pickle

import pytest
import numpy as np
import astropy.units as u

from syotools.spectra.spec_defaults import syn_spectra_library
from syotools.models import Camera, Spectrograph, Telescope, Source, SourcePhotometricExposure, SourceSpectrographicExposure

def _do_calculation(tel, inst, exp, mode=None, source=None, snr=10.0, exptime=100, bandpass=None, target="magnitude", verbose=False):

    if verbose:
        print(target)

    # a setting for spectrographs
    if mode is not None:
        inst.mode = mode

    if target == "magnitude":
        if isinstance(inst, Spectrograph):
            raise NotImplementedError("Spectrographs cannot currently solve for limiting magnitude")
        
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
        exp.unknown = target
        tel.add_camera(inst)
        inst.add_exposure(exp)

        result = exp.snr

    return result



def compute_observation(telescope, instrument="hri", sed="G2V Star", magnitude=20.0, snr=10.0, exptime=100, redshift=0, extinction=0, bandpass="johnson,v", target="magnitude", verbose=False):
    
    if verbose:
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
    result = []
    if instrument in ["camera", "hri"]:
        inst = Camera()
        inst.set_from_sei('HRI')
        exp = SourcePhotometricExposure()

        exp.source = source
        exp.verbose = verbose

        result.append(_do_calculation(tel, inst, exp, source=source, snr=snr, exptime=exptime, bandpass=bandpass, target=target, verbose=verbose))
        
    elif instrument in ["spectroscopy", "uvi"]:
        inst = Spectrograph()
        inst.set_from_sei('UVI')
        inst.bandnames = inst.modes
        exp = SourceSpectrographicExposure() 
        exp.source = source
        exp.verbose = verbose

        for mode in inst.modes:
            result.append(_do_calculation(tel, inst, exp, mode=mode, source=source, snr=snr, exptime=exptime, bandpass=bandpass, target=target, verbose=verbose))

    return result

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

