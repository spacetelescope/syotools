#import syotools.environment
import pytest

from syotools.models import Spectrograph
from syotools.models.coronagraph import Coronagraph
from syotools.models.telescope import Telescope
from syotools.models.camera import Camera

from syotools.spectra.spec_defaults import syn_spectra_library

from syotools.wrappers.camera_wrapper import camera_exptime
from syotools.wrappers.uvspec_wrapper import uvspec_exptime

import pickle

from syotools.models import Camera, Telescope, Source, SourcePhotometricExposure
import numpy as np, astropy.units as u


def compute_observation(telescope, instrument = sed, magnitude=30, snr=10, exptime=100, redshift=0, extinction=0, bandpass="johnson,v", target="magnitude")
    # create a Telescope, Camera, and Exposure 
	tel = Telescope()
    inst = Camera()
	tel.set_from_sei(telescope)
	if instrument in ["camera", "hri"]:
		inst.set_from_sei('HRI')
        tel.add_camera(inst)
	hri.add_exposure(exp)
    elif instrument in ["spectroscopy", "uvi"]:
        inst.set_from_sei('UVI')
        tel.add_spectrograph(inst)

	source = Source()
	
    source.set_sed(sed, magnitude, redshift, extinction, bandpass=bandpass)   

    exp = SourcePhotometricExposure() 
    exp.source = source

    exp.unknown = target
    if target == "magnitude":
        
        exp.exptime = [[exptime, exptime, exptime, exptime, exptime, exptime, exptime, exptime, exptime, exptime], 'hr']
        exp._snr = [snr] * u.Unit('electron(1/2)')  

        inst.add_exposure(exp)

    elif target == "exptime":
        exp._snr = [snr_goal] * u.Unit('electron(1/2)')  

        inst.add_exposure(exp)

    elif target == "snr":

        exp.exptime = [[exptime, exptime, exptime, exptime, exptime, exptime, exptime, exptime, exptime, exptime], 'hr']

        inst.add_exposure(exp)

    return getattr(exp, target), tel

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
        if e == 0:
            # Can't calculate relative difference if expected is zero
            if a != 0:
                all_within_tolerance = False
                differences.append((i, a, e, "inf"))
        else:
            rel_diff = abs(a - e) / abs(e)
            pct_diff = 100 * rel_diff

            if rel_diff > rel_tol:
                all_within_tolerance = False
                differences.append((i, a, e, pct_diff))

    if not all_within_tolerance:
        print("\nValues exceeding relative tolerance:")
        for i, a, e, pct in differences:
            if pct == "inf":
                print(f"  Index {i}: actual={a}, expected={e}, difference=infinite (division by zero)")
            else:
                print(f"  Index {i}: actual={a:.6g}, expected={e:.6g}, difference={pct:.2f}%")

    return all_within_tolerance

    @pytest.mark.parametrize("telescope, sed, magnitude, snr, exptime, redshift, extinction, bandpass, target, expected", test_setups)
def test_etc_calculation(telescope, sed, magnitude, snr, exptime, redshift, extinction, bandpass, target, expected):
    result = compute_observation(telescope, sed, magnitude, snr, exptime, redshift, extinction, bandpass, target)
    assert check_relative_diff(result, expected, 0.0005) #1e-3)

telescopes = ["EAC1", "EAC2", "EAC3"]
sed = syn_spectra_library
redshift = np.logspace(-0.2,5,20)
extinction = np.linspace(0,6,10)
snr= np.logspace(0.01,1e2,20)
exptime = np.logspace(0.1,1e6,20)
magnitude = np.linspace(5, 30, 20)
bandpass = ["johnson,v", "galex,fuv", "2mass,j"]

