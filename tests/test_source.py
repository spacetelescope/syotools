import copy

import numpy as np
from astropy import units as u

from syotools.models.source import Source
import synphot as syn

from syotools.spectra.spec_defaults import syn_spectra_library

def test_source_ranges(verbose=False):
    source = Source()
    MIN_WAVE = 90 * u.nm
    MAX_WAVE = 3000 * u.nm

    passing = False
    num_passing = 0

    for sed in syn_spectra_library:
        source.set_sed(sed, 30, 0, 0)
        wave = source.sed.waveset
        # as a test, just make sure we have SOMETHING that covers
        # the full wavelength range
        if wave[0] < MIN_WAVE and wave[-1] > MAX_WAVE:
            passing = True
            num_passing += 1
            if verbose:
                print("*" , end="")

        if verbose:
            print(f"{sed:34s} ({wave[0]:.3f} -- {wave[-1]:.3f})")
        
    if verbose:
        print(f"Full coverage of {num_passing}/{len(syn_spectra_library)}")
        print("-----------------------")

    assert passing

def test_sed_norm(verbose=False):
    source = Source()
    redshift = 0. # changes to these are not implemented yet 
    extinction = 0.
    magnitude = 10.
    template = "G2V Star"

    sed1 = copy.deepcopy(source.sed)
    source.set_sed(template, magnitude, redshift, extinction, bandpass="2mass,j")
    sed2 = copy.deepcopy(source.sed)

    if verbose:
        print("Initial Median Flux:", np.median(sed1(sed1.waveset)))
        print("Post-Norm Median Flux:", np.median(sed2(sed2.waveset)))
        print("-----------------------")

    assert np.median(sed1(sed1.waveset)) != np.median(sed2(sed2.waveset))

if __name__ == "__main__":
    test_sed_norm(verbose=True)
    test_source_ranges(verbose=True)
