import numpy as np
from astropy import units as u

import synphot as syn

from syotools.models.source_exposure import SourceExposure
from syotools.models.instrument import Instrument

class Mock_Instrument(Instrument):
    def __init__(self):
        self.configuration = {"detector": {"total_qe": 1}}

    def _print_initcon(self, verbose):
        pass

class Mock_Telescope():
    def __init__(self):
        self.effective_area = 10000 * u.cm**2


class Mock_SourceExposure(SourceExposure):

    def __init__(self, snr=10, exptime=30, nexp=1, fsource=10, fsky=1, thermal=0, dark=0, read_noise=0):
        self.snr = snr
        self.n_exp = nexp
        self.exptime = exptime
        self.fsource = fsource
        self.fsky = fsky
        self.thermal = thermal
        self.dark = dark
        self.read_noise = read_noise 
        self.verbose = False
        self.instrument = Mock_Instrument()
        self.telescope = Mock_Telescope()

    @property
    def exptime(self):
        return self._exptime

    @exptime.setter
    def exptime(self, new_exptime):
        self._exptime = new_exptime * u.s

    @property
    def snr(self):
        return self._snr

    @snr.setter
    def snr(self, new_snr):
        self._snr = new_snr * u.dimensionless_unscaled

    def process_observation(self, source, band):
        fsource_countrate = self.fsource * u.ct/u.s
        fsky_countrate = self.fsky * u.ct/u.s
        thermal_countrate = self.thermal * u.ct / u.s
        dark_current = self.dark * u.ct / u.s
        read_noise = self.read_noise * u.ct

        return fsource_countrate, fsky_countrate, thermal_countrate, dark_current, read_noise







def test_exptime(verbose=False):
    source_exposure = Mock_SourceExposure(snr=10)
    exptime_1 = source_exposure._update_exptime(None, None)

    source_exposure.snr = 20
    exptime_2 = source_exposure._update_exptime(None, None)

    if verbose:
        print("Exptime SNR=10:", exptime_1)
        print("Exptime SNR=20:", exptime_2)
        print("Ratio (4 expected):", exptime_1/exptime_2)
        print("-----------------------")

    assert exptime_2 == exptime_1 * 4

def test_snr(verbose=False):
    exptime = 1
    source_exposure = Mock_SourceExposure(exptime=exptime)
    snr_1 = source_exposure._update_snr(None, None)

    source_exposure.exptime = 2
    snr_2 = source_exposure._update_snr(None, None)

    if verbose:
        print("SNR Exptime=1:", snr_1)
        print("SNR Exptime=2:", snr_2)
        print("Ratio (sqrt 2 expected):", snr_2/snr_1)
        print("-----------------------")

    assert snr_2 == snr_1 * np.sqrt(2)

def test_snr_exptime(verbose=False):
    exptime_1 = 30 * u.s
    source_exposure = Mock_SourceExposure(exptime=exptime_1.value)
    snr_1 = source_exposure._update_snr(None, None)

    source_exposure.snr = snr_1
    exptime_2 = source_exposure._update_exptime(None, None)

    if verbose:
        print("Initial Exptime:", exptime_1)
        print("Roundtrip Exptime:", exptime_2)
        print("Ratio (1 expected):", exptime_2/exptime_1)
        print("-----------------------")

    assert np.round(exptime_2, 6) == np.round(exptime_1,6)

def test_flux_snr(verbose=False):
    fsource = 10.0
    fsky = 0
    source_exposure = Mock_SourceExposure(fsource = fsource, fsky= fsky)
    snr_1 = source_exposure._update_snr(None, None)

    source_exposure.fsource = 20.0
    snr_2 = source_exposure._update_snr(None, None)

    if verbose:
        print("SNR Exptime=1:", snr_1)
        print("SNR Exptime=2:", snr_2)
        print("Ratio (sqrt 2 expected):", snr_2/snr_1)
        print("-----------------------")

    assert snr_2 == snr_1 * np.sqrt(2)

def test_mag(verbose=False):
    bandpass = syn.spectrum.SpectralElement(syn.models.Gaussian1D, amplitude=1, mean=5000, stddev=400)
    band = {"bandpass": bandpass}

    snr = 10
    fsky = 0
    source_exposure = Mock_SourceExposure(snr=snr, fsky=fsky)
    mag_1 = source_exposure._update_magnitude(None, band)

    source_exposure.snr = 100
    mag_2 = source_exposure._update_magnitude(None, band)

    if verbose:
        print("Mag SNR=10:", mag_1)
        print("Mag SNR=100:", mag_2)
        print("Difference (5 expected):", mag_1 - mag_2)
        print("-----------------------")

    assert np.round(mag_1, 6).value == np.round(mag_2, 6).value + 5

def test_dark_exptime(verbose=False):
    source_exposure = Mock_SourceExposure(snr=10)
    exptime_1 = source_exposure._update_exptime(None, None)

    source_exposure.dark = 1
    exptime_2 = source_exposure._update_exptime(None, None)

    if verbose:
        print("Exptime Dark=0:", exptime_1)
        print("Exptime Dark=1:", exptime_2)
        print("Ratio >1:", exptime_2/exptime_1)
        print("-----------------------")

    assert exptime_2 > exptime_1

def test_thermal_exptime(verbose=False):
    source_exposure = Mock_SourceExposure(snr=10)
    exptime_1 = source_exposure._update_exptime(None, None)

    source_exposure.thermal = 1
    exptime_2 = source_exposure._update_exptime(None, None)

    if verbose:
        print("Exptime Thermal=0:", exptime_1)
        print("Exptime Thermal=1:", exptime_2)
        print("Ratio >1:", exptime_2/exptime_1)
        print("-----------------------")

    assert exptime_2 > exptime_1

def test_readnoise_exptime(verbose=False):
    source_exposure = Mock_SourceExposure(snr=10)
    exptime_1 = source_exposure._update_exptime(None, None)

    source_exposure.read_noise = 1
    exptime_2 = source_exposure._update_exptime(None, None)

    if verbose:
        print("Exptime Readnoise=0:", exptime_1)
        print("Exptime Readnoise=1:", exptime_2)
        print("Ratio >1:", exptime_2/exptime_1)
        print("-----------------------")

    assert exptime_2 > exptime_1

if __name__ == "__main__":
    test_exptime(verbose=True)
    test_snr(verbose=True)
    test_snr_exptime(verbose=True)
    test_flux_snr(verbose=True)
    test_mag(verbose=True)
    test_dark_exptime(verbose=True)
    test_thermal_exptime(verbose=True)
    test_readnoise_exptime(verbose=True)
