import numpy as np
from astropy import units as u

import synphot as syn

from syotools.models.source import Source
from syotools.models.source_exposure import SourceExposure, SourcePhotometricExposure
from syotools.models.instrument import Instrument
from syotools.models.telescope import Telescope

def test_mag(verbose=False):
    tel = Telescope()
    tel.set_from_hwome("EAC5")
    inst = tel.instruments["HRI_S.HRI_S_UVIS_Imager"]

    inst.band = "HRI_S_UVIS.HRI_Johnson_V"

    source_exposure = SourcePhotometricExposure()
    inst.add_exposure(source_exposure)
    source_exposure.exptime = 1
    source_exposure.snr = 30
    source_exposure.calculate_magnitude()
    mag_1 = source_exposure.magnitude[0]

    source_exposure.exptime = 100
    source_exposure.snr = 30
    source_exposure.calculate_magnitude()
    mag_2 = source_exposure.magnitude[0]

    if verbose:
        print("Mag SNR=10:", mag_1)
        print("Mag SNR=100:", mag_2)
        print("Difference (5 expected):", mag_2 - mag_1)
        print("-----------------------")

    assert np.round(np.abs(mag_2.to_value(u.ABmag) - mag_1.to_value(u.ABmag) - 5), 6) < 0.1

def test_snr_mag(verbose=False):
    # This is a separate test to ensure conversions are being done correctly
    tel = Telescope()
    tel.set_from_hwome("EAC5")
    inst = tel.instruments["HRI_S.HRI_S_UVIS_Imager"]

    source = Source()
    template = "Flat (AB)"
    redshift = 0
    extinction = 0
    source.set_sed("Flat (AB)", 25, redshift, extinction, bandpass="johnson,v")
    exp = SourcePhotometricExposure()
    inst.add_exposure(exp)
    exp.source = source

    # Set this SNR
    snr_1 = 10
    inst.add_exposure(exp)
    exp.exptime = 1 * u.hr
    exp.snr = snr_1 * u.dimensionless_unscaled
    exp.unknown = "magnitude"
    mag_1 = exp.magnitude[0].value

    source.set_sed("Flat (AB)", mag_1, redshift, extinction, bandpass="johnson,v")

    exp.unknown = "snr"
    snr_2 = exp.snr[0]

    if verbose:
        print("Initial SNR:", snr_1)
        print("Roundtrip SNR:", snr_2)
        print("Ratio (1 expected):", snr_2/snr_1)
        print("-----------------------")

    assert np.round(np.abs(snr_2-snr_1), 6) < 0.1

if __name__ == "__main__":
    test_mag(verbose=True)
    test_snr_mag(verbose=True)
