import copy
import numpy as np
from astropy import units as u

import pytest

from syotools.models.instrument import Instrument
from syotools.models.telescope import Telescope
from syotools.models.source import Source
from syotools.models.source_exposure import SourcePhotometricExposure

geometries = [{
                "geometry": "point"
            },{
                "geometry": "gaussian2d",
                "major": 0.3 * u.arcsec,
                "minor": 0.2 * u.arcsec,
                "norm_method": "surf_scale",
                "surf_area_units": "arcsec^2"
            },{
                "geometry": "flat",
                "major": 0.3 * u.deg,
                "minor": 0.2 * u.deg,
                "norm_method": "surf_center",
                "surf_area_units": "sr"
            },{
                "geometry": "sersic",
                "major": 0.3 * u.rad,
                "minor": 0.2 * u.rad,
                "norm_method": "surf_scale",
                "sersic_index": 1.0,
                "surf_area_units": "sr"
            },{
                "geometry": "sersic_scale",
                "major": 0.3,
                "minor": 0.2,
                "norm_method": "integ_infinity",
                "sersic_index": 1.0,
                "surf_area_units": "sr"
            },{
                "geometry": "power",
                "norm_method": "surf_center",
                "power_index": 1,
                "r_core": 0.005 * u.arcsec,
                "surf_area_units": "arcsec^2"
            }]
testnames = [x["geometry"] for x in geometries]

@pytest.mark.parametrize("geometries", geometries, ids=testnames)
def test_shape(geometries, verbose=False):
    snr = None
    telescope = Telescope()
    telescope.set_from_hwome("EAC5")
    instrument = telescope.instruments["HRI_S.HRI_S_UVIS_Imager"]

    # set a gaussian source
    source = Source()

    template = "NGC 1068"
    redshift = 0
    extinction = 0
    source.set_sed(template, 30., redshift, extinction, geometry=geometries)   

    exp = SourcePhotometricExposure()
    exp.source = source

    exp.exptime = 1 * u.h
    instrument.add_exposure(exp)
    exp.unknown = "snr"

    snr = exp.snr
    

    if verbose:
        print("Test geometry")
        #print(config)
        print("Shape:", geometries["geometry"])
        print("Geometry:", geometries)
        print("SNR:", snr)
        print("-----------------------")

    # prove it completed
    assert snr is not None

if __name__ == "__main__":
    for geometry in geometries:
        test_shape(geometry, verbose=True)
