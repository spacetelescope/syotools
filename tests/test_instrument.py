import copy
import numpy as np
from astropy import units as u

import synphot as syn

from syotools.models.instrument import Instrument
from syotools.models.telescope import Telescope

class Mock_Instrument(Instrument):
    def __init__(self, telescope):
        self.configuration = {"detector": {"total_qe": 1, "thermal": 200.0 * u.K}, "diffraction_limit": 500 * u.AA}

        self.telescope = telescope

    def _print_initcon(self, verbose):
        pass

class Mock_Telescope(Telescope):
    def __init__(self, diameter=10):
        self.effective_diameter = diameter


def test_diameter(verbose=False):
    diameter_1 = 10 * u.m
    area_2 = 1767145.8676442 * u.cm**2
    telescope = Mock_Telescope()

    telescope.effective_diameter = diameter_1
    out_diameter_1 = telescope.effective_diameter
    out_area_1 = telescope.effective_area

    telescope.effective_area = area_2
    out_diameter_2 = telescope.effective_diameter
    out_area_2 = telescope.effective_area

    if verbose:
        print("Diameter1:", out_diameter_1)
        print("Diameter2:", out_diameter_2)
        print("Ratio (3/2 expected):", out_diameter_2/out_diameter_1)
        print("Area1:", out_area_1)
        print("Area2:", out_area_2)
        print("Ratio (2.25 expected):", out_area_2/out_area_1)
        print("-----------------------")

    assert np.round(out_area_2/out_area_1,6) == 2.25

def test_diff_limit(verbose=False):
    telescope = Mock_Telescope()
    instrument = Mock_Instrument(telescope)

    diff_limit_1 = instrument.diffraction_limit(5000 * u.AA)

    diff_limit_2 = instrument.diffraction_limit(10000 * u.AA)    

    if verbose:
        print("Diff Limit 1:", diff_limit_1)
        print("Diff Limit 2:", diff_limit_2)
        print("Ratio (2 expected):", diff_limit_2/diff_limit_1)
        print("-----------------------")

    assert diff_limit_2 == diff_limit_1*2

def test_planck(verbose=False):
    telescope = Mock_Telescope()
    instrument = Mock_Instrument(telescope)

    flux_1 = instrument._planck(5000 * u.AA)
    flux_2 = instrument._planck(50000 * u.AA)
    flux_3 = instrument._planck(500000 * u.AA)


    if verbose:
        print("Temperature:", instrument.configuration["detector"]["thermal"])
        print("Flux 1:", flux_1)
        print("Flux 2:", flux_2)
        print("Flux 3:", flux_3)
        print("-----------------------")

    assert flux_1 < flux_2

def test_save_instrument(verbose=False):
    new_pixel_scale = 0.1 * u.arcsec / u.pix
    wavelength = np.arange(5000,10000, 1000) << u.AA
    telescope = Telescope()
    telescope.set_from_hwome("EAC5")
    instrument = telescope.instruments["HRI_S.HRI_S_UVIS_Imager"]

    pixel_scale_1 = instrument.configuration["pixel_scale"]

    thermal_1 = instrument._c_thermal(wavelength)

    # Now let's dump it, and then blank it out.
    config = copy.deepcopy(instrument.save_to_dict())
    #print(config)
    instrument.configuration={}

    # Modify the saved config
    config["pixel_scale"] = new_pixel_scale

    # If it worked, we should now load the modified config:
    instrument.load_from_dict(config)

    pixel_scale_2 = instrument.configuration["pixel_scale"]

    thermal_2 = instrument._c_thermal(wavelength)

    # try putting the pixel scale back (and make sure thermal_3 matches thermal_1)
    instrument.configuration["pixel_scale"] = pixel_scale_1
    
    thermal_3 = instrument._c_thermal(wavelength)

    if verbose:
        print("Save and Load an Instrument")
        #print(config)
        print("Pixel_scale before:", pixel_scale_1)
        print("Pixel_scale after:", pixel_scale_2)
        print("Thermal before:", thermal_1(wavelength))
        print("Thermal after:", thermal_2(wavelength))
        print("Thermal edited:", thermal_3(wavelength))
        print("-----------------------")


    assert pixel_scale_2 > pixel_scale_1

if __name__ == "__main__":

    test_diameter(verbose=True)
    test_diff_limit(verbose=True)
    test_planck(verbose=True)
    test_save_instrument(verbose=True)
