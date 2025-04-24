
from syotools.models.camera import Camera
from syotools.models.telescope import Telescope
from astropy import units as u
import pytest

def test_spatial_resolution():
    tel: Telescope = Telescope()
    tel.set_from_sei('EAC1')
    expected = (1.22 * u.rad * 500 * u.nm / tel.effective_aperture.to(u.nm)).to(u.arcsec)
    actual = tel.diffraction_limit(500 * u.nm).value
    assert actual == pytest.approx(expected.value)
