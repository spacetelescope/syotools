import numpy as np
import astropy.units as u

from photutils.geometry import rectangular_overlap_grid

from syotools.models.multispec import MultiSpec
from syotools.models.source_exposure import SourceIFSExposure

class IFS(MultiSpec):

    def extraction_mask(self, x, y, band, xsamp, ysamp, extraction_aperture):
        """
        Draw an extraction mask.
        For IFUs, this is a spaxel-wide slit

        Parameters
        ----------
        mask : np.ndarray
            a 2D mask that draws the extraction aperture
        """
        wave = band["effective_wavelength"]
        if "image_slicer" in self.configuration:
            if extraction_aperture is None or np.isclose(extraction_aperture, 0*u.arcsec):
                height = 3 * self.fwhm_psf(wave).to_value(u.arcsec)
                width = self.configuration["image_slicer"]["spaxel_angle"].to_value(u.arcsec)
            else:
                width = self.configuration["image_slicer"]["spaxel_angle"].to_value(u.arcsec)
                height = extraction_aperture.to_value(u.arcsec) * 2 # because it's a half-height
        else:
            raise ValueError("Incomplete instrument slit specification.")
        
        mask = rectangular_overlap_grid(np.min(x), np.max(x), np.min(y), np.max(y), x.shape[1], y.shape[0], width, height, 0, 0, 4)

        return mask, height

    def create_exposure(self, source=None):
        new_exposure = SourceIFSExposure()
        if source is not None:
            new_exposure.source = source
        self.add_exposure(new_exposure)
        return new_exposure
