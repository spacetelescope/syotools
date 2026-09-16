from syotools.models.multispec import MultiSpec
from syotools.models.source_exposure import SourceMOSExposure

class MOS(MultiSpec):

    def extraction_mask(self, x, y, band, xsamp, ysamp, extraction_aperture):
        """
        Draw an extraction mask.
        For MOSes, this is the size of the microshutter

        Parameters
        ----------
        mask : np.ndarray
            a 2D mask that draws the extraction aperture
        """
        wave = band["effective_wavelength"]
        if "microshutter" in self.configuration:
            if extraction_aperture is not  None or extraction_aperture > 0*u.pix**2:
                warnings.warn("Ignoring extraction aperture size for microshutter array")
            height = self.configuration["microshutter"]["microshutter_height"].to_value(u.arcsec)
            width = self.configuration["microshutter"]["microshutter_width"].to_value(u.arcsec)
        else:
            raise ValueError("Incomplete instrument slit specification.")
        
        mask = rectangular_overlap_grid(np.min(x), np.max(x), np.min(y), np.max(y), x.shape[1], y.shape[0], width, height, 0, 0, 2)

        return mask, height

    def create_exposure(self, source=None):
        new_exposure = SourceMOSExposure()
        if source is not None:
            new_exposure.source = source
        self.add_exposure(new_exposure)
        return new_exposure
