import copy

from astropy import units as u
from astropy.modeling.functional_models import AiryDisk2D
from photutils.geometry import elliptical_overlap_grid, rectangular_overlap_grid
import scipy.special as sp
import numpy as np

from syotools.models.base import PersistentModel

MIN_CLIP = 1e-10

class Profile(PersistentModel):

    def __init__(self, telescope, instrument, geometry, wavelen):
        self.telescope = telescope
        self.instrument = instrument
        self.geometry = geometry
        self.wavelen = wavelen

        configuration = self.recover("instrument.configuration")
        self.pixel_scale = configuration["pixel_scale"].to_value(u.arcsec/u.pix)
        self.pix_area_sqarcsec = self.pixel_scale**2

    @property
    def pixelscale(self):
        """
        Return a per-pixel normalization factor for the appropriate area unit.

        Returns
        -------
        normfactor: float
            Normalization factor, unitless

        Raises
        ------
        ValueError
            Raised on invalid area unit
        """
        if self.geometry["surf_area_units"] in ['sr']:
            arcsec2 = u.arcsec * u.arcsec
            normfactor = self.pix_area_sqarcsec / u.sr.to(arcsec2)  # convert area in steradians to area in pixels
        elif self.geometry["surf_area_units"] in ['arcsec^2', None]: # 'None' should be an option because integrated flux
                                                        # shouldn't have units (internally, the grid is arcsec)
            normfactor = self.pix_area_sqarcsec
        else:
            msg = f"Unsupported surface area unit: {self.geometry['surf_area_units']}"
            raise ValueError(msg)

        return normfactor


    def generate_profile(self):
        """
        Make a 2D grid to add the profile to
        """


        pa_radians = (self.geometry.get("pa", 0) * u.deg).to(u.rad)

        xval = np.arange(-50,51,1) * self.pixel_scale
        yval = np.arange(-50,51,1) * self.pixel_scale

        xt,yt = np.meshgrid(xval,yval)

        self.x = xt * np.cos(pa_radians) + yt * np.sin(pa_radians)
        self.y = -xt * np.sin(pa_radians) + yt * np.cos(pa_radians)

        xsamp = ysamp = self.pixel_scale

        return self.x, self.y, xt, yt, xsamp, ysamp

    def point_profile(self):
        effective_diameter = self.recover("telescope.effective_diameter")

        Rz = 1.2196698912665045 * u.rad
        radius = (Rz * self.wavelen.to(u.AA)/effective_diameter.to(u.m))
        #print("Radius", radius)
        airymodel = AiryDisk2D(amplitude=1, x_0=0, y_0=0, radius=radius.to_value(u.arcsec))

        profile = airymodel(self.x,self.y)

        # # dist is in arcsec and actually an angle
        # dist = np.sqrt(x**2.0 + y**2.0)
        
        # print(effective_diameter)
        # x = 2*np.pi/wavelen.to_value(u.m) * effective_diameter/2.0 * np.sin(dist * u.arcsec)
        # print(x)
        # x = x.value
        # profile = (2 * sp.j1(x)/x)**2
        # # The Bessel Function of the first kind first order is 0 at r=0, so the middle is inf.
        # profile[50,50] = 1

        norm_method = self.geometry.get("norm_method", "integ_infinity")

        if norm_method == "surf_scale":
            norm_val = 1
        elif norm_method == "surf_center":
            norm_val = 1
        elif norm_method == "integ_infinity":
            #print("Profilesum", np.sum(profile))
            #print("Integration", ((4 * radius.to(u.arcsec)**2)/(np.pi * Rz.to(u.arcsec)**2)).to(u.dimensionless_unscaled)) # from Astropy
            norm_val = 1/np.sum(profile)

        profile = profile * norm_val

        return profile

    def sersic_profile(self):
        major = self.quant_to_val(self.geometry["major"], unit=u.arcsec)
        minor = self.quant_to_val(self.geometry["minor"], unit=u.arcsec)
        index = self.geometry["sersic_index"]

        # the actual value of b. Formula taken from astropy's sersic2d shape.
        b = sp.gammaincinv(2*index,0.5)

        dist = np.sqrt((self.x / major)**2.0 + (self.y / minor)**2.0)
        # This is Equation 1 of Graham & Driver (2005) 2005PASA...22..118G
        profile = np.exp( -b * (dist**(1.0 / index) - 1) )

        # Sersic profiles are highly centralized, so we need to oversample the central pixel
        # to get the appropriate flux. This is the difference between sampling
        # and integrating, and unfortunately we're sampling this function.
        dist = np.sqrt((self.x/101. / major)**2.0 + (self.y/101. / minor)**2.0)
        central_pixel = np.exp( -b * (dist**(1.0 / index) - 1) )

        profile[50,50] = np.sum(central_pixel) / 101**2

        if self.geometry["norm_method"] == "surf_scale":
            norm_val = self.pixelscale
        elif self.geometry["norm_method"] == "surf_center":
            norm_val = self.pixelscale * np.e**(-1*b)
        elif self.geometry["norm_method"] == "integ_infinity":
            # integrate the Sersic profile to get the total flux for normalization, including flux outside the FOV
            # http://ned.ipac.caltech.edu/level5/March05/Graham/Graham2.html
            integral = major * minor * 2 * np.pi * index * np.exp(b)/(b**(2*index))* sp.gamma(2 * index)
            norm_val = self.pixelscale / integral

        profile = profile * norm_val

        return profile

    def gaussian_profile(self):
        # The gaussian profile is actually a scale-sersic of index 0.5
        original_geometry = copy.deepcopy(self.geometry)

        self.geometry["major"] = self.quant_to_val(original_geometry["major"], unit=u.arcsec) * np.sqrt(2.0) # to match the usual definition of a Gaussian
        self.geometry["minor"] = self.quant_to_val(original_geometry["minor"], unit=u.arcsec) * np.sqrt(2.0) # to match the usual definition of a Gaussian
        self.geometry["shape"] = "sersic_scale"
        self.geometry["sersic_index"] = 0.5

        profile = self.sersic_scale_profile()

        self.geometry = original_geometry

        return profile

    def sersic_scale_profile(self):
        major = self.quant_to_val(self.geometry["major"], unit=u.arcsec)
        minor = self.quant_to_val(self.geometry["minor"], unit=u.arcsec)
        index = self.geometry["sersic_index"]

        dist = np.sqrt((self.x / major) ** 2.0 + (self.y / minor) ** 2.0)
        # This is Equation 14 of Graham & Driver (2005) 2005PASA...22..118G
        profile = np.exp(-dist ** (1.0 / index))

        # Sersic profiles are highly centralized, so we need to oversample the central pixel
        # to get the appropriate flux. This is the difference between sampling
        # and integrating, and unfortunately we're sampling this function.
        dist = np.sqrt((self.x/101. / major)**2.0 + (self.y/101. / minor)**2.0)
        central_pixel = np.exp(-dist ** (1.0 / index))

        profile[50,50] = np.sum(central_pixel) / 101**2

        if self.geometry["norm_method"] == "surf_scale":
            norm_val = self.pixelscale * np.e
        elif self.geometry["norm_method"] == "surf_center":
            norm_val = self.pixelscale
        elif self.geometry["norm_method"] == "integ_infinity":
            # integrate the Sersic profile to get the total flux for normalization, including flux outside the FOV
            # http://ned.ipac.caltech.edu/level5/March05/Graham/Graham2.html
            integral = major * minor * 2 * np.pi * index * sp.gamma(2 * index)
            norm_val = self.pixelscale / integral

        profile = profile * norm_val

        return profile

    def flat_profile(self):

        major = self.quant_to_val(self.geometry["major"], unit=u.arcsec)
        minor = self.quant_to_val(self.geometry["minor"], unit=u.arcsec)

        # We are not having elliptical_overlap_grid rotate the ellipse because we 
        # already had generate_profile rotate the coordinate system
        profile = elliptical_overlap_grid(np.min(self.x), np.max(self.x), np.min(self.y), np.max(self.y), self.x.shape[1], self.y.shape[0], major, minor, 0, 1, 1)

        # dist = np.sqrt((x / major) ** 2.0 + (y / minor) ** 2.0)

        # profile[dist < 1] = 1.0
        if self.geometry["norm_method"] in ("surf_scale", "surf_center"):
            norm_val = self.pixelscale
        elif self.geometry["norm_method"] == "integ_infinity":
            norm_val = self.pixelscale / np.sum(profile)

        profile = profile * norm_val

        return profile

    def power_profile(self):
        power_index = self.geometry['power_index']
        r_core = self.quant_to_val(self.geometry["r_core"], unit=u.arcsec)

        if power_index <= 0:
            raise ValueError('Power Law Index must be positive, not {}'.format(power_index))

        dist = np.sqrt((self.x/r_core)**2.0 + (self.y/r_core)**2.0).value
        profile = (dist.clip(MIN_CLIP, np.max(dist)))**(-1*power_index)
        # flatten the central portion. Everything within the core radius is set to 1.
        profile[np.where(dist <= 1.0)] = 1.0

        if self.geometry["norm_method"] in ["surf_scale", "surf_center"]:
            norm_val = self.pixelscale
        elif self.geometry["norm_method"] in ["integ_infinity"]:
            integral = np.pi * r_core**2 + 2 * np.pi * r_core**2/(power_index - 2)
            norm_val = self.pixelscale/integral

        profile = profile * norm_val

        return profile

