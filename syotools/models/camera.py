#!/usr/bin/env python
"""
Created on Fri Oct 14 21:31:18 2016
@author: gkanarek, tumlinson
"""
import numpy as np
import astropy.constants as const
import astropy.units as u

from syotools.models.base import PersistentModel
from syotools.models.source_exposure import SourcePhotometricExposure
from syotools.defaults import default_camera
from syotools.spectra.utils import mag_from_sed
from hwo_sci_eng.utils import read_yaml 

def nice_print(arr):
    """
    Utility to make the verbose output more readable.
    """

    if isinstance(arr, u.Quantity):
        l = ['{:.2f}'.format(i) for i in arr.value]
    else:
        l = ['{:.2f}'.format(i) for i in arr]
    return ', '.join(l)

class Camera(PersistentModel):
    """
    The basic camera class, which provides parameter storage for
    optimization.

    Attributes: #adapted from the original in Telescope.py
        telescope    - the Telescope object associated with the camera
        exposures    - the list of Exposures taken with this camera

        name         - name of the camera (string)
        n_bands      - number of wavelength bands (int)
        n_channels   - number of channels (int)
        pivotwave    - central wavelengths for bands, in nanometers (float array)
        bandnames    - names of bands (string list)
        channels     - grouping of bands into channels [UV, Optical, IR],
                       and indicating the reference band for pixel size (list of tuples)
        fiducials    - fiducial wavelength of the band, for reference (float array)
        total_qe     - total quantum efficiency in each band (float array)
        ap_corr      - magnitude correction factor for aperture size (float array)
        bandpass_r   - resolution in each bandpass (float array)
        dark_current - dark current values in each band (float array)
        detector_rn  - read noise for the detector in each band (float array)
        sky_sigma    - sky background emission (float array)

        _default_model - used by PersistentModel

    The following are attributes I haven't included, and the justification:
        R_effective - this doesn't seem to be used anywhere
    """

    def __init__(self, **kw):

        self.telescope = None
        self.exposures = []
        self.name = ''
        self.pivotwave = np.zeros(1, dtype=float) * u.nm
        self.bandnames = ['']
        self.channels = [([],0)]
        self.fiducials = np.zeros(1, dtype=float) * u.nm
        self.total_qe = np.zeros(1, dtype=float) * u.dimensionless_unscaled
        self.ap_corr = np.zeros(1, dtype=float) * u.dimensionless_unscaled
        self.bandpass_r = np.zeros(1, dtype=float) * u.dimensionless_unscaled
        self.dark_current = np.zeros(1, dtype=float) * (u.electron / u.s / u.pixel)
        self.detector_rn = np.zeros(1, dtype=float) * (u.electron / u.pixel)**0.5
        self.sky_sigma = np.zeros(1, dtype=float) * u.dimensionless_unscaled
        super().__init__(default_camera, **kw)

    @property
    def pixel_size(self):
        """
        Compute the pixel size as a function of pivot wavelength.

        Use the reference band for each channel as the fiducial: U-band for UV
        and optical, J-band for IR.
        """

        pixsize = np.zeros(self.n_bands, dtype=float)
        
        fiducials, effective_aperture = self.recover('fiducials', 'telescope.effective_aperture')
        
        for ref, bands in enumerate(self.channels):
            pxs = (0.5 * fiducials[0][ref] * u.Unit(fiducials[1])* u.rad / effective_aperture).to(u.arcsec).value
            pixsize[bands] = pxs

        return pixsize * u.arcsec / u.pix

    @property
    def n_bands(self):
        return len(self.bandnames)

    @property
    def n_channels(self):
        return len(self.channels)

    @property
    def derived_bandpass(self):
        """
        Calculate the bandpasses.
        """

        #Convert to Quantity for calculations.
        pivotwave, bandpass_r = self.recover('pivotwave','bandpass_r')

        return np.array(pivotwave[0]) / np.array(bandpass_r[0])

    @property
    def ab_zeropoint(self):
        """
        AB-magnitude zero points as per Marc Postman's equation.
        """
        pivotwave = self.pivotwave[0] * u.nm
        abzp = 5509900. * (u.photon / u.s / u.cm**2) / pivotwave

        return abzp


    @property
    def fwhm_psf(self):
        """
        Calculate the FWHM of the camera's PSF.
        """
        #Convert to Quantity for calculations.
        pivotwave, effective_aperture = self.recover('pivotwave', 'telescope.effective_aperture')
        diff_limit, diff_fwhm = self.recover('telescope.diff_limit_wavelength',
                                             'telescope.diff_limit_fwhm')

        #fwhm = (1.22 * u.rad * pivotwave / aperture).to(u.arcsec)
        pivots = pivotwave[0] * u.Unit(pivotwave[1]) 
        fwhm = (1.03 * u.rad * pivots / effective_aperture).to(u.arcsec)
        
        #only use these values where the wavelength is greater than the diffraction limit
        fwhm = np.where(pivots.value > diff_limit[0], fwhm.value, diff_fwhm.value) * u.arcsec

        return fwhm

    def _print_initcon(self, verbose):
        if verbose: #These are our initial conditions
            print('Telescope diffraction limit: {}'.format(self.telescope.diff_limit_wavelength))
            print('Telescope effective_aperture: {}'.format(self.telescope.effective_aperture))
            print('Telescope temperature: {}'.format(self.telescope.temperature))
            print('Pivot waves: {}'.format(nice_print(self.pivotwave[0] * u.Unit(self.pivotwave[1]))))
            print('Pixel sizes: {}'.format(nice_print(self.pixel_size)))
            print('AB mag zero points: {}'.format(nice_print(self.ab_zeropoint)))
            print('Quantum efficiency: {}'.format(nice_print(self.total_qe[0] * u.Unit(self.total_qe[1]))))
            print('Aperture correction: {}'.format(nice_print(self.ap_corr[0] * u.Unit(self.ap_corr[1]))))
            print('Bandpass resolution: {}'.format(nice_print(self.bandpass_r[0] * u.Unit(self.bandpass_r[1]))))
            print('Derived_bandpass: {}'.format(nice_print(self.derived_bandpass)))
            print('Detector read noise: {}'.format(nice_print(self.detector_rn[0] * u.Unit(self.detector_rn[1]))))
            print('Dark rate: {}'.format(nice_print(self.dark_current[0] * u.Unit(self.dark_current[1]))))

    def _fsky(self, verbose=True):
        """
        Calculate the sky flux as per Eq 6 in the SNR equation paper.
        """

        (f0, D, dlam, Phi, fwhm, Sigma) = self.recover('ab_zeropoint',
                'telescope.effective_aperture', 'derived_bandpass',
                'pixel_size', 'fwhm_psf', 'sky_sigma')

        D = D.to(u.cm)
        m = 10.**(-0.4 * np.array(Sigma[0])) / u.arcsec**2
        Npix = self._sn_box(False)

        if verbose:
            print('Sky brightness: {}'.format(nice_print(Sigma[0])))

        fsky = f0 * np.pi / 4. * D**2 * (dlam*u.nm) * m * (Phi**2 * Npix) * u.pix

        return fsky

    def _sn_box(self, verbose):
        """
        Calculate the number of pixels in the SNR computation box.
        """

        (Phi, fwhm_psf) = self.recover('pixel_size', 'fwhm_psf')
        sn_box = np.round(3. * fwhm_psf / Phi)

        if verbose:
            print('PSF width: {}'.format(nice_print(fwhm_psf)))
            print('SN box width: {}'.format(nice_print(sn_box)))

        return sn_box**2 / u.pix #don't want pix**2 units

    def c_thermal(self, verbose=True):
        """
        Calculate the thermal emission counts for the telescope.
        """

        #Convert to Quantities for calculation.
        (bandpass, pivotwave, aperture, ota_emissivity,
         total_qe, pixel_size) = self.recover('derived_bandpass', 'pivotwave',
                'telescope.effective_aperture',  'telescope.ota_emissivity',
                'total_qe', 'pixel_size')

        box = self._sn_box(verbose)

        bandwidth = (bandpass * u.nm).to(u.cm)

        h = const.h.to(u.erg * u.s) # Planck's constant erg s
        c = const.c.to(u.cm / u.s) # speed of light [cm / s]

        pivots = pivotwave[0] * u.Unit(pivotwave[1])
        energy_per_photon = h * c / pivots.to(u.cm) / u.ph

        D = aperture.to(u.cm) # telescope diameter in cm

        Omega = (pixel_size**2 * box * u.pix).to(u.sr)

        planck = self.planck
        QE = total_qe[0] * u.Unit(total_qe[1])
        qepephot = QE * planck / energy_per_photon

        if verbose:
            print('Planck spectrum: {}'.format(nice_print(planck)))
            print('QE * Planck / E_phot: {}'.format(nice_print(qepephot)))
            print('E_phot: {}'.format(nice_print(energy_per_photon)))
            print('Omega: {}'.format(nice_print(Omega)))

        thermal = (ota_emissivity[0] * planck / energy_per_photon *
    			(np.pi / 4. * D**2) * QE * Omega * bandwidth )

        return thermal

    @property
    def planck(self):
        """
        Planck spectrum for the various wave bands.
        """
        #Convert to Quantities for calculation
        pivotwave, temperature = self.recover('pivotwave', 'telescope.temperature')

        pivots = pivotwave[0] * u.Unit(pivotwave[1])
        wave = pivots.to('cm')
        if isinstance(temperature, u.Quantity):
            temp = temperature
        else:
            temps = temperature[0] * u.Unit(temperature[1])
            temp = temps.to('K')
        h = const.h.to(u.erg * u.s) # Planck's constant erg s 
        c = const.c.to(u.cm / u.s) # speed of light [cm / s] 
        k = const.k_B.to(u.erg / u.K) # Boltzmann's constant [erg deg K^-1] 
        x = 2. * h * c**2 / wave**5 
        exponent = (h * c / (wave * k * temp)) 
    
        result = (x / (np.exp(exponent)-1.)).to(u.erg / u.s / u.cm**3) / u.sr

        return result

    def interpolate_at_bands(self, sed):
        """
        Interpolate an SED to obtain magnitudes for the camera's wavebands.
        """
        return mag_from_sed(sed, self)

    def interpolate_source_at_bands(self, source):
        """
        Interpolate an SED to obtain magnitudes for the camera's wavebands.
        """
        return mag_from_source(self, source)

    def create_exposure(self, source=None):
        new_exposure = SourcePhotometricExposure()
        if source is not None:
            new_exposure.source = source
        self.add_exposure(new_exposure)
        return new_exposure

    def add_exposure(self, exposure):
        self.exposures.append(exposure)
        exposure.camera = self
        exposure.telescope = self.telescope
        exposure.calculate()
            
    def set_from_sei(self, name): 

        if ('HRI' in name): hri = read_yaml.hri()
        
        # the "hri" dictionary returned by read_yaml is nested, and therefore awkward  
        # when summoning individual entries. And often, we do not need the individual 
        # components. So, we are going to break this dictionary up and carry the 
        # pieces separately:  

        self.UVIS = hri['UVIS']  

        self.NIR = hri['NIR'] 


