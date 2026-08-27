#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 30 12:31:11 2017
@author: gkanarek, jt
"""
import numpy as np
from scipy.interpolate import interp1d
import astropy.units as u
import astropy.constants as const

import synphot as syn
from synphot.models import Empirical1D

from syotools.models.base import PersistentModel

from syotools.defaults import default_exposure
from syotools.models.source import Source

SPECTRAL_RADIANCE = u.W / (u.m**2 * u.sr * u.um)
PHOTON_SPECTRAL_RADIANCE = u.photon / (u.cm**2 * u.s * u.nm * u.arcsec**2)
SPECTRAL_RADIANCE_CGS = u.erg / (u.s * u.cm**2 * u.arcsec**2 * u.nm)

class SourceExposure(PersistentModel):
    """
    The base Source exposure class, which provides parameter storage for
    optimization, and all exposure-specific calculations. Since the
    Nov 2024 refactor, this class uses the Source object to handle
    astrophysical source information. Also, all JSON encoding has been
    stripped out.

    The SNR, exptime, and limiting magnitude can each be calculated from the
    other two. To trigger such calculations when parameters are updated, we
    will need to create property setters.

    Attributes:
        telescope    - the Telescope model instance associated with this exposure
        camera       - the Camera model instance associated with this exposure
        spectrograph - the Spectrograph model instance (if applicable) associated
                       with this exposure
        ifs          - the IFS model instance (if applicable) associated with this exposure

        exp_id       - a unique exposure ID, used for save/load purposes (string)
                        NOTE: THIS HAS NO DEFAULT, A NEW EXP_ID IS CREATED
                        WHENEVER A NEW CALCULATION IS SAVED.
        source       - Source object that this Exposure will observe
        n_exp        - the desired number of exposures (integer)
        exptime      - the desired exposure time (float array)
        snr          - the desired S/N ratio (float array)
        magnitude    - either the input source magnitude, in which case this is
                       equal to the SED interpolated to the desired wavelengths,
                       or the limiting magnitude of the exposure (float array)
        unknown      - a flag to indicate which variable should be calculated
                       ('snr', 'exptime', or 'magnitude'). this should generally
                       be set by the tool, and not be available to users. (string)
        sources      - list of source objects to be added to this exposure

        _default_model - used by PersistentModel
    """

    def __init__(self, default_model=default_exposure, **kw):

        self.source = Source() # this is the Source object, returns a flat spectrum by default.
                        # currently an Exposure can have only one Source

        self.telescope = None
        self.instrument = None

        self.exp_id = ''
        self.n_exp = 1
        self._exptime = np.ones(1, dtype=float) * u.h
        self._snr = np.zeros(1, dtype=float)
        self._magnitude = np.zeros(1, dtype=float) * u.ABmag
        self._unknown = "" # one of 'snr', 'magnitude', 'exptime'
        self._interp_flux = np.zeros(1, dtype=float) * u.dimensionless_unscaled # the source SED interpolated to the Spectrograph wavelength grid

        self.verbose = False # set this to True for debugging purposes
        self._disable = True #set this to disable recalculating (when updating several attributes at the same time)
        #super().__init__(default_model, **kw)

    def disable(self):
        self._disable = True

    def enable(self):
        self._disable = False
        self.calculate()

    #Property wrappers for the three possible unknowns, so that we can auto-
    #calculate whenever they're set, and to prevent overwriting previous
    #calculations by accident.

    @property
    def unknown(self):
        return self._unknown

    @unknown.setter
    def unknown(self, new_unknown):
        valid_unknowns = ("exptime", "snr", "magnitude")
        if new_unknown in valid_unknowns:
            self._unknown = new_unknown
            self.enable() # once this is set, enable calculation (which immediately runs one)
        else:
            raise KeyError(f"Cannot solve for {new_unknown}, unrecognized unknown.")

    def _ensure_array(self, quant, nb=None):
        """
        Ensure that the given Quantity is an array, propagating if necessary.
        """
        if self.instrument is None:
            nb = 1
        elif nb is None:
            nb = self.recover('instrument.n_bands')
        val = quant 
        if isinstance(quant, (int, float)):
            q = np.full(nb, val)
        if isinstance(quant, (u.Quantity)) and quant.isscalar:
            q = np.full(nb, val)
        elif len(quant) < nb:
            if len(quant) > 1:
                q = np.full(nb, val[0])
            else:
                q = np.full(nb, val)
        elif len(quant) > nb:
            q = quant[0:nb]
        else:
            q = val

        if not isinstance(q, u.Quantity):
            q = q << quant.unit

        return q

    def _ensure_quantity(self, quant, unit, nb=None):
        """
        Ensure given quantity is an astropy unit.Quantity
        of appropriate type
        """
        if isinstance(quant, u.Quantity):
            # just see if this crashes.
            try:
                quant.to(unit)
            except:
                raise ValueError(f"Quantity {quant} unit is not convertible to {unit}.")
        else:
            quant = quant << unit
        quant = self._ensure_array(quant, nb=nb)
        return quant

    @property
    def exptime(self):
        return self._exptime

    @exptime.setter
    def exptime(self, new_exptime):
        if self.unknown == "exptime":
            return
        self._exptime = self._ensure_quantity(new_exptime, u.s)
        self.calculate()

    @property
    def snr(self):
        return self._snr

    @snr.setter
    def snr(self, new_snr):
        if self.unknown == "snr":
            return
        self._snr = self._ensure_quantity(new_snr, u.dimensionless_unscaled)
        self.calculate()

    @property 
    def magnitude(self):
        return self._magnitude

    @magnitude.setter 
    def magnitude(self, new_magnitude):
        if self.unknown == "magnitude":
            return
        self._magnitude = self._ensure_quantity(new_magnitude, u.ABmag)
        self.calculate()

    # Code from pyEDITH, astrophysical_scene.py calc_zodi_flux
    # Courtesy of Eleonora Alei
    def calc_zodi_flux(
        self,
        wave: u.Quantity,
        sn_box: u.Quantity,
        pixel_scale: u.Quantity,
        # starshade: bool = False,
        # ss_elongation: u.Quantity = None,
    ) -> u.Quantity:
        """

        Calculate the zodiacal light flux for given celestial coordinates and wavelengths.

        This function computes the zodiacal light flux based on the target's position in the sky,
        observation wavelengths, and whether a starshade is used. It uses the model from
        Leinert et al. (1998) to calculate the zodiacal light intensity.

        Parameters
        ----------
        dec : Quantity
            Declination of targets in degrees (J2000 equatorial coordinate).
        ra : Quantity
            Right ascension of targets in degrees (J2000 equatorial coordinate).
        wave : Quantity
            Wavelengths in microns (vector of length nlambda).
        F0 : Quantity
            Flux zero points at wavelengths wave (vector of length nlambda).

        Returns
        -------
        np.ndarray
            Zodi surface brightness in units of photons s^-1 cm^-2 arcsec^-2 nm^-1 / F0.
            This is equivalent to 10^(-0.4*magOmega_ZL).
            - Multiply by F0 to get photons s^-1 cm^-2 arcsec^-2 nm^-1
            - Multiply by energy of photons to get erg s^-1 cm^-2 arcsec^-2 nm^-1
            The output array has dimensions (nlambda, nstars).

        Raises
        ------
        ValueError
            If F0 and lambd have different lengths, or if starshade mode is inconsistent with ss_elongation.

        Note
        ----
        - The function uses the zodiacal light model from Leinert et al. (1998).
        - For coronagraph mode, it assumes observations near solar longitude of 135 degrees.
        - Starshade functionality is currently not fully implemented.

        References:

        Leinert, C., et al. (1998). The 1997 reference of diffuse night sky brightness.
        Astronomy and Astrophysics Supplement Series, 127(1), 1-99.
        """

        # if starshade and ss_elongation is None:
        #     raise ValueError(
        #         "ERROR. You have set the STARSHADE flag. Must specify SS_ELONGATION in degrees."
        #     )
        # if not starshade and ss_elongation is not None:
        #     raise ValueError(
        #         "ERROR. You must enable STARSHADE mode if you are setting SS_ELONGATION in degrees."
        #     )

        F0 = 5509900. * (u.photon / u.s / u.cm**2) / wave

        # This code is in nanometers
        wave = wave.to(u.nm)

        # Convert equatorial coordinates to ecliptic coordinates
        coords = self.source.coords
        ecl_coords = coords.barycentrictrueecliptic
        beta = ecl_coords.lat.rad

        # all we need is the sine of the latitude
        # Use absolute value of the sin of beta (symmetry about ecliptic plane)
        sinbeta = np.abs(np.sin(beta))

        # SOURCE: Leinert et al. (1998)
        # Define solar longitude and beta values for interpolation
        beta_vector = np.array([0.0, 5, 10, 15, 20, 25, 30, 45, 60, 75]) * u.deg
        sollong_vector = (
            np.array(
                [
                    0,
                    5,
                    10,
                    15,
                    20,
                    25,
                    30,
                    35,
                    40,
                    45,
                    60,
                    75,
                    90,
                    105,
                    120,
                    135,
                    150,
                    165,
                    180.0,
                ]
            )
            * u.deg
        )

        # Table 17 values (assumed to be in some brightness units)
        table17 = (
            np.array(
                [
                    [-1, -1, -1, 3140, 1610, 985, 640, 275, 150, 100],
                    [-1, -1, -1, 2940, 1540, 945, 625, 271, 150, 100],
                    [-1, -1, 4740, 2470, 1370, 865, 590, 264, 148, 100],
                    [11500, 6780, 3440, 1860, 1110, 755, 525, 251, 146, 100],
                    [6400, 4480, 2410, 1410, 910, 635, 454, 237, 141, 99],
                    [3840, 2830, 1730, 1100, 749, 545, 410, 223, 136, 97],
                    [2480, 1870, 1220, 845, 615, 467, 365, 207, 131, 95],
                    [1650, 1270, 910, 680, 510, 397, 320, 193, 125, 93],
                    [1180, 940, 700, 530, 416, 338, 282, 179, 120, 92],
                    [910, 730, 555, 442, 356, 292, 250, 166, 116, 90],
                    [505, 442, 352, 292, 243, 209, 183, 134, 104, 86],
                    [338, 317, 269, 227, 196, 172, 151, 116, 93, 82],
                    [259, 251, 225, 193, 166, 147, 132, 104, 86, 79],
                    [212, 210, 197, 170, 150, 133, 119, 96, 82, 77],
                    [188, 186, 177, 154, 138, 125, 113, 90, 77, 74],
                    [179, 178, 166, 147, 134, 122, 110, 90, 77, 73],
                    [179, 178, 165, 148, 137, 127, 116, 96, 79, 72],
                    [196, 192, 179, 165, 151, 141, 131, 104, 82, 72],
                    [230, 212, 195, 178, 163, 148, 134, 105, 83, 72],
                ]
            )
            * SPECTRAL_RADIANCE
        )
        # For coronagraph, assume observations near solar longitude of 135 degrees
        j = np.argmin(np.abs(sollong_vector - 135 * u.deg))
        k = np.argmin(np.abs(sollong_vector - 90 * u.deg))

        # Interpolate to get zodi brightness factor


        interp = interp1d(
            np.sin(beta_vector),
            table17[j] / table17[k, 0],
            kind="cubic",
            fill_value="extrapolate",
        )
        # this specifically selects the 135 degree longitude, and interpolates to the chosen ecliptic latitude
        f = interp(sinbeta) * u.dimensionless_unscaled

        # Wavelength dependence (fits to Table 19 in Leinert et al 1998)
        zodi_lambd = (
            np.array(
                [
                    0.2,
                    0.3,
                    0.4,
                    0.5,
                    0.7,
                    0.9,
                    1.0,
                    1.2,
                    2.2,
                    3.5,
                    4.8,
                    12,
                    25,
                    60,
                    100,
                    140,
                ]
            )
            * u.micron
        )
        zodi_blambd = (
            np.array(
                [
                    2.5e-8,
                    5.3e-7,
                    2.2e-6,
                    2.6e-6,
                    2.0e-6,
                    1.3e-6,
                    1.2e-6,
                    8.1e-7,
                    1.7e-7,
                    5.2e-8,
                    1.2e-7,
                    7.5e-7,
                    3.2e-7,
                    1.8e-8,
                    3.2e-9,
                    6.9e-10,
                ]
            )
            * SPECTRAL_RADIANCE
        )

        # Convert to erg s^-1 cm^-2 arcsec^-2 angstrom^-1
        zodi_blambd = zodi_blambd.to(SPECTRAL_RADIANCE_CGS)
        zodi_lambd = zodi_lambd.to(u.nm)

        # SYOTools actually needs the zodi in flux units (syn.units.PHOTLAM)
        # so we do not need to call out for a magnitude calculation here.
        interp = interp1d(zodi_lambd.value, zodi_blambd.value, kind="cubic", bounds_error=False, fill_value=0.0)
        # flux is an array in SPECTRAL_RADIANCE_CGS units
        flux = interp(wave) << SPECTRAL_RADIANCE_CGS
        flux_zodi = f * flux # apply the scaling relative to 90 degrees ecliptic (f)

        # # Interpolate to get zodi brightness at desired wavelengths
        # interp = interp1d(
        #     np.log10(zodi_lambd.value), np.log10(zodi_blambd.value), kind="cubic"
        # )
        # blambd = 10 ** interp(np.log10(wave.to_value(u.nm))) * zodi_blambd.unit

        # # Convert to photon flux
        # # I90fabsfco = blambd / (u.h * u.c / lambd)

        # I90fabsfco = blambd.to(
        #     PHOTON_SPECTRAL_RADIANCE,
        #     equivalencies=u.spectral_density(wave),
        # )
        # # Divide by F0
        # I90fabsfco = I90fabsfco / F0

        # # Calculate final zodi flux
        # nlambda = len(wave)
        # flux_zodi = f * I90fabsfco

        # omega is the size of the extraction box in steradians
        Omega = (pixel_scale**2 * sn_box).to(u.sr)
        flux_zodi *= Omega

        # from matplotlib import pyplot as plt
        # print(flux)
        # print(f)
        # print(flux_zodi)
        # print("Omega", Omega)
        # print(zodi_blambd)
        # print(zodi_lambd.to_value(u.nm))
        # print(wave)
        # plt.plot(wave, flux_zodi)
        # plt.plot(zodi_lambd.to_value(u.nm), zodi_blambd)
        # plt.show()

        # now convert to PHOTLAM
        sky = syn.spectrum.SourceSpectrum(Empirical1D, points=wave, lookup_table=syn.units.convert_flux(wave, flux_zodi, syn.units.PHOTLAM))

        return sky  # 1/arcsec^2 (UNITS OF SPECTRAL RADIANCE) - original, now PHOTLAM

    @property
    def interpolated_sed(self):
        """
        The exposure's (old style) SED interpolated at the camera bandpasses.
        """
        if not self.camera:
            return self.sed
        sed = self.recover('sed')
        return self.camera.interpolate_at_bands(sed)

    def process_observation(self, source, band, verbose=False):
        """
        Process the entire observation, up through the point we compute SNR/Exptime/Mag
        
        The components of flux in the observation are: 
        1. The source (assumed to be a point, but let's give it a size 
           parameter) 
        2. Sky background (assumed uniform across the aperture) 
        3. Thermal self-emission (assumed uniform across the aperture).
        At the moment we only model the heat of the detector itself

        4. Dark current (assumed uniform across the aperture).
        This is the additional current flowing regardless of photons hitting
        the detector. It doesn't care about the detector QE or filter wheel.

        5. Read noise (assumed uniform across the aperture)
        The previous terms were all signal that accumulates with time. Read
        noise is the uncertainty introduced by the detector readout process 
        itself; a fixed value per exposure.

        Once we've computed all of these values, we can proceed to the
        exposure time/SNR/magnitude calculations.

        At that point, the difference between imaging and spectroscopy matter.
        
        For imaging:
        * All non-spatially-uniform components have (size * psf size) 
        compared to (aperture size), light losses adjusted accordingly, and 
        integrated over the bandpass + QE (source, sky) or QE (thermal) to be 
        single value(s)
        * All uniform components processed for the aperture size

        For spectroscopy:
        * All non-spatially-uniform components are convolved with a 
        response function equal to the resolving power of the instrument
        and then have their (size * psf size) compared to slit size 
        (width * height, if applicable), light losses adjusted accordingly, 
        and convolved with the bandpass+QE (source, sky) or QE (thermal), 
        then convolved with a response function equal to the resolving power 
        of the instrument.
        * All uniform components processed for the height of the slit * 
        resolving power.
        """

        configuration, c_thermal, _sn_box, transform_flux = self.recover("instrument.configuration", "instrument._c_thermal", "instrument._sn_box", "instrument.transform_flux")
        pixel_scale = configuration["pixel_scale"]
        for detector in configuration["detector"]:
            dark_current = configuration["detector"]["dark_current"]
            qe = configuration["detector"]["total_qe"]
            read_noise = configuration["detector"]["read_noise"]

        if band["kind"] in ("disperser", "ifs"):
            R = band["resolution"]
            waveunit = band["bandpass"].waveset.unit
            wavepix = np.linspace(band["bandpass"].waveset[0], band["bandpass"].waveset[-1], 1000) # using the bandpass wavelengths leads to weird fringing
            delta_lambda = wavepix/R
            pixel = np.cumsum(1.0 / delta_lambda * np.gradient(wavepix))
            pixel_integer = np.arange(int(pixel[0]), int(pixel[-1]))
            wave = np.interp(pixel_integer, pixel, wavepix) << waveunit
            # Or just use the instrument bandpass
            #wave = band["bandpass"].waveset

            dw = wave[1:] - wave[:-1]
            good = np.where(dw != 0)[0]
            wave = wave[good]
            dw = dw[good]
            #dw = np.append(dw, dw[-1])
        else:
            dw = 1
            wave = source.sed.waveset
        syn.utils.validate_wavelengths(wave)
        self.wave = wave

        # set up an appropriately sized aperture
        sn_box = _sn_box(self.wave, False)

        sn_box = np.median(sn_box)

        # fsource is:
        # shaped
        # goes through the full optical path + QE
        # accumulates over time
        rel_area = 1
        # scale source radius to the aperture size - we get all of the flux if it's smaller than the aperture
        if source.radius > 0 * u.arcsec:
            area = np.pi * (source.radius/pixel_scale)**2
        else:
            area = np.pi * (np.median(self.instrument.fwhm_psf(self.wave))/pixel_scale)**2
        if area > sn_box:
            rel_area = (sn_box/area)
        flux_source = source.sed * rel_area

        sky = self.calc_zodi_flux(wave, sn_box, pixel_scale)

        print(sky)
        print(self.instrument.sky(wave))


        # fsky is:
        # uniform
        # goes through the full optical path QE
        # accumulates over time
        # Synphot doesn't like dividing a spectrum by an area unit. 
        # Rest assured, sky was supposed to be in ABMag/arcsec**2, so 
        # ABMag/arcsec**2 * pixels**2 * arcsec**2/pixel**2 is flux.
        flux_sky = sky * (sn_box * pixel_scale**2).value
        #print("Skyflux", flux_sky(flux_sky.waveset))

        # thermal is:
        # uniform
        # goes through the filter wheel and QE
        # accumulates over time
        thermal = c_thermal(self.wave)

        # print("Source", flux_source.waveset)
        # print("Sky", flux_sky.waveset)
        # print("Thermal", thermal.waveset)
        # total_band = band["bandpass"] * qe
        # total_flux = total_band(total_band.waveset)
        # b1, b2 = total_band.waveset.min(), total_band.waveset.max()
        # a1, a2 = flux_source.waveset.min(), flux_source.waveset.max()
        # print("WAVE EDGES", a1, b1, a2, b2, a2 < b1, b2 < a1)
        # print("Disjoint", total_band.check_overlap(flux_source))
        # print("Valid band", total_band.waveset[total_flux > 0])
        # print("Band", (band["bandpass"]* qe).waveset)
        # print("QE", qe.waveset)



        # apply internal effects within telescope & instrument
        fsource = syn.observation.Observation(flux_source, band["bandpass"] * qe, binset=self.wave, force="taper")
        fsky = syn.observation.Observation(flux_sky, band["bandpass"] * qe, binset=self.wave, force="taper")
        self.thermal = syn.observation.Observation(thermal, band["bandpass"] * qe, binset=self.wave, force="taper")

        # dark is:
        # uniform
        # only within detector
        # accumulates over time
        dark = dark_current * sn_box# * qy

        # readnoise is:
        # uniform
        # only within detector
        # single event at read time
        read_noise = read_noise * np.sqrt(sn_box) # * u.electron**0.5 / u.pix**0.5

        fsource_countrate = transform_flux(fsource, wave) * dw
        fsky_countrate = transform_flux(fsky, wave) * dw
        thermal_countrate = transform_flux(self.thermal, wave) * dw
        if dw == 1:
            self.wave = band["bandpass"].pivot()

        return fsource_countrate, fsky_countrate, thermal_countrate, dark, read_noise

    def calculate(self, custom_band=None):
        """
        Wrapper to calculate the exposure time, SNR, or limiting magnitude,
        based on the other two. The "unknown" attribute controls which of these
        parameters is calculated.
        """
        if self._disable:
            return False
        if self.instrument is None or self.telescope is None:
            return False
        result = {'magnitude': self.calculate_magnitude,
                'exptime': self.calculate_exptime,
                'snr': self.calculate_snr}[self.unknown](custom_band=custom_band)

        return result

    def calculate_exptime(self, custom_band=None):
        """
        Calculate for exposure times. If a band has been passed in, do that. Otherwise, do all of them in the channel.

        Parameters
        ----------
        band : _type_, optional
            _description_, by default None
        """
        configuration, band, all_bands = self.recover("instrument.configuration", "instrument.band", "instrument.bands")
        if custom_band is not None:
            bands = [custom_band]
        else:
            if band is None:
                bands = all_bands
            else:
                bands = [band]
        self._exptime = []
        _initial_band = self.instrument.band
        _snr_temp = self._ensure_quantity(self._snr, u.dimensionless_unscaled, len(bands))
        for idx, band in enumerate(bands):
            # because a multiple-in, multiple-out is a valid use case
            self._snr = _snr_temp[idx]
            self.instrument.band = band
            result = self._update_exptime(self.source, configuration["band"][band])
            self._exptime.append(result)
        self._snr = _snr_temp
        self.instrument.band = _initial_band

        return True

    def calculate_snr(self, custom_band=None):
        """
        Calculate for SNR. If a band has been passed in, do that. Otherwise, do all of them in the channel.

        Parameters
        ----------
        band : _type_, optional
            _description_, by default None
        """
        configuration, band, all_bands = self.recover("instrument.configuration", "instrument.band", "instrument.bands")
        if custom_band is not None:
            bands = [custom_band]
        else:
            if band is None:
                bands = all_bands
            else:
                bands = [band]
        self._snr = []
        _initial_band = self.instrument.band
        _exptime_temp = self._ensure_quantity(self._exptime, u.s, len(bands))
        for idx, band in enumerate(bands):
            # because a multiple-in, multiple-out is a valid use case
            self._exptime = _exptime_temp[idx]
            self.instrument.band = band
            result = self._update_snr(self.source, configuration["band"][band])
            self._snr.append(result)
        self._exptime = _exptime_temp
        self.instrument.band = _initial_band

        return True

    def calculate_magnitude(self, custom_band=None):
        """
        Calculate for magnitudes. If a band has been passed in, do that. Otherwise, do all of them in the channel.

        Parameters
        ----------
        band : _type_, optional
            _description_, by default None
        """
        configuration, band, all_bands = self.recover("instrument.configuration", "instrument.band", "instrument.bands")
        if custom_band is not None:
            bands = [custom_band]
        else:
            if band is None:
                bands = all_bands
            else:
                bands = [band]
        self._magnitude = []
        _initial_band = self.instrument.band
        _exptime_temp = self._ensure_quantity(self._exptime, u.s, len(bands))
        _snr_temp = self._ensure_quantity(self._snr, u.dimensionless_unscaled, len(bands))
        for idx, band in enumerate(bands):
            # because a multiple-in, multiple-out is a valid use case
            self._exptime = _exptime_temp[idx]
            self._snr = _snr_temp[idx]
            self.instrument.band = band
            result = self._update_magnitude(self.source, configuration["band"][band])
            self._magnitude.append(result)

        self._exptime = _exptime_temp
        self._snr = _snr_temp
        self.instrument.band = _initial_band

        return True

    def _update_exptime(self, source, band):
        """
        Calculate the exposure time to achieve the desired S/N for the
        given SED.
        """
        self.instrument._print_initcon(self.verbose)

        (_snr, _nexp) = self.recover('_snr', 'n_exp')

        # all of these are now rates, in the extraction aperture (except read_noise)
        fsource_countrate, fsky_countrate, thermal_countrate, dark_current, read_noise = self.process_observation(source, band)

        snr2 = -(_snr**2)

        a = (fsource_countrate)**2
        b = snr2 * (fsource_countrate + (fsky_countrate + thermal_countrate + dark_current)) * u.ct
        c = snr2 * read_noise**2 * _nexp
        texp = ((-b + np.sqrt(b**2 - 4*a*c)) / (2*a)).to(u.s)

        if self.verbose:
            print("Fstar:", fsource_countrate)
            print("Texp:", texp)


        _exptime = texp

        return _exptime

    def _update_magnitude(self, source, band):
        """
        Calculate the limiting magnitude given the desired S/N and exposure
        time.
        """
        self.instrument._print_initcon(self.verbose)

        (_snr, _exptime, _nexp) = self.recover('snr', 'exptime', 'n_exp')
        effective_area = self.recover("telescope.effective_area")
        configuration, ab_zeropoint = self.recover("instrument.configuration", "instrument.ab_zeropoint")
        qe = configuration["detector"]["total_qe"]

        # all of these are now rates, in the extraction aperture (except read_noise)
        fsource_countrate, fsky_countrate, thermal_countrate, dark_current, read_noise = self.process_observation(source, band)

        read_noise /= u.ct**0.5
        _exptime = _exptime.to(u.s)

        snr2 = (_snr ** 2) * u.ct
        f0 = ab_zeropoint(band)

        fsky_counts = fsource_countrate * _exptime
        thermal_counts = fsource_countrate * _exptime
        dark_counts = dark_current * _exptime

        # Original equation is SNR = Sc / sqrt(Sc + Dc*Npix + Thermal*Npix + Sky*Npix + Rn**2*Nreads*Npix)
        # Rearrange: Sc/SNR = sqrt(Sc + Dc*Npix + Thermal*Npix + Sky*Npix + Rn**2*Nreads*Npix)
        # Square and collect terms of SC: 
        # Rearranged it becomes -Sc**2/SNR**2 + Sc**1 * 1 + Sc**0 * (Dc*Npix + Thermal*Npix + Sky*Npix + Rn**2*Nreads*Npix)
        # 
        # Our outputs from process_observation already have Npix applied, and read_noise is multiplied by the square root of Npix.
        a0 = -1 / snr2
        b0 = 1
        c0 = (dark_counts + thermal_counts + fsky_counts + read_noise**2 * _nexp)

        sc = (-b0 - np.sqrt(b0**2 - 4.0 * a0 * c0))/(2.0 * a0)
        # now get the source flux in counts.
        #sc = sc * _exptime
        # Convert counts back to photons
        phot_energy = const.h.to(u.erg * u.s) * const.c.to(u.cm / u.s) / band["bandpass"].pivot().to(u.cm) / u.ct
        photons = sc * phot_energy / effective_area
        flux = syn.units.convert_flux(band["bandpass"].pivot(), photons, syn.units.PHOTLAM, area=effective_area)
        fnu = syn.units.convert_flux(band["bandpass"].pivot(), photons, syn.units.FNU, area=effective_area)
        mag = -2.5*np.log10(flux.value / band["bandpass"].efficiency()/f0.value)
        print("Mag1", mag)
        fnu = fnu / band["bandpass"].efficiency()
        mag = -2.5*np.log10(fnu.value) + 8.90
        #mag = flux.to_value(u.ABmag)
        # Remove the impact of the bandpass
        print(sc, photons, flux, fnu)
        print("SNR", _snr)
        print("A0:", a0)
        print("B0:", b0)
        print("C0:", c0)
        print("Mag:", mag)

        # Convert to AB Magnitudes

        a0 = (_exptime)**2
        b0 = snr2 * _exptime
        c0 = snr2 * ((fsky_countrate + thermal_countrate + dark_current) * _exptime + (read_noise**2 * _nexp)) / u.ct
        k = (-b0 + np.sqrt(b0**2 - 4. * a0 * c0)) / (2. * a0)

        flux = (4. * k) / (f0 * effective_area)# * (band["bandpass"]*qe).equivwidth().to(u.nm))

        flux /= band["bandpass"].tlambda()

        _magnitude = -2.5 * np.log10(np.array(flux)) * u.mag('AB')

        return _magnitude

    def _update_snr(self, source, band):
        """
        Calculate the SNR for the given exposure time and source SED.
        """

        self.instrument._print_initcon(self.verbose)

        (_exptime, _nexp) = self.recover('_exptime', 'n_exp')

        # all of these are now rates, in the extraction aperture (except read_noise)
        fsource_countrate, fsky_countrate, thermal_countrate, dark_current, read_noise = self.process_observation(source, band)

        time_per_exposure = _exptime / _nexp

        signal_counts = (fsource_countrate * _exptime).to(u.ct)
        shot_noise_in_signal = np.sqrt(signal_counts)

        sky_counts = (fsky_countrate * _exptime).to(u.ct)
        shot_noise_in_sky = np.sqrt(sky_counts)

        read_counts = (read_noise**2 * _nexp / u.ct).to(u.ct)

        dark_counts = (dark_current * _exptime).to(u.ct)

        thermal_counts = (thermal_countrate * _exptime).to(u.ct)

        snr = signal_counts / np.sqrt(signal_counts + sky_counts + read_counts
                                      + dark_counts + thermal_counts)
        _snr = snr.value * u.dimensionless_unscaled

        if self.verbose:
            print('# of exposures: {}'.format(_nexp))
            print('Time per exposure: {}'.format(time_per_exposure))
            print('Signal counts: {}'.format(self.nice_print(signal_counts)))
            print('Signal shot noise: {}'.format(self.nice_print(shot_noise_in_signal)))
            print('Sky counts: {}'.format(self.nice_print(sky_counts)))
            print('Sky shot noise: {}'.format(self.nice_print(shot_noise_in_sky)))
            print('Total read noise: {}'.format(self.nice_print(read_counts)))
            print('Dark current noise: {}'.format(self.nice_print(dark_counts)))
            print('Thermal counts: {}'.format(self.nice_print(thermal_counts)))
            print('SNR: {}'.format(snr))
            
        return _snr

    def add_source(self, new_source):
        self.source = new_source

class SourcePhotometricExposure(SourceExposure):
    """ A subclass of the base Exposure model, for photometric ETC calculations """
    pass



class SourceSpectrographicExposure(SourceExposure):
    """
    A subclass of the base Exposure model, for spectroscopic ETC calculations.
    """

    def calculate_magnitude(self, custom_band=None):
        """
        Not supported, make this an error
        """
        raise ValueError("Magnitude calculation not supported for Spectroscopy")

class SourceIFSExposure(SourceExposure):
    """ 
    This is currently a subclass of Spectrographic exposure that accepts multiple
    sources and produces multiple returns. 
    """
    def __init__(self, default_model=default_exposure, **kw):

        # need this before so super().__init__ has somewhere to put the default source
        self.sources = []
        self._snrs = []
        self._exptimes = []
        self._wavelength = []

        super().__init__(default_model, **kw)
        # Do this after, because by default super().__init__ loads a default source
        self.sources = []

    def add_source(self, source):
        # and now the magic: create a master wavelength array from all of the sources.
        self.sources.append(source)
        for source in self.sources:
            self._wavelength = syn.utils.merge_wavelengths(self._wavelength, syn.models.get_waveset(source.sed.model))

    @property
    def num_sources(self):
        return len(self.sources)

    @property
    def source(self):
        return self.sources[-1]

    @source.setter
    def source(self, new_source):
        self.add_source(new_source)

    @property
    def exptimes(self):
        return self._exptimes

    @exptimes.setter
    def exptimes(self, new_exptime):
        print("Did not set exptimes")

    @property
    def snrs(self):
        return self._snrs

    @snrs.setter
    def snrs(self, new_snr):
        print("Did not set snrs")

    def calculate_exptime(self, custom_band=None):
        """
        Calculate for exposure times. If a band has been passed in, do that. 
        Otherwise, do all of the bands in the channel.


        Parameters
        ----------
        band : _type_, optional
            _description_, by default None
        """
        configuration, band, all_bands = self.recover("instrument.configuration", "instrument.band", "instrument.bands")
        if custom_band is not None:
            bands = [custom_band]
        else:
            if band is None:
                bands = all_bands
            else:
                bands = [band]
        self._exptime = []
        self._exptimes = []
        _snr_temp = self._ensure_array(self._snr, len(bands))
        # The unique thing about IFS is it has multiple sources
        for source in self.sources:
            _single_exptime = []
            for idx,band in enumerate(bands):
                # because a multiple-in, multiple-out is a valid use case
                self._snr = _snr_temp[idx]
                result = self._update_exptime(source, configuration["band"][band])
                _single_exptime.append(result)
            self._exptimes.append(_single_exptime)
        # find the highest exposure time amongst the set of sources
        self._exptime = np.max(self._exptimes,axis=0)
        
        self._snr = _snr_temp

        return True

    def calculate_snr(self, custom_band=None):
        """
        Calculate for SNR. If a band has been passed in, do that. 
        Otherwise, do all of the bands in the channel.

        Parameters
        ----------
        band : _type_, optional
            _description_, by default None
        """
        configuration, band, all_bands = self.recover("instrument.configuration", "instrument.band", "instrument.bands")
        if custom_band is not None:
            bands = [custom_band]
        else:
            if band is None:
                bands = all_bands
            else:
                bands = [band]
        self._snr = []
        self._snrs = []
        _exptime_temp =  self._ensure_array(self._exptime, len(bands))
        # The unique thing about IFS is it has multiple sources
        for source in self.sources:
            _single_snr = []
            for idx, band in enumerate(bands):
                # because a multiple-in, multiple-out is a valid use case
                self._exptime = _exptime_temp[idx]
                result = self._update_snr(source, configuration["band"][band])
                _single_snr.append(result)
            self._snrs.append(_single_snr)
        # find the highest exposure time amongst the set of sources
        self._snr = np.max(self._snrs,axis=0)

        self._exptime = _exptime_temp

        return True

    def calculate_magnitude(self, custom_band=None):
        """
        Not supported, make this an error
        """
        raise ValueError("Magnitude calculation not supported for IFS Spectroscopy")

class SourceCoronagraphicExposure(SourceExposure):
    """
    A subclass of the base Exposure model, for coronagraphic imaging calculations.
    """

    def calculate(self):
        """
        Wrapper to calculate the exposure time, SNR, or limiting magnitude,
        based on the other two. The "unknown" attribute controls which of these
        parameters is calculated.
        JT - THIS PART DOESNT WORK FOR CORON YET
        """
        if self._disable:
            return False
        if self.camera is None or self.telescope is None:
            return False
        status = {'magnitude': self._update_magnitude,
                  'exptime': self._update_exptime,
                  'snr': self._update_snr}[self.unknown]()
        return status

    #Calculation methods
    def _update_exptime(self):
        """
        Calculate the exposure time to achieve the desired S/N for the
        given SED.
        """
        print("Doesn't exist yet, pull it from camera class")

        return False #completed successfully

    def _update_magnitude(self):
        """
        Calculate the limiting magnitude given the desired S/N and exposure
        time.
        """

        print("Doesn't exist yet, pull it from camera class")

        return False #completed successfully

    def _update_snr(self):
        """
        Calculate the SNR for the given exposure time and planet properties.
        Follows Mennesson et al. 2024
        """

        self.instrument._print_initcon(self.verbose)

        print(' telescope inside the Coron exposure object ',
         self.telescope.effective_aperture)

        #serialize with JsonUnit for transportation
        self._snr = 10.

        return True #completed successfully
