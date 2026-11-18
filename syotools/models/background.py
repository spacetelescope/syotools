import warnings
import datetime

import numpy as np
from scipy.interpolate import CubicSpline, RegularGridInterpolator
import astropy.units as u
import astropy.constants as const
import synphot as syn
from synphot.models import Empirical1D, ConstFlux1D

from astropy.utils.exceptions import AstropyUserWarning

SPECTRAL_RADIANCE = u.W / (u.m**2 * u.sr * u.um)
PHOTON_SPECTRAL_RADIANCE = u.photon / (u.cm**2 * u.s * u.nm * u.arcsec**2)
SPECTRAL_RADIANCE_CGS = u.erg / (u.s * u.cm**2 * u.arcsec**2 * u.nm)

ZODI_INTERP = None
FLUX_INTERP = None

def load_leinert():

    """
    Load the wavelength array from Leinert et al. (1998)
    """
    global ZODI_INTERP
    global FLUX_INTERP

    #####
    # 1. The zodi strength
    #####

    # SOURCE: Leinert et al. (1998)
    # Define solar longitude and beta values for interpolation
    lat_vector = np.array([0.0, 5, 10, 15, 20, 25, 30, 45, 60, 75]) * u.deg
    long_vector = (
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
                [-1, -1, -1, 3140, 1610, 985, 640, 275, 150, 100], # 0
                [-1, -1, -1, 2940, 1540, 945, 625, 271, 150, 100], # 5
                [-1, -1, 4740, 2470, 1370, 865, 590, 264, 148, 100], # 10
                [11500, 6780, 3440, 1860, 1110, 755, 525, 251, 146, 100], # 15
                [6400, 4480, 2410, 1410, 910, 635, 454, 237, 141, 99], # 20
                [3840, 2830, 1730, 1100, 749, 545, 410, 223, 136, 97], # 25
                [2480, 1870, 1220, 845, 615, 467, 365, 207, 131, 95], # 30
                [1650, 1270, 910, 680, 510, 397, 320, 193, 125, 93], # 35
                [1180, 940, 700, 530, 416, 338, 282, 179, 120, 92], # 40
                [910, 730, 555, 442, 356, 292, 250, 166, 116, 90], # 45
                [505, 442, 352, 292, 243, 209, 183, 134, 104, 86], # 60
                [338, 317, 269, 227, 196, 172, 151, 116, 93, 82], # 75
                [259, 251, 225, 193, 166, 147, 132, 104, 86, 79], # 90
                [212, 210, 197, 170, 150, 133, 119, 96, 82, 77], # 105
                [188, 186, 177, 154, 138, 125, 113, 90, 77, 74], # 120
                [179, 178, 166, 147, 134, 122, 110, 90, 77, 73], # 135
                [179, 178, 165, 148, 137, 127, 116, 96, 79, 72], # 150
                [196, 192, 179, 165, 151, 141, 131, 104, 82, 72], # 165
                [230, 212, 195, 178, 163, 148, 134, 105, 83, 72], # 180
            ]
        )
        * SPECTRAL_RADIANCE
    )

    ZODI_INTERP = RegularGridInterpolator((long_vector, lat_vector), table17, method="cubic")

    #####
    # 2. The wavelength dependence
    #####

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

    FLUX_INTERP = CubicSpline(zodi_lambd, zodi_blambd)

load_leinert()


# Code from pyEDITH, astrophysical_scene.py calc_zodi_flux
# Courtesy of Eleonora Alei
def calc_zodi_flux(
    source,
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
    wave : Quantity
        Wavelengths in microns (vector of length nlambda).
    sn_box : Quantity
        Size of the extraction aperture in pixels squared
    pixel_scale : QUANTITY
        Dimensions of a single pixel (assumed to be square) in arcseconds.

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

    # This code is in nanometers
    wave = wave.to(u.nm)

    # Convert equatorial coordinates to ecliptic coordinates
    coords = source.coords
    ecl_coords = coords.barycentrictrueecliptic
    date = source.date
    vernal_equinox = datetime.date(source.date.year, 3, 20)
    sun_delta = (date - vernal_equinox).days
    if sun_delta < 0:
        sun_delta = 365.25+sun_delta
    sun_angle = sun_delta/365.25 * 360.0


    target_sun_angle = ecl_coords.lon.deg - sun_angle
    if target_sun_angle < 0:
        target_sun_angle = 360+target_sun_angle

    if target_sun_angle > 360:
        # not sure how this would happen...
        target_sun_angle = target_sun_angle - 360

    # now correct for the fact that it's symmetric about 180 degrees:
    # (and there is no data from 180-360)
    if target_sun_angle > 180:
        target_sun_angle = 360 - target_sun_angle

    #print(ecl_coords.lon, ecl_coords.lat, sun_angle, sun_delta, target_sun_angle, ecl_coords.lat.deg)

    # we haven't mirrored the latitude, so take the absolute value
    factor = ZODI_INTERP((target_sun_angle, np.abs(ecl_coords.lat.deg)))
    # # For coronagraph, assume observations near solar longitude of 135 degrees
    # # For all other instruments, use the actual ecliptic longitude
    # j = np.argmin(np.abs(sollong_vector - ecl_coords.lon.deg))
    # k = np.argmin(np.abs(sollong_vector - 90 * u.deg))

    # # Interpolate to get zodi brightness factor


    # interp = interp1d(
    #     np.sin(beta_vector),
    #     table17[j] / table17[k, 0],
    #     kind="cubic",
    #     fill_value="extrapolate",
    # )
    # # this specifically selects the 135 degree longitude, and interpolates to the chosen ecliptic latitude
    # f = interp(sinbeta) * u.dimensionless_unscaled


    # SYOTools actually needs the zodi in flux units (syn.units.PHOTLAM)
    # so we do not need to call out for a magnitude calculation here.
    flux = FLUX_INTERP(wave) << SPECTRAL_RADIANCE_CGS
    flux_zodi = factor * flux # apply the scaling relative to 90 degrees ecliptic (f)

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
    with warnings.catch_warnings():
        warnings.filterwarnings(
            'ignore', message=r'.*contained negative flux or throughput.*',
            category=AstropyUserWarning)
        sky = syn.spectrum.SourceSpectrum(Empirical1D, points=wave, lookup_table=syn.units.convert_flux(wave, flux_zodi, syn.units.PHOTLAM))

    return sky  # 1/arcsec^2 (UNITS OF SPECTRAL RADIANCE) - original, now PHOTLAM
