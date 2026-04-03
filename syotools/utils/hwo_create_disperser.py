import numpy as np
from scipy import interpolate as interp
from astropy.io import fits

dispersers = [{"name": "UVPASS", "min": 900, "max": 3600, "R": 20000, "diff_limit": 3000},
              {"name": "OPTPASS", "min": 3200, "max": 9000, "R": 20000, "diff_limit": 5000},
              {"name": "IRPASS", "min": 8000, "max": 18000, "R": 20000, "diff_limit": 10000},
              {"name": "HIRES", "min": 900, "max": 18000, "R": 100000, "diff_limit": 3000}]

def load_readnoise(wavepix, filename=None):
    return np.ones_like(wavepix) * 0

def load_dark(wavepix, filename=None):
    return 4.2e-6 * np.ones_like(wavepix)

def load_sky(wavepix, filename=None):
    return np.zeros_like(wavepix)

def load_aeff(wavepix, filename=None, wave_lo=None, wave_hi=None):
    # relative to a 15m diameter mirror, in cm^2
    wave_temp = np.arange(wave_lo, wave_hi+33,33)
    flux_temp = np.ones_like(wave_temp)
    flux_temp[0] = 0
    flux_temp[1] = 0
    flux_temp[2] = 0
    flux_temp[-1] = 0
    flux_temp[-2] = 0
    flux_temp[-3] = 0
    interpolator = interp.Akima1DInterpolator(wave_temp, flux_temp)
    flux = interpolator(wavepix)
    return flux * 7500 * np.pi ** 2

def load_bef(wavepix, filename=None):
    # in erg/cm2/pixel/s
    return np.ones_like(wavepix) * 1.75e-14

def load_disp_width(wavepix, filename=None, diff_limit=None):
    disp_width = np.ones_like(wavepix)
    if isinstance(diff_limit, float):

        scaled_diffraction = 1.22 * wavepix # this is just for scaling purposes so I omit the / D part of lambda / D
        scaled_diffraction /= 1.22 * diff_limit # to scale this diffraction to pixels, assuming it's set up for 
        scaled_diffraction[np.where(scaled_diffraction < 1)] = 1

        disp_width = scaled_diffraction * 2
    return disp_width

def load_xdisp_width(wavepix, filename=None, diff_limit=None):
    xdisp_width = np.ones_like(wavepix)
    if isinstance(diff_limit, float):

        scaled_diffraction = 1.22 * wavepix # this is just for scaling purposes so I omit the / D part of lambda / D
        scaled_diffraction /= 1.22 * diff_limit # to scale this diffraction to pixels, assuming it's set up for 
        scaled_diffraction[np.where(scaled_diffraction < 1)] = 1

        xdisp_width = scaled_diffraction
    return xdisp_width

def make_file(dispersers):
    hdulist = []

    hdulist.append(fits.PrimaryHDU())

    head = hdulist[0].header

    for disperser in dispersers:
        name = disperser["name"]
        wave_lo = disperser["min"]
        wave_hi = disperser["max"]
        r = disperser["R"]
        diff_limit = disperser["diff_limit"]
        wave = np.arange(wave_lo, wave_hi, 1)

        dlds = wave/r

        # from pandeia.engine Instrument.get_wave_pix()
        pixel = np.cumsum(1.0 / dlds * np.gradient(wave))
        pixel_integer = np.arange(int(pixel[0]), int(pixel[-1]))
        wavepix = np.interp(pixel_integer, pixel, wave)

        readnoise = load_readnoise(wavepix, None) # not used
        sky = load_sky(wavepix, None) # complicated - Lorentzian at ~1210 Angstroms (Ly-alpha?), not used
        dark = load_dark(wavepix, None) # not used
        Aeff = load_aeff(wavepix, None, wave_lo=wave_lo, wave_hi=wave_hi) # Throughput as effective area in cm^2
        BEF = load_bef(wavepix, None) # Single peak of 1.75e-13 at 1000 Angstroms, Background Emission Function
        Disp_Width = load_disp_width(wavepix, None, diff_limit=diff_limit) # Complex function with minimum at 1075-1175A, not used
        XDisp_Width = load_xdisp_width(wavepix, None, diff_limit=diff_limit) # Fourth order with troughs at 1025-1075A and 1290-1325A, not used

        wavelength_col = fits.Column(name="Wavelength", format="D", unit="Angstrom", array=wavepix)
        readnoise_col = fits.Column(name="ReadNoise", format="D", unit="count", array=readnoise)
        dark_col = fits.Column(name="Dark", format="D", unit="count pix-1 s-1", array=dark)
        sky_col = fits.Column(name="Sky", format="D", unit="count pix-1 s-1", array=sky)
        aeff_col = fits.Column(name="Aeff", format="D", unit="cm2", array=Aeff)
        bef_col = fits.Column(name="BEF", format="D", unit="cm-2 erg pix-1 s-1", array=BEF)
        disp_col = fits.Column(name="Disp_Width", format="D", unit="pix", array=Disp_Width)
        xdisp_col = fits.Column(name="XDisp_Width", format="D", unit="pix", array=XDisp_Width)

        hdu = fits.BinTableHDU.from_columns([wavelength_col, readnoise_col, dark_col, sky_col, aeff_col, bef_col, disp_col, xdisp_col], name=name.upper())
        hdu.header["NAME"] = name
        hdu.header["R"] = r
        hdu.header["WAVE_LO"] = wave_lo
        hdu.header["WAVE_HI"] = wave_hi
        hdu.header["COMMENT"] = "(Angstroms) (Cnt/Exposure) (Cnt/s/res)   (Cnt/s/res)      (cm^2)  (Erg/cm^2/s/res)   (N_resels)     (N_resels)"
        hdulist.append(hdu)
        head[name.upper()] = f"{name} (R = {r})"

    hdulist[0].header = head

    fitsobj = fits.HDUList(hdulist)
    fitsobj.writeto("output.fits", overwrite=True)

make_file(dispersers)
