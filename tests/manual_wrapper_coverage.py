from syotools.wrappers.camera_wrapper import camera_exptime, camera_snr, camera_magnitude
from syotools.wrappers.uvspec_wrapper import uvspec_snr, uvspec_exptime
from syotools.wrappers.common import compute_observation

camera_exptime("EAC5", "Flat (AB)", 30, 10, silent=False)

camera_snr("EAC5", "Flat (AB)", 30, 10, silent=False)

camera_magnitude("EAC5", "Flat (AB)", 30, 10, silent=False)

uvspec_exptime("EAC5", "G150M", "Flat (AB)", 30, 10, silent=False)

uvspec_snr("EAC5", "G150M", "Flat (AB)", 30, 10, silent=False)

compute_observation("EAC5", instrument="camera", sed="Flat (AB)", magnitude=25.0, snr=10, redshift=0, extinction=0, bandpass="galex,fuv", target="exptime", verbose=True)
compute_observation("EAC5", instrument="camera", sed="Flat (AB)", magnitude=25.0, exptime=1, redshift=0, extinction=0, bandpass="galex,fuv", target="snr", verbose=True)
compute_observation("EAC5", instrument="camera", sed="Flat (AB)", exptime=1, snr=10, redshift=0, extinction=0, bandpass="galex,fuv", target="magnitude", verbose=True)
compute_observation("EAC5", instrument="spectroscopy", sed="Flat (AB)", magnitude=25.0, snr=10, redshift=0, extinction=0, bandpass="galex,fuv", target="exptime", verbose=True)
compute_observation("EAC5", instrument="spectroscopy", sed="Flat (AB)", magnitude=25.0, exptime=1, redshift=0, extinction=0, bandpass="galex,fuv", target="snr", verbose=True)
compute_observation("EAC5", instrument="ifs", sed="Flat (AB)", magnitude=25.0, snr=10, redshift=0, extinction=0, bandpass="galex,fuv", target="exptime", verbose=True)
compute_observation("EAC5", instrument="ifs", sed="Flat (AB)", magnitude=25.0, exptime=1, redshift=0, extinction=0, bandpass="galex,fuv", target="snr", verbose=True)
