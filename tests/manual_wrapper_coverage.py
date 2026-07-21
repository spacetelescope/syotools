from matplotlib import pyplot as plt
from syotools.wrappers.camera_wrapper import camera_exptime, camera_snr, camera_magnitude
from syotools.wrappers.uvspec_wrapper import uvspec_snr, uvspec_exptime
from syotools.wrappers.common import compute_observation

# print("---- Camera Wrapper EXPtime -----")
# result,hri = camera_exptime("EAC5", "Flat (AB)", 30, 10, silent=True)
# print(result)
# print()

# print("---- Camera Wrapper SNR -----")
# result,hri = camera_snr("EAC5", "Flat (AB)", 30, 10, silent=True)
# print(result)
# print()

# print("---- Camera Wrapper Magnitude -----")
# result,hri = camera_magnitude("EAC5", "Flat (AB)", 30, 10, silent=True)
# print(result)
# print()

# print("---- UVSpec Wrapper EXPtime -----")
# wave,result,uvi = uvspec_exptime("EAC5", "G150M", "Flat (AB)", 30, 10, silent=True)
# print(result)
# plt.plot(wave, result)
# plt.ylim(1e1,1e22)
# plt.yscale("log")
# plt.show()
# print()

print("---- UVSpec Wrapper SNR -----")
wave,result,uvi = uvspec_snr("EAC5", "G150M", "Flat (AB)", 30, 10, silent=True)
print(result)
print()

# print("---- Camera Common Exptime -----")
# result = compute_observation("EAC5", instrument="camera", sed="Flat (AB)", magnitude=25.0, snr=10, redshift=0, extinction=0, bandpass="galex,fuv", target="exptime", verbose=False)
# print(result)
# result = compute_observation("EAC5", instrument="camera", sed="Flat (AB)", magnitude=25.0, exptime=1, redshift=0, extinction=0, bandpass="galex,fuv", target="snr", verbose=False)
# print(result)
# result = compute_observation("EAC5", instrument="camera", sed="Flat (AB)", exptime=1, snr=10, redshift=0, extinction=0, bandpass="galex,fuv", target="magnitude", verbose=False)
# print(result)
# result = compute_observation("EAC5", instrument="spectroscopy", sed="Flat (AB)", magnitude=25.0, snr=10, redshift=0, extinction=0, bandpass="galex,fuv", target="exptime", verbose=False)
# print(result)
# result = compute_observation("EAC5", instrument="spectroscopy", sed="Flat (AB)", magnitude=25.0, exptime=1, redshift=0, extinction=0, bandpass="galex,fuv", target="snr", verbose=False)
# print(result)
# print("---- IFS Common Exptime -----")
# result = compute_observation("EAC5", instrument="ifs", sed="Flat (AB)", magnitude=25.0, snr=10, redshift=0, extinction=0, bandpass="galex,fuv", target="exptime", verbose=False)
# # plt.plot(wave, result)
# # plt.ylim(1e1,1e22)
# # plt.yscale("log")
# # plt.show()
#print(result)
result = compute_observation("EAC5", instrument="ifs", sed="Flat (AB)", magnitude=25.0, exptime=1, redshift=0, extinction=0, bandpass="galex,fuv", target="snr", verbose=False)
print(result)
