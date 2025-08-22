import syotools.environment
import pytest

from syotools.models import Spectrograph
from syotools.models.coronagraph import Coronagraph
from syotools.models.telescope import Telescope
from syotools.models.camera import Camera

from syotools.wrappers.camera_wrapper import camera_exptime
from syotools.wrappers.uvspec_wrapper import uvspec_exptime

import pickle

UVSPEC_BASELINE_PICKLE = "tests/baselines/uvspec_exptime.pickle"
CAMERA_BASELINE_PICKLE = "tests/baselines/camera_exptime.pickle"

def check_relative_diff(actual, expected, rel_tol=0.25):
    """
    Simple function to check if two lists are approximately equal within a relative tolerance.
    Reports percentage differences for values that exceed the tolerance.

    Args:
        actual: List of actual values
        expected: List of expected values
        rel_tol: Relative tolerance (default: 0.25 or 25%)

    Returns:
        True if all values are within tolerance, False otherwise
    """
    if len(actual) != len(expected):
        print(f"Lists have different lengths: actual={len(actual)}, expected={len(expected)}")
        return False

    all_within_tolerance = True
    differences = []

    for i, (a, e) in enumerate(zip(actual, expected)):
        if e == 0:
            # Can't calculate relative difference if expected is zero
            if a != 0:
                all_within_tolerance = False
                differences.append((i, a, e, "inf"))
        else:
            rel_diff = abs(a - e) / abs(e)
            pct_diff = 100 * rel_diff

            if rel_diff > rel_tol:
                all_within_tolerance = False
                differences.append((i, a, e, pct_diff))

    if not all_within_tolerance:
        print("\nValues exceeding relative tolerance:")
        for i, a, e, pct in differences:
            if pct == "inf":
                print(f"  Index {i}: actual={a}, expected={e}, difference=infinite (division by zero)")
            else:
                print(f"  Index {i}: actual={a:.6g}, expected={e:.6g}, difference={pct:.2f}%")

    return all_within_tolerance

def test_telescope_json():
    tel1, tel2 = Telescope(), Telescope()
    tel1.add_camera(Camera())
    print(tel1.cameras)
    tel2.add_camera(Camera())
    print(tel1.cameras)
    tel1.set_from_sei("EAC1")
    tel2.set_from_sei("EAC2")
    tel1.cameras[0].set_from_sei('HRI')
    tel2.cameras[0].set_from_sei('HRI')

    assert tel1.name != tel2.name
    assert tel1.effective_aperture != tel2.effective_aperture
    assert all([c.telescope == tel1 for c in tel1.cameras])
    assert tel1.cameras != tel2.cameras


camera_exptime_baseline = pickle.load(open(CAMERA_BASELINE_PICKLE, "rb"))
uvspec_exptime_baseline = pickle.load(open(UVSPEC_BASELINE_PICKLE, "rb"))
BASELINE_SIZE = 50

@pytest.mark.parametrize("magnitude, snr_goal, expected", camera_exptime_baseline)
def test_camera_exptime_calculation(magnitude, snr_goal, expected):
    exp_times, camera = camera_exptime('EAC1', 'G2V Star', magnitude, snr_goal, True)
    assert check_relative_diff([q.value for q in exp_times], expected, 0.25) #1e-3)

@pytest.mark.parametrize("magnitude, snr_goal, expected", uvspec_exptime_baseline)
def test_uvspec_exptime_calculation(magnitude, snr_goal, expected):
    wave, exp_times, uvi = uvspec_exptime('EAC1', 'G180M', 'G2V Star', magnitude, snr_goal, True)
    expected = [q.value for q in expected]
    assert check_relative_diff([q.value for q in exp_times], expected, 0.25) #1e-3)

# def test_generate_camera_exptime_baseline():
#     baseline = [(mag, snr, [q.value for q in camera_exptime('EAC1', 'G2V Star', mag, snr, True)[0]]) for mag in range(4, 35) for snr in range(3, 10)]
#     with open(CAMERA_BASELINE_PICKLE, "wb") as f:
#         random.seed(7)
#         pickle.dump(random.sample(baseline, BASELINE_SIZE, ), f)

# def test_print_things():
#     print(uvspec_exptime('EAC1', 'G180M', 'G2V Star', 20, 5.0, True)[1][0:10])

# def test_generate_uvspec_exptime_baseline():
#     # telescope, mode, template, uvmag, snr_goal
#     baseline = [(mag, snr, [exp for exp in uvspec_exptime('EAC1', 'G180M', 'G2V Star', mag, snr, True)[1]])
#                   for mag in range(4, 35)
#                   for snr in range(3, 10)]

#     with open(UVSPEC_BASELINE_PICKLE, "wb") as f:
#         random.seed(7)
#         pickle.dump(random.sample(baseline, BASELINE_SIZE), f)

# Okay, coronagraph models are not actually implemented in the structure
# that camera and uvspec are.  See cg_example.ipynb for how different
# things are right now.  The methods in CoronagraphicExposureSpec are
# all stubs which do nothing.

# def coronagraph_exptime(eta_p, raw_contrast, sigma_DeltaC):
#     rad_file = 'syotools/coronagraph/planets/F2V_5.e-1fCO2_1.e6H2Volc_1.e10BIF.out_toa.rad'

#     cac, tel = Coronagraph(), Telescope()
#     tel.set_from_json('EAC1')
#     tel.add_coronagraph(cac)

#     cac.eta_p = eta_p
#     cac.raw_contrast = raw_contrast
#     cac.sigma_DeltaC = sigma_DeltaC

#     cac_exp = cac.create_exposure()
#     cac._calc_count_rates_imaging()
#     cac._signal_to_noise
#     return cac_exp.exptime

# def test_generate_coronagraph_baseline():
#     baseline = [(eta_p, raw_contrast, sigma_DeltaC, coronagraph_exptime(eta_p, raw_contrast, sigma_DeltaC))
#                 for eta_p in [0.2, 0.4, 0.6]
#                 for raw_contrast in [1e-10, 1e-9, 1e-8]
#                 for sigma_DeltaC in [0.0, 5e-12, ]]
#     for result in baseline:
#         print(result)

def create_eac(name: str):
    telescope = Telescope()
    telescope.set_from_sei(name)
    camera = Camera()
    telescope.add_camera(camera)
    spec = Spectrograph()
    telescope.add_spectrograph(spec)
    camera.set_from_sei('HRI')
    return telescope

def test_telescope_pickle():
  """
  This test verifies that the telescope can be pickled with a camera and
  spectrograph attached. Pickling is required to support the DISRA notebook
  framework. Document gotchas related to pickling found by this test here
  so that if it starts failing people have some clues about where to look.
  - Some astropy quantities related to FITS file data retain weakrefs
    (not pickle) to the table data. Using copy() removes these.
    The spectrograph class had this issue, resolved with:
     self.wave = table['Wavelength'].copy()
  """
  import pickle
  telescope = create_eac("EAC1")
  exp = telescope.cameras[0].create_exposure()
  tele2 = pickle.loads(pickle.dumps(telescope))
  assert telescope.name == tele2.name
