import pytest

from syotools.models.telescope import Telescope
from syotools.models.camera import Camera

from syotools.wrappers.camera_wrapper import camera_exptime
from syotools.wrappers.uvspec_wrapper import uvspec_exptime

import pickle

UVSPEC_BASELINE_PICKLE = "tests/baselines/uvspec_exptime.pickle"
CAMERA_BASELINE_PICKLE = "tests/baselines/camera_exptime.pickle"

def test_telescope():
    tel1, tel2 = Telescope(), Telescope()
    tel1.add_camera(Camera())
    print(tel1.cameras)
    tel2.add_camera(Camera())
    print(tel1.cameras)
    tel1.set_from_json("EAC1")
    tel2.set_from_json("EAC2")

    assert tel1.name != tel2.name
    assert tel1.aperture != tel2.aperture
    assert all([c.telescope == tel1 for c in tel1.cameras])
    assert tel1.cameras != tel2.cameras


camera_exptime_baseline = pickle.load(open(CAMERA_BASELINE_PICKLE, "rb"))
uvspec_exptime_baseline = pickle.load(open(UVSPEC_BASELINE_PICKLE, "rb"))

@pytest.mark.parametrize("magnitude, snr_goal, expected", camera_exptime_baseline)
def test_camera_exptime_calculation(magnitude, snr_goal, expected):
    exp_times, camera = camera_exptime('EAC1', 'G2V Star', magnitude, snr_goal, True)
    assert [q.value for q in exp_times] == pytest.approx(expected)

@pytest.mark.parametrize("magnitude, snr_goal, expected", uvspec_exptime_baseline)
def test_uvspec_exptime_calculation(magnitude, snr_goal, expected):
    wave, exp_times, uvi = uvspec_exptime('EAC1', 'G180M', 'G2V Star', magnitude, snr_goal, True)
    assert [q.value for q in exp_times] == pytest.approx(expected)

# def test_generate_camera_exptime_baseline():
#     baseline = [(mag, snr, [q.value for q in camera_exptime('EAC1', 'G2V Star', mag, snr, True)[0]]) for mag in range(4, 35) for snr in range(3, 10)]
#     with open(CAMERA_BASELINE_PICKLE, "wb") as f:
#         pickle.dump(baseline, f)

# def test_generate_uvspec_exptime_baseline():
#     # telescope, mode, template, uvmag, snr_goal
#     baseline = [(mag, snr, [q.value for q in uvspec_exptime('EAC1', 'G180M', 'G2V Star', mag, snr, True)[1]])
#                   for mag in range(4, 35)
#                   for snr in range(3, 10)]

#     with open(UVSPEC_BASELINE_PICKLE, "wb") as f:
#         pickle.dump(baseline, f)


def create_eac(name: str):
    telescope = Telescope()
    telescope.set_from_json(name)
    telescope.add_camera(Camera())
    return telescope
