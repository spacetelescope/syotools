from syotools.models import Telescope

tel = Telescope()
tel.set_from_hwome("EAC5")

def test_generic():
    suitable_instruments, suitable_filters = tel.find_instrument_with()
    #print(suitable_instruments)
    #print(suitable_filters)
    assert len(suitable_filters) > 30

def test_filter():
    suitable_instruments, suitable_filters = tel.find_instrument_with(kind="filter")
    assert len(suitable_filters) > 15

def test_disperser():
    suitable_instruments, suitable_filters = tel.find_instrument_with(kind="disperser")
        # for x in suitable_filters:
        #     print(x)
    assert len(suitable_filters) > 10

def test_wavelength():
    suitable_instruments, suitable_filters = tel.find_instrument_with(wavelength=4500)
    # for x in suitable_filters:
    #     print(x)
    assert len(suitable_filters) > 2

def test_disperser_wavelength():
    suitable_instruments, suitable_filters = tel.find_instrument_with(kind="disperser", wavelength=4500.0)
    # for x in suitable_filters:
    #     print(x)
    assert len(suitable_filters) > 1

def test_disperser_wrange():
    suitable_instruments, suitable_filters = tel.find_instrument_with(kind="disperser", wavelength=[5500, 6500])
    # for x in suitable_filters:
    #     print(x)
    assert len(suitable_filters) > 1

def test_resolution():
    suitable_instruments, suitable_filters = tel.find_instrument_with(resolution = 2000.0)
    # for x in suitable_filters:
    #     print(x, suitable_filters[x])
    assert len(suitable_filters) > 5

def test_resrange():
    suitable_instruments, suitable_filters = tel.find_instrument_with(resolution = {"min":5000.0, "max": 20000.0})
    # for x in suitable_filters:
    #     print(x)
    assert len(suitable_filters) > 2

def test_waveres():
    suitable_instruments, suitable_filters = tel.find_instrument_with(wavelength={"wave_min": 5000, "wave_max": 6000}, resolution = {"min":5000.0, "max": 20000.0})
    # for x in suitable_filters:
    #     print(x)
    assert len(suitable_filters) > 0

def test_instrument():
    suitable_instruments, suitable_filters = tel.find_instrument_with(instrument="IFS")
    assert len(suitable_filters) > 3
