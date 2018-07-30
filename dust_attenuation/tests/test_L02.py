import numpy as np
import pytest

import astropy.units as u
from astropy.modeling import InputParameterError

from ..averages import L02
from .helpers import _invalid_x_range


@pytest.mark.parametrize("Av_invalid", [-1.0, -0.00001, -10])
def test_invalid_Av_input(Av_invalid):
    with pytest.raises(InputParameterError) as exc:
        tmodel = L02(Av=Av_invalid)
    assert exc.value.args[0] == 'parameter Av must be positive'


@pytest.mark.parametrize("x_invalid", [-1.0, 0.09, 11, 100.])
def test_invalid_wavenumbers(x_invalid):
    _invalid_x_range(x_invalid, L02(), 'L02')


@pytest.mark.parametrize("x_invalid_wavenumber",
                         [-1.0, 0.09, 11, 100.]/u.micron)
def test_invalid_wavenumbers_imicron(x_invalid_wavenumber):
    _invalid_x_range(x_invalid_wavenumber, L02(), 'L02')


@pytest.mark.parametrize("x_invalid_micron",
                         u.micron/[-1.0, 0.08, 11, 100.])
def test_invalid_micron(x_invalid_micron):
    _invalid_x_range(x_invalid_micron, L02(), 'L02')


@pytest.mark.parametrize("x_invalid_angstrom",
                         u.angstrom*1e4/[-1.0, 0.09, 11, 100.])
def test_invalid_angstrom(x_invalid_angstrom):
    _invalid_x_range(x_invalid_angstrom, L02(), 'L02')


def get_axav_cor_vals(Av):
    # testing wavenumbers. Validity range 0.097 - 0.18 microns, equation (14)
    # in Leitherer (2002, ApJS, Volume 140, Issue 2, pp. 303-329).
    # The reference values were computed using this code.
    x = np.array([0.097, 0.10136842, 0.10573684, 0.11010526, 0.11447368,
                  0.11884211, 0.12321053, 0.12757895, 0.13194737, 0.13631579,
                  0.14068421, 0.14505263, 0.14942105, 0.15378947, 0.15815789,
                  0.16252632, 0.16689474, 0.17126316, 0.17563158, 0.18])

    # add units
    x = x*u.micron

    # correct values generated using this code
    if Av == 0.2:
        cor_vals = np.array([0.70521057, 0.67702016, 0.65233204, 0.63054867,
                             0.6111968, 0.59389768, 0.57834525, 0.56428989,
                             0.55152634, 0.53988449, 0.52922231, 0.51942042,
                             0.51037781, 0.50200855, 0.49423907, 0.48700606,
                             0.4802548, 0.47393773, 0.46801335, 0.46244529])
    elif Av == 1.0:
        cor_vals = np.array([3.52605287, 3.38510078, 3.26166022, 3.15274335,
                             3.05598398, 2.96948842, 2.89172624, 2.82144944,
                             2.75763171, 2.69942246, 2.64611156, 2.59710209,
                             2.55188907, 2.51004276, 2.47119533, 2.4350303,
                             2.401274, 2.36968867, 2.34006677, 2.31222646])
    elif Av == 2.4:
        cor_vals = np.array([8.46252688, 8.12424187, 7.82798452, 7.56658403,
                             7.33436156, 7.1267722, 6.94014298, 6.77147866,
                             6.61831611, 6.47861391, 6.35066774, 6.23304501,
                             6.12453378, 6.02410262, 5.93086878, 5.84407272,
                             5.76305761, 5.6872528, 5.61616025, 5.54934349])
    elif Av == 5.0:
        cor_vals = np.array([17.63026434, 16.92550389, 16.30830109,
                             15.76371673, 15.27991992, 14.84744209,
                             14.45863121, 14.1072472, 13.78815856,
                             13.49711231, 13.2305578, 12.98551044,
                             12.75944537, 12.55021379, 12.35597663,
                             12.17515149, 12.00637002, 11.84844333,
                             11.70033386, 11.56113228])
    elif Av == 10.0:
        cor_vals = np.array([35.26052867, 33.85100778, 32.61660217,
                             31.52743346, 30.55983983, 29.69488417,
                             28.91726242, 28.21449441, 27.57631713,
                             26.99422463, 26.4611156, 25.97102088,
                             25.51889074, 25.10042759, 24.71195326,
                             24.35030298, 24.01274003, 23.69688666,
                             23.40066773, 23.12226456])
    else:
        cor_vals = np.array([0.0])

    return (x, cor_vals)


@pytest.mark.parametrize("Av", [0.2, 1.0, 2.4, 5.0, 10.0])
def test_attenuation_L02_values(Av):
    # get the correct values
    x, cor_vals = get_axav_cor_vals(Av)

    # initialize extinction model
    tmodel = L02(Av=Av)

    # test.
    np.testing.assert_allclose(tmodel(x), cor_vals, atol=1e-7)


@pytest.mark.parametrize("Av", [0.2, 1.0, 2.4, 5.0, 10.0])
def test_attenuation_L02_attenuate_values_Av(Av):
    # get the correct values
    x, cor_vals = get_axav_cor_vals(Av)

    # calculate the cor_vals in fractional units
    cor_vals = np.power(10.0, -0.4*(cor_vals))

    # initialize extinction model
    tmodel = L02(Av=Av)

    # test
    np.testing.assert_allclose(tmodel.attenuate(x), cor_vals, atol=1e-6)
