import numpy as np
import pytest

import astropy.units as u
from astropy.modeling import InputParameterError

from ..averages import C00
from .helpers import _invalid_x_range


@pytest.mark.parametrize("Av_invalid", [-1.0, -0.00001, -10])
def test_invalid_Av_input(Av_invalid):
    with pytest.raises(InputParameterError) as exc:
        tmodel = C00(Av=Av_invalid)
    assert exc.value.args[0] == 'parameter Av must be positive'


@pytest.mark.parametrize("x_invalid", [-1.0, 0.1, 10.1, 100.])
def test_invalid_wavenumbers(x_invalid):
    _invalid_x_range(x_invalid, C00(Av=1), 'C00')


@pytest.mark.parametrize("x_invalid_wavenumber",
                         [-1.0, 0.1, 10.1, 100.]/u.micron)
def test_invalid_wavenumbers_imicron(x_invalid_wavenumber):
    _invalid_x_range(x_invalid_wavenumber, C00(Av=1), 'C00')


@pytest.mark.parametrize("x_invalid_micron",
                         u.micron/[-1.0, 0.1, 10.1, 100.])
def test_invalid_micron(x_invalid_micron):
    _invalid_x_range(x_invalid_micron, C00(Av=1), 'C00')


@pytest.mark.parametrize("x_invalid_angstrom",
                         u.angstrom*1e4/[-1.0, 0.1, 10.1, 100.])
def test_invalid_micron(x_invalid_angstrom):
    _invalid_x_range(x_invalid_angstrom, C00(Av=1), 'C00')


def get_axav_cor_vals(Av):
    # testing wavenumbers. Validity range 0.12 - 2.2 microns, equation (4) 
    # in Calzetti (2000, ApJ, Volume 533, Issue 2, pp. 682-695)
    x = np.array([0.12, 0.22947368, 0.33894737, 0.44842105, 0.55789474,
                  0.66736842, 0.77684211, 0.88631579, 0.99578947,
                  1.10526316, 1.21473684, 1.32421053, 1.43368421, 1.54315789,
                  1.65263158, 1.76210526,  1.87157895,  1.98105263, 
                  2.09052632,  2.2])

    # add units
    x = x*u.micron

    # correct values generated using this code
    if Av == 0.2:
        cor_vals = np.array([0.59848769, 0.40617258, 0.31227495, 0.24549295,
                             0.19684966, 0.16078594, 0.13194974, 0.11023698,
                             0.09329826, 0.07971503, 0.06858008, 0.05928619,
                             0.05141164, 0.04465435, 0.0387923, 0.03365862,
                             0.02912551, 0.0250934, 0.02148359, 0.])
    elif Av == 1.0:
        cor_vals = np.array([2.99243843, 2.03086288, 1.56137474, 1.22746473,
                             0.98424828, 0.80392969, 0.65974871, 0.55118488,
                             0.46649132, 0.39857516, 0.34290038, 0.29643097,
                             0.25705821, 0.22327176, 0.19396148, 0.1682931,
                             0.14562754, 0.125467, 0.10741794, 0.])
    elif Av == 2.4:
        cor_vals = np.array([7.18185222, 4.8740709, 3.74729937, 2.94591535,
                             2.36219588, 1.92943125, 1.58339689, 1.32284371,
                             1.11957917, 0.95658037, 0.82296091, 0.71143433,
                             0.6169397, 0.53585223, 0.46550756, 0.40390344,
                             0.34950611, 0.3011208, 0.25780305, 0.])
    elif Av == 5.0:
        cor_vals = np.array([14.96219214, 10.15431438, 7.80687369, 6.13732365,
                             4.92124141, 4.01964844, 3.29874353, 2.75592439,
                             2.33245661, 1.99287578, 1.71450189, 1.48215485,
                             1.28529105, 1.11635882, 0.96980742, 0.84146551,
                             0.72813772, 0.62733501, 0.53708968, 0.])
    elif Av == 10.0:
        cor_vals = np.array([29.92438427, 20.30862876, 15.61374738,
                             12.27464731, 9.84248281, 8.03929687, 6.59748706,
                             5.51184878, 4.66491322, 3.98575156, 3.42900378,
                             2.96430969, 2.5705821, 2.23271764, 1.93961483,
                             1.68293101, 1.45627545, 1.25467002, 1.07417936,
                             0.])
    else:
        cor_vals = np.array([0.0])

    return (x, cor_vals)


@pytest.mark.parametrize("Av", [0.2, 1.0, 2.4, 5.0, 10.0])
def test_attenuation_C00_values(Av):
    # get the correct values
    x, cor_vals = get_axav_cor_vals(Av)

    # initialize extinction model
    tmodel = C00(Av=Av)
    
    # test. Needed to decreased atol to 1e-7 because of Av=0.2 case
    np.testing.assert_allclose(tmodel(x), cor_vals,atol=1e-7)


@pytest.mark.parametrize("Av", [0.2, 1.0, 2.4, 5.0, 10.0])
def test_attenuation_C00_attenuate_values_Av(Av):
    # get the correct values
    x, cor_vals = get_axav_cor_vals(Av)

    # calculate the cor_vals in fractional units
    cor_vals = np.power(10.0, -0.4*(cor_vals))

    # initialize extinction model
    tmodel = C00(Av=Av)

    # test
    np.testing.assert_allclose(tmodel.attenuate(x), cor_vals,atol=1e-10)
