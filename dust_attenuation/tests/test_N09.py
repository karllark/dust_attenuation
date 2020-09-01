import numpy as np
import pytest

import astropy.units as u
from astropy.modeling import InputParameterError

from ..shapes import N09
from .helpers import _invalid_x_range


@pytest.mark.parametrize("x0_invalid", [-1.0, -0.00001, -10])
def test_invalid_x0_input(x0_invalid):
    with pytest.raises(InputParameterError) as exc:
        N09(x0=x0_invalid)
    assert exc.value.args[0] == "parameter x0 must be positive"


@pytest.mark.parametrize("gamma_invalid", [-1.0, -0.00001, -10])
def test_invalid_gamma_input(gamma_invalid):
    with pytest.raises(InputParameterError) as exc:
        N09(gamma=gamma_invalid)
    assert exc.value.args[0] == "parameter gamma must be positive"


@pytest.mark.parametrize("ampl_invalid", [-1.0, -0.00001, -10])
def test_invalid_ampl_input(ampl_invalid):
    with pytest.raises(InputParameterError) as exc:
        N09(ampl=ampl_invalid)
    assert exc.value.args[0] == "parameter ampl must be positive"


@pytest.mark.parametrize("slope_invalid", [-4.0, -3.00001, 10])
def test_invalid_slope_input(slope_invalid):
    with pytest.raises(InputParameterError) as exc:
        N09(slope=slope_invalid)
    assert exc.value.args[0] == "parameter slope must be between -3.0 and 3.0"


@pytest.mark.parametrize("Av_invalid", [-1.0, -0.00001, -10])
def test_invalid_Av_input(Av_invalid):
    with pytest.raises(InputParameterError) as exc:
        N09(Av=Av_invalid)
    assert exc.value.args[0] == "parameter Av must be positive"


@pytest.mark.parametrize("x_invalid", [-1.0, 0.09, 11, 100.0])
def test_invalid_wavenumbers(x_invalid):
    _invalid_x_range(x_invalid, N09(), "N09")


@pytest.mark.parametrize("x_invalid_wavenumber", [-1.0, 0.09, 11, 100.0] / u.micron)
def test_invalid_wavenumbers_imicron(x_invalid_wavenumber):
    _invalid_x_range(x_invalid_wavenumber, N09(), "N09")


@pytest.mark.parametrize("x_invalid_micron", u.micron / [-1.0, 0.09, 11, 100.0])
def test_invalid_micron(x_invalid_micron):
    _invalid_x_range(x_invalid_micron, N09(), "N09")


@pytest.mark.parametrize(
    "x_invalid_angstrom", u.angstrom * 1e4 / [-1.0, 0.09, 11, 100.0]
)
def test_invalid_angstrom(x_invalid_angstrom):
    _invalid_x_range(x_invalid_angstrom, N09(), "N09")


def get_axav_cor_vals(x0, gamma, ampl, slope, Av):
    # testing wavenumbers. Validity range 0.097 - 2.2 microns
    # correct values generated using this code
    x = np.array(
        [
            0.097,
            0.1021387,
            0.10785232,
            0.11424305,
            0.12143884,
            0.12960205,
            0.13894181,
            0.14973226,
            0.16233985,
            0.17726577,
            0.19521425,
            0.21720684,
            0.24478387,
            0.28038172,
            0.32809516,
            0.39537786,
            0.49737488,
            0.67029261,
            1.02752154,
            2.2,
        ]
    )

    # add units
    x = x * u.micron
    if x0 == 0.2175:
        if gamma == 0.035:
            if ampl == 0.0:
                if slope == -1.0:
                    if Av == 0.2:
                        cor_vals = np.array(
                            [
                                0.0,
                                0.0476727,
                                0.13119835,
                                0.2454022,
                                0.38269147,
                                0.53675915,
                                0.70305229,
                                0.87839772,
                                1.06100214,
                                1.25045205,
                                1.44771379,
                                1.65513354,
                                1.87246302,
                                2.11132361,
                                2.36899125,
                                2.64703323,
                                2.94714046,
                                3.27112745,
                                3.62093228,
                                3.99861665,
                            ]
                        )
                    elif Av == 1.0:
                        cor_vals = np.array(
                            [
                                0.0,
                                0.23836349,
                                0.65599175,
                                1.22701102,
                                1.91345735,
                                2.68379577,
                                3.51526143,
                                4.3919886,
                                5.30501068,
                                6.25226024,
                                7.23856897,
                                8.2756677,
                                9.36231512,
                                10.55661806,
                                11.84495623,
                                13.23516615,
                                14.73570232,
                                16.35563723,
                                18.10466138,
                                19.99308327,
                            ]
                        )
                elif slope == 0.0:
                    if Av == 0.2:
                        cor_vals = np.array(
                            [
                                0.0,
                                0.08906314,
                                0.15989324,
                                0.22192162,
                                0.27510497,
                                0.32019651,
                                0.35840547,
                                0.39094108,
                                0.41901258,
                                0.4438292,
                                0.46660019,
                                0.48853477,
                                0.50976024,
                                0.53336569,
                                0.55822929,
                                0.58445936,
                                0.6121642,
                                0.64145213,
                                0.67243149,
                                0.70521057,
                            ]
                        )
                    elif Av == 1.0:
                        cor_vals = np.array(
                            [
                                0.0,
                                0.44531568,
                                0.79946622,
                                1.1096081,
                                1.37552487,
                                1.60098256,
                                1.79202733,
                                1.95470538,
                                2.09506288,
                                2.219146,
                                2.33300094,
                                2.44267386,
                                2.54880118,
                                2.66682847,
                                2.79114647,
                                2.92229678,
                                3.06082098,
                                3.20726067,
                                3.36215744,
                                3.52605287,
                            ]
                        )
                elif slope == 1.0:
                    if Av == 0.2:
                        cor_vals = np.array(
                            [
                                0.0,
                                0.16638962,
                                0.19486411,
                                0.20068771,
                                0.1977644,
                                0.19100896,
                                0.18270971,
                                0.17399285,
                                0.16547708,
                                0.15753052,
                                0.1503859,
                                0.14419756,
                                0.13877737,
                                0.13473963,
                                0.1315412,
                                0.12904739,
                                0.12715546,
                                0.12578563,
                                0.12487505,
                                0.1243735,
                            ]
                        )
                    elif Av == 1.0:
                        cor_vals = np.array(
                            [
                                0.0,
                                0.8319481,
                                0.97432055,
                                1.00343853,
                                0.98882198,
                                0.95504478,
                                0.91354854,
                                0.86996426,
                                0.82738541,
                                0.78765259,
                                0.75192948,
                                0.72098782,
                                0.69388686,
                                0.67369815,
                                0.65770599,
                                0.64523697,
                                0.6357773,
                                0.62892817,
                                0.62437526,
                                0.62186751,
                            ]
                        )
            elif ampl == 5.0:
                if slope == -1.0:
                    if Av == 0.2:
                        cor_vals = np.array(
                            [
                                1.59292377e-05,
                                4.78405545e-02,
                                1.31886067e-01,
                                2.47453670e-01,
                                3.88138032e-01,
                                5.51226180e-01,
                                7.46382511e-01,
                                1.05355562e00,
                                1.68604838e00,
                                1.49782806e00,
                                1.54913409e00,
                                1.71243926e00,
                                1.91097429e00,
                                2.13993978e00,
                                2.39165273e00,
                                2.66577127e00,
                                2.96311962e00,
                                3.28506963e00,
                                3.63331261e00,
                                4.00976290e00,
                            ]
                        )
                    elif Av == 1.0:
                        cor_vals = np.array(
                            [
                                7.96461884e-05,
                                2.39202772e-01,
                                6.59430336e-01,
                                1.23726835e00,
                                1.94069016e00,
                                2.75613090e00,
                                3.73191256e00,
                                5.26777812e00,
                                8.43024191e00,
                                7.48914031e00,
                                7.74567043e00,
                                8.56219628e00,
                                9.55487144e00,
                                1.06996989e01,
                                1.19582636e01,
                                1.33288564e01,
                                1.48155981e01,
                                1.64253482e01,
                                1.81665630e01,
                                2.00488145e01,
                            ]
                        )
                elif slope == 0.0:
                    if Av == 0.2:
                        cor_vals = np.array(
                            [
                                6.37169507e-05,
                                8.93767277e-02,
                                1.60731375e-01,
                                2.23776797e-01,
                                2.79020339e-01,
                                3.28826623e-01,
                                3.80494563e-01,
                                4.68897130e-01,
                                6.65856788e-01,
                                5.31631607e-01,
                                4.99288090e-01,
                                5.05449319e-01,
                                5.20244560e-01,
                                5.40594753e-01,
                                5.63569247e-01,
                                5.88596677e-01,
                                6.15483300e-01,
                                6.44186129e-01,
                                6.74730598e-01,
                                7.07176366e-01,
                            ]
                        )
                    elif Av == 1.0:
                        cor_vals = np.array(
                            [
                                3.18584754e-04,
                                4.46883639e-01,
                                8.03656874e-01,
                                1.11888399e00,
                                1.39510169e00,
                                1.64413311e00,
                                1.90247281e00,
                                2.34448565e00,
                                3.32928394e00,
                                2.65815804e00,
                                2.49644045e00,
                                2.52724660e00,
                                2.60122280e00,
                                2.70297377e00,
                                2.81784624e00,
                                2.94298338e00,
                                3.07741650e00,
                                3.22093064e00,
                                3.37365299e00,
                                3.53588183e00,
                            ]
                        )
                elif slope == 1.0:
                    if Av == 0.2:
                        cor_vals = np.array(
                            [
                                2.54867803e-04,
                                1.66975478e-01,
                                1.95885551e-01,
                                2.02365377e-01,
                                2.00579029e-01,
                                1.96157134e-01,
                                1.93970397e-01,
                                2.08688097e-01,
                                2.62961174e-01,
                                1.88694666e-01,
                                1.60921252e-01,
                                1.49190118e-01,
                                1.41631630e-01,
                                1.36565846e-01,
                                1.32799504e-01,
                                1.29960905e-01,
                                1.27844887e-01,
                                1.26321757e-01,
                                1.25302012e-01,
                                1.24720195e-01,
                            ]
                        )
                    elif Av == 1.0:
                        cor_vals = np.array(
                            [
                                1.27433901e-03,
                                8.34877391e-01,
                                9.79427753e-01,
                                1.01182688e00,
                                1.00289514e00,
                                9.80785670e-01,
                                9.69851987e-01,
                                1.04344049e00,
                                1.31480587e00,
                                9.43473329e-01,
                                8.04606261e-01,
                                7.45950589e-01,
                                7.08158148e-01,
                                6.82829230e-01,
                                6.63997521e-01,
                                6.49804526e-01,
                                6.39224435e-01,
                                6.31608786e-01,
                                6.26510060e-01,
                                6.23600977e-01,
                            ]
                        )

    else:
        cor_vals = np.array([0.0])

    return (x, cor_vals)


@pytest.mark.parametrize("x0", [0.2175])
@pytest.mark.parametrize("gamma", [0.035])
@pytest.mark.parametrize("ampl", [0.0, 5.0])
@pytest.mark.parametrize("slope", [-1.0, 0.0, 1.0])
@pytest.mark.parametrize("Av", [0.2, 1.0])
def test_attenuation_N09_values(x0, gamma, ampl, slope, Av):
    # get the correct values
    x, cor_vals = get_axav_cor_vals(x0, gamma, ampl, slope, Av)

    # initialize model
    tmodel = N09(x0=x0, gamma=gamma, ampl=ampl, slope=slope, Av=Av)

    # test. Needed to decreased atol to 1e-7 because of Av=0.2 case
    np.testing.assert_allclose(tmodel(x), cor_vals[::-1], atol=1e-7)


@pytest.mark.parametrize("x0", [0.2175])
@pytest.mark.parametrize("gamma", [0.035])
@pytest.mark.parametrize("ampl", [0.0, 5.0])
@pytest.mark.parametrize("slope", [-1.0, 0.0, 1.0])
@pytest.mark.parametrize("Av", [0.2, 1.0])
def test_attenuation_N09_attenuate_values(x0, gamma, ampl, slope, Av):
    # get the correct values
    x, cor_vals = get_axav_cor_vals(x0, gamma, ampl, slope, Av)

    # calculate the cor_vals in fractional units
    cor_vals = np.power(10.0, -0.4 * (cor_vals))

    # initialize model
    tmodel = N09(x0=x0, gamma=gamma, ampl=ampl, slope=slope, Av=Av)

    # test
    np.testing.assert_allclose(tmodel.attenuate(x), cor_vals[::-1], atol=1e-6)
