# -*- coding: utf-8 -*-

import numpy as np
import astropy.units as u

from .baseclasses import BaseAttAvModel
from .helpers import _test_valid_x_range

from .averages import C00, L02
from astropy.modeling import Parameter, InputParameterError

__all__ = ['N09', 'SBL18']

x_range_N09 = [0.097, 2.2]
x_range_SBL18 = [0.097, 2.2]


class N09(BaseAttAvModel):
    """
    Attenuation curve using a modified version of the Calzetti law
    allowing for a varying UV slope and the presence of a UV bump.

    Parameters
    ----------
    x0: float
        Central wavelength of the UV bump (in microns).

    gamma: float
        Width (FWHM) of thhe UV bump (in microns).

    ampl: float
        Amplitude of the UV bump.

    slope: float
        Slope of the power law.

    Av: float
        attenuation in V band.

    Raises
    ------
    InputParameterError
       Input Av values outside of defined range

    Notes
    -----

    The original formalism is from Noll et al, A&A 507, 1793–1813 (2009).

    Example:

    .. plot::
        :include-source:

        import matplotlib.pyplot as plt
        import numpy as np
        import astropy.units as u

        from dust_attenuation.shapes import N09

        fig, ax = plt.subplots()

        # generate the curves and plot them
        x = np.arange(0.5,10,0.1)/u.micron

        slopes = [-1, -0.5, 0, 0.5, 1]
        for slope in slopes:
            att_model = N09(Av=1,ampl=3.5,slope=slope)
            ax.plot(x,att_model(1/x),label=r'$\delta$ = %.2f' % (slope))

        ax.set_xlabel('$x$ [$\mu m^{-1}$]')
        ax.set_ylabel('A(x) [mag]')

        ax.legend(loc='best')
        plt.title("Noll09 with varying slopes")
        plt.show()
    """

    x_range = x_range_N09

    # Did not want to create a new base class only for this particular model.
    # So parameters are defined here.

    x0 = Parameter(description="bump: centroid", default=0.2175, min=0)

    gamma = Parameter(description="bump: width (FWHM)", default=0.035, min=0)

    ampl = Parameter(description="bump: amplitude ", default=0, min=0)

    slope = Parameter(description="slope: slope of the power law",
                            default=0., min=-3., max=3.)

    # Had to redefine the Av parameter though it is defined in the parent class.
    # May be beause new Parameters are defined here?...

    Av = Parameter(description="Av: attenuation in V band ",
                   default=1.0, min=0.0)

    @x0.validator
    def x0(self, value):
        """
        Check that x0 is in the valid range

        Parameters
        ----------
        value: float
            x0 value to check

        Raises
        ------
        InputParameterError
           Input x0 values outside of defined range
        """

        if (value < 0.0):
            raise InputParameterError("parameter x0 must be positive")

    @gamma.validator
    def gamma(self, value):
        """
        Check that gamma is in the valid range

        Parameters
        ----------
        value: float
            gamma value to check

        Raises
        ------
        InputParameterError
           Input gamma values outside of defined range
        """

        if (value < 0.0):
            raise InputParameterError("parameter gamma must be positive")

    @ampl.validator
    def ampl(self, value):
        """
        Check that ampl is in the valid range

        Parameters
        ----------
        value: float
            ampl value to check

        Raises
        ------
        InputParameterError
           Input ampl values outside of defined range
        """

        if (value < 0.0):
            raise InputParameterError("parameter ampl must be positive")

    @slope.validator
    def slope(self, value):
        """
        Check that the slope is in the valid range

        Parameters
        ----------
        value: float
            slope value to check

        Raises
        ------
        InputParameterError
           Input slope values outside of defined range
        """

        if (value < -3.0) or (value > 3.0):
            raise InputParameterError("parameter slope must be between "
                                      "-3.0 and 3.0")

    @Av.validator
    def Av(self, value):
        """
        Check that Av is in the valid range

        Parameters
        ----------
        value: float
            Av value to check

        Raises
        ------
        InputParameterError
           Input Av values outside of defined range
        """

        if (value < 0.0):
            raise InputParameterError("parameter Av must be positive")


    # Rv from Calzetti 2000
    Rv_C00 = 4.05


    def uv_bump(self, x, x0, gamma, ampl):
        """
        Drude profile for computing the UV bump.

        Parameters
        ----------
        x: np array (float)
           expects wavelengths in [micron]

        x0: float
           Central wavelength of the UV bump (in microns).

        gamma: float
           Width (FWHM) of the UV bump (in microns).

        ampl: float
           Amplitude of the UV bump.

        Returns
        -------
        np array (float)
           lorentzian-like Drude profile

        Raises
        ------
        ValueError
           Input x values outside of defined range

        """
        return ampl * (x**2 * gamma**2 /
                            ((x**2 - x0**2)**2 + x**2 * gamma**2))


    def power_law(self, x, slope):
        """ Power law normalised at 0.55 microns (V band).

        Parameters
        ----------
        x: np array (float)
           expects wavelengths in [micron]

        slope: float
           slope of the power law

        Returns
        -------
        powlaw: np array (float)
           power law
        """

        return (x / 0.55)**slope


    def k_lambda(self, x, x0, gamma, ampl, slope):
        """ Compute the starburst reddening curve k'(λ)=A(λ)/E(B-V)
        using recipe of Calzetti 2000 and Leitherer 2002

        Parameters
        ----------
        in_x: np array (float)
           expects either x in units of wavelengths or frequency
           or assumes wavelengths in [micron]
           internally microns are used

        x0: float
           Central wavelength of the UV bump (in microns).

        gamma: float
           Width (FWHM) of thhe UV bump (in microns).

        ampl: float
           Amplitude of the UV bump.

        slope: float
           Slope of the power law.

        Returns
        -------
        k_lambda: np array (float)
           k_lambda(x) reddening curve

        Raises
        ------
        ValueError
           Input x values outside of defined range

        """
        # convert to wavenumbers (1/micron) if x input in units
        # otherwise, assume x in appropriate wavenumber units
        with u.add_enabled_equivalencies(u.spectral()):
            x_quant = u.Quantity(x, u.micron, dtype=np.float64)

        # strip the quantity to avoid needing to add units to all the
        #    polynomical coefficients
        x = x_quant.value

        # check that the wavenumbers are within the defined range
        _test_valid_x_range(x, x_range_N09, 'N09')

        # setup the axEbv vectors
        axEbv = np.zeros(len(x))

        # Compute reddening curve using Calzetti 2000
        mask_C00 = x > 0.15
        axEbv[mask_C00] = C00().k_lambda(x[mask_C00])

        # Use recipe of Leitherer 2002 below 0.15 microns
        mask_L02 = x <= 0.15
        axEbv[mask_L02] = L02().k_lambda(x[mask_L02])

        # Add the UV bump using the Drude profile
        axEbv += self.uv_bump(x, x0, gamma, ampl)

        # Multiply the reddening curve with a power law with varying slope
        axEbv *= self.power_law(x, slope)

        return axEbv


    def evaluate(self, x, x0, gamma, ampl, slope, Av):
        """
        C00 function

        Parameters
        ----------
        x: np array (float)
           expects either x in units of wavelengths or frequency
           or assumes wavelengths in [micron]
           internally microns are used

        x0: float
           Central wavelength of the UV bump (in microns).

        gamma: float
           Width (FWHM) of thhe UV bump (in microns).

        ampl: float
           Amplitude of the UV bump.

        slope: float
           Slope of the power law.

        Av: float
           attenuation in V band.

        Returns
        -------
        att: np array (float)
            Att(x) attenuation curve [mag]

        Raises
        ------
        ValueError
           Input x values outside of defined range
        """

        axEbv = self.k_lambda(x, x0, gamma, ampl, slope)
        ax = axEbv / self.Rv_C00 * Av

        return ax



class SBL18(N09):
    """
    Attenuation curve using a modified version of the Calzetti law
    allowing for a varying UV slope and the presence of a UV bump.

    Parameters
    ----------
    x0: float
        Central wavelength of the UV bump (in microns).

    gamma: float
        Width (FWHM) of thhe UV bump (in microns).

    ampl: float
        Amplitude of the UV bump.

    slope: float
        Slope of the power law.

    Av: float
        attenuation in V band.

    Raises
    ------
    InputParameterError
       Input Av values outside of defined range

    Notes
    -----

    Modification of N09 formalism: the UV bump was added before applying
    the power law correction in N09, in this class the UV bump is now added
    after the power law correction.
    This modification is first mentionned in in ApJ, Volume 859, Issue 1,
    article id. 11, 17 pp. (2018)

    Example:

    .. plot::
        :include-source:

        import matplotlib.pyplot as plt
        import numpy as np
        import astropy.units as u

        from dust_attenuation.shapes import SBL18

        fig, ax = plt.subplots()

        # generate the curves and plot them
        x = np.arange(0.5,10,0.1)/u.micron

        slopes = [-1, -0.5, 0, 0.5, 1]
        for slope in slopes:
            att_model = SBL18(Av=1,ampl=3.5,slope=slope)
            ax.plot(x,att_model(1/x),label=r'$\delta$ = %.2f' % (slope))

        ax.set_xlabel('$x$ [$\mu m^{-1}$]')
        ax.set_ylabel('A(x) [mag]')

        ax.legend(loc='best')
        plt.title("SBL18 with varying slopes")
        plt.show()
    """

    def k_lambda(self, x, x0, gamma, ampl, slope):
        """ Compute the starburst reddening curve k'(λ)=A(λ)/E(B-V)
        using recipe of Calzetti 2000 and Leitherer 2002

        Parameters
        ----------
        in_x: np array (float)
           expects either x in units of wavelengths or frequency
           or assumes wavelengths in [micron]
           internally microns are used

        x0: float
           Central wavelength of the UV bump (in microns).

        gamma: float
           Width (FWHM) of thhe UV bump (in microns).

        ampl: float
           Amplitude of the UV bump.

        slope: float
           Slope of the power law.

        Returns
        -------
        k_lambda: np array (float)
           k_lambda(x) reddening curve

        Raises
        ------
        ValueError
           Input x values outside of defined range

        """
        # convert to wavenumbers (1/micron) if x input in units
        # otherwise, assume x in appropriate wavenumber units
        with u.add_enabled_equivalencies(u.spectral()):
            x_quant = u.Quantity(x, u.micron, dtype=np.float64)

        # strip the quantity to avoid needing to add units to all the
        #    polynomical coefficients
        x = x_quant.value

        # check that the wavenumbers are within the defined range
        _test_valid_x_range(x, x_range_SBL18, 'SBL18')

        # setup the axEbv vectors
        axEbv = np.zeros(len(x))

        # Compute reddening curve using Calzetti 2000
        mask_C00 = x > 0.15
        axEbv[mask_C00] = C00().k_lambda(x[mask_C00])

        # Use recipe of Leitherer 2002 below 0.15 microns
        mask_L02 = x <= 0.15
        axEbv[mask_L02] = L02().k_lambda(x[mask_L02])

        # Multiply the reddening curve with a power law with varying slope
        axEbv *= self.power_law(x, slope)

        # Add the UV bump using the Drude profile
        axEbv += self.uv_bump(x, x0, gamma, ampl)

        return axEbv
