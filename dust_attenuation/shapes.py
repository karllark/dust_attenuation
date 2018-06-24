# -*- coding: utf-8 -*-

import numpy as np
import astropy.units as u

from .baseclasses import BaseAttAvModel
from .helpers import _test_valid_x_range

from .averages import C00, Leitherer02
from astropy.modeling import Parameter

__all__ = ['Calzmod']

x_range_Calzmod = [0.097, 2.2]


class Calzmod(BaseAttAvModel):
    """
    Attenuation curve using a modified version of the Calzett law
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

    The original formalism was done by Noll et al, A&A 507, 1793–1813 (2009).
    However the UV bump was added before applying the power law correction.
    In this class we use the version of Salim et al 2018 (The Astrophysical
    Journal, Volume 859, Issue 1, article id. 11, 17 pp.) where the UV bump is
    added after the power law correction.

    Example:

    .. plot::
        :include-source:

        import matplotlib.pyplot as plt
        import numpy as np
        import astropy.units as u

        from dust_attenuation.shapes import Calzmod

        fig, ax = plt.subplots()

        # generate the curves and plot them
        x = np.arange(0.5,10,0.1)/u.micron

        slopes = [-1, -0.5, 0, 0.5, 1]
        for slope in slopes:
            att_model = Calzmod(Av=1,ampl=3.5,slope=slope)
            ax.plot(x,att_model(1/x),label=r'$\delta$ = %.2f' % (slope))

        ax.set_xlabel('$x$ [$\mu m^{-1}$]')
        ax.set_ylabel('A(x) [mag]')

        ax.legend(loc='best')
        plt.title("Calzmod with varying slopes")
        plt.show()
    """
 
    x_range = x_range_Calzmod
     
    # Did not want to create a new base class only for this particular model.
    # So parameters are defined here.

    x0 = Parameter(description="bump: centroid", default=0.2175, min=0)

    gamma = Parameter(description="bump: width (FWHM)", default=0.035, min=0)

    ampl = Parameter(description="bump: amplitude ", default=0, min=0)

    slope = Parameter(description="slope: slope of the power law",
                            default=0., min=-2., max=2.)

    # Had to redefine the Av parameter though it is defined in the parent class.
    # May be beause new Parameters are defined here?...

    Av = Parameter(description="Av: attenuation in V band ",
                   default=1.0, min=0.0)

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


    """
    def __init__(self):

        self.Rv = 4.05

        # In Python 2: super(C00, self) 
        # In Python 3: super() but super(C00, self) still works
        super(Calzmod, self).__init__()

    """
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
        _test_valid_x_range(x, x_range_Calzmod, 'Calzmod')

        # setup the axEbv vectors
        axEbv = np.zeros(len(x))

        # Compute reddening curve using Calzetti 2000
        # Wavelength range taken from Salim+18
        mask_C00 = x > 0.15
        axEbv[mask_C00] = C00().k_lambda(x[mask_C00])

        # Use recipe of Leitherer 2002 below 0.15 microns
        # Wavelength range taken from Salim+18
        mask_leit02 = x <= 0.15
        axEbv[mask_leit02] = Leitherer02().k_lambda(x[mask_leit02])

        # Compute the new Rv depending on the value of the power law slope
        #Rv_mod = self.Rv_C00 / ((self.Rv_C00 + 1) * (0.44/0.55)**slope - self.Rv_C00)
        # Correct the reddening curve so that E(B-V) remains unchanged when slope !=0
        #axEbv *= Rv_mod / self.Rv_C00
        # Note: this does not seem to work well. Need to check later.

        # Multiply the reddening curve with a power law with varying slope
        axEbv *= self.power_law(x, slope)

        # Add the UV bump using the Drude profile
        axEbv += self.uv_bump(x, x0, gamma, ampl)

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
        # convert to wavenumbers (1/micron) if x input in units
        # otherwise, assume x in appropriate wavenumber units
        with u.add_enabled_equivalencies(u.spectral()):
            x_quant = u.Quantity(x, u.micron, dtype=np.float64)

        # strip the quantity to avoid needing to add units to all the
        #    polynomical coefficients
        x = x_quant.value

        # check that the wavenumbers are within the defined range
        _test_valid_x_range(x, x_range_Calzmod, 'Calzmod')

        axEbv = self.k_lambda(x, x0, gamma, ampl, slope)
        ax = axEbv / self.Rv_C00 * Av

        return ax
