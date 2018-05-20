# -*- coding: utf-8 -*-

import numpy as np
import astropy.units as u

from .base_classes import (BaseAttAvModel, _test_valid_x_range)

__all__ = ['C00']

x_range_C00 = [0.12, 2.2]


class C00(BaseAttAvModel):
    """
    Attenuation curve of Calzetti et al. (2000)

    Parameters
    ----------
    Av: float
        attenuation in V band

    Raises
    ------
    InputParameterError
       Input Av values outside of defined range

    Notes
    -----

    From Calzetti (2000, ApJ, Volume 533, Issue 2, pp. 682-695)

    Example:

    .. plot::
        :include-source:

        import numpy as np
        import matplotlib.pyplot as plt
        import astropy.units as u

        from dust_attenuation.averages import C00

        fig, ax = plt.subplots()

        # generate the curves and plot them
        x = np.arange(0.12,2.2,0.1)*u.micron

        Avs = [0.1,0.5,1.0,2.0,5.0]
        for cur_Av in Avs:
           att_model = C00(Av=cur_Av)
           ax.plot(1/x,att_model(x),label=r'A$_V$ = %.2f mag' % (cur_Av))

        ax.set_xlabel('$x$ [$\mu m^{-1}$]')
        ax.set_ylabel('$A(x)$ [mag]')

        ax.legend(loc='best')
        plt.show()
    """

    x_range = x_range_C00

    @staticmethod
    def evaluate(x, Av):
        """
        C00 function

        Parameters
        ----------
        in_x: float
           expects either x in units of wavelengths or frequency
           or assumes wavelengths in [micron]

           internally microns are used

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
        _test_valid_x_range(x, x_range_C00, 'C00')

        # setup the ax vectors
        n_x = len(x)
        axEbv = np.zeros(n_x)

        # Rv is fixed to 4.05
        Rv = 4.05

        # define the ranges
        uv2vis_indxs = np.where(np.logical_and(0.12 <= x, x < 0.63))
        nir_indxs = np.where(np.logical_and(0.63 <= x, x < 2.2))

        axEbv[uv2vis_indxs] = (2.659 * (-2.156 +
                                        1.509 * 1 / x[uv2vis_indxs] -
                                        0.198 * 1 / x[uv2vis_indxs] ** 2 +
                                        0.011 * 1 / x[uv2vis_indxs] ** 3) + Rv)

        axEbv[nir_indxs] = 2.659 * (-1.857 + 1.040 * 1 / x[nir_indxs]) + Rv

        ax = axEbv / Rv * Av

        return ax
