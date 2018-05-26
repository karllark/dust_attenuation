# -*- coding: utf-8 -*-
import numpy as np

from astropy.modeling import (Fittable1DModel,
                              Parameter,
                              InputParameterError)

__all__ = ['BaseAttModel', 'BaseAttAvModel', 'BaseAtttauVModel']


class BaseAttModel(Fittable1DModel):
    """
    Base Attenuation Model.  Do not use.
    """
    inputs = ('x',)
    outputs = ('ax',)

    def attenuate(self, x):
        """
        Calculate the attenuation as a fraction

        Parameters
        ----------
        x: float
           expects either x in units of wavelengths or frequency
           or assumes wavelengths in [micron]

           internally microns are used

        Returns
        -------
        frac_att: np array (float)
           fractional attenuation as a function of x
        """
        # get the attenuation curve
        ax = self(x)

        # return fractional attenuation
        return np.power(10.0, -0.4*ax)


class BaseAttAvModel(BaseAttModel):
    """
    Base attenuation Av Model.  Do not use.
    """
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


class BaseAtttauVModel(BaseAttModel):
    """
    Base attenuation tau_V Model.  Do not use.
    """
    tau_V = Parameter(description="tau_V: optical depth in V band ",
                      default=1.0, min=0.25, max=50.00)

    @tau_V.validator
    def tau_V(self, value):
        """
        Check that tau_V is in the valid range

        Parameters
        ----------
        value: float
            tau_V value to check

        Raises
        ------
        InputParameterError
           Input tau_V values outside of defined range
        """
        if not (self.tau_V_range[0] <= value <= self.tau_V_range[1]):
            raise InputParameterError("parameter tau_V must be between "
                                      + str(self.tau_V_range[0])
                                      + " and "
                                      + str(self.tau_V_range[1]))
