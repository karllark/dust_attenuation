# -*- coding: utf-8 -*-
# Main module for dust_attenuation


import numpy as np
import astropy.units as u

from astropy.modeling import (Model, Fittable1DModel,
                              Parameter, InputParameterError)


__all__ = ['C00']

x_range_C00 = [0.12, 2.2]

def _test_valid_x_range(x, x_range, outname):
    """
    Test if any of the x values are outside of the valid range

    Parameters
    ----------
    x : float array
       wavelength in microns

    x_range: 2 floats
       allowed min/max of x

    outname: str
       name of curve for error message
    """
    if np.logical_or(np.any(x < x_range[0]),
                     np.any(x > x_range[1])):
        raise ValueError('Input x outside of range defined for ' + outname
                         + ' ['
                         + str(x_range[0])
                         + ' <= x <= '
                         + str(x_range[1])
                         + ', x has units micron]')

class BaseAttModel(Fittable1DModel):
    """
    Base Attenuation Model.  Do not use.
    """
    inputs = ('x',)
    outputs = ('axav',)

    def attenuated(self, x, Av=None, Ebv=None):
        """
        Calculate the attenuation as a fraction

        Parameters
        ----------
        x: float
           expects either x in units of wavelengths or frequency
           or assumes wavelengths in [micron]

           internally microns are used

        Av: float
           A(V) value of dust column
           Av or Ebv must be set

        Ebv: float
           E(B-V) value of dust column
           Av or Ebv must be set

        Returns
        -------
        frac_att: np array (float)
           fractional attenuation as a function of x
        """
        # get the attenuation curve
        axav = self(x)

        # check that av or ebv is set
        if (Av is None) and (Ebv is None):
            raise InputParameterError('neither Av or Ebv passed, one required')

        # if Av is not set and Ebv set, convert to Av
        if Av is None:
            Av = self.Rv*Ebv

        # return fractional attenuation
        return np.power(10.0, -0.4*axav*Av)


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



