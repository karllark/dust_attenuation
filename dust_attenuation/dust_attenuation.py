# -*- coding: utf-8 -*-
# Main module for dust_attenuation


import numpy as np
import astropy.units as u

from astropy.modeling import (Model, Fittable1DModel,
                              Parameter, InputParameterError)


__all__ = ['C00']

x_range_C00 = [1.0/2.2, 1.0/0.12]

def _test_valid_x_range(x, x_range, outname):
    """
    Test if any of the x values are outside of the valid range

    Parameters
    ----------
    x : float array
       wavenumbers in inverse microns

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
                         + ', x has units 1/micron]')



class C00(object):
    """ Class for computing the attenuation curve of Calzetti et al. (2000)

    """
  
    def __init__(self,wavelength):
        """
        wavelength: np.array of floats
             Wavelength in microns.
        """
        self.wvl = np.array(wavelength)
         

    def ax_Ebv(self):
        """Compute the starburst reddenig curve of Calzetti et al. (2000)
            A(λ)/E(B-V)_continuum between 0.12 and 2.2 microns.
     
        Parameters
        ----------

        Returns
        -------
        axEbv: np.array of floats
             reddening curve

        """
        axEbv = np.zeros(len(self.wvl))

        # Attenuation between 0.12 and 0.63 microns
        mask = (self.wvl >= 0.12) & (self.wvl < 0.63)
        axEbv[mask] = 2.659 * (-2.156 + 
                               1.509 / self.wvl[mask] -
                               0.198 / self.wvl[mask] ** 2 +
                               0.011 / self.wvl[mask] ** 3) + 4.05

        # Attenuation between 0.63 and 2.2 microns
        mask = (self.wvl >= 0.63) & (self.wvl < 2.2)
        axEbv[mask] = 2.659 * (-1.857 + 1.040 / self.wvl[mask]) + 4.05

        return axEbv


    def ax_Av(self):
        """ Compute the normalised attenuation curve of Calzetti et al. (2000)

        Parameters
        ---------- 

        Returns
        -------
        axAv: np.array of floats
             normalised attenuation curve

        """

        axAv = self.ax_Ebv() / 4.05

        return axAv

    def ax(self, Ebv_s):
        """ Compute the normalised attenuation curve of Calzetti et al. (2000)

        Parameters
        ---------- 
        Ebv_s: float 
           color excess of the stellar continuum (E(B-V))

        Returns
        -------
        axAv: np.array of floats
             normalised attenuation curve

        """

        ax = self.ax_Ebv() * Ebv_s

        return ax

