# -*- coding: utf-8 -*-
# Main module for dust_attenuation


import numpy as np
import astropy.units as u

from .base_classes import BaseAttAvModel, _test_valid_x_range

from astropy.modeling import (Model, Fittable1DModel,
                              Parameter, InputParameterError)


__all__ = ['S07']

x_range_S07 = [1.0/2.5, 1.0/30.]

class S07(Fittable1DModel):
    """
    S07 kvt extinction model calculation

    Parameters
    ----------
    kvt_amp : float
      amplitude

    Notes
    -----
    S07 extinction model

    From Kemper, Vriend, & Tielens (2004)

    Applicable for Mid-Infrared

    A custom extinction curve constructed by two components:
    silicate profile & exponent 1.7 power-law.

    Example showing a S07 curve with components identified.

    .. plot::
        :include-source:

        import numpy as np
        import matplotlib.pyplot as plt
        import astropy.units as u

        from dust_extinction.dust_extinction import S07

    """
    inputs = ('x',)
    outputs = ('tau')#('axav',)

    kvt_wav = [8.0, 8.2, 8.4, 8.6, 8.8, 9.0, 9.2, 9.4, 9.6, 9.7, 
               9.75, 9.8, 10.0, 10.2, 10.4, 10.6,10.8, 11.0, 11.2, 11.4, 
               11.6, 11.8, 12.0, 12.2, 12.4, 12.6, 12.7]

    kvt_int = [.06, .09, .16, .275, .415, .575, .755, .895, .98, .99, 
            1.0, .99, .94, .83, .745, .655, .58, .525, .43, .35, 
            .27, .20, .13, .09, .06, .045, .04314]

    amp = Parameter(description="kvt term: amplitude", default=1, min=0.0)

    @staticmethod
    def drude(in_lam, cen_wave, inten, frac_fwhm):
        return inten*frac_fwhm**2 / ((in_lam/cen_wave-cen_wave/in_lam)**2 + frac_fwhm**2)

    # Extend kvt profile to shorter wavelengths
    kvt_wav_short = in_x[in_x < min(kvt_wav)]
    kvt_int_short = min(kvt_int) * np.exp(2.03*(wav_short-min(kvt_wav)))

    # Extend kvt profile to longer wavelengths
    kvt_wav_long = in_x[in_x > max(kvt_wav)]
    kvt_int_long = np.zeros(len(kvt_wav_long)) # Need to input the value instead of using zeros

    spline_x = np.concatenate([kvt_wav_short, kvt_wav, kvt_wav_long])
    spline_y = np.concatenate([kvt_int_short, kvt_int, kvt_int_long])

    spline_rep = interpolate.splrep(spline_x, spline_y)
    new_spline_y = interpolate.splev(in_x, spline_rep, der=0)

    def evaluate(self, in_x, amp):
        ext = self.drude(in_x, 18, 0.4, 0.247) + new_spline_y
        
        # assuing beta is 0.1 
        beta = 0.1
        tau = amp * ((1.0-beta)*ext + beta*(9.7/in_x)**1.7)
        return tau