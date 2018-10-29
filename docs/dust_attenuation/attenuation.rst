.. _AttvsExt:

#############################
Attenuation versus Extinction
#############################

Attenuation
===========

Dust attenuation refers to the general impact on the spectrum of an object due
to the presence of dust.  In general, attenuation is used to indicate that the
geometry of the sources and dust in a system is more complex than a single star
with a foreground screen of dust.  Examples of such systems include dusty
galaxies (composed of many stars) and stars with circumstellar dust.

Measuring attenuation for an object is nominally as straightforward as measuring
extinction.  The spectrum of the object with dust is compared to an equivalent
object with no dust.  In reality, this is quite hard as the usual objects for
such work are galaxies and no two galaxies are similar enough.  This can and has
been done for samples of galaxies (e.g., starbursts: Calzetti et al. 2000) or by
using a model for the stellar populations in that galaxy.

Theoretically calculating the attenuation curves can be done using dust
radiative transfer models.  As a result, the effects of multiple sources  and
scattering into the observed beam are sometimes referred to as radiative
transfer effects.  Such radiative transfer calculations require specifying the
geometry of the photon emitters (stars) and dust grains as well as the grain
details (size, shape, composition).

TBD: add in 2 star example from SFChapter illustrating the two effects
in the simplest case.

Attenuation includes the effects of having multiple sources extinguished by
different columns of dust and the scattering of photons into the observation
beam.  Unlike extinction, these two effects mean that the attenuation is not
directly proportional to the amount dust in the system. Hence the ratio of
attenuations at two different wavelengths *varies* with the amount total system
dust. This is illustrated below with `WG00` shell, clumpy, mw models where the
left plot shows the total attenuation as a function of wavelength and the right
plot shows the same curves normalized by Att(V).

.. plot::

      import numpy as np
      import matplotlib.pyplot as plt
      import astropy.units as u

      from dust_attenuation.radiative_transfer import WG00

      fig, ax = plt.subplots(ncols=2)

      # generate the curves and plot them
      ix = np.arange(1.0/3.0,1.0/0.1,0.1)/u.micron
      x = 1./ix

      # defined for normalization
      x_Vband = 0.55

      att_model = WG00(tau_V = 0.5, geometry = 'shell',
                       dust_type = 'mw', dust_distribution = 'clumpy')
      ax[0].plot(x,att_model(x), label = r'$\tau(V) = 0.5$')
      ax[1].plot(x,att_model(x)/att_model(x_Vband), label = r'$\tau(V) = 0.5$')

      att_model = WG00(tau_V = 1.5, geometry = 'shell',
                       dust_type = 'mw', dust_distribution = 'clumpy')
      ax[0].plot(x,att_model(x), label = r'$\tau(V) = 1.5$')
      ax[1].plot(x,att_model(x)/att_model(x_Vband), label = r'$\tau(V) = 1.5$')

      att_model = WG00(tau_V = 2.5, geometry = 'shell',
                       dust_type = 'mw', dust_distribution = 'clumpy')
      ax[0].plot(x,att_model(x), label = r'$\tau(V) = 2.5$')
      ax[1].plot(x,att_model(x)/att_model(x_Vband), label = r'$\tau(V) = 2.5$')

      ax[0].set_title('Total Attenuation')
      ax[1].set_title('Normalized Attenuation')
      ax[0].set_xlabel(r'$\lambda$ [$\mu m$]')
      ax[1].set_xlabel(r'$\lambda$ [$\mu m$]')
      ax[0].set_ylabel(r'$Att(\lambda)$')
      ax[1].set_ylabel(r'$Att(\lambda)/Att(V)$')

      ax[0].set_xscale('log')
      ax[1].set_xscale('log')
      ax[0].set_xlim(0.09,4.0)
      ax[1].set_xlim(0.09,4.0)

      ax[0].legend(loc='best')
      ax[1].legend(loc='best')
      plt.tight_layout()
      plt.show()

Extinction
==========

Interstellar dust extinction is the result of photons being absorbed or
scattered *out* of the line-of-sight by dust grains.  Extinction is
explicitly linked to the specific geometry of a single star observed
through a column of dust.

Both dust absorption and scattering out of the line-of-sight are processes
that are directly proportional to the amount of dust along the line-of-sight.
As a result, the ratio of dust extinctions at two different wavelengths
does not vary with different amounts of dust.  This makes the measurement
and/or theoretical calculation of extinction much simpler than the more
general case of attenuation.

The separate package `dust_extinction package
<http://dust-extinction.readthedocs.io/>`_ exists to provide extinction
models.

Note: all extinction curves are attenuation curves, but not all attenuation
curves are extinction curves.
