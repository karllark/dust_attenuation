#############
Model Flavors
#############

There are three different types of models: averages, radiative transfer based,
and shape fitting.

Average models
==============

These models provide averages from the literature with the ability to
interpolate between the observed data points.  In general, these average
models have shapes that are not dependent on the amount of dust.

The `COO` average attenuation model is based on a small number of
starburst galaxies observed
in the ultraviolet with the International Ultraviolet Explorer (IUE)
supplemented with ground-based optical spectroscopy,
near-infrared photometry, and
Infrared Space Observatory (ISO) far-infrared photometry
(Calzetti et al. 2000).

.. plot::

      import math
      import numpy as np
      import matplotlib.pyplot as plt
      import astropy.units as u

      from dust_attenuation.averages import C00

      fig, ax = plt.subplots()

      # generate the curves and plot them
      ix = np.arange(1.0/2.2, 1.0/0.12 , 0.1)/u.micron
      x = 1/ix
      att_model = C00(Av = 1.0)
      ax.plot(x,att_model(x), label = 'C00')

      ax.set_xlabel('$\lambda$ [$\mu m$]')
      ax.set_ylabel('$Att(\lambda)/Att(V)$')

      ax.set_xscale('log')
      ax.set_xlim(0.1, 3.0)

      ax.legend(loc='best')
      plt.tight_layout()
      plt.show()


Radiative Transfer Based Models
===============================

These models provide attenuation predictions based on dust radiative transfer
calculations.  The attenuation curve strength and wavelength dependent shape
are based on the amount of dust, star/dust geometry, and other
parameters.

Witt & Gordon 2000 (WG00)
-------------------------

The `WG00` attenuation models are based on DIRTY radiative transfer
calculations for spherical galactic environments (shell, dusty, cloudy)
with homogeneous or clumpy local dust distributions using
empirical Milky Way (MW) and Small Magellanic Cloud (SMC)
dust grain properties (Witt & Gordon 2000).
The `WG00` models were chosen to span the range of
possible star/dust geometries and types of dust grains.

Example `WG00` models showing variation in shape with amount of dust.

.. plot::

      import numpy as np
      import matplotlib.pyplot as plt
      import astropy.units as u

      from dust_attenuation.radiative_transfer import WG00

      fig, ax = plt.subplots()

      # generate the curves and plot them
      ix = np.arange(1.0/3.0,1.0/0.1,0.1)/u.micron
      x = 1./ix

      # defined for normalization
      x_Vband = 0.55

      att_model = WG00(tau_V = 0.25, geometry = 'shell',
                       dust_type = 'mw', dust_distribution = 'clumpy')
      ax.plot(x,att_model(x)/att_model(x_Vband), label = r'$\tau(V) = 0.25$')

      att_model = WG00(tau_V = 1.0, geometry = 'shell',
                       dust_type = 'mw', dust_distribution = 'clumpy')
      ax.plot(x,att_model(x)/att_model(x_Vband), label = r'$\tau(V) = 1.0$')

      att_model = WG00(tau_V = 5.0, geometry = 'shell',
                       dust_type = 'mw', dust_distribution = 'clumpy')
      ax.plot(x,att_model(x)/att_model(x_Vband), label = r'$\tau(V) = 5.0$')

      att_model = WG00(tau_V = 50.0, geometry = 'shell',
                       dust_type = 'mw', dust_distribution = 'clumpy')
      ax.plot(x,att_model(x)/att_model(x_Vband), label = r'$\tau(V) = 50.0$')

      ax.set_xlabel(r'$\lambda$ [$\mu m$]')
      ax.set_ylabel(r'$Att(\lambda)/Att(V)$')

      ax.set_xscale('log')
      ax.set_xlim(0.09,4.0)

      ax.set_title('WG00 Shell, clumpy, MW')

      ax.legend(loc='best')
      plt.tight_layout()
      plt.show()

Example `WG00` models showing shape variation with different types of
dust grains.

.. plot::

      import numpy as np
      import matplotlib.pyplot as plt
      import astropy.units as u

      from dust_attenuation.radiative_transfer import WG00

      fig, ax = plt.subplots()

      # generate the curves and plot them
      ix = np.arange(1.0/3.0,1.0/0.1,0.1)/u.micron
      x = 1./ix

      # defined for normalization
      x_Vband = 0.55

      att_model = WG00(tau_V = 1.0, geometry = 'shell',
                       dust_type = 'mw', dust_distribution = 'clumpy')
      ax.plot(x,att_model(x)/att_model(x_Vband), label = 'MW')

      att_model = WG00(tau_V = 1.0, geometry = 'shell',
                       dust_type = 'smc', dust_distribution = 'clumpy')
      ax.plot(x,att_model(x)/att_model(x_Vband), label = 'SMC')

      ax.set_xlabel(r'$\lambda$ [$\mu m$]')
      ax.set_ylabel(r'$Att(\lambda)/Att(V)$')

      ax.set_xscale('log')
      ax.set_xlim(0.09,4.0)

      ax.set_title(r'WG00 Shell, clumpy, $\tau(V) = 1.0$')

      ax.legend(loc='best')
      plt.tight_layout()
      plt.show()


Example `WG00` models showing shape variation with different spherical galactic
environments.

.. plot::

      import numpy as np
      import matplotlib.pyplot as plt
      import astropy.units as u

      from dust_attenuation.radiative_transfer import WG00

      fig, ax = plt.subplots()

      # generate the curves and plot them
      ix = np.arange(1.0/3.0,1.0/0.1,0.1)/u.micron
      x = 1./ix

      # defined for normalization
      x_Vband = 0.55

      att_model = WG00(tau_V = 1.0, geometry = 'shell',
                       dust_type = 'mw', dust_distribution = 'clumpy')
      ax.plot(x,att_model(x)/att_model(x_Vband), label = 'Shell')

      att_model = WG00(tau_V = 1.0, geometry = 'dusty',
                       dust_type = 'mw', dust_distribution = 'clumpy')
      ax.plot(x,att_model(x)/att_model(x_Vband), label = 'Dusty')

      att_model = WG00(tau_V = 1.0, geometry = 'cloudy',
                       dust_type = 'mw', dust_distribution = 'clumpy')
      ax.plot(x,att_model(x)/att_model(x_Vband), label = 'Cloudy')

      ax.set_xlabel(r'$\lambda$ [$\mu m$]')
      ax.set_ylabel(r'$Att(\lambda)/Att(V)$')

      ax.set_xscale('log')
      ax.set_xlim(0.09,4.0)

      ax.set_title(r'WG00, clumpy, $\tau(V) = 1.0$')

      ax.legend(loc='best')
      plt.tight_layout()
      plt.show()

Example `WG00` models showing shape variation with local dust distributions.

.. plot::

      import numpy as np
      import matplotlib.pyplot as plt
      import astropy.units as u

      from dust_attenuation.radiative_transfer import WG00

      fig, ax = plt.subplots()

      # generate the curves and plot them
      ix = np.arange(1.0/3.0,1.0/0.1,0.1)/u.micron
      x = 1./ix

      # defined for normalization
      x_Vband = 0.55

      att_model = WG00(tau_V = 1.0, geometry = 'shell',
                       dust_type = 'mw', dust_distribution = 'homogeneous')
      ax.plot(x,att_model(x)/att_model(x_Vband),label = 'homogeneous')

      att_model = WG00(tau_V = 1.0, geometry = 'dusty',
                       dust_type = 'mw', dust_distribution = 'clumpy')
      ax.plot(x,att_model(x)/att_model(x_Vband),label = 'clumpy')

      ax.set_xlabel(r'$\lambda$ [$\mu m$]')
      ax.set_ylabel(r'$Att(\lambda)/Att(V)$')

      ax.set_xscale('log')
      ax.set_xlim(0.09,4.0)

      ax.set_title(r'WG00, Shell, $\tau(V) = 1.0$')

      ax.legend(loc='best')
      plt.tight_layout()
      plt.show()

Shape fitting models
====================

These models allow for more arbitrary shapes to be modeled than the
other model flavors.

Calzmod: modified Calzetti law of Noll09
----------------------------------------

Noll+09 first introduced a modified version of the Calzetti 2000 law, allowing
for a varying slope and the presence of a UV bump.

Example `Calzmod` models showing variation in slopes.
A UV bump with an amplitude of 3.5 is added to the C00 law.

.. plot::

      import matplotlib.pyplot as plt
      import numpy as np
      import astropy.units as u

      from dust_attenuation.averages import C00
      from dust_attenuation.shapes import Calzmod

      fig, ax = plt.subplots()

      # generate the curves and plot them
      x = np.arange(1/2.2, 1/0.12,0.1)/u.micron

      # Original Calzetti law
      C00_model = C00(Av=1)
      ax.plot(x, C00_model(1/x), label='C00', color='black', lw=2.5, ls='--')

      slopes = [-1, -0.5, 0, 0.5, 1]
      for slope in slopes:
          att_model = Calzmod(Av=1, ampl=3.5, slope=slope)
          ax.plot(x, att_model(1/x), label=r'$\delta$ = %.2f' % (slope))

      ax.set_xlabel('$x$ [$\mu m^{-1}$]')
      ax.set_ylabel('A(x) [mag]')

      ax.legend(loc='best')
      plt.title("Calzmod with varying slopes")
      plt.show()

Example `Calzmod` models showing variation in UV bump amplitude.
The central wavelength of the UV bump and its width are kept fixed 
to 0.2175 and 0.035 microns respectively.

.. plot::

      import matplotlib.pyplot as plt
      import numpy as np
      import astropy.units as u

      from dust_attenuation.averages import C00
      from dust_attenuation.shapes import Calzmod

      fig, ax = plt.subplots()

      # generate the curves and plot them
      x = np.arange(1/2.2, 1/0.12,0.1)/u.micron

      # Original Calzetti law
      C00_model = C00(Av=1)
      ax.plot(x, C00_model(1/x), label='C00', color='black', lw=2.5, ls='--')

      amplitudes = [0, 1, 3.5, 7, 10]
      for ampl in amplitudes:
          att_model = Calzmod(Av=1, ampl=ampl, slope=0)
          ax.plot(x, att_model(1/x), label = 'ampl = %.2f' % (ampl))

      ax.set_xlabel('$x$ [$\mu m^{-1}$]')
      ax.set_ylabel('A(x) [mag]')

      ax.legend(loc='best')
      plt.title("Calzmod with varying UV bump amplitude")
      plt.show()


Charlot & Fall.
---------------

Others
------
