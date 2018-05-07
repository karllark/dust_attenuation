#####################
Fit Attenuation Curves
#####################

The ``dust_attenuation`` package is built on the `astropy.modeling
<http://docs.astropy.org/en/stable/modeling/>`_ package.  Fitting is
done in the standard way for this package where the model is initialized
with a starting point (either the default or user input), the fitter
is chosen, and the fit performed.

Example: C00 Fit
=================

In this example, we create some artificial curve with the C00 model
and explain how to fit an attenuation curve with it.

.. plot::
   :include-source:

   import matplotlib.pyplot as plt
   import numpy as np

   from astropy.modeling.fitting import LevMarLSQFitter
   import astropy.units as u

   from dust_attenuation.C00 import C00
   
   # Create artificial attenuation curve with C00 and Av = 1.3 mag
   # Better sampling using wavenumbers
   x = np.arange(0.5,8.0,0.1)/u.micron
   # Convert to microns
   x=1/x

   att_model = C00(Av=1.3)
   y = att_model(x)
   # add some noise 
   noise = np.random.normal(0, 0.2, y.shape)
   y+=noise

   # initialize the model
   c00_init = C00()

   # pick the fitter
   fit = LevMarLSQFitter()

   # fit the data to the FM90 model using the fitter
   #   use the initialized model as the starting point
   c00_fit = fit(c00_init, x.value, y)

   print ('Fit results:\n', c00_fit)
   # plot the observed data, initial guess, and final fit
   fig, ax = plt.subplots()

   ax.plot(1/x, y, 'ko', label='Observed Curve')
   ax.plot(1/x.value, c00_init(x.value), label='Initial guess')
   ax.plot(1/x.value, c00_fit(x.value), label='Fitted model')

   ax.set_xlabel('$x$ [$\mu m^{-1}$]')
   ax.set_ylabel('$Ax $')

   ax.set_title('Example C00 Fit ')

   ax.legend(loc='best')
   plt.tight_layout()
   plt.show()



Example: Use WG00 to fit C00
=================

In this example, we are using the WG00 attenuation curves to 
fit the original Calzetti attenuation curve with Av = 1 mag.
The 2 configurations best fittind the C00 curves are for a SMC 
in either a SHELL geometry with clumpy dust distribution or a
DUSTY geometry with homogeneous dust distribution.


.. plot::
   :include-source:

   import numpy as np
   import matplotlib.pyplot as plt
   from astropy.modeling.fitting import LevMarLSQFitter
   import astropy.units as u

   from dust_attenuation.C00 import C00
   from dust_attenuation.WG00 import WG00

   # Generate the C00 curve with Av = 1mag and add some noise
   x = np.arange(1/2,1/0.15,0.1)/u.micron
   x=1/x
   att_model = C00(Av=1)
   y = att_model(x)
   noise = np.random.normal(0, 0.05, y.shape)
   y+=noise

   # Convert A_lambda to tau_lambda
   y /= 1.086

   # Wavelength of V band
   x_Vband = 0.55

   geometries = ['shell', 'cloudy', 'dusty']
   dust_types = ['MW', 'SMC']
   dust_distribs = ['homogeneous', 'clumpy']

   # initialize the model
   WG00_init = WG00(tau_V=2)

   # pick the fitter
   fit = LevMarLSQFitter()

   # plot the observed data, initial guess, and final fit
   plt.figure(figsize=(15,9))

   plt.plot(1/x, y, 'ko', label='C00')

   # Loop over the different configurations
   for geo in geometries:
      for dust in dust_types:
          for distrib in dust_distribs:
            
              label = geo + '_' + dust + '_' + distrib[0]
            
              if geo == 'cloudy': color = 'red'
              elif geo == 'dusty': color = 'blue'
              elif geo == 'shell': color = 'green'
            
              if dust == 'MW': marker = 'o'
              elif dust == 'SMC': marker = '^'
                
              if distrib == 'homogeneous': ls = '--'
              if distrib == 'clumpy':  ls = '-'
            
                
              WG00_init.get_model(geometry = geo,
                                  dust_type = dust,
                                  dust_distribution = distrib)


              # fit the data to the FM90 model using the fitter
              #   use the initialized model as the starting point
              WG00_fit = fit(WG00_init, x.value, y)

              plt.plot(1/x.value, WG00_fit(x.value) / WG00_fit(x_Vband),
                       label = label, ls = ls, lw = 2, color = color,
                       marker = marker, markevery = 10, markersize = 8 )

            
   plt.xlabel('$x$ [$\mu m^{-1}$]',size=16)
   plt.ylabel(r'$\tau / \tau_V $',size=16)

   plt.title('Example: fit C00 with WG00', size =20)
   plt.tick_params(labelsize=15)
   plt.legend(loc='upper left',fontsize=18)
   plt.tight_layout()
   plt.show()
