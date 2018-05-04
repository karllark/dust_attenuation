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

   from dust_extinction.dust_extinction import C00
   
   # Create artificial attenuation curve with C00 and Av = 1.3 mag
   x = np.arange(0.5,8.0,0.1)/u.micron
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

   # plot the observed data, initial guess, and final fit
   fig, ax = plt.subplots()

   ax.plot(x, y, 'ko', label='Observed Curve')
   ax.plot(x.value, c00_init(x.value), label='Initial guess')
   ax.plot(x.value, c00_fit(x.value), label='Fitted model')

   ax.set_xlabel('$x$ [$\mu m^{-1}$]')
   ax.set_ylabel('$Ax $')

   ax.set_title('Example C00 Fit ')

   ax.legend(loc='best')
   plt.tight_layout()
   plt.show()
