# -*- coding: utf-8 -*-

import numpy as np
import astropy.units as u
import pkg_resources

from astropy.io import ascii
from astropy.modeling.tabular import tabular_model

from .baseclasses import BaseAtttauVModel
from .helpers import _test_valid_x_range


__all__ = ["WG00"]

x_range_WG00 = [0.1, 3.0001]


class WG00(BaseAtttauVModel):
    r"""
    Attenuation curve of Witt & Gordon (2000)

    Parameters
    ----------
    tau_v: float
        optical depth in V band

    Raises
    ------
    InputParameterError
       Input Av values outside of defined range

    Notes
    -----
    From Witt & Gordon (2000, ApJ, Volume 528, pp. 799-816)

    Example:

    .. plot::
        :include-source:

        import numpy as np
        import matplotlib.pyplot as plt
        import astropy.units as u

        from dust_attenuation.radiative_transfer import WG00

        fig, ax = plt.subplots(1,2, figsize=(10,6))

        # generate the curves and plot them
        # Use 1/microns for a better sampling
        x = np.arange(0.35,10.0,0.1)/u.micron

        x_Vband = 0.55 # microns

        tau_Vs = [0.25,0.4,1.1,17.0,46.0]
        for tau_V in tau_Vs[::-1]:
           att_model = WG00(tau_V = tau_V, geometry = 'cloudy',
                            dust_type = 'mw', dust_distribution = 'clumpy')
           ax[0].plot(x,att_model(1/x),label=r'$\tau_V$ = %.2f mag' % (tau_V))
           ax[1].plot(x,att_model(1/x)/att_model(x_Vband),
                      label=r'$\tau_V$ = %.2f mag' % (tau_V))

        ax[0].set_xlabel(r'$x$ [$\mu m^{-1}$]')
        ax[0].set_ylabel(r'$Att(x)$ [mag]')
        ax[1].set_xlabel(r'$x$ [$\mu m^{-1}$]')
        ax[1].set_ylabel(r'$Att(x)/Att_V$')


        ax[0].legend(loc='best')
        ax[1].legend(loc='best')
        fig.suptitle(r'CLOUDY / MW / clumpy model',size=15)
        plt.tight_layout()
        fig.subplots_adjust(top=0.88)

        plt.show()

    """

    tau_V_range = [0.25, 50.0]
    x_range = x_range_WG00

    def __init__(
        self, tau_V, geometry="dusty", dust_type="mw", dust_distribution="clumpy"
    ):
        """
        Load the attenuation curves for a given geometry, dust type and
        dust distribution.

        Parameters
        ----------
        tau_V: float
           optical depth in V band

        geometry: string
           'shell', 'cloudy' or 'dusty'

        dust_type: string
           'mw' or 'smc'

        dust_distribution: string
           'homogeneous' or 'clumpy'

        Returns
        -------
        Attx: np array (float)
            Att(x) attenuation curve [mag]

        """
        # Ensure strings are lower cases
        self.geometry = geometry.lower()
        self.dust_type = dust_type.lower()
        self.dust_distribution = dust_distribution.lower()

        data_path = pkg_resources.resource_filename("dust_attenuation", "data/WG00/")

        data = ascii.read(data_path + self.geometry + ".txt", header_start=0)

        if self.dust_type == "mw":
            start = 0
        elif self.dust_type == "smc":
            start = 25

        # Column names
        tau_colname = "tau"
        tau_att_colname = "tau_att"
        fsca_colname = "f(sca)"
        fdir_colname = "f(dir)"
        fesc_colname = "f(esc)"

        if self.dust_distribution == "clumpy":
            tau_att_colname += "_c"
            fsca_colname += "_c"
            fdir_colname += "_c"
            fesc_colname += "_c"

        elif self.dust_distribution == "homogeneous":
            tau_att_colname += "_h"
            fsca_colname += "_h"
            fdir_colname += "_h"
            fesc_colname += "_h"

        tau_att_list = []
        tau_list = []
        fsca_list = []
        fdir_list = []
        fesc_list = []

        len_data = len(data["lambda"])
        # number of lines between 2 models
        steps = 25

        counter = start
        while counter < len_data:
            tau_att_list.append(
                np.array(data[tau_att_colname][counter : counter + steps])
            )
            tau_list.append(np.array(data[tau_colname][counter : counter + steps]))
            fsca_list.append(np.array(data[fsca_colname][counter : counter + steps]))
            fdir_list.append(np.array(data[fdir_colname][counter : counter + steps]))
            fesc_list.append(np.array(data[fesc_colname][counter : counter + steps]))

            counter += int(2 * steps)

        # Convert to np.array and take transpose to have (wvl, tau_V)
        tau_att_table = np.array(tau_att_list).T
        tau_table = np.array(tau_list).T
        fsca_table = np.array(fsca_list).T
        fdir_table = np.array(fdir_list).T
        fesc_table = np.array(fesc_list).T

        # wavelength grid. It is the same for all the models
        wvl = np.array(data["lambda"][0:25])
        self.wvl_grid = wvl

        # Grid for the optical depth
        tau_V_grid = np.array(
            [
                0.25,
                0.5,
                0.75,
                1.0,
                1.5,
                2.0,
                2.5,
                3.0,
                3.5,
                4.0,
                4.5,
                5.0,
                5.5,
                6.0,
                7.0,
                8.0,
                9.0,
                10.0,
                15.0,
                20.0,
                25.0,
                30.0,
                35.0,
                40.0,
                45.0,
                50.0,
            ]
        )

        # Create a 2D tabular model for tau_att and all flux fraction
        tab = tabular_model(2, name="2D_table")

        # Values corresponding to the x and y grid points
        gridpoints = (wvl, tau_V_grid)

        self.model = tab(
            gridpoints,
            lookup_table=tau_att_table,
            name="tau_att_WG00",
            bounds_error=False,
            fill_value=None,
            method="linear",
        )

        self.tau = tab(
            gridpoints,
            lookup_table=tau_table,
            name="tau_WG00",
            bounds_error=False,
            fill_value=None,
            method="linear",
        )

        self.fsca = tab(
            gridpoints,
            lookup_table=fsca_table,
            name="fsca_WG00",
            bounds_error=False,
            fill_value=None,
            method="linear",
        )

        self.fdir = tab(
            gridpoints,
            lookup_table=fdir_table,
            name="fdir_WG00",
            bounds_error=False,
            fill_value=None,
            method="linear",
        )

        self.fesc = tab(
            gridpoints,
            lookup_table=fesc_table,
            name="fesc_WG00",
            bounds_error=False,
            fill_value=None,
            method="linear",
        )

        # In Python 2: super(WG00, self)
        # In Python 3: super() but super(WG00, self) still works
        super(WG00, self).__init__(tau_V=tau_V)

    def evaluate(self, x, tau_V):
        """
        WG00 function

        Parameters
        ----------
        x: float
           expects either x in units of wavelengths or frequency
           or assumes wavelengths in [micron]

           internally microns are used

        tau_V: float
           optical depth in V band

        Returns
        -------
        Attx: np array (float)
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
        _test_valid_x_range(x, self.x_range, "WG00")

        # setup the ax vectors
        n_x = len(x)

        xinterp = 1e4 * x
        yinterp = tau_V * np.ones(n_x)

        taux = self.model(xinterp, yinterp)

        # Convert optical depth to attenuation
        Attx = 1.086 * taux

        return Attx

    def get_extinction(self, x, tau_V):
        """
        Return the extinction at a given wavelength and
        V-band optical depth.

        Parameters
        ----------
        x: float
           expects either x in units of wavelengths or frequency
           or assumes wavelengths in [micron]

           internally microns are used

        tau_V: float
           optical depth in V band

        Returns
        -------
        ext: np array (float)
            ext(x) extinction curve [mag]

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
        _test_valid_x_range(x, self.x_range, "WG00")

        # setup the ax vectors
        x = np.atleast_1d(x)
        n_x = len(x)

        xinterp = 1e4 * x
        yinterp = tau_V * np.ones(n_x)

        return self.tau(xinterp, yinterp) * 1.086

    def get_fsca(self, x, tau_V):
        """
        Return the scattered flux fraction  at a given wavelength and
        V-band optical depth.

        Parameters
        ----------
        x: float
           expects either x in units of wavelengths or frequency
           or assumes wavelengths in [micron]

           internally microns are used

        tau_V: float
           optical depth in V band

        Returns
        -------
        fsca: np array (float)
            fsca(x) scattered flux fraction

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
        _test_valid_x_range(x, self.x_range, "WG00")

        # setup the ax vectors
        x = np.atleast_1d(x)
        n_x = len(x)

        xinterp = 1e4 * x
        yinterp = tau_V * np.ones(n_x)

        return self.fsca(xinterp, yinterp)

    def get_fdir(self, x, tau_V):
        """
        Return the direct attenuated stellar flux fraction  at a given
        wavelength and V-band optical depth.

        Parameters
        ----------
        x: float
           expects either x in units of wavelengths or frequency
           or assumes wavelengths in [micron]

           internally microns are used

        tau_V: float
           optical depth in V band

        Returns
        -------
        fsca: np array (float)
            fsca(x) scattered flux fraction

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
        _test_valid_x_range(x, self.x_range, "WG00")

        # setup the ax vectors
        x = np.atleast_1d(x)
        n_x = len(x)

        xinterp = 1e4 * x
        yinterp = tau_V * np.ones(n_x)

        return self.fdir(xinterp, yinterp)

    def get_fesc(self, x, tau_V):
        """
        Return the total escaping flux fraction  at a given wavelength and
        V-band optical depth.

        Parameters
        ----------
        x: float
           expects either x in units of wavelengths or frequency
           or assumes wavelengths in [micron]

           internally microns are used

        tau_V: float
           optical depth in V band

        Returns
        -------
        fsca: np array (float)
            fsca(x) scattered flux fraction

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
        _test_valid_x_range(x, self.x_range, "WG00")

        # setup the ax vectors
        x = np.atleast_1d(x)
        n_x = len(x)

        xinterp = 1e4 * x
        yinterp = tau_V * np.ones(n_x)

        return self.fesc(xinterp, yinterp)

    def get_albedo(self, x):
        """
        Return the albedo in function of wavelength for the corresponding
        dust type (SMC or MW). The albedo gives the probability a photon
        is scattered from a dust grain.

        Parameters
        ----------
        x: float
           expects either x in units of wavelengths or frequency
           or assumes wavelengths in [micron]

           internally microns are used

        Returns
        -------
        albedo: np array (float)
            alb(x) albedo

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
        _test_valid_x_range(x, self.x_range, "WG00")

        # setup the ax vectors
        x = np.atleast_1d(x)

        alb_MW = np.array(
            [
                0.320,
                0.409,
                0.481,
                0.526,
                0.542,
                0.536,
                0.503,
                0.432,
                0.371,
                0.389,
                0.437,
                0.470,
                0.486,
                0.499,
                0.506,
                0.498,
                0.502,
                0.491,
                0.481,
                0.500,
                0.473,
                0.457,
                0.448,
                0.424,
                0.400,
            ]
        )

        alb_SMC = np.array(
            [
                0.400,
                0.449,
                0.473,
                0.494,
                0.508,
                0.524,
                0.529,
                0.528,
                0.523,
                0.520,
                0.516,
                0.511,
                0.505,
                0.513,
                0.515,
                0.498,
                0.494,
                0.489,
                0.484,
                0.493,
                0.475,
                0.465,
                0.439,
                0.417,
                0.400,
            ]
        )

        if self.dust_type == "smc":
            albedo = alb_SMC
        elif self.dust_type == "mw":
            albedo = alb_MW

        tab = tabular_model(1, name="Tabular1D")
        alb_fit = tab(
            self.wvl_grid,
            lookup_table=albedo,
            name="albedo",
            bounds_error=False,
            fill_value=None,
            method="linear",
        )

        xinterp = 1e4 * x

        return alb_fit(xinterp)

    def get_scattering_phase_function(self, x):
        """
        Return the scattering phase function in function of wavelength for the
        corresponding dust type (SMC or MW). The scattering phase
        function gives the angle at which the photon scatters.

        Parameters
        ----------
        x: float
           expects either x in units of wavelengths or frequency
           or assumes wavelengths in [micron]

           internally microns are used

        Returns
        -------
        g: np array (float)
            g(x) scattering phase function

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
        _test_valid_x_range(x, self.x_range, "WG00")

        # setup the ax vectors
        x = np.atleast_1d(x)

        g_MW = np.array(
            [
                0.800,
                0.783,
                0.767,
                0.756,
                0.745,
                0.736,
                0.727,
                0.720,
                0.712,
                0.707,
                0.702,
                0.697,
                0.691,
                0.685,
                0.678,
                0.646,
                0.624,
                0.597,
                0.563,
                0.545,
                0.533,
                0.511,
                0.480,
                0.445,
                0.420,
            ]
        )

        g_SMC = np.array(
            [
                0.800,
                0.783,
                0.767,
                0.756,
                0.745,
                0.736,
                0.727,
                0.720,
                0.712,
                0.707,
                0.702,
                0.697,
                0.691,
                0.685,
                0.678,
                0.646,
                0.624,
                0.597,
                0.563,
                0.545,
                0.533,
                0.511,
                0.480,
                0.445,
                0.420,
            ]
        )

        if self.dust_type == "smc":
            g = g_SMC
        elif self.dust_type == "mw":
            g = g_MW

        tab = tabular_model(1, name="Tabular1D")
        g_fit = tab(
            self.wvl_grid,
            lookup_table=g,
            name="albedo",
            bounds_error=False,
            fill_value=None,
            method="linear",
        )

        xinterp = 1e4 * x

        return g_fit(xinterp)
