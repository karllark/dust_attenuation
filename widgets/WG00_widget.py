#! /usr/bin/python

# This widget aims at illustrating the effect of the geometry
# on the attenaution curve for the radiative transfer model of
# Witt & Gordon 2000.
# To execute it, just type "python WG00_widget.py" in the terminal.
# A new window will pop up.

import matplotlib.pyplot as plt
import matplotlib.widgets as wgt
from matplotlib.patches import Circle, PathPatch, Wedge
import numpy as np
from dust_attenuation.radiative_transfer import WG00
import astropy.units as u
import matplotlib.gridspec as gridspec


class WG00_widget:

    def __init__(self):

        self.param = {'dust_type': 'MW',
                      'tau_V':1,
                      'geometry':'SHELL',
                      'dust_distrib':'homogeneous'}

        # generate the curves and plot them
        # Use 1/microns for a better sampling
        self.x = np.arange(0.35,10.0,0.1)/u.micron

        self.x_Vband = 0.55

        self.fig, self.ax = plt.subplots(figsize=(15,10))
        # set up subplot grid
        gs1 = gridspec.GridSpec(4,2)

        self.axatt = plt.subplot2grid((4,2), (0,0), colspan=2, rowspan=3)
        self.axFF = plt.subplot2grid((4,2), (3,0), colspan=1, rowspan=1)
        self.axalb = plt.subplot2grid((4,2), (3,1), colspan=1, rowspan=1)


        self.fig.canvas.set_window_title('widget for WG00 RT model')
        plt.subplots_adjust(bottom=0.15)
        plt.subplots_adjust(left=0.35)
        plt.subplots_adjust(hspace=0.8)


        axcolor = 'lightgoldenrodyellow'

        self.rax_sketch = plt.axes([.03,.68,.25,.32])
        self.rax_sketch.axis('off')


        rax = plt.axes([0.05, 0.52,0.15,0.15], facecolor=axcolor)
        dust_type = wgt.RadioButtons(rax,
                                    ('MW', 'SMC'))
        dust_type.on_clicked(self.update_dust_type)

        rax = plt.axes([0.05, 0.335, 0.15, 0.15], facecolor=axcolor)
        geometry = wgt.RadioButtons(rax,
                                    ('SHELL', 'CLOUDY', 'DUSTY'))
        geometry.on_clicked(self.update_geometry)


        rax = plt.axes([0.05, 0.15, 0.15, 0.15], facecolor=axcolor)
        distrib = wgt.RadioButtons(rax,
                                    ('Homogeneous', 'Clumpy'))
        distrib.on_clicked(self.update_dust_distrib)

        rax = plt.axes([0.05, 0.03, 0.12, 0.08], facecolor=axcolor)
        norm = wgt.CheckButtons(rax, (r'Normalised to A$_V$',), (True,))
        norm.on_clicked(self.update_norm)
        self.norm = True

        rax = plt.axes([.3,.05,.6,.0275])
        tau_V = wgt.Slider(rax, r'$\tau_V$', 0.5, 50,
                           valinit=self.param['tau_V'])
        tau_V.on_changed(self.update_tau_V)

        # Initialise WG00 model
        self.att_model = WG00(tau_V = self.param['tau_V'],
                              geometry = self.param['geometry'],
                              dust_type = self.param['dust_type'],
                              dust_distribution = self.param['dust_distrib'])

        self.update_sketch()
        self.update_att_curve()
        plt.show()

    def update_tau_V(self, val):
        self.param['tau_V'] = val
        self.att_model = WG00(tau_V = self.param['tau_V'],
                              geometry = self.param['geometry'],
                              dust_type = self.param['dust_type'],
                              dust_distribution = self.param['dust_distrib'])

        self.update_att_curve()


    def update_geometry(self, val):
        self.param['geometry'] = val
        self.att_model = WG00(tau_V = self.param['tau_V'],
                              geometry = self.param['geometry'],
                              dust_type = self.param['dust_type'],
                              dust_distribution = self.param['dust_distrib'])

        self.update_att_curve()

    def update_dust_type(self, val):
        self.param['dust_type'] = val
        self.att_model = WG00(tau_V = self.param['tau_V'],
                              geometry = self.param['geometry'],
                              dust_type = self.param['dust_type'],
                              dust_distribution = self.param['dust_distrib'])

        self.update_att_curve()

    def update_dust_distrib(self, val):
        self.param['dust_distrib'] = val
        self.att_model = WG00(tau_V = self.param['tau_V'],
                              geometry = self.param['geometry'],
                              dust_type = self.param['dust_type'],
                              dust_distribution = self.param['dust_distrib'])

        self.update_att_curve()

    def update_norm(self, val):
        self.norm = not self.norm
        self.update_att_curve()

    def update_att_curve(self):
        self.att = self.att_model(1/self.x)
        self.att_V = self.att_model(self.x_Vband)
        self.ext = self.att_model.get_extinction(1/self.x, self.param['tau_V'])
        self.ext_V = self.param['tau_V'] * 1.086
        self.fsca = self.att_model.get_fsca(1/self.x, self.param['tau_V'])
        self.fdir = self.att_model.get_fdir(1/self.x, self.param['tau_V'])
        self.fesc = self.att_model.get_fesc(1/self.x, self.param['tau_V'])
        self.alb = self.att_model.get_albedo(1/self.x)
        self.g = self.att_model.get_scattering_phase_function(1/self.x)
        self.update_plot()


    def update_sketch(self):
        # clear previous plot
        self.rax_sketch.clear()

        # Fixing random state for reproducibility
        np.random.seed(1234567890)

        if self.param['geometry'] == 'SHELL':
            Rs = [0, 0.3]
            Rd = [0.3, 1]
        elif self.param['geometry'] == 'DUSTY':
            Rs = [0, 1]
            Rd = [0, 1]
        # in WG00 Rd=0.69, Rs=1, here we normalised Rd=1
        elif self.param['geometry'] == 'CLOUDY':
            Rs = [0, 1.45]
            Rd = [0, 1]

        rad_max=max(Rs[1],Rd[1])

        if self.param['dust_distrib'] == 'Clumpy':
            clumpy = True
        else:
            clumpy = False

        # set transparency
        tauV_lim = [0.5, 75]
        alpha = self.param['tau_V'] / (tauV_lim[1]-tauV_lim[0])
        if alpha>1: alpha=1

        #plot stars
        x = np.arange(-Rs[1],Rs[1],0.2*rad_max)
        y=np.arange(-Rs[1],Rs[1],0.2*rad_max)
        X,Y=np.meshgrid(x,y)
        mask=X**2+Y**2 < Rs[1]**2 * 0.99
        self.rax_sketch.plot(X[mask], Y[mask], "*", color='orange',
                             ms=10, zorder=1)

        if clumpy:
            # Plot clumpiness
            marker = 'o'#(15,1,60)
            num = 150
            size=np.random.rand(num)*500
            t = np.random.uniform(0.0, 2.0*np.pi, num)
            r = np.sqrt(np.random.uniform(Rd[0]**2, Rd[1]**2, num))
            x = r * np.cos(t)
            y = r * np.sin(t)
            self.rax_sketch.scatter(x, y, marker=marker, color='black',
                                    s=size, alpha=alpha, zorder=2)

        #Plot star ring
        star_ring=Wedge(0, Rs[1], 0, 360, width=Rs[1], alpha=alpha,
                        color='lightyellow', zorder=0)
        self.rax_sketch.add_patch(star_ring)

        # Plot dust
        if not clumpy:
            dust_ring=Wedge(0, Rd[1], 0, 360, width=1-Rd[0], alpha=alpha,
                            color='black', zorder=2)
        else:

            dust_ring=Wedge(0, Rd[1], 0, 360, width=Rd[1], alpha=alpha/2,
                            color='grey', zorder=2)

        self.rax_sketch.add_patch(dust_ring)

        self.rax_sketch.set_xlim(-1.1*rad_max, 1.1*rad_max)
        self.rax_sketch.set_ylim(-1.1*rad_max, 1.1*rad_max)
        self.rax_sketch.axis('off')


    def update_plot(self):
        if len(self.axatt.get_lines()) > 0:
            if self.param['dust_type'] == 'MW':
                color='C1'
            elif self.param['dust_type'] == 'SMC':
                color = 'C0'
            if self.norm:
                self.plot_att.set_ydata(self.att/self.att_V)
                self.plot_ext.set_ydata(self.ext/self.ext_V)
                self.axatt.set_ylabel(r'A$_{\lambda}$ / A$_V$', size=16 )
                self.axatt.set_ylim(0,10)
            else:
                self.plot_att.set_ydata(self.att)
                self.plot_ext.set_ydata(self.ext)
                self.axatt.set_ylabel(r'A$_{\lambda}$', size=16 )
                self.axatt.set_ylim(0, max(np.max(self.ext),10))

            self.plot_att.set_color(color)
            self.plot_ext.set_color(color)
            self.axatt.legend(prop={'size': 16})

            self.plot_fesc.set_ydata(self.fesc)
            self.plot_fdir.set_ydata(self.fdir)
            self.plot_fsca.set_ydata(self.fsca)
            #ffmin = np.nanmin([self.fsca, self.fdir,self.fesc])
            #self.axFF.set_ylim(ffmin, 1)

            self.plot_alb.set_ydata(self.alb)
            self.plot_g.set_ydata(self.g)

            self.update_sketch()
            self.fig.canvas.draw()
        else:
            if self.param['dust_type'] == 'MW':
                color='C1'
            elif self.param['dust_type'] == 'SMC':
                color = 'C0'

            self.plot_ext, = self.axatt.plot(self.x, self.ext, ls='--', lw=3,
                                          color=color, label='Extinction')
            self.plot_att, = self.axatt.plot(self.x, self.att/self.att_V,
                                          ls='-', lw=3, color=color,
                                          label='Attenuation')

            self.axatt.set_xlim(0.,10.5)
            self.axatt.set_ylim(0,10)
            self.axatt.legend(prop={'size': 16})
            self.axatt.set_xlabel(r'1/$\lambda$ [$\mu m^{-1}$]', size=16)
            self.axatt.set_ylabel(r'A$_{\lambda}$ / A$_V$', size=16 )
            self.axatt.tick_params(labelsize=14)
            self.axatt.set_title('Attenuation curves from radiative transfer'
                              +' model of\n  Witt & Gordon (2000)', size=20)



            self.plot_fesc, = self.axFF.plot(self.x, self.fesc,
                                          ls='-', lw=1, color='black',
                                          label='esc')

            self.plot_fsca, = self.axFF.plot(self.x, self.fsca,
                                          ls=':', lw=1, color='black',
                                          label='sca')


            self.plot_fdir, = self.axFF.plot(self.x, self.fdir,
                                          ls='--', lw=1, color='black',
                                          label='dir')


            self.axFF.set_xlim(0.,10.5)
            self.axFF.set_yscale('log')
            self.axFF.set_ylim(1e-2,1)
            self.axFF.legend(prop={'size': 10})
            self.axFF.set_xlabel(r'1/$\lambda$ [$\mu m^{-1}$]', size=14)
            #self.axFF.set_ylabel('Flux fraction', size=14 )
            self.axFF.set_title('Flux fraction', size=14 )
            self.axFF.tick_params(labelsize=12)
            #self.axatt.set_title('Attenuation curves from radiative transfer'
            #                  +' model of\n  Witt & Gordon (2000)', size=20)


            self.plot_g, = self.axalb.plot(self.x, self.g,
                                          ls='--', lw=1, color='black',
                                          label='g')

            self.plot_alb, = self.axalb.plot(self.x, self.alb,
                                          ls='-', lw=1, color='black',
                                          label='albedo')

            self.axalb.set_xlim(0.,10.5)
            self.axalb.set_ylim(0,1)
            self.axalb.legend(loc='lower right', prop={'size': 10})
            self.axalb.set_xlabel(r'1/$\lambda$ [$\mu m^{-1}$]', size=14)
            #self.axalb.set_ylabel('Albedo / phase function', size=14 )
            self.axalb.set_title('Albedo / phase function', size=14 )
            self.axalb.tick_params(labelsize=12)


            self.fig.canvas.draw()


commander = WG00_widget()
