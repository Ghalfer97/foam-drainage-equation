import numpy as np
import math as mt
import sys, libconf, io
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cmx
from datetime import datetime


class PlotClass:
    '''
    Ensures that the correct plots are chosen
    based on the conditions specified in the config file
    '''
    def __init__(self, Col):
        self.set_what_to_plot(Col)
        self.Input_Parameters(Col)

        self.Import_Arrays(Col)

        self.Bubbling_Crit_Vel(Col)

        if Col.input_lf == 0.:
            self.Free_Drainage_SS_Height_min(Col)
            self.Free_Drainage_SS_Height_max(Col)
            self.Free_Drainage_SS_Prof(Col)
        self.Critical_Alpha_Line()
        if Col.bubble_velocity > 0.:
            if Col.dimless_bubble_velocity < self.crit_lf_min:
                self.Bubbling_SS_Height(Col)
            self.Bubbling_SS_Prof(Col)
        if self.dimensionless_plots:
            self.Plot_Dimensionless(Col)
        if self.dimensional_plots:
            self.Plot_Dimensional(Col)

    def set_what_to_plot(self, Col):
        self.dimensionless_plots = Col.dimensionless_plots
        self.dimensional_plots = Col.dimensional_plots

    def Input_Parameters(self, Col):
        '''
        Input parameters and variables
        '''

        phifactor_inverse = 1.0 / Col.phifactor  # Prevent multiple floating point divisions
        self.input_lf = Col.input_lf * phifactor_inverse
        self.initial_lf = Col.initial_lf * phifactor_inverse
        self.base_lf = Col.base_lf * phifactor_inverse
        self.crit_lf_min = Col.crit_lf_min * phifactor_inverse
        self.crit_lf_max = Col.crit_lf_max * phifactor_inverse
        self.Total_Liquid_In_Foam = (Col.col_height /
                                     Col.length_scale) * (Col.initial_lf)
        self.crit_lf = Col.crit_lf
        self.Equilibrium_Avg_lf = (((2. / np.sqrt(3.)) *
                                    (Col.gamma /
                                     (Col.water_density * Col.gravity)) *
                                    (np.sqrt(Col.base_lf))) /
                                   (Col.bubble_diameter * Col.col_height)) * (
                                       1. - 1. /
                                       (1. + (np.sqrt(3.) / 2.) *
                                        (Col.bubble_diameter /
                                         (Col.gamma /
                                          (Col.water_density * Col.gravity))) *
                                        Col.col_height * np.sqrt(Col.base_lf)))

    def Import_Arrays(self, Col):
        '''
        Imports the various arrays
        Gives the relevant heights
        '''

        self.num_time_samples = Col.plot_height_array.size
        self.num_time_samples2 = Col.plot_height_increase_array.size
        self.num_space_samples = Col.plot_alpha_array[:, 0].size
        self.col_height = self.num_space_samples * Col.dz
        self.maxtime = int(Col.runtime / Col.dt)
        self.plot_times = Col.dt * np.linspace(0., self.maxtime,
                                               self.num_time_samples)
        self.plot_times2 = Col.dt * np.linspace(0., self.maxtime,
                                                self.num_time_samples2)
        self.plot_cells = np.arange(0., self.col_height, Col.dz)
        self.highest_point = max(Col.plot_height_array)

    def Bubbling_Crit_Vel(self, Col):
        self.dimless_crit_velocity = self.crit_lf_min
        self.dimensional_crit_velocity = (Col.length_scale /
                                          Col.time_scale) * self.crit_lf_min

    def Free_Drainage_SS_Height_min(self, Col):
        self.free_drainage_eq_height_min = 2 * Col.gamma / (
            np.sqrt(3) * Col.bubble_diameter * Col.water_density *
            Col.gravity) * (1 / np.sqrt(Col.phifactor * self.crit_lf_min) -
                            1 / np.sqrt(Col.phifactor * self.base_lf))
        self.dimless_free_drainage_eq_height_min = self.free_drainage_eq_height_min / Col.length_scale

    def Free_Drainage_SS_Height_max(self, Col):
        self.free_drainage_eq_height_max = 2 * Col.gamma / (
            np.sqrt(3) * Col.bubble_diameter * Col.water_density *
            Col.gravity) * (1 / np.sqrt(Col.phifactor * self.crit_lf_max) -
                            1 / np.sqrt(Col.phifactor * self.base_lf))
        self.dimless_free_drainage_eq_height_max = self.free_drainage_eq_height_max / Col.length_scale

    def Free_Drainage_SS_Prof(self, Col):  # free drainage equilibrium profile
        self.free_eq_prof = np.power(
            1 / (mt.sqrt(Col.base_lf)) +
            (mt.sqrt(3) * Col.bubble_diameter * Col.length_scale *
             np.linspace(0.0, self.col_height, self.num_space_samples)) /
            (2 * (Col.gamma / (Col.water_density * Col.gravity))),
            -2)  #Eqns 10.7 - 10.9 H&W
        self.dimless_free_eq_prof = self.free_eq_prof / Col.phifactor

    def Critical_Alpha_Line(self):
        self.crit_lf_line_lf_min = np.array(
            [self.crit_lf_min, self.crit_lf_min])
        self.crit_lf_line_lf_max = np.array(
            [self.crit_lf_max, self.crit_lf_max])
        self.crit_lf_line_z = np.array([0., self.col_height])

    def Bubbling_SS_Height(self, Col):
        '''The expected dependency of the Steady State height on the inputted velocity'''

        numerator = np.sqrt(
            self.crit_lf * self.base_lf) - Col.dimless_bubble_velocity
        denominator = np.sqrt(Col.dimless_bubble_velocity) * (
            np.sqrt(self.base_lf) - np.sqrt(self.crit_lf))
        arg_steadystate = numerator / denominator
        coefficient = 1. / np.sqrt(Col.dimless_bubble_velocity)
        arcoth = .5 * mt.log((1 + arg_steadystate) / (arg_steadystate - 1))
        self.ss_height = coefficient * arcoth

    def Bubbling_SS_Prof(self, Col):
        #calculate steady state equilibrium profile
        '''
        Below is inverse of Eqn 8 H et al. evaluated at alpha_0 (base) to find constant ksi_a
        '''
        arg_epsilon = mt.sqrt(self.base_lf / Col.dimless_bubble_velocity)
        epsilon_arcoth = mt.log((arg_epsilon + 1) / (arg_epsilon - 1)) / 2
        epsilon = -epsilon_arcoth / mt.sqrt(
            Col.dimless_bubble_velocity)  #xi_a (since xi = 0 at base)
        '''
        Eqn 8 H et al., giving equilibrium profile
        '''
        self.prof_points = np.linspace(0., self.col_height,
                                       self.num_space_samples)
        arg_prof = np.sqrt(
            Col.dimless_bubble_velocity) * (self.prof_points - epsilon)
        prof_cotanh = np.cosh(arg_prof) / np.sinh(arg_prof)
        self.eq_prof = Col.dimless_bubble_velocity * prof_cotanh * prof_cotanh  #(Dimensionless) Equilibrium Profile (array)

    '''
    Defines the colour gradient scheme
    This makes it such that at profiles plotted have fainter colours if they are older profiles
    '''

    def Plot_Colors(self):
        self.gradient = 5.  #adjust gradient of colour scale and transparency of first line
        Blues = plt.get_cmap('Blues')  #retrieves 'blues' sequential colour map
        cNorm = colors.Normalize(
            vmin=0, vmax=self.gradient +
            len(self.plot_times))  #normalizes plotting range for color-coding
        self.scalarMap = cmx.ScalarMappable(
            norm=cNorm, cmap=Blues
        )  #creates scalar to colour map using normalized scale and selected colour map

    def Color_Gradient(self, is_dimensional, Col):

        for i in range(0, self.num_time_samples):
            colorVal = self.scalarMap.to_rgba(self.gradient)
            if is_dimensional:
                plt.plot(Col.phifactor * Col.plot_alpha_array[:, i],
                         Col.length_scale * self.plot_cells,
                         color=colorVal)
            else:
                plt.plot(Col.plot_alpha_array[:, i],
                         self.plot_cells,
                         color=colorVal)
            self.gradient += 1

    def Dimless_Profile_Plot(self, Col):

        if Col.bubble_velocity > 0.:
            plt.plot(
                self.eq_prof,
                self.prof_points,
                color='r',
                label='Bubbling Eq. Profile',
                linestyle='--')  #plots theoretical equilibrium profile in red
            conv_title = 'Alpha Profiles for a Rising Column \n' r'of Gas Velocity, $\nu=$' + '{:.6f}'.format(
                Col.dimless_bubble_velocity)
        if self.input_lf == 0.:
            plt.plot(self.dimless_free_eq_prof,
                     np.linspace(0.0, self.col_height, self.num_space_samples),
                     linestyle='--',
                     color='g',
                     label='Free Drainage Eq. Profile')
            if Col.bubble_velocity == 0.:
                conv_title = 'Alpha Profiles for a Column under Free Drainage'
        else:
            conv_title = 'Alpha Profiles for a Column under Forced Drainage \n' r'of Amplitude, $\alpha$ = ' + '{:.4f}'.format(
                self.input_lf)

        plt.plot(self.crit_lf_line_lf_min,
                 self.crit_lf_line_z,
                 color='b',
                 label='Rupture Liquid Fraction (Minimum)',
                 linestyle='--')  #plots critical liquid fraction in blue
        plt.plot(self.crit_lf_line_lf_max,
                 self.crit_lf_line_z,
                 color='yellow',
                 label='Rupture Liquid Fraction (Maximum)',
                 linestyle='--')  #plots critical liquid fraction in blue

        plt.xlabel(r'$\alpha$')
        plt.ylabel('Height, $z$')
        plt.grid(b=True, which='major', linestyle='--')
        plt.title(conv_title)
        plt.legend(loc='best')
        plt.figure()

    def Dimless_Height_Plot(self, Col):

        plt.title('Column Height against Time (Dimensionless)')
        plt.ylabel(r'Foam Height, $h$')
        plt.xlabel(r'Time Elapsed, $\tau$')
        plt.xlim(xmax=Col.runtime, xmin=0)
        plt.grid(b=True, which='major', linestyle='--')
        plt.plot(self.plot_times,
                 Col.plot_height_array,
                 label='Simulated Heights'
                 )  #height against time (for fixed velocity)

        if self.input_lf == 0. and Col.dimless_bubble_velocity == 0.:
            '''
            plt.plot(
                np.linspace(0, Col.runtime, 100),
                np.ones(100) * self.dimless_free_drainage_eq_height_min,
                color='r',
                linestyle='--',
                label='Theoretical Steady State Height (Min Const Crit Lf)'
            )  #theoretical steady state height
            '''
            plt.plot(
                np.linspace(0, Col.runtime, 100),
                np.ones(100) * self.dimless_free_drainage_eq_height_max,
                color='yellow',
                linestyle='--',
                label='Theoretical Steady State Height (Max Const Crit Lf)')
            plt.ylim(ymax=1.3 * max(self.dimless_free_drainage_eq_height_max,
                                    self.highest_point),
                     ymin=0)
            plt.legend(loc='best')
        elif Col.bubble_velocity != 0. and Col.dimless_bubble_velocity < self.crit_lf_min:
            plt.plot(np.linspace(0, Col.runtime, 100),
                     np.ones(100) * self.ss_height,
                     color='r',
                     linestyle='--',
                     label='Theoretical Steady State Height'
                     )  #theoretical steady state height
            plt.ylim(ymax=1.3 * max(self.ss_height, self.highest_point),
                     ymin=0)
            plt.legend(loc='best')
        else:
            plt.ylim(ymax=1.3 * self.highest_point, ymin=0)
        print("Dimensionless Bubble Velocity: ", Col.dimless_bubble_velocity)
        print("Steady State Height: ", np.mean(Col.plot_height_array[-10:]))
        print("Dimensionless Final Height: ", Col.plot_height_array[-1])
        print("Dimensionless Initial Height: ", Col.plot_height_array[0])
        print("Slope Velocity: ",
              ((Col.plot_height_array[-1] - Col.plot_height_array[0]) /
               Col.runtime))
        if self.dimensional_plots:
            plt.figure()

    def Dimless_Liquid_Height_Plot(self, Col):

        plt.title('Liquid Drained against Time (Dimensionless and Normalised)')
        plt.ylabel(r'Normalised Liquid Drained, $h$')
        plt.xlabel(r'Time Elapsed, $\tau$')
        plt.xlim(xmax=Col.runtime, xmin=0)
        plt.ylim(ymin=-0.00001)
        plt.grid(b=True, which='major', linestyle='--')
        plt.plot(self.plot_times, (Col.plot_height_increase_array) /
                 (self.Total_Liquid_In_Foam),
                 label='Normalised Liquid Height'
                 )  #height against time (for fixed velocity)
        #plt.plot(np.linspace(0, Col.runtime, 100), np.ones(100)*self.Equilibrium_Avg_lf, color = 'r', linestyle = '--', label = 'Theoretical End behaviour')
        plt.legend(loc='best')
        if self.dimensional_plots:
            plt.figure()

    def Plot_Dimensionless(self, Col):
        self.Plot_Colors()
        self.Color_Gradient(False, Col)
        self.Dimless_Profile_Plot(Col)
        self.Dimless_Height_Plot(Col)
        self.Dimless_Liquid_Height_Plot(Col)

    def Dimful_Profile_Plot(self, Col):

        if Col.bubble_velocity != 0.:
            plt.plot(
                Col.phifactor * self.eq_prof,
                Col.length_scale * self.prof_points,
                color='r',
                label='Bubbling Eq. Profile',
                linestyle='--')  #plots theoretical equilibrium profile in red
            conv_title = 'Liquid Fraction Profiles for a Rising Column \n' r'of Gas Velocity, $V=$' + '{:.6f}'.format(
                Col.bubble_velocity) + r'$ms^{-1}$'
        if self.input_lf == 0.:
            plt.plot(self.free_eq_prof,
                     np.linspace(0.0, Col.length_scale * self.col_height,
                                 self.num_space_samples),
                     linestyle='--',
                     color='g',
                     label='Free Drainage Eq. Profile')
            if Col.bubble_velocity == 0.:
                conv_title = 'Liquid Fraction Profiles for a Column under Free Drainage'
        else:
            conv_title = 'Liquid Fraction Profiles for a Column under Forced Drainage \n' r'of Amplitude, $\phi_l$ = ' + '{:.4f}'.format(
                Col.input_lf)

        plt.plot(Col.phifactor * self.crit_lf_line_lf_min,
                 Col.length_scale * self.crit_lf_line_z,
                 color='b',
                 label='Rupture Liquid Fraction (Minimum)',
                 linestyle='--')  #plots critical liquid fraction in blue
        plt.plot(Col.phifactor * self.crit_lf_line_lf_max,
                 Col.length_scale * self.crit_lf_line_z,
                 color='yellow',
                 label='Rupture Liquid Fraction (Maximum)',
                 linestyle='--')

        plt.xlabel(r'Liquid Fraction, $\phi_l$')
        plt.ylabel(r'Height, $Z$ ($m$)')
        plt.grid(b=True, which='major', linestyle='--')
        plt.title(conv_title)
        plt.legend(loc='best')
        plt.figure()

    def Dimful_Height_Plot(self, Col):

        plt.title('Column Height against Time')
        plt.ylabel('Foam Height, $H$ ($m$)')
        plt.xlabel('Time Elapsed, $T$ ($s$)')
        plt.xlim(xmax=Col.time_scale * Col.runtime, xmin=0)
        plt.grid(b=True, which='major', linestyle='--')
        plt.plot(Col.time_scale * self.plot_times,
                 Col.length_scale * Col.plot_height_array,
                 label='Simulated Heights'
                 )  #height against time (for fixed velocity)

        if self.input_lf == 0. and Col.dimless_bubble_velocity == 0.:
            '''
            plt.plot(
                np.linspace(0, Col.runtime, 100),
                np.ones(100) * self.free_drainage_eq_height_min,
                color='r',
                linestyle='--',
                label='Theoretical Steady State Height (Min Const Crit lf)'
            )  #theoretical steady state height
            '''
            plt.plot(
                np.linspace(0, Col.runtime, 100),
                np.ones(100) * self.free_drainage_eq_height_max,
                color='yellow',
                linestyle='--',
                label='Theoretical Steady State Height (Max Const Crit lf)')
            plt.ylim(ymax=1.3 * max(self.free_drainage_eq_height_min,
                                    Col.length_scale * self.highest_point),
                     ymin=0)
            plt.legend(loc='best')
        elif Col.bubble_velocity != 0. and Col.dimless_bubble_velocity < self.crit_lf_min:
            plt.plot(np.linspace(0, Col.time_scale * Col.runtime, 100),
                     np.ones(100) * Col.length_scale * self.ss_height,
                     color='r',
                     linestyle='--',
                     label='Theoretical Steady State Height'
                     )  #theoretical steady state height
            plt.ylim(ymax=1.3 * Col.length_scale *
                     max(self.ss_height, self.highest_point),
                     ymin=0)
            plt.legend(loc='best')
        else:
            plt.ylim(ymax=1.3 * Col.length_scale * self.highest_point, ymin=0)
        print("Dimensional Bubbling Velocity: ", Col.bubble_velocity)
        print("Steady State Height: ",
              Col.length_scale * np.mean(Col.plot_height_array[-10:]))
        print("Dimensional Final Height: ",
              Col.length_scale * Col.plot_height_array[-1])
        print("Dimensional Initial Height: ",
              Col.length_scale * Col.plot_height_array[0])
        print("Slope Velocity: ",
              ((Col.plot_height_array[-1] - Col.plot_height_array[0]) /
               Col.runtime) * Col.length_scale / Col.time_scale)
        plt.figure()

    def Dimful_Liquid_Height_Plot(self, Col):

        plt.title(
            'Liquid Drained against Time (Dimensional Normalised Height)')
        plt.ylabel(r' Normalised Liquid Drained, $h$')
        plt.xlabel(r'Time Elapsed, $t$ ($s$)')
        plt.xlim(xmax=Col.time_scale * Col.runtime, xmin=0)
        plt.ylim(ymin=-0.00001)
        plt.grid(b=True, which='major', linestyle='--')
        plt.plot(Col.time_scale * self.plot_times2,
                 (Col.plot_height_increase_array) /
                 (self.Total_Liquid_In_Foam),
                 label='Normalised Liquid Height'
                 )  #height against time (for fixed velocity)
        #plt.plot(np.linspace(0, Col.runtime*Col.time_scale, 100), np.ones(100)*self.Equilibrium_Avg_lf, color = 'r', linestyle = '--', label = 'Theoretical End behaviour')
        plt.legend(loc='best')

    def Theoretical_SSH(self):
        self.theoretical_height_array = []
        self.bubble_vel = 0.0
        #self.bubble_vel = np.linspace(0, self.dimensional_crit_velocity, 1001)
        for i in range(0, 1000):
            if (self.bubble_vel == 0.0):
                self.theoretical_height = 1 / np.sqrt(
                    self.crit_lf_min) - 1 / np.sqrt(self.base_lf)
                #self.theoretical_height2 = self.length_scale * self.theoretical_height
            else:
                numerator1 = np.sqrt(
                    self.crit_lf_min * self.base_lf) - self.bubble_vel
                denominator1 = np.sqrt(self.bubble_vel) * (
                    np.sqrt(self.base_lf) - np.sqrt(self.crit_lf_min))
                arg_steadystate1 = numerator1 / denominator1
                coefficient1 = 1. / np.sqrt(self.bubble_vel)
                #arcoth1 = np.sinh(arg_steadystate1) / np.cosh(arg_steadystate1)
                arcoth1 = .5 * mt.log(
                    (1 + arg_steadystate1) / (arg_steadystate1 - 1))
                self.theoretical_height = coefficient1 * arcoth1
                #self.theoretical_height2 = self.length_scale * self.theoretical_height
            i += 1
            self.bubble_vel += self.dimless_crit_velocity / 1000.
            self.theoretical_height_array.append(self.theoretical_height)

    def Steady_State_Height_vs_Bubble_Velocity(self):
        self.Theoretical_SSH()
        self.bubble_vel = np.linspace(0, self.dimless_crit_velocity, 1000)
        self.critical_height = np.linspace(0, 100, 1000)
        self.simulated_bubble_velocity = [
            0.0000090500085, 0.00009050085, 0.00036200034, 0.0004525004,
            0.0005882505530445465, 0.0006787506381283229,
            0.0008145007657539876, 0.000905000850837764, 0.0011312510635472048,
            0.0013575012762566457, 0.0014932514038823103, 0.001810001701675528,
            0.0019005017867593044, 0.001991001871843081, 0.0020000518803514585,
            0.002009101888859836, 0.0020136268931140247
        ]
        self.simulated_steady_state_heights = [
            20.6, 20.9, 22.1, 22.5, 23.2, 23.8, 24.6, 25.3, 27.2, 29.8, 31.7,
            39.7, 45.2, 58, 61.2, 66.3, 71.3
        ]
        plt.title('Steady State Height vs Bubble Velocity')
        plt.xlabel('Bubble Velocity')
        plt.ylabel('Steady State Height')
        plt.ylim(ymax=90)
        plt.xlim(xmax=self.dimless_crit_velocity + 0.00015)
        plt.grid(b=True, which='major', linestyle='--')
        plt.plot(np.ones(1000) * self.dimless_crit_velocity,
                 self.critical_height,
                 'g',
                 linestyle='--',
                 label='Theoretical Critical Velocity')
        plt.plot(self.bubble_vel,
                 self.theoretical_height_array,
                 'r',
                 linestyle='--',
                 label='Theoretical Prediction')
        #plt.plot(self.dimless_bubble_velocity, np.mean(self.heights[-10:]), 'X', label = 'Calculated Height')
        plt.plot(self.simulated_bubble_velocity,
                 self.simulated_steady_state_heights,
                 'bX',
                 label='Calculated Heights')
        plt.legend(loc='best')

    def Plot_Dimensional(self, Col):
        self.Plot_Colors()
        self.Color_Gradient(True, Col)
        self.Dimful_Profile_Plot(Col)
        self.Dimful_Height_Plot(Col)
        self.Dimful_Liquid_Height_Plot(Col)
        self.Steady_State_Height_vs_Bubble_Velocity()

    def create_plots(self):
        plt.show()
