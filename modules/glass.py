import numpy as np
import math as mt
import sys, libconf, io
from copy import deepcopy


#Class should construct a column of foam
class column:
    '''
    Constructor Initialises the most changed vars
    and sets up the 'alpha array' 
    Important to note that crit lf min is used as the crit liquid fraction in this program
    '''
    def __init__(self, filename):
        with io.open(filename) as Conf_file:  # Config file is read in
            Con = libconf.load(Conf_file)

        self.pour = True

        #Physical Constants
        self.gamma = Con.gamma
        C_factor = Con.C_factor
        f_factor = Con.f_factor
        self.gravity = Con.gravity
        self.water_density = Con.water_density
        water_viscosity = Con.water_viscosity
        self.bubble_diameter = Con.bubdiam  #m
        self.col_height = Con.col_height
        self.dt = Con.dt
        self.dz = Con.dz
        self.crit_lf_min = Con.crit_lf_min
        self.crit_lf_max = Con.crit_lf_max
        self.dimensionless_plots = Con.dimensionless_plots
        self.dimensional_plots = Con.dimensional_plots

        #Calculate parameters
        self.bubble_velocity = Con.bubble_vel
        self.runtime = Con.Dimless_runtime
        self.time_separation = Con.time_separation
        self.bubble_volume = mt.pi * self.bubble_diameter**3 / 6  #TODO unused?
        effective_viscosity = 3. * f_factor * water_viscosity
        self.length_scale = mt.sqrt(
            (C_factor * self.gamma) / (self.water_density * self.gravity))
        self.time_scale = effective_viscosity * (1. / mt.sqrt(
            C_factor * self.gamma * self.water_density * self.gravity))
        self.phifactor = (5.35 * self.length_scale * self.length_scale) / pow(
            ((4. / 3.) * mt.pi * (self.bubble_diameter * 0.5)**3), 2. / 3.)
        self.dimless_bubble_velocity = self.bubble_velocity * self.time_scale / self.length_scale
        maxtime = int(self.runtime / self.dt)
        self.dz1 = 1.0 / self.dz
        self.rz1 = self.dt * self.dz1  #d(tau)/d(xi)

        print("Dimensional Runtime: ", self.runtime * self.time_scale, "s")

        self.iterations = int(maxtime / self.time_separation)
        self.test = (np.linspace(0, self.iterations - 1, 100)).astype(int)

        self.z_height = 1 + int(self.col_height /
                                (self.length_scale * self.dz))
        self.plot_alpha_array = np.zeros([self.z_height, self.iterations + 1])
        self.plot_height_array = np.zeros([self.iterations + 1])
        self.plot_height_increase_array = np.zeros(self.iterations + 1)

        self.alpha_v_time = np.zeros([self.iterations + 1, self.z_height])
        phifactor_inverse = 1.0 / self.phifactor  # Avoid multiple floating point divisions
        self.input_lf = Con.input_lf * phifactor_inverse
        self.initial_lf = Con.initial_lf * phifactor_inverse
        self.base_lf = Con.base_lf * phifactor_inverse
        self.crit_lf = Con.crit_lf_min * phifactor_inverse  # liquid fraction below which bubbles rupture

        self.initialise_evolving_variables()

    def initialise_evolving_variables(self):
        self.dimless_base_liquid_height = 0.
        self.foam_height = 0
        self.alpha = np.zeros(
            self.z_height)  #matrix dim = fractional height x 2
        self.q = np.zeros(self.alpha.size)

    #TODO Function relating to the rate of liquid input, x is time, enter a function of time pref. with range [0,1]
    def apply_liquid_input_rate(self):
        return 1

    #Initialise the Alpha Array
    def initialise_alpha_array(self):
        self.alpha[:2] = self.base_lf
        self.alpha[2:-1] = self.initial_lf

    #Find height of foam at each instance in time
    def update_foam_height(self):
        self.foam_height = len(np.where(self.alpha[:-1] >= self.crit_lf)[0])

    def update_dimless_base_liquid_height(self):
        self.dimless_base_liquid_height += self.q[0] * self.dt
        if self.dimless_base_liquid_height <= 0.:
            self.dimless_base_liquid_height = 0.

    def update_height_variables(self):
        self.update_foam_height()
        self.update_dimless_base_liquid_height()
        #self.critical_liquid_fraction((self.foam_height - 1)*self.dz, self.dimless_base_liquid_height * self.phifactor)

        if self.dimless_bubble_velocity == 0.:
            self.rupture()

    #Rupture of Bubbles scheme goes in here
    def rupture(self):
        self.alpha[self.foam_height:] = 0.

    #Humidity factor determined here: critical liquid fraction is updated as a linear function of the height
    '''
    def critical_liquid_fraction(self, height_of_foam, height_of_liquid):
        self.crit_lf = self.crit_lf_max / self.phifactor
        self.crit_lf_minimum = self.crit_lf_min / self.phifactor
        self.column_height = self.col_height / self.length_scale
        self.crit_lf = self.crit_lf_max + (self.crit_lf_minimum - self.crit_lf_max) * ((self.column_height - height_of_foam * (1. - height_of_liquid / self.column_height) - height_of_liquid) / (self.column_height))
    '''

    def append_to_animation_array(self, i):
        self.alpha_v_time[i, :] = self.alpha
        self.plot_height_increase_array[i] = self.dimless_base_liquid_height

    def append_to_alpha_plot_array(self, i):
        self.plot_alpha_array[:, i] = self.alpha
        self.plot_height_array[i] = (self.foam_height - 1) * self.dz

    #Integration Scheme for Free and Forced Drainage
    def integrate_step_down(self):
        self.q[:self.foam_height -
               1] = self.alpha[1:self.foam_height]**2 + np.sqrt(
                   self.alpha[1:self.foam_height]) * 0.5 * np.ediff1d(
                       self.alpha[:self.foam_height]) * self.dz1

        if (self.pour):
            self.q[self.foam_height - 1] = (self.input_lf *
                                            self.apply_liquid_input_rate())**2
        self.q[self.foam_height:] = 0.
        self.alpha[1:] += self.rz1 * (np.ediff1d(self.q))

    #Integration Scheme for Gas Inflow at Bottom
    def integrate_step_up(self):
        self.q[0] = -self.base_lf * self.dimless_bubble_velocity + self.alpha[
            0] * self.alpha[0] + np.sqrt(self.alpha[0]) * 0.5 * np.ediff1d(
                self.alpha[0:2]) * self.dz1
        self.q[1:self.foam_height - 1] = -self.alpha[
            0:self.foam_height -
            2] * self.dimless_bubble_velocity + self.alpha[
                1:self.foam_height -
                1] * self.alpha[1:self.foam_height - 1] + np.sqrt(
                    self.alpha[1:self.foam_height - 1]) * 0.5 * np.ediff1d(
                        self.alpha[1:self.foam_height]) * self.dz1

        diff = np.ediff1d(self.alpha[self.foam_height - 1:self.foam_height +
                                     1])
        self.q[self.foam_height -
               1] = self.input_lf * self.input_lf - self.alpha[
                   self.foam_height -
                   2] * self.dimless_bubble_velocity + self.alpha[
                       self.foam_height -
                       1] * self.alpha[self.foam_height - 1] + np.sqrt(
                           self.alpha[self.foam_height -
                                      1]) * 0.5 * diff * self.dz1
        self.q[self.foam_height:] = 0.

        self.alpha[1:] += self.rz1 * np.ediff1d(self.q)

    def integrate(self):

        if self.input_lf != 0.:
            if self.dimless_bubble_velocity == 0.:
                self.integrate_step_down()
            else:
                print('Sorry, that set-up is currently invalid.')
                sys.exit()
        else:
            if self.dimless_bubble_velocity != 0.:
                self.integrate_step_up()

            if self.dimless_bubble_velocity == 0.:
                self.integrate_step_down()

    #Calls Integration scheme in loop and functions which save the array for later use
    def simulate(self):
        self.initialise_alpha_array()
        self.append_to_alpha_plot_array(0)
        self.append_to_animation_array(0)

        for i in range(self.iterations):  #time-evolve the system
            for j in range(self.time_separation):
                self.update_height_variables()
                self.integrate()

            self.append_to_alpha_plot_array(i + 1)
            self.append_to_animation_array(i + 1)

            if i in self.test:  # Says what percent of the way through the simulation it is at
                print('\r' + "{:0.2f} %".format(100 * i /
                                                (self.iterations - 1)),
                      flush=True,
                      end='')

        print('')
