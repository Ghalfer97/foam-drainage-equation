import numpy as np
import math as mt  # Very Sloppy
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.patches as mpatches
import matplotlib.animation as animation
from matplotlib.widgets import Slider, Button
import matplotlib, io, libconf


#Class animates a previously determined column instance
class animator:
    font = {
        'family': 'serif',
        'weight': 'normal',
        'size': 12
    }  # Dictates what the text style on plots is
    matplotlib.rc('font', **font)

    #Constructor sets up the animation and slider animation
    def __init__(self, Col):
        self.col = Col
        self.x = np.linspace(0., self.col.z_height,
                             self.col.z_height) * self.col.length_scale * self.col.dz
        self.y = np.linspace(0., .36, 2)
        self.interval = 200
        self.delay = 200
        self.pause = False
        self.background_colour = '#bdf8a3'
        colours = [self.background_colour, '#fffe7a',
                   '#b04e0f']  #Setting up custom colour map
        self.cm = colors.LinearSegmentedColormap.from_list('My', colours, 1000)
        self.dimless_col_height = self.col.col_height / self.col.length_scale

    def set_initial_conditions(self):
        self.col.initialise_evolving_variables()
        self.col.initialise_alpha_array()
        self.col.update_height_variables()
        self.col.integrate()  #TODO maybe this gets t=1, not t=0

        #TODO don't need temp?
        self.temp = np.zeros((1, self.col.z_height))
        self.temp[0, :] = self.col.phifactor * self.col.alpha[0]
        self.image_one = plt.pcolor(self.x,
                                    self.y,
                                    self.temp,
                                    norm=colors.Normalize(vmin=0.,
                                                          vmax=self.col.base_lf),
                                    cmap=self.cm)
        self.image_two, = plt.plot(self.x, self.col.phifactor * self.col.alpha, 'r')

        #Relevant for the animation's attributes
        self.glass = self.ax.add_patch(
            mpatches.Rectangle((-.09, 0),
                               width=.035,
                               height=0.35,
                               figure=self.fig,
                               edgecolor='lightgoldenrodyellow',
                               fill=False,
                               alpha=1))
        self.liquid = self.ax.add_patch(
            mpatches.Rectangle(
                (-.09, 0),
                width=.035,
                height=(0.35 / (self.col.col_height)) *
                self.col.dimless_base_liquid_height * self.col.length_scale,
                figure=self.fig,
                facecolor='k',
                alpha=1))
        self.head = self.ax.add_patch(
            mpatches.Rectangle(
                (-.09, 0),
                width=.035,
                height=(0.35 / (self.col.col_height)) * (self.col.foam_height - 1) *
                self.col.dz * self.col.length_scale,
                figure=self.fig,
                facecolor='#fffe7a',
                alpha=1))

        return self.head, self.liquid

    def set_up_canvas(self):
        self.fig = plt.figure()

        self.ax = plt.axes(xlim=(-.1, self.x[-1]), ylim=(-0.01, .5))

        self.ax.xaxis.set_major_locator(
            plt.NullLocator())  #This couple of lines makes partial x-ticks
        self.ax.set_xticks(
            np.append(
                self.ax.get_xticks(),
                np.linspace(0., (self.col.z_height) * self.col.length_scale * self.col.dz,
                            5)))

        self.ax.yaxis.set_major_locator(
            plt.NullLocator())  #This couple of lines makes partial y-ticks
        self.ax.set_yticks(
            np.append(self.ax.get_yticks(), np.linspace(0., .36, 7)))

        self.ax.spines['left'].set_position(
            'zero')  #Puts the axis in the middle
        self.ax.spines['right'].set_color(None)

        plt.xlabel('Height (m)', x=.7)  # Labelling the axes
        plt.ylabel('Liquid Fraction', x=.1)

        self.ax.axvline(-0.03, color='k')

        self.ax.set_facecolor('#bdf8a3')  # Colour of background
        self.ax.patch.set_alpha(0.8)

        self.fig.patch.set_facecolor('#02d8e9')  # Colour of Border
        self.fig.patch.set_alpha(0.7)

        self.time_template = 'Time = %.2f s'  # Timer Value and position
        self.time_text = self.ax.text(0.65,
                                      0.9,
                                      '',
                                      transform=self.ax.transAxes)

        self.points_glass = [[-.065, 0], [-0.045, 0], [-0.04, self.col.col_height],
                             [-.07, self.col.col_height]]
        self.points_liquid = [[-.065, 0], [-0.045, 0], [-.065, 0], [-0.045, 0]]

        self.fig.canvas.mpl_connect('button_press_event', self.click_pause)

    #Pauses/Unpauses the Animation on a screen click
    def click_pause(self, event):
        if self.pause:
            self.window[0].event_source.stop()

        else:
            self.window[0].event_source.start()

        self.pause ^= True

    #Gives the timer a time value
    def timer(self, time):
        self.time_text.set_text(
            self.time_template %
            (time * self.col.dt * self.col.time_scale * self.col.time_separation))

    #Makes the foam head which represents the foam
    def get_head_height(self, time):
        self.head.set_height(
            (0.35 / (self.col.col_height)) *
            (1. - ((self.col.dimless_base_liquid_height * self.col.phifactor) /
                   self.dimless_col_height)) * (self.col.foam_height - 1) * self.col.dz *
            self.col.length_scale)

    #Makes the foam rise to correct location
    def get_head_length(self, time):
        self.head.set_xy(
            (-.09, (0.35 / (self.col.col_height)) * self.col.dimless_base_liquid_height *
             self.col.phifactor * self.col.length_scale))

    #Makes the black Rectangle which represents the fluid
    def get_liquid_height(self, time):
        self.liquid.set_height(
            (0.35 / (self.col.col_height)) * self.col.dimless_base_liquid_height *
            self.col.phifactor * self.col.length_scale)

    #Calls alpha array and update array info
    def update_images(self, time):
        temp = self.col.phifactor * self.col.alpha
        self.temp[0:1, :] = self.col.alpha

        self.image_one = plt.pcolor(self.x,
                                    self.y,
                                    self.temp,
                                    norm=colors.Normalize(vmin=0.,
                                                          vmax=self.col.base_lf),
                                    cmap=self.cm)
        self.image_two.set_data(self.x, temp)

    def get_next_frame(self, time):
        for i in range(self.col.time_separation):
            self.col.update_height_variables()
            self.col.integrate()

        self.timer(time)
        self.get_head_height(time)
        self.get_head_length(time)
        self.get_liquid_height(time)
        self.update_images(time)
        return self.head, self.liquid, self.image_one, self.image_two, self.time_text

    def generate_animation(self):
        self.window = animation.FuncAnimation(self.fig,
                                    self.get_next_frame,
                                    np.arange(1, self.col.iterations),
                                    init_func=self.set_initial_conditions,
                                    interval=self.interval,
                                    repeat=True,
                                    blit=True,
                                    repeat_delay=self.delay,
        )

    def animate(self):
        self.set_up_canvas()
        self.generate_animation()
        plt.show()
