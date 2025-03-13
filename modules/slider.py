import numpy as np
import math as mt  # Very Sloppy
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.patches as mpatches
import matplotlib.animation as animation
from matplotlib.widgets import Slider, Button
import matplotlib, io, libconf


#Class animates a previously determined column instance
class slider:
    font = {
        'family': 'serif',
        'weight': 'normal',
        'size': 12
    }  # Dictates what the text style on plots is
    matplotlib.rc('font', **font)

    #Constructor sets up the animation and slider animation
    def __init__(self, Col):
        self.first_temp = np.zeros((1, Col.z_height))
        self.x = np.linspace(0., Col.z_height,
                             Col.z_height) * Col.length_scale * Col.dz
        self.y = np.linspace(0., .36, 2)
        self.background_colour = '#bdf8a3'
        colours = [self.background_colour, '#fffe7a',
                   '#b04e0f']  #Setting up custom colour map
        self.cm = colors.LinearSegmentedColormap.from_list('My', colours, 1000)
        self.dimless_col_height = Col.col_height / Col.length_scale

    #The slider portion of the code

    #Set up for the slider canvas stuff
    #TODO split up into functions
    def set_up_canvas(self, Col):
        fig = plt.figure()  # Make figure

        ax = plt.axes(xlim=(-.1, self.x[-1]), ylim=(0.0, .5))
        ax.xaxis.set_major_locator(
            plt.NullLocator())  #This couple of lines makes partial x-ticks
        ax.set_xticks(
            np.append(
                ax.get_xticks(),
                np.linspace(0., (Col.z_height) * Col.length_scale * Col.dz,
                            5)))

        ax.yaxis.set_major_locator(
            plt.NullLocator())  #This couple of lines makes partial y-ticks
        ax.set_yticks(np.append(ax.get_yticks(), np.linspace(0., .36, 7)))

        ax.spines['left'].set_position('zero')  #Puts the axis in the middle
        ax.spines['right'].set_color(None)

        plt.xlabel('Height (m)', x=.7)  # Labelling the axes
        plt.ylabel('Liquid Fraction', x=.1)
        ax.axvline(-0.03, color='k')

        ax.set_facecolor(
            self.background_colour)  # Setting colour of background
        ax.patch.set_alpha(0.8)

        fig.patch.set_facecolor('#02d8e9')  # Colour of Border of window
        fig.patch.set_alpha(0.7)

        plt.subplots_adjust(
            bottom=0.25
        )  # Raises the bottom of the picture so that the slider can fit

        self.first_temp[0, :] = Col.phifactor * Col.alpha_v_time[
            0, :]  # Plots the initial conditions
        im1, = ax.plot(self.x, Col.phifactor * Col.alpha_v_time[0, :], 'b')

        if (Col.input_lf == 0.0 and Col.dimless_bubble_velocity == 0.0):
            '''Analytic Solution for Free Drainage: Steady State flow solution'''
            lamda2 = Col.gamma / (Col.water_density * Col.gravity)
            z = np.arange(0, self.dimless_col_height, 0.001)
            phi_array = []
            for i, x in enumerate(z):
                if (x == 0):
                    phi2 = .36
                else:
                    phi = (mt.sqrt(3) / 2.0) * (Col.bubble_diameter /
                                                lamda2) * x
                    phi += 1.0 / mt.sqrt(Col.base_lf)
                    phi2 = 1.0 / (phi**2)
                phi_array.append(phi2)
            im3, = ax.plot(z, phi_array, '--r', label='Analytic Solution')
            plt.legend()

        self.im2 = ax.pcolormesh(self.x,
                                 self.y,
                                 self.first_temp,
                                 norm=colors.Normalize(vmin=0.,
                                                       vmax=Col.base_lf),
                                 cmap=self.cm)

        glass = ax.add_patch(
            mpatches.Rectangle((-.09, 0),
                               width=.035,
                               height=0.35,
                               figure=fig,
                               edgecolor='w',
                               fill=False,
                               alpha=1))  #'lightgoldenrodyellow'
        liquid = ax.add_patch(
            mpatches.Rectangle(
                (-.09, 0),
                width=.035,
                height=(0.35 / (Col.col_height)) *
                Col.plot_height_increase_array[0] * Col.length_scale,
                figure=fig,
                facecolor='k',
                alpha=1))
        head = ax.add_patch(
            mpatches.Rectangle((-.09, 0),
                               width=.035,
                               height=(0.35 / (Col.col_height)) *
                               Col.plot_height_array[0] * Col.length_scale,
                               figure=fig,
                               facecolor='#fffe7a',
                               ec='w',
                               alpha=1))  #facecolor='#fffe7a'

        axtime = plt.axes([.2, 0.1, 0.65, 0.03], facecolor='#d648d7'
                          )  #Gives slider position values and slider colours
        self.stime = Slider(axtime,
                            'Time (s)',
                            0.0,
                            Col.dt * Col.time_scale * Col.time_separation *
                            (Col.alpha_v_time[:, 0].size - 1),
                            valinit=0.0)

        def reset(event):
            self.stime.reset(
            )  # Function allows reset button to make everything zero again

        def update(val):
            time = int(
                self.stime.val * 1 /
                (Col.dt * Col.time_scale * Col.time_separation)
            )  #Function chooses image number to display using the slider

            im1.set_xdata(self.x[:])
            im1.set_ydata(Col.phifactor * Col.alpha_v_time[time, :].ravel())

            self.im2.set_array(Col.phifactor *
                               Col.alpha_v_time[time, :].ravel())
            col = self.im2.get_facecolor()

            liquid.set_height((0.35 / (Col.col_height)) *
                              Col.plot_height_increase_array[time] *
                              Col.phifactor * Col.length_scale)
            head.set_xy((-.09, (0.35 / (Col.col_height)) *
                         Col.plot_height_increase_array[time] * Col.phifactor *
                         Col.length_scale))
            head.set_height(
                (0.35 / (Col.col_height)) *
                (1. - ((Col.plot_height_increase_array[time] * Col.phifactor) /
                       self.dimless_col_height)) *
                Col.plot_height_array[time] * Col.length_scale)
            head.set_facecolor(col[int(Col.z_height / 2)])
            # Accounts for the liquid draining out; otherwise artificially raising the height

            if (Col.input_lf > 0.0 and Col.initial_lf < 0.01):
                '''Analytic solution with varying times for a solitary wave: wetting of a dry foam'''
                tao = time / Col.phifactor * (Col.time_separation * Col.dt)
                alpha_array = []
                x1 = 1 + self.x / Col.length_scale
                x0 = Col.z_height * Col.dz - 1.
                #x0 = self.dimless_col_height
                input_lf = Col.input_lf * Col.phifactor
                for i, x in enumerate(x1):
                    if ((x < (x0 - input_lf * tao))):
                        alpha = 0.
                    else:
                        alpha = input_lf * pow(
                            np.tanh(
                                mt.sqrt(input_lf) *
                                ((x - x0) + input_lf * tao)), 2)
                    alpha_array.append(alpha)
                im4, = ax.plot(self.x,
                               alpha_array,
                               '--r',
                               label='Analytic Solution')
                print(self.x / Col.length_scale)
                print(Col.z_height * Col.dz)
                print(self.dimless_col_height * Col.length_scale)
            fig.canvas.draw()
            if (Col.input_lf > 0.0 and Col.initial_lf < 0.01):
                im4.remove()

        self.stime.on_changed(update)

        resetax = plt.axes([0.8, 0.025, 0.1, 0.04])
        self.button = Button(resetax,
                             'Reset',
                             color='lightgoldenrodyellow',
                             hovercolor='0.975')
        self.button.on_clicked(reset)

    def create_slider(self, Col):
        self.set_up_canvas(Col)
        plt.show()
