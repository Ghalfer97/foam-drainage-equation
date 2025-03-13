#!/usr/bin/python

import sys
sys.path.insert(1, './modules')

import numpy as np
import math as mt  # Very Sloppy
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.patches as mpatches
import matplotlib.animation as animation
from matplotlib.widgets import Slider, Button
import matplotlib, io, libconf
from datetime import datetime
from modules.glass import column
from modules.slider import slider
from modules.animate import animator
from modules.FDE_Plots import PlotClass

simulated = False

while (True):
    # Reflect any updates to config file
    action = input('''
Enter one of the following keys to perform an action:
    s = interactive slider style animation
    a = automatic animation of the column
    p = output plots
    r = rerun the simulation
    q = quit
    
Please enter choice:''')

    if action == 's':
        if (simulated == False):
            print('Running simulation...')
            col = column("config.cfg")
            col.simulate()
            simulated = True

        start_time = datetime.now()
        print('Working on it boss...')
        slide = slider(col)
        slide.create_slider(col)
        end_time = datetime.now()
        print('Duration of slider: {}'.format(end_time - start_time))

    elif action == 'a':
        start_time = datetime.now()
        print('Working on it boss...')
        col = column("config.cfg")
        anime = animator(col)
        anime.animate()
        end_time = datetime.now()
        print('Duration of animation: {}'.format(end_time - start_time))

    elif action == 'p':
        if (simulated == False):
            print('Running simulation...')
            col = column("config.cfg")
            col.simulate()
            simulated = True

        plots = PlotClass(col)
        plots.create_plots()

    elif action == 'r':
        print('Running simulation...')
        col = column("config.cfg")
        col.simulate()

    else:
        break
