"""
CoRNy Figures
    Functions for Plotting and output
    REDUCED VERSION
"""


import pandas as pd
import numpy as np
import os
import re
import matplotlib
import matplotlib.pyplot as plt


def color_ramp( ramp_colors ): 
    from colour import Color
    from matplotlib.colors import LinearSegmentedColormap

    color_ramp = LinearSegmentedColormap.from_list( 'my_list', [ Color( c1 ).rgb for c1 in ramp_colors ] )
    plt.figure( figsize = (15,3))
    plt.imshow( [list(np.arange(0, len( ramp_colors ) , 0.1)) ] , interpolation='nearest', origin='lower', cmap= color_ramp )
    plt.xticks([])
    plt.yticks([])
    return color_ramp

##

def label_abc(axes, fontsize=10):
    abc = 'abcdefghijklmnopqrstuvwxyz'
    i = 0
    for ax in axes.flatten():
        ax.set_title('%s)' % abc[i], loc='left', fontsize=str(fontsize))
        i += 1

# Very Short Fig version, just to arrange the subplot environment
class Figure:

    def __init__(self, title='', layout=(1,1), size=(11.69,8.27), abc=False,
                 gridspec_kw={}, save=None, save_dpi=250, fig=None):
        
        self.layout = layout
        self.size   = size
        self.title  = title
        self.abc    = abc
        self.gridspec_kw = gridspec_kw
        self.save   = save
        self.save_dpi = save_dpi
        
    # Entering `with` statement
    def __enter__(self):
        self.fig, self.axes = plt.subplots(self.layout[0], self.layout[1],
                                           figsize=self.size, gridspec_kw=self.gridspec_kw)
        self.fig.suptitle(self.title)
        return(self.axes) # makes possibe: with Fig() as ax: ax.change

    # Exiting `with` statement
    def __exit__(self, type, value, traceback):
        
        if self.abc == True:
            label_abc(self.axes)
        
        if self.save:
            # Check and create folder
            parent_folders = os.path.dirname(self.save)
            if parent_folders and not os.path.exists(parent_folders):
                os.makedirs(parent_folders)
            # Save and close single plot
            self.fig.savefig(self.save, bbox_inches='tight',
                frameon=False, dpi=self.save_dpi, transparent=True)
            plt.show()
            plt.close()

