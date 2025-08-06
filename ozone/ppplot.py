import matplotlib.pyplot as plt
import matplotlib as mpl
def prepare_plot(reset=False, defaults=False, close_all=True, fontsize=12.,
        font={'family':'serif', 'serif': ['DejaVu Sans']}, linewidth=1.5,
        markersize=4.0, axis=None, cmap='cividis', equal_AR=False):
    '''
    This function sets parameters for plotting.

    Inputs:
    -------
        reset: if True, will reset to default parameters before setting
            input arguments
        defaults: if True, will reset to default parameters and then
            immediately return (input arguments not set)
        close_all: if True, will close all current figures
        fontsize: font size
        font: font
        linewidth: line width
        markersize: size of markers
        axis: axis limits [2*ndims]
        cmap: colormap
        equal_AR: if True, will set equal aspect ratio (only affects 2D)
    '''
    if reset or defaults:
        mpl.rcdefaults() # return to default settings
        if defaults:
            return

    # Use tex syntax
    plt.rc('text', usetex=True)

    # Set parameters
    if close_all:
        plt.close("all")
    mpl.rcParams['font.size'] = fontsize
    plt.rc('font',**font)
    # font={'family':'serif', 'serif': ['computer modern roman']}
    mpl.rcParams['lines.linewidth'] = linewidth
    mpl.rcParams['lines.markersize'] = markersize
    # Note: cividis is more colorblind-friendly than viridis (which is already pretty colorblind-friendly)
    mpl.rcParams['image.cmap'] = cmap
    # Color cycle
    # plt.style.use('tableau-colorblind10')
    try:
        plt.style.use('seaborn-colorblind')
    except OSError:
        plt.style.use('seaborn-v0_8-colorblind')

    if axis is not None:
        plt.axis(axis)
    if equal_AR:
        plt.gca().set_aspect('equal', adjustable='box')
