import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.patches as patches

from pathlib import Path
from matplotlib.ticker import AutoMinorLocator


def cm_to_inches(cm: float) -> float:
    """
    Function Details
    ================
    Returns centimeters as inches.

    Uses the conversion rate to convert a value given in centimeters to inches.
    Useful for matplotlib plotting.

    Parameters
    ----------
    cm : float
        Value of the desired figure size in centimeters.

    Returns
    -------
    inches : float
        Value of the desired figure size in inches.

    See Also
    --------
    None

    Notes
    -----
    Conversion rate given to 6 decimal places, but inches rounded to 2 decimal
    places.

    Examples
    --------
    >>> cm = 15
    >>> inches = cm_to_inches(cm=cm)
    >>> inches
    5.91

    ----------------------------------------------------------------------------
    Update History
    ==============

    25/05/2024
    ----------
    Documentation updated.

    """
    return round(cm * 0.393701, 2)


def S4_plot(wavelength : list,
            transmission : list,
            reflection : list,
            out_path : str,
            plot_dict : dict) -> None:
    """
    Function Details
    ================
    Plot transmission and reflection as a standard plot.

    Plot S4 transmission and reflection as a function of wavelength.

    Parameters
    ----------
    wavelength, transmission, reflection : list
        Wavelength, transmission, and reflection arrays from S4 spectral output.
    out_path : string
        Path to save.
    plot_dict : dictionary
        Plot settings dictionary containing:
            {
                "width": plot width,\n
                "height": plot height,\n
                "dpi": dots per square inch,\n
                "grid": True/False,\n
                "legend_loc": legend location,\n
                "legend_col": legend column number,\n
                "legend_size": size of legend text,\n
                "axis_fontsize": font size for axis labels,\n
                "label_size": size for tick labels
            }

    Returns
    -------
    None

    See Also
    --------
    matplotlib library

    Notes
    -----
    Plots a standard x, y graph for transmission ad reflection from the spectrum
    output of S4. The inputs are x, y1, y2 and requires a plot dictionary that
    governs the plot settings.

    Example
    -------
    None

    ----------------------------------------------------------------------------
    Update History
    ==============

    25/05/2024
    ----------
    Documentation updated.

    """
    fig, (ax1, ax2) = plt.subplots(
        nrows=2,
        ncols=1,
        figsize=[
            cm_to_inches(cm=plot_dict["width"]),
            cm_to_inches(cm=plot_dict["height"]) * 2],
        dpi=plot_dict["dpi"])
    line1 = ax1.plot(
        wavelength,
        transmission,
        'blue',
        lw=2,
        label='Transmission')
    line2 = ax2.plot(
        wavelength,
        reflection,
        'red',
        lw=2,
        label='Reflection')
    ax1.grid(
        visible=True,
        alpha=0.5)
    ax1.legend(
        frameon=True,
        loc=plot_dict["legend_loc"],
        ncol=plot_dict["legend_col"],
        prop={"size": plot_dict["legend_size"]})
    ax1.set_xlabel(
        'Wavelength [nm]',
        fontsize=plot_dict["axis_fontsize"],
        fontweight='bold')
    ax1.set_ylabel(
        'Transmission [au]',
        fontsize=plot_dict["axis_fontsize"],
        fontweight='bold')
    ax1.tick_params(
        axis='both',
        which='major',
        labelsize=plot_dict["label_size"])
    ax1.set_ylim(0, 1)
    ax1.xaxis.set_minor_locator(AutoMinorLocator())
    ax1.yaxis.set_minor_locator(AutoMinorLocator())
    ax2.grid(
        visible=True,
        alpha=0.5)
    ax2.legend(
        frameon=True,
        loc=plot_dict["legend_loc"],
        ncol=plot_dict["legend_col"],
        prop={"size": plot_dict["legend_size"]})
    ax2.set_xlabel(
        'Wavelength [nm]',
        fontsize=plot_dict["axis_fontsize"],
        fontweight='bold')
    ax2.set_ylabel(
        'Reflection [au]',
        fontsize=plot_dict["axis_fontsize"],
        fontweight='bold')
    ax2.tick_params(
        axis='both',
        which='major',
        labelsize=plot_dict["label_size"])
    ax2.set_ylim(0, 1)
    ax2.xaxis.set_minor_locator(AutoMinorLocator())
    ax2.yaxis.set_minor_locator(AutoMinorLocator())
    plt.savefig(
        out_path,
        bbox_inches='tight')
    fig.clf()
    plt.cla()
    plt.close(fig)


def Bragg_plots(wavelength : list,
                transmission : list,
                reflection : list,
                phase : list,
                out_path : str,
                plot_dict : dict) -> None:
    """
    Function Details
    ================
    Plot transmission, reflection, and phase as a standard plot.

    Parameters
    ----------
    wavelength, transmission, reflection, phase: list
        Wavelength, transmission, reflection, and phase arrays.
    out_path: string
        Path to save.
    plot_dict : dictionary
        Plot settings dictionary containing:
            {
                "width": plot width,\n
                "height": plot height,\n
                "dpi": dots per square inch,\n
                "grid": True/False,\n
                "legend_loc": legend location,\n
                "legend_col": legend column number,\n
                "legend_size": size of legend text,\n
                "axis_fontsize": font size for axis labels,\n
                "label_size": size for tick labels
            }

    Returns
    -------
    None.

    See Also
    --------
    matplotlib library

    Notes
    -----
    None.

    Example
    -------
    None.

    ----------------------------------------------------------------------------
    Update History
    ==============

    25/05/2024
    ----------
    Created.

    """
    fig, (ax1, ax2, ax3) = plt.subplots(
        nrows=3,
        ncols=1,
        figsize=[
            cm_to_inches(cm=plot_dict["width"]),
            cm_to_inches(cm=plot_dict["height"]) * 3],
        dpi=plot_dict["dpi"])
    ax1.plot(
        wavelength,
        transmission,
        'blue',
        lw=2,
        label='Transmission')
    ax2.plot(
        wavelength,
        reflection,
        'red',
        lw=2,
        label='Transmission')
    ax3.plot(
        wavelength,
        phase,
        'green',
        lw=2,
        label='Transmission')
    ax1.grid(
        visible=True,
        alpha=0.5)
    ax1.legend(
        frameon=True,
        loc=plot_dict["legend_loc"],
        ncol=plot_dict["legend_col"],
        prop={"size": plot_dict["legend_size"]})
    ax1.set_ylabel(
        'Transmission [au]',
        fontsize=plot_dict["axis_fontsize"],
        fontweight='bold')
    ax1.tick_params(
        axis='both',
        which='major',
        labelsize=plot_dict["label_size"])
    ax1.set_xlim(min(wavelength), max(wavelength))
    ax1.set_ylim(0, 1)
    ax1.xaxis.set_minor_locator(AutoMinorLocator())
    ax1.yaxis.set_minor_locator(AutoMinorLocator())
    ax2.grid(
        visible=True,
        alpha=0.5)
    ax2.legend(
        frameon=True,
        loc=plot_dict["legend_loc"],
        ncol=plot_dict["legend_col"],
        prop={"size": plot_dict["legend_size"]})
    ax2.set_ylabel(
        'Reflection [au]',
        fontsize=plot_dict["axis_fontsize"],
        fontweight='bold')
    ax2.tick_params(
        axis='both',
        which='major',
        labelsize=plot_dict["label_size"])
    ax2.set_xlim(min(wavelength), max(wavelength))
    ax2.set_ylim(0, 1)
    ax2.xaxis.set_minor_locator(AutoMinorLocator())
    ax2.yaxis.set_minor_locator(AutoMinorLocator())
    ax3.grid(
        visible=True,
        alpha=0.5)
    ax3.legend(
        frameon=True,
        loc=plot_dict["legend_loc"],
        ncol=plot_dict["legend_col"],
        prop={"size": plot_dict["legend_size"]})
    ax3.set_xlabel(
        'Wavelength [nm]',
        fontsize=plot_dict["axis_fontsize"],
        fontweight='bold')
    ax3.set_ylabel(
        'Phase [pi]',
        fontsize=plot_dict["axis_fontsize"],
        fontweight='bold')
    ax3.tick_params(
        axis='both',
        which='major',
        labelsize=plot_dict["label_size"])
    ax3.set_xlim(min(wavelength), max(wavelength))
    ax3.set_ylim(-2, 2)
    ax1.xaxis.set_minor_locator(AutoMinorLocator())
    ax1.yaxis.set_minor_locator(AutoMinorLocator())
    plt.savefig(
        out_path,
        bbox_inches='tight')
    fig.clf()
    plt.cla()
    plt.close(fig)


def stack_plots(stack_n : list,
                stack_t : list,
                plot_dict : dict,
                out_path : str) -> None:
    """
    Function Details
    ================
    Plot layer indices and thickness values as a function of number of layers.

    Parameters
    ----------
    stack_n, stack_t: list
        List of the refractive indices and layer thicknesses (nm) for the Bragg
        stack.
    out_path: string
        Path to save.
    plot_dict : dictionary
        Plot settings dictionary containing:
            {
                "width": plot width,\n
                "height": plot height,\n
                "dpi": dots per square inch,\n
                "grid": True/False,\n
                "legend_loc": legend location,\n
                "legend_col": legend column number,\n
                "legend_size": size of legend text,\n
                "axis_fontsize": font size for axis labels,\n
                "label_size": size for tick labels
            }

    Returns
    -------
    None.

    See Also
    --------
    None.

    Notes
    -----
    None.

    Example
    -------
    None.

    ----------------------------------------------------------------------------
    Update History
    ==============

    29/05/2024
    ----------
    Created.

    """
    fig, (ax1, ax2) = plt.subplots(
        nrows=1,
        ncols=2,
        figsize=[
            cm_to_inches(cm=plot_dict["width"]),
            cm_to_inches(cm=plot_dict["height"])],
        dpi=plot_dict["dpi"])
    ax1.plot(
        [i for i in range(len(stack_n))],
        stack_n,
        'blue',
        lw=2,
        label='Index Profile')
    ax1.set_xlabel(
        'Layer [#]',
        fontsize=plot_dict["axis_fontsize"],
        fontweight='bold')
    ax1.set_ylabel(
        'Refractive Index [RIU]',
        fontsize=plot_dict["axis_fontsize"],
        fontweight='bold')
    ax1.xaxis.set_minor_locator(AutoMinorLocator())
    ax1.yaxis.set_minor_locator(AutoMinorLocator())
    ax2.plot(
        [i for i in range(len(stack_t))],
        stack_t,
        'red',
        lw=2,
        label='Thickness Profile')
    ax2.set_xlabel(
        'Layer [#]',
        fontsize=plot_dict["axis_fontsize"],
        fontweight='bold')
    ax2.set_ylabel(
        'Layer Thickness [nm]',
        fontsize=plot_dict["axis_fontsize"],
        fontweight='bold')
    ax2.xaxis.set_minor_locator(AutoMinorLocator())
    ax2.yaxis.set_minor_locator(AutoMinorLocator())
    fig.tight_layout()
    plt.savefig(
        out_path,
        bbox_inches='tight')
    fig.clf()
    plt.cla()
    plt.close(fig)
