#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 17 20:23:10 2020

@author: aguimera
"""

import matplotlib.pyplot as plt
from scipy import signal
import numpy as np


class DebugFilterPlot():
    """
    Debug visualization tool for filter frequency response.

    Provides utilities to plot the magnitude and phase response of filter
    designs to aid in filter development and debugging.

    Attributes
    ----------
    nFFT : int
        FFT size for frequency response computation (default: 2**15).
    Fig : matplotlib.figure.Figure or None
        Matplotlib figure object, None until InitFigure is called.
    AxM : matplotlib.axes.Axes
        Magnitude response axes.
    AxP : matplotlib.axes.Axes
        Phase response axes.
    """
    nFFT = 2**15

    def __init__(self):
        """
        Initialize the debug filter plot tool.

        Parameters
        ----------
        None

        Returns
        -------
        None
        """
        self.Fig = None

    def InitFigure(self):
        """
        Create figure and axes for plotting filter response.

        Parameters
        ----------
        None

        Returns
        -------
        None
            Creates self.Fig, self.AxM, and self.AxP.
        """
        self.Fig, self.AxM = plt.subplots()
        self.AxP = plt.twinx(self.AxM)

    def PlotResponse(self, sos, Fs):
        """
        Plot the magnitude and phase response of a filter.

        Parameters
        ----------
        sos : array-like
            Second-order sections representation of the filter.
        Fs : float
            Sampling frequency of the filter.

        Returns
        -------
        None
            Plots magnitude and phase response on the figure.
        """
        if self.Fig is None:
            self.InitFigure()

        # w, h = signal.freqz(b, a,
        w, h = signal.sosfreqz(sos,
                               worN=self.nFFT,
                               fs=2*np.pi*Fs)
        ff = (w/(2*np.pi))[1:]
        self.AxM.semilogx(ff, np.abs(h[1:]))
        self.AxP.semilogx(ff, np.unwrap(np.angle(h)[1:]), '--')
        plt.show()


class ColorBarPlot():
    """
    Utility for generating and displaying colorbars for image plots.

    Manages a collection of image objects and displays their colorbars
    in a convenient grid layout.

    Attributes
    ----------
    ImgDicts : dict
        Dictionary of image objects where key is the image name/id and value
        is the image object from matplotlib.
    Fig : matplotlib.figure.Figure or None
        Matplotlib figure for colorbars, None until GenColorBars is called.
    """
    ImgDicts = {}

    def __init__(self):
        """
        Initialize the colorbar plotter.

        Parameters
        ----------
        None

        Returns
        -------
        None
        """
        self.Fig = None

    def GenColorBars(self):
        """
        Generate colorbars for all images in ImgDicts.

        Creates a new figure with one axis per image and displays
        the colorbar for each.

        Parameters
        ----------
        None

        Returns
        -------
        None
            Creates self.Fig with colorbars.
        """
        self.Fig, Ax = plt.subplots(1, len(self.ImgDicts))
        try:
            ax = Ax[0]
        except:
            Ax = [Ax, ]

        for iax, (name, img) in enumerate(self.ImgDicts.items()):
            plt.colorbar(img, cax=Ax[iax])
            Ax[iax].set_title(name)


def DrawBarScale(Ax, Location='Bottom Left',
                 xsize=None, ysize=None, xoff=0.1, yoff=0.1,
                 xlabelpad=-0.04, ylabelpad=-0.04,
                 xunit='sec', yunit='mV', LineWidth=5, Color='k',
                 FontSize=None, ylabel=None, xlabel=None):
    """
    Draw a scale bar on the specified axes for indicating data magnitudes.

    Adds a scale bar with optional labels to a matplotlib axes, useful for
    publications and presentations to indicate the scale of plotted data.

    Parameters
    ----------
    Ax : matplotlib.axes.Axes
        Axes object to draw the scale bar on.
    Location : str, optional
        Position of scale bar: 'Bottom Left' (default), 'Bottom Rigth', 'Top Left', 'Top Rigth'.
    xsize : float or int, optional
        Horizontal size of the scale bar. If None, defaults to 1/5 of axis width.
    ysize : float or int, optional
        Vertical size of the scale bar. If None, defaults to 1/5 of axis height.
    xoff : float, optional
        Horizontal offset from the specified location (default: 0.1).
    yoff : float, optional
        Vertical offset from the specified location (default: 0.1).
    xlabelpad : float, optional
        Padding for horizontal label relative to bar (default: -0.04).
    ylabelpad : float, optional
        Padding for vertical label relative to bar (default: -0.04).
    xunit : str, optional
        Unit label for horizontal scale (default: 'sec').
    yunit : str, optional
        Unit label for vertical scale (default: 'mV').
    LineWidth : int, optional
        Width of scale bar lines (default: 5).
    Color : str, optional
        Color of scale bar (default: 'k' for black).
    FontSize : int or str, optional
        Font size for labels. If None, uses default.
    ylabel : str, optional
        Custom label for vertical scale. If None, 'ysize yunit' is used.
    xlabel : str, optional
        Custom label for horizontal scale. If None, 'xsize xunit' is used.

    Returns
    -------
    None
        Modifies the axes in-place by adding scale bar graphics and labels.
    """

    # calculate length of the bars
    xmin, xmax, ymin, ymax = Ax.axis()
    AxTrans = Ax.transAxes
    if xsize is None:
        xsize = (xmax - xmin)/5
        xsize = int(np.round(xsize, 0))
    if ysize is None:
        ysize = (ymax - ymin)/5
        ysize = int(np.round(ysize, 0))
    xlen = 1/((xmax - xmin)/xsize)  # length in axes coord
    ylen = 1/((ymax - ymin)/ysize)

    # calculate locations
    if Location == 'Bottom Rigth':
        xoff = 1 - xoff
        ylabelpad = - ylabelpad
        xlen = - xlen
    elif Location == 'Top Left':
        yoff = 1 - yoff
        ylen = - ylen
        xlabelpad = -xlabelpad
    elif Location == 'Top Rigth':
        xoff = 1 - xoff
        ylabelpad = - ylabelpad
        xlen = - xlen
        yoff = 1 - yoff
        ylen = - ylen
        xlabelpad = -xlabelpad
    xdraw = xoff + xlen
    ydraw = yoff + ylen

    # Draw lines
    Ax.hlines(yoff, xoff, xdraw,
              Color,
              linewidth=LineWidth,
              transform=AxTrans,
              clip_on=False)

    if xlabel is None:
        xlabel = str(xsize) + ' ' + xunit

    Ax.text(xoff + xlen/2,
            yoff + xlabelpad,
            xlabel,
            horizontalalignment='center',
            verticalalignment='center',
            fontsize=FontSize,
            transform=AxTrans)

    Ax.vlines(xoff, yoff, ydraw,
              Color,
              linewidth=LineWidth,
              transform=AxTrans,
              clip_on=False)

    if ylabel is None:
        ylabel = str(ysize) + ' ' + yunit

    Ax.text(xoff + ylabelpad,
            yoff + ylen/2,
            ylabel,
            horizontalalignment='center',
            verticalalignment='center',
            rotation='vertical',
            fontsize=FontSize,
            transform=AxTrans)
