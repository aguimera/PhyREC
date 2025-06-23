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
    nFFT = 2**15

    def __init__(self):
        self.Fig = None

    def InitFigure(self):
        self.Fig, self.AxM = plt.subplots()
        self.AxP = plt.twinx(self.AxM)

    def PlotResponse(self, sos, Fs):
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
    ImgDicts = {}

    def __init__(self):
        self.Fig = None

    def GenColorBars(self):
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

