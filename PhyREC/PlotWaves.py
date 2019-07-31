#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 26 12:11:53 2017

@author: aguimera
"""

import numpy as np
import matplotlib.pyplot as plt
import quantities as pq
from scipy import signal
import matplotlib.colors as colors
from collections import OrderedDict
from scipy.interpolate import interp2d
from matplotlib.widgets import Slider, Button, TextBox
from matplotlib.artist import ArtistInspector

#from NeoInterface import NeoTrain

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


class SpecSlot():

    AxKwargs = {      
                    'xaxis': {'visible': False,
                      },
                    'yaxis': {'visible': True,
                      },
                    'ylabel':'Freq [Hz]',
                    }
    def UpdateLineKwargs(self, LineKwargs):
        pass
#        self.LineKwargs.update(LineKwargs)
#        UpdateTreeDictProp(self.Line, self.LineKwargs)

    def UpdateAxKwargs(self, AxKwargs):
        pass
#        self.AxKwargs.update(AxKwargs)
#        UpdateTreeDictProp(self.Ax, self.AxKwargs)
        
    def __init__(self, Signal, Units=None, Position=None, DispName=None):
        self.Fres = 5.0
        self.TimeRes = 0.2
        self.Fmin = 1
        self.Fmax = 150
        self.MinPSD = None
        self.MaxPSD = None
        self.MaxPSDrange = None
        self.LogNormalize = True
        self.Cmap = 'jet'
        self.LogScale = False

        self.Ax = None
        self.CAx = None
        self.Fig = None

        self.Signal = Signal

        self.units = Units
        self.Position = Position
        if DispName is None:
            if Signal is not None:
                self.DispName = self.Signal.name
        else:
            self.DispName = DispName

    def GetSignal(self, Time, Units=None):
        if Units is None:
            _Units = self.units
        else:
            _Units = Units
        sig = self.Signal.GetSignal(Time, _Units)
        self.units = sig.units
        return sig

    def PlotSignal(self, Time, Units=None):
        sig = self.GetSignal(Time, Units)

        nFFT = int(2**(np.around(np.log2(sig.sampling_rate.magnitude/self.Fres))+1))
        Ts = sig.sampling_period.magnitude
        noverlap = int((Ts*nFFT - self.TimeRes)/Ts)

        f, t, Sxx = signal.spectrogram(sig,
                                       fs=sig.sampling_rate,
                                       window='hanning',
                                       nperseg=nFFT,
                                       noverlap=noverlap,
                                       scaling='density',
                                       axis=0)

        finds = np.where((self.Fmin < f) & (f < self.Fmax))[0][1:]
        r, g, c = Sxx.shape
        data = Sxx.reshape((r, c))[finds][:]

        if self.MaxPSD is None:
            MaxPSD = data.max()
        else:
            MaxPSD = self.MaxPSD

        if self.MinPSD is None:
            if self.MaxPSDrange is None:
                MinPSD = 10**(np.log10(MaxPSD)-5)
            else:
                MinPSD = 10**(np.log10(MaxPSD)-self.MaxPSDrange)
        else:
            MinPSD = self.MinPSD

        if self.LogNormalize:
            Norm = colors.LogNorm(MinPSD, MaxPSD)
        else:
            Norm = colors.Normalize(MinPSD, MaxPSD)

        if self.Ax is None:
            self.Fig, self.Ax = plt.subplots()

        x = t + sig.t_start.magnitude
        y = f[finds].magnitude
        img = self.Ax.imshow(data,
                             cmap='jet',
                             norm=Norm,
                             interpolation='bilinear',
                             origin='lower',
                             aspect='auto',
                             extent=(np.min(x), np.max(x), np.min(y), np.max(y)))

        cbar = plt.colorbar(img, ax=self.CAx, fraction=0.8)
        cbar.ax.tick_params(length=1, labelsize='xx-small')

        su = str(sig.units).split(' ')[-1]
        label = "[{}^2]".format(su)
        cbar.set_label(label, fontsize='xx-small')

        if self.LogScale:
            self.Ax.set_yscale('log')


def UpdateTreeDictProp(obj, prop):
    ains = ArtistInspector(obj)
    validp = ains.get_setters()
    for p in prop.keys():
        if p in validp:
            obj.set(**{p: prop[p]})
        else:
            obj2 = getattr(obj, 'get_' + p)()
            UpdateTreeDictProp(obj2, prop[p])


class SpikeSlot():
    DefLineKwargs = {'color': 'r',
                     'linestyle': '-.',
                     'alpha': 0.5,
                     'linewidth': 0.5,
                     }
    DefAxKwargs = {}

    def UpdateLineKwargs(self, LineKwargs):
        self.LineKwargs.update(LineKwargs)
        UpdateTreeDictProp(self.Line, self.LineKwargs)

    def UpdateAxKwargs(self, AxKwargs):
        self.AxKwargs.update(AxKwargs)
        UpdateTreeDictProp(self.Ax, self.AxKwargs)

    def __init__(self, Signal, Units='s',
                 Position=None, Ax=None, AxKwargs=None,
                 **LineKwargs):

        self.LineKwargs = self.DefLineKwargs.copy()
        self.AxKwargs = self.DefAxKwargs.copy()
        self.Signal = Signal
#        self.Signal.__class__ = NeoTrain
        self.name = self.Signal.name
        self.Position = Position
        self.Ax = Ax
        self.units = Units

        if AxKwargs is not None:
            self.AxKwargs.update(AxKwargs)

        if self.Ax is not None:
            UpdateTreeDictProp(self.Ax, self.AxKwargs)
        self.LineKwargs.update(LineKwargs)

    def GetSignal(self, Time, Units=None):
        if Units is None:
            _Units = self.units
        else:
            _Units = Units
        sig = self.Signal.GetSignal(Time, _Units)
        self.units = sig.units
        return sig

    def PlotSignal(self, Time, Units=None):
        if self.Ax is None:
            self.Fig, self.Ax = plt.subplots()

        sig = self.GetSignal(Time, Units)

        xmin, xmax, ymin, ymax = self.Ax.axis()

        self.Lines = self.Ax.vlines(sig,
                                    ymin, ymax,
                                    **self.LineKwargs
                                    )


class WaveSlot():

    DefLineKwargs = {'color': 'k',
                     'linestyle': '-',
                     'alpha': 1,
                     'linewidth': 0.5,
                     'clip_on': True,
                     }
    DefAxKwargs = {}

    def UpdateLineKwargs(self, LineKwargs):
        self.LineKwargs.update(LineKwargs)
        UpdateTreeDictProp(self.Line, self.LineKwargs)

    def UpdateAxKwargs(self, AxKwargs):
        self.AxKwargs.update(AxKwargs)
        UpdateTreeDictProp(self.Ax, self.AxKwargs)

    def __init__(self, Signal, Units=None, UnitsInLabel=False,
                 Position=None, Ax=None, AxKwargs=None,
                 **LineKwargs):
        self.LineKwargs = self.DefLineKwargs.copy()
        self.AxKwargs = self.DefAxKwargs.copy()

        self.Signal = Signal
        self.name = self.Signal.name

        self.units = Units

        self.Position = Position
        self.UnitsInLabel = UnitsInLabel

        self.Ax = Ax
        if AxKwargs is not None:
            self.AxKwargs.update(AxKwargs)

        if self.Ax is not None:
            UpdateTreeDictProp(self.Ax, self.AxKwargs)

        self.LineKwargs.update(LineKwargs)
        if 'label' not in self.LineKwargs:
            self.LineKwargs.update({'label': self.name})
        else:
            self.name = self.LineKwargs['label']

    def GetSignal(self, Time, Units=None):
        if Units is None:
            _Units = self.units
        else:
            _Units = Units
        sig = self.Signal.GetSignal(Time, _Units)
        self.units = sig.units
        return sig

    def PlotSignal(self, Time, Units=None):
        if self.Ax is None:
            self.Fig, self.Ax = plt.subplots()

        sig = self.GetSignal(Time, Units)

        if self.UnitsInLabel is True:
            su = str(sig.units).split(' ')[-1]
            label = "{} [{}]".format(self.name, su)
        else:
            label = self.name
        self.LineKwargs.update({'label': label})

        self._PlotSignal(sig)

    def _PlotSignal(self, sig):
        self.Lines = self.Ax.plot(sig.times.rescale('s'),
                                  sig,
                                  **self.LineKwargs
                                  )
        self.Ax.set_xlim(left=sig.t_start.rescale('s').magnitude,
                         right=sig.t_stop.rescale('s').magnitude)
        self.Line = self.Lines[0]

    def CalcAvarage(self, TimeAvg, TimesEvent, Units=None,
                    PltStd=False, StdAlpha=0.2,
                    PlotTrials=False, TrialsColor='k', TrialsAlpha=0.01):
        avsig = self.GetSignal(None, Units)
        avg = np.array([])

        Ts = avsig.sampling_period
        nSamps = int((TimeAvg[1]-TimeAvg[0])/Ts)
        t = np.arange(nSamps)*Ts + TimeAvg[0]

        for et in TimesEvent:
            start = et+TimeAvg[0]
            stop = et+TimeAvg[1]

            st = np.array(avsig.GetSignal((start, stop))[:nSamps])
            try:
                avg = np.hstack([avg, st]) if avg.size else st
                if PlotTrials:
                    self.Ax.plot(t, st,
                                 color=TrialsColor,
                                 alpha=TrialsAlpha)
#                                 clip_on=self.clip_on)
            except:
                print ('Error', nSamps, et, avg.shape, st.shape)

        MeanT = np.mean(avg, axis=1)

        MeanTsig = avsig.duplicate_with_new_array(signal=MeanT*avsig.units)
        MeanTsig.t_start = TimeAvg[0]
        MeanTsig.name = MeanTsig.name
        MeanTsig.annotate(Process='LED averaging')
#        self._PlotSignal(MeanTsig, label=self.DispName + ' Avg')
        self._PlotSignal(MeanTsig)


        if PltStd:
            StdT = np.std(avg, axis=1)
            self.Ax.fill_between(t, MeanT+StdT, MeanT-StdT,
                                 alpha=StdAlpha,
                                 facecolor=TrialsColor,
                                 edgecolor=None)

        ylim = self.Ax.get_ylim()
        self.Ax.vlines((0,), ylim[0], ylim[1], 'r', 'dashdot', alpha=0.5)
        return MeanTsig


class ControlFigure():

    def __init__(self, pltSL, figsize=(20*0.394, 5*0.394)):

        self.pltSL = pltSL

        TMax = np.max([sl.Signal.t_stop.rescale('s') for sl in pltSL.Slots])
        TMin = np.min([sl.Signal.t_start.rescale('s') for sl in pltSL.Slots])

        self.Fig, ax = plt.subplots(4, 1, figsize=figsize)
        self.sTstart = Slider(ax[0],
                              label='TStart [s]',
                              valmax=TMax,
                              valmin=TMin,
                              valinit=TMin)

        self.sTshow = Slider(ax[1],
                             label='TShow [s]',
                             valmax=TMax-TMin,
                             valmin=0,
                             valinit=(TMax-TMin)/10)

        self.sTshow.on_changed(self.Update)
        self.sTstart.on_changed(self.Update)

        self.TextStart = TextBox(ax[2],
                                 'Start time [s]',
                                 initial='0')
        self.TextStart.on_submit(self.submit_start)

        self.TextStop = TextBox(ax[3],
                                 'Stop time [s]',
                                 initial='10')
        self.TextStop.on_submit(self.submit_stop)
        
        self.Refresh = True
        self.OldStart = 0
        self.OldStop = 0

    def Update(self, val):
        twind = (self.sTstart.val * pq.s,
                 self.sTstart.val * pq.s + self.sTshow.val * pq.s)

        if self.Refresh:
            self.UpdateGraph(twind)

    def UpdateGraph(self, twind):
        self.pltSL.PlotChannels(twind)
        self.pltSL.Fig.canvas.draw()

    def submit_start(self, text):
        try:
            val = float(text)
        except:
            print('bad number')
            return

        if self.OldStart == val:
            return

        if val > self.OldStop:
            self.TextStart.set_val(str(self.sTstart.val))
            return
    
        self.OldStart = val
        self.sTstart.set_val(float(text))

    def submit_stop(self, text):
        try:
            val = float(text)
        except:
            print('bad number')
            return

        if self.OldStop == val:
            return
        
        if val < self.sTstart.val:
            self.TextStop.set_val(str(self.sTstart.val + self.sTshow.val))
            return
        
        self.OldStop = val
    
        show = val - self.sTstart.val 
        self.sTshow.set_val(show)

    def SetTimes(self, twind):
        self.Refresh = False
        self.TextStop.set_val(str(twind[1]))
        self.TextStart.set_val(str(twind[0]))
        self.TextStop.set_val(str(twind[1]))
        self.Refresh = True


class PlotSlots():
    ScaleBarKwargs = {'Location': 'Bottom Left',
                      'xsize': None,
                      'ysize': None,
                      'xoff': 0.1,
                      'yoff': 0.1,
                      'xlabelpad': -0.04,
                      'ylabelpad': -0.04,
                      'xunit': 'sec',
                      'yunit': None,
                      'LineWidth': 5,
                      'Color': 'k',
                      'FontSize': None}

    RcGeneralParams = {
#                       'axes.spines.left': False,
#                       'axes.spines.bottom': False,
#                       'axes.spines.top': False,
#                       'axes.spines.right': False,
                       }

    FigKwargs = {}

    gridspec_Kwargs = {'width_ratios': (15, 1)}

    TimeAxisProp = {'xaxis': {'visible': True, },
                    'xlabel': 'Time [s]',
                    }

    LegendKwargs = {'fontsize': 'xx-small',
                    'ncol': 5,
                    'loc': 'upper right',
                    'frameon': False}

    def UpdateFigKwargs(self, FigKwargs):
        self.FigKwargs.update(FigKwargs)
        UpdateTreeDictProp(self.Fig, self.FigKwargs)

    def _GenerateFigure(self):

        Pos = []
        for isl, sl in enumerate(self.Slots):
            if sl.Position is None:
                sl.Position = isl
            Pos.append(sl.Position)

        self.Fig, A = plt.subplots(max(Pos) + 1, 2,
                                   sharex=True,
                                   gridspec_kw=self.gridspec_Kwargs,
                                   )

        if len(A.shape) == 1:
            A = A[:, None].transpose()
        self.Axs = [a[0] for a in A]
        self.CAxs = [a[1] for a in A]

        for ca in self.CAxs:
            ca.axis('off')

        for sl in self.Slots:
            if isinstance(sl, SpecSlot):
                sl.CAx = self.CAxs[sl.Position]
            sl.Ax = self.Axs[sl.Position]
            sl.Fig = self.Fig
            UpdateTreeDictProp(sl.Ax, sl.AxKwargs)


    def __init__(self, Slots, Fig=None, FigKwargs=None, RcGeneralParams=None,
                 AxKwargs=None, TimeAxis=-1, CalcSignal=True,
                 ScaleBarAx=None, LiveControl=False):

        if RcGeneralParams is not None:
            self.RcGeneralParams.update(RcGeneralParams)
        plt.rcParams.update(self.RcGeneralParams)

        if FigKwargs is not None:
            self.FigKwargs.update(FigKwargs)

        self.Slots = Slots
        if CalcSignal:
            for sl in self.Slots:
                sig = sl.Signal
                sl.Signal = sig.GetSignal(None)

        self.ScaleBarAx = ScaleBarAx

        if LiveControl:
            self.CtrFig = ControlFigure(self)
        else:
            self.CtrFig = None

        if Fig is None:
            self._GenerateFigure()
        else:
            self.Fig = Fig
            self.Axs = []
            self.CAxs = []
            for sl in self.Slots:
                self.Axs.append(sl.Ax)

        for sl in self.Slots:
            if AxKwargs is not None:
                sl.UpdateAxKwargs(AxKwargs)

        self.TimeAxis = TimeAxis
        if self.TimeAxis is not None:
            sl = self.Slots[TimeAxis]
            sl.UpdateAxKwargs(self.TimeAxisProp)

        UpdateTreeDictProp(self.Fig, self.FigKwargs)
        self.SortSlotsAx()

    def SortSlotsAx(self):
        self.SlotsInAxs = {}
        for ax in self.Axs:
            sll = []
            for sl in self.Slots:
                if sl.Ax == ax:
                    sll.append(sl)
            self.SlotsInAxs.update({ax: sll})

    def ClearAxes(self):
        for sl in self.Slots:
            while sl.Ax.lines:
                sl.Ax.lines[0].remove()

    def FormatFigure(self):

        if self.ScaleBarAx is not None:
            if self.ScaleBarKwargs['yunit'] is None:
                sl = self.SlotsInAxs[self.Axs[self.ScaleBarAx]][0]
                su = str(sl.units).split(' ')[-1]
                self.ScaleBarKwargs['yunit'] = su
            DrawBarScale(self.Axs[self.ScaleBarAx], **self.ScaleBarKwargs)

    def AddLegend(self, **LegendKwargs):
        self.LegendKwargs.update(LegendKwargs)
        for Ax in self.Axs:
            Ax.legend(**self.LegendKwargs)

    def PlotChannels(self, Time, Units=None, FormatFigure=True):
        self.ClearAxes()
        print('plot channels')
        for sl in self.Slots:
            sl.PlotSignal(Time, Units=Units)
#        if Time is not None:
#            if Time[0] is not None:
#                sl.Ax.set_xlim(left=Time[0].magnitude)
#            if Time[1] is not None:
#                sl.Ax.set_xlim(right=Time[1].magnitude)

        self.current_time = sl.Ax.get_xlim()
        
        if self.CtrFig is not None:
            self.CtrFig.SetTimes(self.current_time)
            
    def PlotEvents(self, Times, color='r', alpha=0.5,
                   Labels=None, lAx=0, fontsize='xx-small', LabPosition='top'):

        self.Texts = []
        if Labels is not None:
            for ilbl, lbl in enumerate(Labels):
                for ax in self.Axs:
                    ylim = ax.get_ylim()
                    ax.vlines(Times[ilbl], ylim[0], ylim[1],
                              color=color,
                              alpha=alpha)
                lax = self.Axs[lAx]
                if LabPosition == 'top':
                    ylim = lax.get_ylim()[1]
                else:
                    ylim = lax.get_ylim()[0]
                txt = lax.text(Times[ilbl], ylim, lbl, fontsize=fontsize)
                self.Texts.append(txt)
            return

        EventLines = []
        for ax in self.Axs:
            ylim = ax.get_ylim()
            lines = ax.vlines(Times, ylim[0], ylim[1],
                              color=color,
                              alpha=alpha)
#            EventLines.append(lines[0])

        return EventLines

    def PlotEventAvarage(self, TimeAvg, TimesEvent, Units=None, PltStd=False,
                         StdAlpha=0.2,
                         PlotTrials=False, TrialsColor='k', TrialsAlpha=0.01,
                         ClearAxes=True, AvgColor=None):

        if ClearAxes:
            self.ClearAxes()

        MeanSigs = []
        for sl in self.Slots:
            if AvgColor is not None:
                sl.Color = AvgColor
            MeanSig = sl.CalcAvarage(TimeAvg, TimesEvent,
                                     Units=Units,
                                     PlotTrials=PlotTrials,
                                     TrialsColor=TrialsColor,
                                     TrialsAlpha=TrialsAlpha,
                                     PltStd=PltStd,
                                     StdAlpha=StdAlpha)
            MeanSigs.append(MeanSig)

        sl.Ax.set_xlim(left=TimeAvg[0].magnitude)
        sl.Ax.set_xlim(right=TimeAvg[1].magnitude)

        self.FormatFigure()
        return MeanSigs
