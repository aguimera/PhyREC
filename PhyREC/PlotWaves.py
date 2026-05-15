#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
PlotWaves module - Comprehensive visualization tools for electrophysiology signals.

This module provides a suite of classes for plotting and visualizing various types
of electrophysiological data including raw waveforms, spectrograms, spike trains,
and 2D image data. It supports interactive visualization with live controls and
animation capabilities.

Classes:
    SlotBase: Base class for signal plotting slots
    WavesColorSlot: Slot for color-mapped multi-channel waveforms
    SpecSlot: Slot for spectrogram visualization
    SpikeSlot: Slot for spike train visualization
    WaveSlot: Slot for time-series waveform plotting
    ImgSlot: Slot for 2D image data
    ControlFigure: Interactive control panel for visualization navigation
    PlotSlots: Main container for managing multiple plot slots

Created on Wed Jul 26 12:11:53 2017

@author: aguimera
"""

import numpy as np
import matplotlib.pyplot as plt
import quantities as pq
import matplotlib.colors as colors
from matplotlib.widgets import Slider, Button, TextBox
from matplotlib.artist import ArtistInspector
import PhyREC.SignalProcess as Spro
from PhyREC import SpectColBars
from PhyREC.SupportFunctions import DrawBarScale
import datetime
from tsdownsample import MinMaxLTTBDownsampler as DownSampler


def UpdateTreeDictProp(obj, prop):
    """
    Recursively update matplotlib artist properties from a nested dictionary.

    This function traverses nested property dictionaries and applies them to
    matplotlib artists. For each property, it checks if it's a direct setter
    on the object, and if not, recursively applies it to the appropriate sub-object.

    Parameters
    ----------
    obj : matplotlib artist object
        A matplotlib artist object to update.
    prop : dict
        Nested dictionary of properties where keys are property names
        and values are either property values or dicts for nested updates.

    Returns
    -------
    None
        Modifies the object in-place.
    """
    ains = ArtistInspector(obj)
    validp = ains.get_setters()
    for p in prop.keys():
        if p in validp:
            obj.set(**{p: prop[p]})
        else:
            obj2 = getattr(obj, 'get_' + p)()
            UpdateTreeDictProp(obj2, prop[p])


class SlotBase():
    """
    Base class for signal plotting slots.

    Provides common functionality for managing and plotting signals with time-based
    slicing, unit conversion, and axis property updates. Serves as the parent class
    for specific signal visualization slot types.
    """

    def CheckTime(self, Time):
        """
        Validate and normalize time window specifications.

        Ensures provided time bounds are valid and within the signal duration.
        If a single time value is provided, generates a window starting at that time.

        Parameters
        ----------
        Time : tuple or None
            (start_time, stop_time) specification. Can be None
            to use full signal range. Single-element sequences
            are converted to windows of size sampling_period.

        Returns
        -------
        tuple
            (Tstart, Tstop) - validated time window within signal bounds.
        """
        if Time is None:
            return (self.Signal.t_start, self.Signal.t_stop)

        if len(Time) == 1:
            Time = (Time[0], Time[0] + self.Signal.sampling_period)

        if Time[0] is None or Time[0] < self.Signal.t_start:
            Tstart = self.Signal.t_start
        else:
            Tstart = Time[0]

        if Time[1] is None or Time[1] > self.Signal.t_stop:
            Tstop = self.Signal.t_stop
        else:
            Tstop = Time[1]

        return (Tstart, Tstop)

    def GetSignal(self, Time, Units=None):
        """
        Extract and optionally rescale signal data for a specified time window.

        Retrieves the signal for the given time window and converts to the specified
        units if provided. Updates the stored units to match the returned signal.

        Parameters
        ----------
        Time : tuple or None
            Time window (start_time, stop_time). See CheckTime().
        Units : quantities.Quantity, optional
            Target units for rescaling. Uses self.units if None.

        Returns
        -------
        neo.AnalogSignal
            Signal data for the time window in requested units.
        """
        if Units is None:
            _Units = self.units
        else:
            _Units = Units
        Time = self.CheckTime(Time)
        sig = self.Signal.time_slice(Time[0], Time[1])
        if _Units is not None:
            sig = sig.rescale(_Units)
        self.units = sig.units
        return sig

    def UpdateAxKwargs(self, AxKwargs):
        """
        Update axis properties from a keyword arguments dictionary.

        Applies nested property updates to the matplotlib axis object using
        the UpdateTreeDictProp function for recursive property setting.

        Parameters
        ----------
        AxKwargs : dict
            Dictionary of axis properties to update, can be nested.

        Returns
        -------
        None
            Modifies the axis in-place.
        """
        # pass
        self.AxKwargs.update(AxKwargs)
        UpdateTreeDictProp(self.Ax, self.AxKwargs)

#
# class WavesColorSlot(SlotBase):
#     """
#     Visualization slot for color-mapped multi-channel waveforms.
#
#     Displays multi-channel signals as a 2D color-mapped image where channels
#     are arranged vertically and time runs horizontally. Useful for visualizing
#     high-density electrode array data.
#
#     Attributes:
#         DefImKwargs (dict): Default imshow() keyword arguments.
#         DefAxKwargs (dict): Default axis property settings.
#     """
#     DefImKwargs = {'cmap': 'viridis',
#                    'interpolation': 'none',
#                    }
#
#     DefAxKwargs = {'ylabel': 'Channels',
#                    'xaxis': {'visible': False,
#                              },
#                    'yaxis': {'visible': True,
#                              },
#                    }
#
#     def __init__(self, Signal, Units=None, Position=None,
#                  imKwargs=None, AxKwargs=None, Ax=None, MaxPoints=None):
#         """
#         Initialize a color-mapped waveform slot.
#
#         Parameters:
#             Signal (neo.AnalogSignal): Multi-channel signal to plot.
#             Units (quantities.Quantity, optional): Units for signal scaling.
#             Position (int, optional): Subplot position for this slot.
#             imKwargs (dict, optional): Additional imshow() keyword arguments.
#             AxKwargs (dict, optional): Additional axis property settings.
#             Ax (matplotlib.axes.Axes, optional): Axis to plot into. Creates new if None.
#             MaxPoints (int, optional): Maximum number of time points to display
#                                        (triggers downsampling if exceeded).
#         """
#
#         self.AxKwargs = self.DefAxKwargs.copy()
#         self.imKwargs = self.DefImKwargs.copy()
#         self.Position = Position
#         self.MaxPoints = MaxPoints
#
#         self.Signal = Signal
#         self.name = self.Signal.name
#         self.units = Units
#
#         self.Ax = Ax
#         if AxKwargs is not None:
#             self.AxKwargs.update(AxKwargs)
#         if imKwargs is not None:
#             self.imKwargs.update(imKwargs)
#
#         if self.Ax is not None:
#             UpdateTreeDictProp(self.Ax, self.AxKwargs)
#
#     def PlotSignal(self, Time, Units=None):
#         """
#         Plot multi-channel signal as a color-mapped image.
#
#         Renders the signal data as a 2D heatmap with optional downsampling
#         for performance optimization.
#
#         Args:
#             Time (tuple or None): Time window to plot. See CheckTime() for details.
#             Units (quantities.Quantity, optional): Units for signal rescaling.
#
#         Returns:
#             None: Updates self.img and self.current_time attributes.
#         """
#         sig = self.GetSignal(Time, Units)
#
#         if self.MaxPoints is not None:
#             sig = Spro.Resample(sig, MaxPoints=self.MaxPoints)
#
#         img = self.Ax.imshow(np.array(sig).astype(float).transpose(),
#                              aspect='auto',
#                              extent=(sig.t_start, sig.t_stop,
#                                      0, sig.shape[1]),
#                              **self.imKwargs,
#                              )
#         self.img = img
#         self.current_time = (sig.t_start.rescale('s'),
#                              sig.t_stop.rescale('s'))


class SpecSlot(SlotBase):
    """
    Visualization slot for spectrogram data.

    Displays time-frequency spectrograms as color-mapped images where time runs
    horizontally and frequency runs vertically. Supports automatic spectrogram
    computation from raw signals or direct display of pre-computed spectrograms.
    Useful for analyzing frequency content changes over time in electrophysiological
    data.

    Attributes:
        DefspecKwargs (dict): Default spectrogram computation parameters.
        DefAvgSpectKwargs (dict): Default averaged spectrogram parameters.
        DefAxKwargs (dict): Default axis property settings.
        DefImKwargs (dict): Default imshow() display parameters.
    """
    DefspecKwargs = {'Fmax': 100 * pq.Hz,
                     'Fmin': 0.5 * pq.Hz,
                     'Fres': 0.5 * pq.Hz,
                     'TimeRes': 0.1 * pq.s,
                     'Zscored': True,
                     }

    DefAvgSpectKwargs = {'SpecArgs': {'Fmax': 100 * pq.Hz,
                                      'Fmin': 0.5 * pq.Hz,
                                      'Fres': 0.5 * pq.Hz,
                                      'TimeRes': 0.1 * pq.s,
                                      'Zscored': False,
                                      },
                         'AvgSpectNorm': 'Zscore',
                         'AvgSpectNormTime': None,
                         }

    DefAxKwargs = {'ylabel': 'Freq [Hz]',
                   'xaxis': {'visible': False,
                             },
                   'yaxis': {'visible': True,
                             },
                   }

    DefImKwargs = {
        'norm': colors.Normalize(-3, 3),
        'cmap': 'seismic',
        'interpolation': 'bilinear',
    }

    def UpdateLineKwargs(self, LineKwargs):
        """
        Update line drawing properties.

        Currently not implemented.

        Parameters
        ----------
        LineKwargs : dict
            Line style parameters.

        Returns
        -------
        None
        """
        pass

    def __init__(self, Signal, Units=None, Position=None, imKwargs=None,
                 specKwargs=None, AxKwargs=None, Ax=None, MaxPoints=10000,
                 AvgSpectKwargs=None):
        """
        Initialize a spectrogram visualization slot.

        Parameters
        ----------
        Signal : neo.AnalogSignal
            Signal to compute spectrogram from.
        Units : quantities.Quantity, optional
            Units for signal rescaling.
        Position : int, optional
            Subplot position for this slot.
        imKwargs : dict, optional
            Additional imshow() keyword arguments.
        specKwargs : dict, optional
            Spectrogram computation parameters.
        AxKwargs : dict, optional
            Additional axis property settings.
        Ax : matplotlib.axes.Axes, optional
            Axis to plot into. Creates new if None.
        MaxPoints : int, optional
            Maximum time points for spectrogram display (default: 10000).
        AvgSpectKwargs : dict, optional
            Parameters for averaged spectrogram computation.
        """

        self.MaxPoints = MaxPoints
        self.specKwargs = self.DefspecKwargs.copy()
        self.AxKwargs = self.DefAxKwargs.copy()
        self.imKwargs = self.DefImKwargs.copy()
        self.AvgSpectKwargs = self.DefAvgSpectKwargs.copy()

        if specKwargs is not None:
            self.specKwargs.update(specKwargs)

        self.Signal = Signal
        self.name = self.Signal.name

        self.units = Units
        self.Position = Position

        self.Ax = Ax
        if AxKwargs is not None:
            self.AxKwargs.update(AxKwargs)
        if imKwargs is not None:
            self.imKwargs.update(imKwargs)
        if AvgSpectKwargs is not None:
            self.AvgSpectKwargs.update(AvgSpectKwargs)
        if self.Ax is not None:
            UpdateTreeDictProp(self.Ax, self.AxKwargs)

    def PlotSignal(self, Time, Units=None):
        """
        Plot spectrogram for the specified time window.

        Computes or retrieves spectrogram data and displays it as a color-mapped
        image with time on x-axis and frequency on y-axis. Supports automatic
        spectrogram computation from raw signals or direct display of pre-computed
        spectrograms with optional downsampling.

        Parameters
        ----------
        Time : tuple or None
            Time window to plot. See CheckTime() for details.
        Units : quantities.Quantity, optional
            Units for signal rescaling.

        Returns
        -------
        None
            Updates self.img and self.current_time attributes.
        """
        sig = self.GetSignal(Time, Units)

        if 'spec' in sig.annotations:
            spec = sig
        else:
            spec = Spro.Spectrogram(sig, **self.specKwargs)

        f = spec.annotations['Freq']
        data = Spro.Resample(spec, MaxPoints=self.MaxPoints)
        data_plain = np.asarray(data.magnitude).transpose()
        img = self.Ax.imshow(data_plain,
                             origin='lower',
                             aspect='auto',
                             extent=(spec.t_start.magnitude, spec.t_stop.magnitude,
                                     np.min(f.magnitude), np.max(f.magnitude)),
                             **self.imKwargs,
                             )
        self.img = img
        self.current_time = (sig.t_start.rescale('s'),
                             sig.t_stop.rescale('s'))

    def CalcAvarage(self, TimeAvg, TimesEvent, Units=None, **Kwargs):
        """
        Calculate and plot averaged spectrogram over triggered events.

        Extracts signal windows around trigger times, computes spectrograms
        for each trial, and averages them to produce a single representative
        spectrogram with optional normalization.

        Parameters
        ----------
        TimeAvg : tuple
            Time window relative to trigger events (before, after).
        TimesEvent : array-like
            Times of trigger events.
        Units : quantities.Quantity, optional
            Units for signal rescaling.
        **Kwargs
            Additional keyword arguments passed to averaging function.

        Returns
        -------
        neo.AnalogSignal
            The averaged spectrogram data.
        """

        sig = self.GetSignal((None, None), Units)

        spect = Spro.AvgSpectrogram(sig,
                                    TimesEvent=TimesEvent,
                                    TimeAvg=TimeAvg,
                                    **self.AvgSpectKwargs)

        f = spect.annotations['Freq']
        img = self.Ax.imshow(np.array(spect).transpose(),
                             origin='lower',
                             aspect='auto',
                             extent=(spect.t_start, spect.t_stop,
                                     np.min(f), np.max(f)),
                             **self.imKwargs,
                             )
        self.img = img
        return spect


class SpikeSlot(SlotBase):
    """
    Visualization slot for spike train data.

    Displays neural spike times as vertical lines on a time axis. Useful for
    visualizing action potentials, threshold crossings, and other discrete
    event times in electrophysiology data.

    Attributes:
        DefLineKwargs (dict): Default vertical line style parameters.
        DefAxKwargs (dict): Default axis property settings.
    """
    DefLineKwargs = {'color': 'r',
                     'linestyle': '-.',
                     'alpha': 0.5,
                     'linewidth': 0.5,
                     }
    DefAxKwargs = {}

    def UpdateLineKwargs(self, LineKwargs):
        """
        Update spike line drawing properties.

        Parameters
        ----------
        LineKwargs : dict
            Line style parameters (color, linestyle, linewidth, etc.).

        Returns
        -------
        None
            Modifies the plotted lines in-place.
        """
        self.LineKwargs.update(LineKwargs)
        UpdateTreeDictProp(self.Line, self.LineKwargs)

    def __init__(self, Signal, Units='s',
                 Position=None, Ax=None, AxKwargs=None,
                 **LineKwargs):
        """
        Initialize a spike train visualization slot.

        Parameters
        ----------
        Signal : neo.SpikeTrain
            Spike times to visualize.
        Units : str, optional
            Units for spike time display (default: 's').
        Position : int, optional
            Subplot position for this slot.
        Ax : matplotlib.axes.Axes, optional
            Axis to plot into. Creates new if None.
        AxKwargs : dict, optional
            Additional axis property settings.
        **LineKwargs
            Additional vertical line style keyword arguments.
        """

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

    def PlotSignal(self, Time, Units=None):
        """
        Plot spike train as vertical lines.

        Displays spike times as vertical lines spanning the current y-axis limits,
        allowing visualization of when spikes occurred in relation to other signals.

        Parameters
        ----------
        Time : tuple or None
            Time window to plot. See CheckTime() for details.
        Units : quantities.Quantity, optional
            Units for spike time display.

        Returns
        -------
        None
            Updates self.Lines attribute.
        """
        if self.Ax is None:
            self.Fig, self.Ax = plt.subplots()

        sig = self.GetSignal(Time, Units)

        xmin, xmax, ymin, ymax = self.Ax.axis()

        self.Lines = self.Ax.vlines(sig,
                                    ymin, ymax,
                                    **self.LineKwargs
                                    )


class WaveSlot(SlotBase):
    """
    Visualization slot for continuous time-series waveforms.

    Displays raw analog signals as line plots with support for trial-by-trial
    visualization, averaged waveforms, standard deviation bands, and automatic
    downsampling for large datasets. The primary class for visualizing continuous
    electrophysiology recordings.

    Attributes:
        DefTrialLineKwargs (dict): Default line style for individual trial plotting.
        DefLineKwargs (dict): Default line style for main signal plotting.
        DefAxKwargs (dict): Default axis property settings.
        MaxPlotPoints (int): Maximum data points to plot before downsampling (100000).
    """
    DefTrialLineKwargs = {'color': 'k',
                          'linestyle': '-',
                          'alpha': 0.05,
                          'linewidth': 0.5,
                          'clip_on': True,
                          }

    DefLineKwargs = {'color': 'k',
                     'linestyle': '-',
                     'alpha': 1,
                     'linewidth': 0.5,
                     'clip_on': True,
                     }
    DefAxKwargs = {}
    MaxPlotPoints = int(1e5)

    def UpdateLineKwargs(self, LineKwargs):
        """
        Update waveform line drawing properties.

        Parameters
        ----------
        LineKwargs : dict
            Line style parameters (color, linestyle, linewidth, etc.).

        Returns
        -------
        None
            Modifies the plotted line in-place.
        """
        self.LineKwargs.update(LineKwargs)
        UpdateTreeDictProp(self.Line, self.LineKwargs)

    def __init__(self, Signal, Units=None, UnitsInLabel=False,
                 Position=None, Ax=None, AxKwargs=None, TrialProcessChain=None,
                 DownSampling=True, Sampler=None,
                 **LineKwargs):
        """
        Initialize a time-series waveform visualization slot.

        Parameters
        ----------
        Signal : neo.AnalogSignal
            Continuous signal to plot.
        Units : quantities.Quantity, optional
            Units for signal scaling.
        UnitsInLabel : bool, optional
            If True, include units in plot label (default: False).
        Position : int, optional
            Subplot position for this slot.
        Ax : matplotlib.axes.Axes, optional
            Axis to plot into. Creates new if None.
        AxKwargs : dict, optional
            Additional axis property settings.
        TrialProcessChain : list, optional
            Processing functions to apply to trials during averaged waveform calculation.
        DownSampling : bool, optional
            Enable automatic downsampling (default: True).
        Sampler : object, optional
            Custom downsampler instance. Uses MinMaxLTTBDownsampler if None.
        **LineKwargs
            Additional line style keyword arguments.
        """

        self.DownSampling = DownSampling
        if Sampler is None:
            self.DownSampler = DownSampler()
        else:
            self.DownSampler = Sampler

        self.TrialLineKwargs = self.DefTrialLineKwargs.copy()
        self.LineKwargs = self.DefLineKwargs.copy()
        self.AxKwargs = self.DefAxKwargs.copy()

        self.TrialProcessChain = TrialProcessChain

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

        self.TrialLineKwargs['color'] = self.LineKwargs['color']

        self.current_time = None

    def PlotSignal(self, Time, Units=None):
        """
        Plot signal waveform as a line plot.

        Renders the signal with automatic downsampling if needed, formats labels
        with units if requested, and updates axis limits to match the signal window.

        Parameters
        ----------
        Time : tuple or None
            Time window to plot. See CheckTime() for details.
        Units : quantities.Quantity, optional
            Units for signal rescaling.

        Returns
        -------
        None
            Updates self.Lines, self.current_time, and axis properties.
        """
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
        """
        Internal method to render signal data on the axis.

        Handles the actual matplotlib plotting with automatic downsampling for
        large datasets using the configured sampler. Updates current_time and
        axis limits to match the signal window.

        Parameters
        ----------
        sig : neo.AnalogSignal
            Signal data to plot.

        Returns
        -------
        None
            Updates self.Lines, self.current_time, and axis limits.
        """
        if sig.size > self.MaxPlotPoints:
            idx = self.DownSampler.downsample(np.asarray(sig)[:, 0], n_out=self.MaxPlotPoints)
            self.Lines = self.Ax.plot(sig.times[idx].rescale('s'),
                                      np.asarray(sig)[idx, 0],
                                      **self.LineKwargs
                                      )
        else:
            self.Lines = self.Ax.plot(sig.times.rescale('s'),
                                      sig,
                                      **self.LineKwargs
                                      )
        self.current_time = (sig.t_start.rescale('s'),
                             sig.t_stop.rescale('s'))
        self.Ax.set_xlim(left=self.current_time[0],
                         right=self.current_time[1])
        self.Line = self.Lines[0]

    def CalcAvarage(self, TimeAvg, TimesEvent, Units=None,
                    PlotMean=True, PlotStd=False, PlotTrials=False,
                    TrialLineKwargs=None, StdAlpha=0.2, TrialProcessChain=None,
                    **kwargs):
        """
        Calculate and plot averaged waveform over triggered events.

        Extracts signal windows around trigger times, averages them to compute
        the mean response, and optionally displays individual trials and/or
        standard deviation as filled regions. Supports preprocessing of individual
        trials before averaging through a processing chain.

        Parameters
        ----------
        TimeAvg : tuple
            Time window relative to trigger events (before, after).
        TimesEvent : array-like
            Times of trigger events.
        Units : quantities.Quantity, optional
            Units for signal rescaling.
        PlotMean : bool, optional
            Plot the averaged waveform (default: True).
        PlotStd : bool, optional
            Show standard deviation as shaded region (default: False).
        PlotTrials : bool, optional
            Overlay individual trial traces (default: False).
        TrialLineKwargs : dict, optional
            Line style for trial plotting.
        StdAlpha : float, optional
            Transparency of std dev region (default: 0.2).
        TrialProcessChain : list, optional
            Processing functions for trials before averaging.
        **kwargs
            Additional keyword arguments.

        Returns
        -------
        neo.AnalogSignal
            The averaged waveform with std annotation.
        """

        if TrialLineKwargs is not None:
            self.TrialLineKwargs.update(TrialLineKwargs)

        if TrialProcessChain is not None:
            self.TrialProcessChain = TrialProcessChain

        sig = self.GetSignal((None, None), Units)
        avg = Spro.TrigAveraging(sig,
                                 TimesEvent=TimesEvent,
                                 TimeAvg=TimeAvg,
                                 TrialProcessChain=self.TrialProcessChain)

        if PlotMean:
            self._PlotSignal(avg)

        if PlotTrials:
            acc = avg.annotations['acc']
            self.Ax.plot(acc.times,
                         acc,
                         **self.TrialLineKwargs)

        if PlotStd:
            std = avg.annotations['std']
            self.Ax.fill_between(std.times,
                                 np.array(avg + std).flatten(),
                                 np.array(avg - std).flatten(),
                                 alpha=StdAlpha,
                                 facecolor=self.LineKwargs['color'],
                                 edgecolor=None,
                                 clip_on=False)
        return avg


class ImgSlot(SlotBase):
    """
    Visualization slot for 2D image data (e.g., spatial imaging maps).

    Displays 2D image arrays as color-mapped heatmaps. Useful for visualizing
    spatially-resolved measurements such as voltage maps from imaging systems
    or electrode array maps.

    Attributes:
        DefAxKwargs (dict): Default axis property settings.
        DefImKwargs (dict): Default imshow() display parameters.
    """
    DefAxKwargs = {}
    DefImKwargs = {
        'vmin': -10,
        'vmax': 10,
        'cmap': 'seismic',
        'interpolation': 'bicubic',
    }

    def __init__(self, Signal, Ax=None,
                 AxKwargs=None, Units=None, imKwargs=None, ):
        """
        Initialize a 2D image data visualization slot.

        Parameters
        ----------
        Signal : array-like
            2D/3D array-like data where first dimension is time/frames.
        Ax : matplotlib.axes.Axes, optional
            Axis to plot into.
        AxKwargs : dict, optional
            Additional axis property settings.
        Units : quantities.Quantity, optional
            Units for data scaling.
        imKwargs : dict, optional
            Additional imshow() keyword arguments.
        """
        self.Signal = Signal
        self.Ax = Ax

        self.Map = True

        self.current_time = None
        self.units = Units

        self.imKwargs = self.DefImKwargs.copy()
        if imKwargs is not None:
            self.imKwargs.update(imKwargs)

        self.Img = self.Ax.imshow(np.array(self.Signal[0, :, :]),
                                  **self.imKwargs,
                                   )

    def PlotSignal(self, Time, Units=None):
        """
        Update 2D image display for a new time frame.

        Updates the displayed image with data from the specified time window,
        refreshes the title with frame information, and triggers canvas redraw.

        Parameters
        ----------
        Time : tuple or None
            Time window to plot. See CheckTime() for details.
        Units : quantities.Quantity, optional
            Units for data rescaling.

        Returns
        -------
        None
            Updates self.Img and redraws the canvas.
        """
        sig = self.GetSignal(Time, Units)

        self.Img.set_array(np.array(sig[0, :, :]))
        self.Ax.set_title('{0}\n{1:.3f}'.format(sig.name, sig.t_start))
        self.current_time = (sig.t_start.rescale('s'),
                             sig.t_stop.rescale('s'))
        self.Ax.figure.canvas.draw()

    def UpdateAxKwargs(self, AxKwargs):
        """
        Update axis properties (not implemented for ImgSlot).

        Parameters
        ----------
        AxKwargs : dict
            Axis property parameters (unused).

        Returns
        -------
        None
            Method not implemented for this slot type.
        """
        pass


class ControlFigure():
    """
    Interactive control panel for visualization navigation and animation.

    Provides sliders and controls for navigating through signal data, setting
    time windows, and triggering animations of map sequences. Integrates with
    PlotSlots to update visualizations in real-time.

    Attributes:
        TimeLineKwargs (dict): Default style for time indicator lines during animation.
    """
    TimeLineKwargs = {'color': 'g',
                      'linewidth': 2,
                      'linestyle': '-.',
                      }

    def __init__(self, pltSL, AxsAnimationLines=None, figsize=(20 * 0.394, 5 * 0.394)):
        """
        Initialize interactive control panel for visualization.

        Creates sliders and text input controls for navigating time windows,
        setting animation parameters, and controlling visualization updates.

        Parameters
        ----------
            pltSL (PlotSlots): PlotSlots instance to control.
            AxsAnimationLines (list, optional): Axes where to draw time indicator lines
                                               during animation.
            figsize (tuple, optional): Figure size in inches (default: (7.87, 1.97)).

        """
        self.pltSL = pltSL

        self.pltSLFigs = set([ax.get_figure() for ax in pltSL.Axs])

        self.MapSlots = []
        for sl in pltSL.Slots:
            if hasattr(sl, 'Map'):
                self.MapSlots.append(sl)

        TMax = []
        TMin = []
        for sl in pltSL.Slots:
            if not hasattr(sl, 'Map'):
                TMax.append(sl.Signal.t_stop.rescale('s'))
                TMin.append(sl.Signal.t_start.rescale('s'))

        TMax = np.max(TMax)
        TMin = np.min(TMin)

        self.Fig, ax = plt.subplots(10, 1, figsize=figsize)
        self.sTstart = Slider(ax[0],
                              label='TStart [s]',
                              valmax=TMax,
                              valmin=TMin,
                              valinit=TMin)

        self.sTshow = Slider(ax[1],
                             label='TShow [s]',
                             valmax=TMax - TMin,
                             valmin=0,
                             valinit=(TMax - TMin) / 10)

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

        # self.bStart = Button(ax[4],
        #                      label='Start')
        # self.bStart.on_clicked(self.StartAnimation)
        # self.Timer = None
        # self.TextInterval = TextBox(ax[5],
        #                             'Interval [ms]',
        #                             initial='2000')

        self.bStartMapAni = Button(ax[6],
                                   label='Annimate Maps START')
        self.TextMapTimeInterval = TextBox(ax[7],
                                           'Frame Interval [Seg]',
                                           initial='0.01')
        self.TextMapRefreshTime = TextBox(ax[8],
                                          'Frame Refresh [Seg]',
                                          initial='0.2')
        self.bStartMapAni.on_clicked(self.StartMapAnimation)
        self.TimerMap = None

        self.bStartSetZero = Button(ax[9],
                                    label='Set Zero at Start Time')
        self.bStartSetZero.on_clicked(self.BtSetZero)

        self.AxsAnimationLines = AxsAnimationLines
        self.TimeLines = None

    def BtSetZero(self, val):
        """
        Set zero baseline at the start time window.

        Applies zero-baseline correction to all signals marked with 'LiveZero'
        annotation using the mean value over the start time window.

        Parameters
        ----------
        val
            Event value (unused, callback parameter).

        Returns
        -------
        None
            Modifies signals in-place and triggers visualization update.
        """
        Twind = (self.sTstart.val * pq.s,
                 (self.sTstart.val + 5) * pq.s)

        for sl in self.pltSL.Slots:
            if 'LiveZero' in sl.Signal.annotations:
                if sl.Signal.annotations['LiveZero']:
                    # print(sl.Signal.name)
                    sl.Signal = Spro.SetZero(sl.Signal, TWind=Twind)
        self.Update(None)

    def StartMapAnimation(self, val):
        """
        Start/stop animation of 2D map sequences.

        Toggles animation mode where map visualizations cycle through a time
        window in discrete steps, with optional time indicator lines on other axes.

        Parameters
        ----------
        val
            Event value (unused, callback parameter).

        Returns
        -------
        None
            Creates or stops animation timer and updates button label.
        """
        if self.TimerMap is not None:
            self.bStartMapAni.label.set_label('Annimate Maps START')
            self.TimerMap.stop()
            self.TimerMap = None
            return

        twind = (self.sTstart.val * pq.s,
                 self.sTstart.val * pq.s + self.sTshow.val * pq.s)

        interval = float(self.TextMapTimeInterval.text) * pq.s
        refresh = float(self.TextMapRefreshTime.text) * pq.s
        self.MapCount = 0
        self.MapTimes = np.arange(twind[0], twind[1], interval) * pq.s
        # print(interval)
        # print(self.MapTimes)

        if self.AxsAnimationLines is not None:
            self.TimeLines = []
            for ax in self.AxsAnimationLines:
                ymin, ymax = ax.get_ylim()
                self.TimeLines.append(
                    ax.plot([self.MapTimes[0], self.MapTimes[0]], [ymin, ymax], **self.TimeLineKwargs)[0])
        else:
            self.TimeLines = None

        self.bStartMapAni.label.set_label('Annimate Maps STOP')
        self.TimerMap = self.Fig.canvas.new_timer(
            interval=refresh.rescale('ms'))
        self.TimerMap.add_callback(self.UpdateMapAnimation)
        print('Start', interval)
        self.OldTime = datetime.datetime.now()
        self.Fig.canvas.draw()
        self.TimerMap.start()

    def UpdateMapAnimation(self):
        """
        Update map display for the next frame in animation sequence.

        Advances the animation counter, updates map visualizations, and moves
        time indicator lines on synchronized axes. Loops back to start when
        the end of the time window is reached.

        Returns
        -------
        None
            Updates map displays and time lines.
        """
        t = datetime.datetime.now()
        # print(t-self.OldTime)
        # print(self.MapCount, self.MapTimes[self.MapCount])
        self.OldTime = t
        for sl in self.MapSlots:
            sl.PlotSignal((self.MapTimes[self.MapCount],
                           self.MapTimes[self.MapCount] + 1 * pq.s))

        if self.TimeLines is not None:
            for tline in self.TimeLines:
                tline.set_xdata([self.MapTimes[self.MapCount],
                                 self.MapTimes[self.MapCount]])

        if self.MapCount >= (len(self.MapTimes) - 1):
            self.MapCount = 0
        else:
            self.MapCount += 1

    def StartAnimation(self, val):
        """
        Start/stop continuous animation through signal time.

        Toggles automatic scrolling mode where the display continuously advances
        through the signal, useful for scanning through large datasets.

        Parameters
        ----------
        val
            Event value (unused, callback parameter).

        Returns
        -------
        None
            Creates or stops animation timer and updates button label.
        """
        if self.Timer is not None:
            self.bStart.label.set_label('Start')
            self.Timer.stop()
            self.Timer = None
            if self.TimeLines is not None:
                for l in self.TimeLines:
                    del l
                self.TimeLines = None
            return

        try:
            interval = float(self.TextInterval.text)
        except:
            return

        self.Timer = self.Fig.canvas.new_timer(interval=interval)
        self.Timer.add_callback(self.UpdateAnimation)
        self.Timer.start()
        self.bStart.label.set_label('Stop')

    def UpdateAnimation(self):
        """
        Advance continuous animation by one time step.

        Called by the animation timer to update the display window position
        during continuous animation playback.

        Returns
        -------
        None
            Updates time slider positions.
        """
        self.sTstart.set_val(self.sTstart.val + self.sTshow.val / 2)

    def Update(self, val):
        """
        Update visualization in response to slider changes.

        Called when time sliders are adjusted to update the displayed time window.
        Refreshes all plot slots if refresh is enabled.

        Parameters
        ----------
        val
            Event value (unused, callback parameter).

        Returns
        -------
        None
            Triggers visualization update.
        """
        twind = (self.sTstart.val * pq.s,
                 self.sTstart.val * pq.s + self.sTshow.val * pq.s)

        if self.Refresh:
            self.UpdateGraph(twind)

    def UpdateGraph(self, twind):
        """
        Render plots for the specified time window.

        Updates all plot slots with data from the given time window and refreshes
        the matplotlib canvases.

        Parameters
        ----------
        twind : tuple
            Time window (start_time, stop_time).

        Returns
        -------
        None
            Redraws all figures.
        """
        self.pltSL.PlotChannels(twind)
        for f in self.pltSLFigs:
            f.canvas.draw()

    def submit_start(self, text):
        """
        Handle start time text box submission.

        Validates the entered start time, ensures it doesn't exceed stop time,
        and updates the time sliders accordingly.

        Parameters
        ----------
        text : str
            Text value entered in start time box.

        Returns
        -------
        None
            Updates time sliders and internal state.
        """
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
        """
        Handle stop time text box submission.

        Validates the entered stop time, ensures it's greater than start time,
        and updates the time window sliders accordingly.

        Parameters
        ----------
        text : str
            Text value entered in stop time box.

        Returns
        -------
        None
            Updates time sliders and internal state.
        """
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
        """
        Set time window controls to specified range.

        Updates text boxes and sliders to reflect a specific time window,
        with refresh disabled during updates to avoid multiple redraws.

        Parameters
        ----------
        twind : tuple
            Time window (start_time, stop_time).

        Returns
        -------
        None
            Updates control positions.
        """
        self.Refresh = False
        self.TextStop.set_val(str(np.array(twind[1])))
        self.TextStart.set_val(str(np.array(twind[0])))
        self.TextStop.set_val(str(np.array(twind[1])))
        self.Refresh = True


class PlotSlots():
    """
    Main container for managing multiple plot slots and figure layout.

    Orchestrates the creation of matplotlib figures and axes, manages multiple
    visualization slots (WaveSlot, SpecSlot, etc.), and provides methods for
    batch operations like plotting, event overlay, and formatting. Supports
    interactive controls and animation capabilities.

    Attributes:
        ScaleBarKwargs (dict): Default parameters for scale bar overlays.
        RcGeneralParams (dict): Default matplotlib rc parameter overrides.
        FigKwargs (dict): Default figure-level property settings.
        gridspec_Kwargs (dict): Default gridspec layout parameters.
        TimeAxisProp (dict): Default properties for time axis display.
        LegendKwargs (dict): Default legend styling parameters.
    """
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
        """
        Update figure-level properties.

        Parameters
        ----------
        FigKwargs : dict
            Figure property updates to apply.

        Returns
        -------
        None
            Modifies self.Fig in-place.
        """
        self.FigKwargs.update(FigKwargs)
        UpdateTreeDictProp(self.Fig, self.FigKwargs)

    def _GenerateFigure(self):
        """
        Generate matplotlib figure and axes for all slots.

        Creates a subplot layout based on slot positions and assigns axes
        to each slot. Handles both single and multiple-row layouts.

        Returns
        -------
        None
            Populates self.Fig, self.Axs, and self.CAxs.
        """

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
                 AxKwargs=None, TimeAxis=-1,
                 ScaleBarAx=None, LiveControl=False, AxsAnimationLines=None):
        """
        Initialize a PlotSlots container.

        Sets up figure and axes layout, applies styling, enables optional
        interactive controls, and configures signal visualization slots.

        Parameters
        ----------
        Slots : list
            List of slot objects (WaveSlot, SpecSlot, etc.).
        Fig : matplotlib.figure.Figure, optional
            Existing figure to use. Creates new if None.
        FigKwargs : dict, optional
            Figure-level property overrides.
        RcGeneralParams : dict, optional
            Matplotlib rc parameter overrides.
        AxKwargs : dict, optional
            Common axis property updates for all slots.
        TimeAxis : int or iterable, optional
            Which subplot(s) show the time axis (default: -1, last subplot).
        ScaleBarAx : int, optional
            Subplot index for scale bar overlay.
        LiveControl : bool, optional
            Enable interactive control panel (default: False).
        AxsAnimationLines : list, optional
            Axes for animation time indicator lines.
        """

        if RcGeneralParams is not None:
            self.RcGeneralParams.update(RcGeneralParams)
        plt.rcParams.update(self.RcGeneralParams)

        if FigKwargs is not None:
            self.FigKwargs.update(FigKwargs)

        self.Slots = Slots

        self.ScaleBarAx = ScaleBarAx

        if Fig is None:
            self._GenerateFigure()
        else:
            self.Fig = Fig
            self.Axs = []
            for sl in self.Slots:
                self.Axs.append(sl.Ax)

        for sl in self.Slots:
            if AxKwargs is not None:
                sl.UpdateAxKwargs(AxKwargs)

        self.TimeAxis = TimeAxis
        if self.TimeAxis is not None:
            if hasattr(TimeAxis, '__iter__'):
                for ti in TimeAxis:
                    sl = self.Slots[TimeAxis]
                    sl.UpdateAxKwargs(self.TimeAxisProp)
            else:
                sl = self.Slots[TimeAxis]
                sl.UpdateAxKwargs(self.TimeAxisProp)

        UpdateTreeDictProp(self.Fig, self.FigKwargs)
        self.SortSlotsAx()

        if LiveControl:
            self.CtrFig = ControlFigure(self, AxsAnimationLines)
        else:
            self.CtrFig = None

    def SortSlotsAx(self):
        """
        Create mapping of axes to their associated slots.

        Builds a dictionary associating each subplot axis with all slots
        that share that axis, useful for batch operations on related slots.

        Returns
        -------
        None
            Populates self.SlotsInAxs dictionary.
        """
        self.SlotsInAxs = {}
        for ax in self.Axs:
            sll = []
            for sl in self.Slots:
                if sl.Ax == ax:
                    sll.append(sl)
            self.SlotsInAxs.update({ax: sll})

    def ClearAxes(self):
        """
        Remove all plotted lines from all axes.

        Clears previous plots to prepare for new visualization, useful before
        updating displays with new time windows or data.

        Returns
        -------
        None
            Removes all line objects from all axes.
        """
        for sl in self.Slots:
            while sl.Ax.lines:
                sl.Ax.lines[0].remove()

    def FormatFigure(self):
        """
        Apply formatting to the figure (e.g., scale bars, colorbars).

        Adds decorative and informational elements like scale bars to specified
        axes. Scale bar unit is automatically derived from signal units if not
        explicitly specified.

        Returns
        -------
        None
            Modifies figure in-place.
        """

        if self.ScaleBarAx is not None:
            if self.ScaleBarKwargs['yunit'] is None:
                sl = self.SlotsInAxs[self.Axs[self.ScaleBarAx]][0]
                su = str(sl.units).split(' ')[-1]
                self.ScaleBarKwargs['yunit'] = su
            DrawBarScale(self.Axs[self.ScaleBarAx], **self.ScaleBarKwargs)

    def AddLegend(self, **LegendKwargs):
        """
        Add legend to all subplots.

        Updates legend style settings and applies them to all axes, displaying
        labels for all plotted signals.

        Parameters
        ----------
        **LegendKwargs
            Legend configuration parameters (fontsize, ncol, loc, etc.).

        Returns
        -------
        None
            Adds legends to all axes.
        """
        self.LegendKwargs.update(LegendKwargs)
        for Ax in self.Axs:
            Ax.legend(**self.LegendKwargs)

    def PlotChannels(self, Time, Units=None, FormatFigure=True):
        """
        Plot all signals for the specified time window.

        Updates all slot visualizations with data from the given time window,
        refreshes image colorbars, and updates the internal time tracking.

        Parameters
        ----------
        Time : tuple or None
            Time window to plot. See CheckTime() for details.
        Units : quantities.Quantity, optional
            Units for signal rescaling.
        FormatFigure : bool, optional
            Apply figure formatting (default: True).

        Returns
        -------
        None
            Updates all slot visualizations and refreshes canvases.
        """
        self.ClearAxes()
        print('plot channels')
        SpectColBars.ImgDicts = {}
        for isl, sl in enumerate(self.Slots):
            sl.PlotSignal(Time, Units=Units)
            if hasattr(sl, 'img'):
                SpectColBars.ImgDicts.update({isl: sl.img})
            #        if Time is not None:
            #            if Time[0] is not None:
            #                sl.Ax.set_xlim(left=Time[0].magnitude)
            #            if Time[1] is not None:
            #                sl.Ax.set_xlim(right=Time[1].magnitude)
            if sl.current_time is not None:
                if not hasattr(sl, 'Map'):
                    self.current_time = sl.current_time

        if self.CtrFig is not None:
            self.CtrFig.SetTimes(self.current_time)

    def PlotEvents(self, Times, Labels=None, lAx=0, fontsize='xx-small',
                   LabPosition='top', duration=None, **kwargs):
        """
        Overlay event markers on all subplots.

        Draws vertical lines (or shaded regions for duration-based events) at
        specified times across all axes, with optional text labels at the top
        or bottom of a selected axis.

        Parameters
        ----------
        Times : array-like
            Times of events to mark.
        Labels : array-like, optional
            Text labels for events.
        lAx : int, optional
            Axis index for placing text labels (default: 0).
        fontsize : str, optional
            Font size for text labels (default: 'xx-small').
        LabPosition : str, optional
            Vertical position of labels: 'top' or 'bottom' (default: 'top').
        duration : quantities.Quantity, optional
            Event duration for shading regions. If None, draws vertical lines.
        **kwargs
            Additional matplotlib vlines/vspan keyword arguments.

        Returns
        -------
        None
            Adds event markers and labels to plot.
        """

        xlim = self.Axs[0].get_xlim()

        Times = Times.rescale('s')

        self.Texts = []

        if Labels is not None:
            for ilbl, lbl in enumerate(Labels):
                for ax in self.Axs:
                    ylim = ax.get_ylim()
                    if duration is not None:
                        ax.vspan(Times[ilbl], Times[ilbl] +
                                 duration, ylim[0], ylim[1], **kwargs)
                    else:
                        ax.vlines(Times[ilbl], ylim[0], ylim[1], **kwargs)
                lax = self.Axs[lAx]
                if LabPosition == 'top':
                    ylim = lax.get_ylim()[1]
                else:
                    ylim = lax.get_ylim()[0]
                txt = lax.text(Times[ilbl], ylim, lbl, fontsize=fontsize)
                self.Texts.append(txt)
            return

        # EventLines = []
        for ax in self.Axs:
            ylim = ax.get_ylim()
            if duration is not None:
                lines = ax.axvspan(Times, Times + duration,
                                   ylim[0], ylim[1], **kwargs)
            else:
                lines = ax.vlines(Times, ylim[0], ylim[1], **kwargs)
        #            EventLines.append(lines[0])
        # return EventLines
        self.Axs[0].set_xlim(xlim)

    def PlotEventAvarage(self, TimeAvg, TimesEvent, Units=None, ClearAxes=True,
                         **Avgkwargs):
        """
        Plot averaged signals aligned to trigger events.

        Extracts signal windows around each event time, computes averages
        and standard deviations, and displays the results. Supports visualization
        of event-locked responses across all signal slots.

        Parameters
        ----------
        TimeAvg : tuple
            Time window relative to events (before, after).
        TimesEvent : array-like
            Times of trigger events.
        Units : quantities.Quantity, optional
            Units for signal rescaling.
        ClearAxes : bool, optional
            Clear previous plots before drawing (default: True).
        **Avgkwargs
            Additional averaging parameters (PlotMean, PlotStd, PlotTrials, etc.).

        Returns
        -------
        list
            Averaged signal objects, one per slot.
        """

        if ClearAxes:
            self.ClearAxes()

        MeanSigs = []
        for isl, sl in enumerate(self.Slots):
            print('Calculating Avg ', sl.name, isl)
            MeanSig = sl.CalcAvarage(TimeAvg, TimesEvent, Units=Units,
                                     **Avgkwargs)
            MeanSigs.append(MeanSig)
            if hasattr(sl, 'img'):
                SpectColBars.ImgDicts.update({isl: sl.img})

        sl.Ax.set_xlim(left=TimeAvg[0].magnitude)
        sl.Ax.set_xlim(right=TimeAvg[1].magnitude)

        self.FormatFigure()
        return MeanSigs
