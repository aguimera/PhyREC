#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
SignalProcess module - Signal processing and analysis tools for electrophysiology data.

This module provides comprehensive signal processing functions for electrophysiological
recordings including filtering, resampling, spectral analysis, averaging, and statistical
measures. Supports both univariate (single channel) and multivariate (multichannel) analysis.

Functions cover:
    - Basic operations: filtering, resampling, downsampling, detrending
    - Spectral analysis: spectrograms, averaged spectrograms
    - Data transformation: derivative, absolute value, power, z-score
    - Triggering and averaging: event-locked averaging, triggered spike detection
    - Correlation analysis: cross-correlation, Pearson correlation with sliding windows
    - Statistical measures: RMS, power, sliding window statistics

Created on Wed Apr 11 10:27:23 2018

@author: aguimera
"""
import numpy as np
from scipy import signal
from fractions import Fraction
from neo.core import AnalogSignal, SpikeTrain
import PhyREC.SignalAnalysis as Ran
from PhyREC.ImageSequence import ImageSequence
import quantities as pq
import matplotlib.mlab as mlab
import scipy.stats as stats
from scipy import interpolate
import sys
from scipy.stats.mstats import zscore
from PhyREC import DbgFplt
from scipy.interpolate import UnivariateSpline
from scipy.signal import medfilt
from numpy.lib.stride_tricks import sliding_window_view


def ApplyProcessChain(sig, ProcessChain):
    """
    Apply a sequence of processing functions to a signal.

    Applies a chain of processing operations (filtering, downsampling, etc.)
    to a signal in sequence. Useful for building complex processing pipelines.

    Parameters
    ----------
    sig : neo.AnalogSignal
        Input signal to process.
    ProcessChain : list or None
        List of dicts with 'function' and 'args' keys,
        or None to skip processing.

    Returns
    -------
    neo.AnalogSignal
        Processed signal after applying all operations in chain.
    """
    if ProcessChain is None:
        return sig

    sl = sig.copy()
    for Proc in ProcessChain:
        sl = Proc['function'](sl, **Proc['args'])

    return sl


def Spectrogram(sig, Fres=2 * pq.Hz, TimeRes=0.01 * pq.s,
                Fmin=1 * pq.Hz, Fmax=200 * pq.Hz, Zscored=True, NormTime=None,
                dtype=float,
                **specKwarg):
    """
    Compute time-frequency spectrogram of a signal.

    Calculates a spectrogram using short-time Fourier transform with configurable
    frequency and time resolution. Supports z-score normalization and time-windowed
    normalization for baseline subtraction.

    Parameters
    ----------
    sig : neo.AnalogSignal
        Input signal to analyze.
    Fres : quantities.Quantity, optional
        Frequency resolution (default: 2 Hz).
    TimeRes : quantities.Quantity, optional
        Time resolution window (default: 0.01 s).
    Fmin : quantities.Quantity, optional
        Minimum frequency to include (default: 1 Hz).
    Fmax : quantities.Quantity, optional
        Maximum frequency to include (default: 200 Hz).
    Zscored : bool, optional
        Apply z-score normalization (default: True).
    NormTime : tuple, optional
        Time window for baseline normalization. If None,
        uses full signal for z-score normalization.
    dtype : type, optional
        Output data type (default: float).
    **specKwarg
        Additional scipy.signal.spectrogram keyword arguments.

    Returns
    -------
    neo.AnalogSignal
        Spectrogram data with frequency array in annotations.
        Shape is (time_bins, freq_bins_subset).
    """
    nFFT = int(2 ** (np.around(np.log2(sig.sampling_rate / Fres)) + 1))
    Ts = sig.sampling_period
    noverlap = int((Ts * nFFT - TimeRes) / Ts)

    f, t, Sxx = signal.spectrogram(sig,
                                   fs=sig.sampling_rate,
                                   nperseg=nFFT,
                                   noverlap=noverlap,
                                   axis=0,
                                   **specKwarg)

    finds = np.where((Fmin < f) & (f < Fmax))[0][1:]
    r, g, c = Sxx.shape
    data = Sxx.reshape((r, c))[finds][:]

    if Zscored and (NormTime is None):
        data = zscore(data, axis=1)

    s = sig.duplicate_with_new_data(data.astype(dtype).transpose())
    s.annotate(Freq=f[finds])
    s.annotate(spec=True)
    s.annotate(nFFT=nFFT)
    s.annotate(WindowTime=nFFT * Ts)
    s.sampling_period = np.mean(t[1:] - t[:-1]) * pq.s
    s.t_start = s.t_start + (nFFT * Ts) / 2

    if Zscored and (NormTime is None):
        return s

    if Zscored:
        NormSig = s.time_slice(NormTime[0], NormTime[1])
        mean = np.mean(NormSig, axis=0)
        std = np.std(NormSig, axis=0)
        return ((s - mean) / std)
    else:
        if NormTime is None:
            return s
        else:
            NormSig = s.time_slice(NormTime[0], NormTime[1])
            mean = np.mean(NormSig, axis=0)
            return (s / mean)


def AvgSpectrogram(sig, TimesEvent, TimeAvg, SpecArgs,
                   AvgSpectNorm='Zscore', AvgSpectNormTime=None,
                   TrialProcessChain=None, **kwargs):
    """
    Compute averaged spectrogram aligned to trigger events.

    Extracts spectrogram segments around trigger times, averages them, and applies
    optional baseline normalization. Useful for analyzing frequency content of
    event-locked responses.

    Parameters
    ----------
    sig : neo.AnalogSignal
        Input signal to analyze.
    TimesEvent : array-like
        Event trigger times.
    TimeAvg : tuple
        Time window relative to events (before, after).
    SpecArgs : dict
        Spectrogram parameters (Fmin, Fmax, Fres, TimeRes, Zscored).
    AvgSpectNorm : str, optional
        Normalization type: 'Zscore' or other (default: 'Zscore').
    AvgSpectNormTime : tuple, optional
        Time window for baseline normalization.
        If None, uses entire averaged spectrogram.
    TrialProcessChain : list, optional
        Processing functions applied to each trial
        before averaging.
    **kwargs
        Additional keyword arguments.

    Returns
    -------
    neo.AnalogSignal
        Averaged and normalized spectrogram.
    """
    Acc = np.array([])
    Trials = 0
    for et in TimesEvent:
        if et + TimeAvg[0] < sig.t_start:
            print(et, ' not valid')
            continue
        elif et + TimeAvg[1] > sig.t_stop:
            print(et, ' not valid')
            continue

        Trials += 1
        s = sig.time_slice(et + TimeAvg[0], et + TimeAvg[1])
        if TrialProcessChain is not None:
            st = ApplyProcessChain(s, TrialProcessChain)
        else:
            st = s

        spect = Spectrogram(st, **SpecArgs)
        Acc = Acc + np.array(spect) if Acc.size else np.array(spect)

    AvgSpect = spect.duplicate_with_new_data(Acc / Trials)
    AvgSpect.t_start = TimeAvg[0] + AvgSpect.annotations['WindowTime'] / 2

    if AvgSpectNorm is None:
        return AvgSpect

    if AvgSpectNormTime is None:
        NormTime = (AvgSpect.t_start, -0 * pq.s)
    else:
        NormTime = AvgSpectNormTime
        if NormTime[0] is None:
            NormTime[0] = AvgSpect.t_start
        if NormTime[1] is None:
            NormTime[1] = AvgSpect.t_stop

    NormSig = AvgSpect.time_slice(NormTime[0], NormTime[1])
    mean = np.mean(NormSig, axis=0)
    std = np.std(NormSig, axis=0)

    if AvgSpectNorm == 'Zscore':
        AvgNormSpect = (AvgSpect - mean) / std
    else:
        AvgNormSpect = AvgSpect / mean

    return AvgNormSpect


def TrigAveraging(sig, TimesEvent, TimeAvg, TrialProcessChain=None):
    """
    Compute trial-averaged signal aligned to trigger events.

    Extracts signal segments around trigger times, optionally applies per-trial
    processing, and computes mean and standard deviation of the ensemble.

    Parameters
    ----------
    sig : neo.AnalogSignal
        Input signal to average.
    TimesEvent : array-like
        Event trigger times.
    TimeAvg : tuple
        Time window relative to events (before, after).
    TrialProcessChain : list, optional
        Processing functions applied to each trial
        before averaging.

    Returns
    -------
    neo.AnalogSignal
        Averaged signal with annotations containing 'std' (std dev signal)
        and 'acc' (all stacked trials for ensemble statistics).
    """
    Ts = sig.sampling_period
    nSamps = int((TimeAvg[1] - TimeAvg[0]) / Ts)
    acc = None
    for et in TimesEvent:
        if et + TimeAvg[0] < sig.t_start:
            print(et, ' not valid')
            continue
        elif et + TimeAvg[1] > sig.t_stop:
            print(et, ' not valid')
            continue

        Samp1 = sig.time_index(et + TimeAvg[0])
        s = sig[Samp1:Samp1 + nSamps, :]

        if TrialProcessChain is not None:
            st = ApplyProcessChain(s, TrialProcessChain)
        else:
            st = s

        st.t_start = TimeAvg[0]
        st.array_annotations = {}
        st.annotations = {}
        if acc is None:
            acc = st
        else:
            acc = acc.merge(st)

    avg = acc.duplicate_with_new_data(np.mean(acc, axis=1))
    std = acc.duplicate_with_new_data(np.std(acc, axis=1))
    avg.annotate(std=std)
    avg.annotate(acc=acc)
    return avg


def Derivative(sig):
    """
    Compute first derivative of signal.

    Calculates the discrete first derivative (dV/dt) of the signal, returning
    a new signal with adjusted start time to align sample positions.

    Parameters
    ----------
    sig : neo.AnalogSignal
        Input signal.

    Returns
    -------
    neo.AnalogSignal
        First derivative with units of [original_units / time].
    """
    derivative_sig = AnalogSignal(
        np.diff(sig.as_quantity(), axis=0) / sig.sampling_period,
        t_start=sig.t_start + sig.sampling_period / 2,
        sampling_period=sig.sampling_period,
        name=sig.name, **sig.annotations)

    return derivative_sig


def DownSampling(sig, Fact, zero_phase=True):
    """
    Downsample signal by integer factor using IIR decimation.

    Reduces sampling rate by an integer factor using scipy's decimate function,
    optionally applying zero-phase filtering to avoid phase distortion.

    Parameters
    ----------
    sig : neo.AnalogSignal
        Input signal to downsample.
    Fact : int
        Decimation factor (output rate = input rate / Fact).
    zero_phase : bool, optional
        Apply zero-phase filtering (default: True).

    Returns
    -------
    neo.AnalogSignal
        Downsampled signal with reduced sampling rate.
    """
    print(sig.sampling_rate, sig.sampling_rate / Fact)
    rs = signal.decimate(np.array(sig),
                         q=Fact,
                         zero_phase=zero_phase,
                         ftype='iir',
                         axis=0)
    sig = sig.duplicate_with_new_data(signal=rs * sig.units)
    sig.sampling_rate = sig.sampling_rate / Fact
    return sig


def RemoveDC(sig, Type='constant'):
    """
    Remove DC offset and low-frequency trends from signal.

    Applies scipy's detrend function to remove constant offset or polynomial trends
    from the signal.

    Parameters
    ----------
    sig : neo.AnalogSignal
        Input signal.
    Type : str, optional
        Detrending type: 'constant', 'linear', etc.
        (default: 'constant' removes DC offset).

    Returns
    -------
    neo.AnalogSignal
        Detrended signal with same units.
    """
    st = np.array(sig)
    st = signal.detrend(st, type=Type, axis=0)
    return sig.duplicate_with_new_data(signal=st * sig.units)


def SetZero(sig, TWind=None):
    """
    Subtract baseline to set zero reference voltage.

    Removes the mean voltage over a specified time window, useful for establishing
    a common baseline reference across recordings.

    Parameters
    ----------
    sig : neo.AnalogSignal
        Input signal.
    TWind : tuple, optional
        Time window for baseline estimation (start, stop).
        If None, uses first 30 seconds of signal.

    Returns
    -------
    neo.AnalogSignal
        Signal with baseline subtracted.
    """
    if TWind is None:
        TWind = (sig.t_start, sig.t_start + 30 * pq.s)
    st = np.array(sig)
    # offset = np.mean(sig.GetSignal(TWind))
    offset = np.mean(sig.time_slice(TWind[0], TWind[1]), axis=0)
    print(sig.name, offset)
    st_corrected = st - offset.magnitude
    return sig.duplicate_with_new_data(signal=st_corrected * sig.units)


def Gain(sig, Gain):
    """
    Apply multiplicative gain to signal.

    Scales signal amplitude by a constant gain factor.

    Parameters
    ----------
    sig : neo.AnalogSignal
        Input signal.
    Gain : float
        Multiplication factor.

    Returns
    -------
    neo.AnalogSignal
        Scaled signal.
    """
    return sig * Gain


def Resample(sig, Fs=None, MaxPoints=None):
    """
    Resample signal to new sampling rate using polyphase filtering.

    Resamples signal using rational resampling with automatic downsampling factor
    selection based on target frequency or maximum point count.

    Parameters
    ----------
    sig : neo.AnalogSignal
        Input signal.
    Fs : quantities.Quantity, optional
        Target sampling rate.
    MaxPoints : int, optional
        Target maximum number of samples. Takes precedence
        over Fs if both provided.

    Returns
    -------
    neo.AnalogSignal
        Resampled signal with new sampling rate.
    """
    if MaxPoints is None:
        f = Fs / sig.sampling_rate
        fact = Fraction(float(f)).limit_denominator()
        dowrate = fact.denominator
        uprate = fact.numerator
    else:
        dowrate = int(sig.times.shape[0] / MaxPoints)
        if dowrate > 0:
            f = float(1 / float(dowrate))
            uprate = 1

    if dowrate > 0:
        print('Down sampling', sig.sampling_rate * f, f, uprate, dowrate)
        rs = signal.resample_poly(np.array(sig), uprate, dowrate)
        sig = sig.duplicate_with_new_data(signal=rs * sig.units)
        sig.sampling_rate = sig.sampling_rate * f
        return sig
    else:
        return sig


def Abs(sig):
    """
    Compute absolute value of signal.

    Takes the element-wise absolute value while preserving signal structure and units.

    Parameters
    ----------
    sig : neo.AnalogSignal
        Input signal.

    Returns
    -------
    neo.AnalogSignal
        Absolute value of signal with same units.
    """
    st = np.array(sig)
    st = np.abs(st)

    return sig.duplicate_with_new_data(signal=st * sig.units)


def power(sig):
    """
    Compute instantaneous power (squared amplitude) of signal.

    Squares the signal values element-wise. Useful for power spectral analysis
    and RMS calculations.

    Parameters
    ----------
    sig : neo.AnalogSignal
        Input signal.

    Returns
    -------
    neo.AnalogSignal
        Squared signal with units of [original_units]^2.
    """
    st = np.array(sig) ** 2
    #    st = st**2
    return sig.duplicate_with_new_data(signal=st * sig.units)


def Filter(sig, Type, Order, Freqs):
    """
    Apply Butterworth IIR filter to signal.

    Implements zero-phase digital filtering using second-order sections (SOS)
    format for numerical stability.

    Parameters
    ----------
    sig : neo.AnalogSignal
        Input signal to filter.
    Type : str
        Filter type: 'lowpass', 'highpass', 'bandpass', 'bandstop'.
    Order : int
        Filter order (pole count).
    Freqs : quantities.Quantity or array-like
        Critical frequency/frequencies.
        For bandpass/bandstop, array of 2 values.

    Returns
    -------
    neo.AnalogSignal
        Filtered signal with same units and structure.
    """

    st = np.array(sig)
    Fs = sig.sampling_rate.magnitude

    # Normalize frequencies to Nyquist
    freqs = np.asarray(Freqs, dtype=float) / (0.5 * Fs)

    # ---- Normalize shape (critical fix) ----
    if freqs.ndim == 0:
        freqs = float(freqs)
    elif freqs.size == 1:
        freqs = float(freqs[0])   # convert [x] → x
    elif freqs.size == 2:
        freqs = np.sort(freqs)    # ensure low < high
    else:
        raise ValueError("Freqs must be scalar or length-2 for band filters")

    # ---- Validate by filter type ----
    Type = Type.lower()

    if Type in ["lowpass", "highpass"]:
        if not np.isscalar(freqs):
            raise ValueError(f"{Type} requires a scalar cutoff frequency")

    elif Type in ["bandpass", "bandstop"]:
        if not (isinstance(freqs, (list, tuple, np.ndarray)) and len(freqs) == 2):
            raise ValueError(f"{Type} requires two cutoff frequencies")

    else:
        raise ValueError(f"Invalid filter type: {Type}")

    # ---- Validate frequency range ----
    if np.any(np.asarray(freqs) <= 0) or np.any(np.asarray(freqs) >= 1):
        raise ValueError("Frequencies must be in (0, Nyquist)")

    # ---- Design + apply filter ----
    sos = signal.butter(Order, freqs, btype=Type, output='sos')
    st = signal.sosfiltfilt(sos, st, axis=0)

    return sig.duplicate_with_new_data(signal=st * sig.units)


def ThresholdTrianGen(sig, RelaxTime=0.4 * pq.s, threshold=None, sign='below'):
    """
    Detect threshold crossings and generate spike train.

    Identifies times when signal crosses a threshold, with minimum time interval
    between detected events to avoid multiple detections of slow threshold crossings.

    Parameters
    ----------
    sig : neo.AnalogSignal
        Input signal to analyze.
    RelaxTime : quantities.Quantity, optional
        Minimum interval between detections
        (default: 0.4 s).
    threshold : float, optional
        Threshold value. If None, uses mean + std of signal.
    sign : str, optional
        Crossing direction: 'below', 'above', or 'both'
        (default: 'below').

    Returns
    -------
    neo.SpikeTrain
        Detected spike times with same t_start and t_stop as input.
    """
    if threshold is None:
        threshold = np.mean(sig) + np.std(sig)
    inttimes = Ran.threshold_detection(signal=sig,
                                       threshold=threshold,
                                       sign=sign,
                                       RelaxTime=RelaxTime)
    inttimes = np.array(inttimes)

    return SpikeTrain(times=inttimes,
                      units='s',
                      t_start=sig.t_start,
                      t_stop=sig.t_stop,
                      )


def SplineSmooth(sig, sFact=2, **kwargs):
    """
    Smooth signal using univariate spline interpolation.

    Fits a smoothing spline to the signal and evaluates at original sample points.
    Reduces noise while preserving signal features.

    Parameters
    ----------
    sig : neo.AnalogSignal
        Input signal to smooth.
    sFact : int, optional
        Smoothing factor: output has sFact times fewer knots
        than input samples (default: 2).
    **kwargs
        Additional scipy UnivariateSpline keyword arguments.

    Returns
    -------
    neo.AnalogSignal
        Smoothed signal.
    """
    s = sig.shape[0] / sFact
    spl = UnivariateSpline(sig.times, sig, s=s)
    return sig.duplicate_with_new_data(signal=spl(sig.times) * sig.units)


def MedianFilt(sig, window_size=None, **kwargs):
    """
    Apply median filter to remove noise and outliers.

    Applies a nonlinear median filter with specified kernel size. Kernel size
    is automatically adjusted to be odd if necessary.

    Parameters
    ----------
    sig : neo.AnalogSignal
        Input signal to filter.
    window_size : quantities.Quantity, optional
        Filter window duration in seconds.
        If None, uses 1/10th of signal duration.
    **kwargs
        Additional keyword arguments.

    Returns
    -------
    neo.AnalogSignal
        Filtered signal.
    """
    if window_size == None:
        kernel_size = int(sig.shape[0] / 10)
    else:
        kernel_size = int(sig.sampling_rate * window_size)
    if bool(kernel_size & 1):  # kernel_size has to be an odd number
        kernel_size = kernel_size
    else:
        kernel_size = kernel_size - 1
    SigFilt = medfilt(np.array(sig).reshape(len(sig), ),
                      kernel_size=kernel_size)  # TODO change by axis
    return sig.duplicate_with_new_data(signal=np.array(SigFilt) * sig.units)


def dbrms(x, axis=-1):
    """
    Compute RMS value in decibels.

    Calculates root mean square in 20*log10 scale (dB).

    Parameters
    ----------
    x : array-like
        Input data.
    axis : int, optional
        Axis along which to compute (default: -1, last axis).

    Returns
    -------
    float or ndarray
        RMS value in dB.
    """
    return 20 * np.log10(np.sqrt(np.mean(x ** 2, axis=axis)))


def rms(x, axis=-1):
    """
    Compute RMS value (root mean square).

    Calculates sqrt(mean(x^2)) along specified axis.

    Parameters
    ----------
    x : array-like
        Input data.
    axis : int, optional
        Axis along which to compute (default: -1, last axis).

    Returns
    -------
    float or ndarray
        RMS value.
    """
    return np.sqrt(np.mean(x ** 2, axis=axis))


def power_sliding(x, axis=-1):
    """
    Compute instantaneous power (mean of squared values).

    Calculates mean(x^2) along specified axis.

    Parameters
    ----------
    x : array-like
        Input data.
    axis : int, optional
        Axis along which to compute (default: -1, last axis).

    Returns
    -------
    float or ndarray
        Power value.
    """
    return np.mean(x ** 2, axis=axis)


def strides_signal(sig, timewidth, steptime=None):
    """
    Create sliding window views of signal data.

    Generates overlapping windows of specified width and step size, useful for
    sliding window analysis and feature extraction.

    Parameters
    ----------
    sig : neo.AnalogSignal
        Input signal (2D or 3D).
    timewidth : quantities.Quantity
        Window duration.
    steptime : quantities.Quantity, optional
        Step size between windows.
        If None, uses timewidth/10.

    Returns
    -------
    tuple
        (strides, time_width, step_time) where:
        - strides: Array of windowed data views
        - time_width: Window duration in seconds
        - step_time: Step duration in seconds
    """
    if steptime is None:
        steptime = timewidth / 10

    window_size = int(timewidth.rescale('s') / sig.sampling_period.rescale('s'))
    time_width = sig.sampling_period.rescale('s') * window_size
    step_size = int(steptime.rescale('s') / sig.sampling_period.rescale('s'))
    step_time = sig.sampling_period.rescale('s') * step_size

    if len(sig.shape) == 3:
        strides = sliding_window_view(sig, window_size, 0)[::step_size, :, :, :]
    elif len(sig.shape) == 2:
        strides = sliding_window_view(sig, window_size, 0)[::step_size, :, :]
    else:
        print('Error in dimension')
        return None

    return strides, time_width, step_time


def sliding_window(sig, timewidth, steptime=None, func=None, **fkwargs):
    """
    Apply function to sliding windows of signal.

    Extracts sliding windows and applies a user-defined function to compute
    statistics or features within each window, returning a new signal with one
    sample per window.

    Parameters
    ----------
    sig : neo.AnalogSignal or neo.ImageSequence
        Input signal.
    timewidth : quantities.Quantity
        Window duration.
    steptime : quantities.Quantity, optional
        Step size between windows.
        If None, uses timewidth/10.
    func : callable
        Function to apply to each window. Should accept array-like
        and return scalar or 1D array.
    **fkwargs
        Additional keyword arguments passed to func.

    Returns
    -------
    neo.AnalogSignal or ImageSequence
        Result with sampling rate corresponding
        to window step time.
    """
    strides, time_width, step_time = strides_signal(sig, timewidth, steptime)

    st = func(strides, **fkwargs)

    if len(sig.shape) == 3:
        ret = ImageSequence(st,
                            units=sig.units,
                            t_start=sig.t_start + time_width / 2,
                            name=sig.name,
                            sampling_rate=(1 / step_time).rescale('Hz'),
                            **sig.annotations)
    else:
        ret = AnalogSignal(signal=st,
                           units=sig.units,
                           t_start=sig.t_start + time_width / 2,
                           name=sig.name,
                           sampling_rate= (1 / step_time).rescale('Hz'),
                           **sig.annotations)
    return ret


def CrossCorr(x1, x2, fs, **fkwargs):
    """
    Compute cross-correlation between signal pairs.

    Calculates normalized cross-correlation and corresponding lags for pairs of
    signals. Returns maximum correlation and lag for each channel pair.

    Parameters
    ----------
    x1 : ndarray
        First signal array (channels, samples).
    x2 : ndarray
        Second signal array (channels, samples).
    fs : quantities.Quantity
        Sampling frequency for lag conversion.
    **fkwargs
        Additional keyword arguments.

    Returns
    -------
    tuple
        (max_corr, lags) where:
        - max_corr: Maximum correlation value per channel pair (dimensionless)
        - lags: Lag times corresponding to maximum correlation
    """
    max_corr = np.ones((x1.shape[0],))
    lags_idx = np.ones((x1.shape[0],))
    for ic, (s1, s2) in enumerate(zip(x1, x2)):
        s1 = s1.ravel()
        s2 = s2.ravel()
        correl = np.correlate(s1, s2, mode='full')
        correl = correl / (np.linalg.norm(s1) * np.linalg.norm(s2))
        idx = np.argmax(correl)
        max_corr[ic] = correl[idx]
        lags_idx[ic] = idx - (s1.shape[0] - 1)

    return max_corr * pq.dimensionless, lags_idx / fs


def PearsonCorr(x1, x2, fs, **fkwargs):
    """
    Compute Pearson correlation coefficient between signal pairs.

    Calculates linear correlation between paired signals. Returns correlation
    coefficient for each channel pair.

    Parameters
    ----------
    x1 : ndarray
        First signal array (channels, samples).
    x2 : ndarray
        Second signal array (channels, samples).
    fs : quantities.Quantity
        Sampling frequency (for compatibility, not used).
    **fkwargs
        Additional keyword arguments passed to scipy.stats.pearsonr.

    Returns
    -------
    quantities.Quantity
        Pearson correlation coefficient per channel pair
        (dimensionless, range -1 to 1).
    """
    pearson_corr = np.ones((x1.shape[0],))
    for ic, (s1, s2) in enumerate(zip(x1, x2)):
        s1 = s1.ravel()
        s2 = s2.ravel()
        res = stats.pearsonr(s1, s2, **fkwargs)
        pearson_corr[ic] = res[0]

    return pearson_corr * pq.dimensionless


def sliding_window_2sigs(sig1, sig2, timewidth, steptime=None, func=CrossCorr, **fkwargs):
    """
    Apply function to sliding windows of two signals.

    Extracts sliding windows from two signals and applies a user-defined function
    to compute correlations or other bivariate statistics within each window pair,
    returning results as a new signal with one sample per window.

    Parameters
    ----------
    sig1 : neo.AnalogSignal
        First input signal.
    sig2 : neo.AnalogSignal
        Second input signal.
    timewidth : quantities.Quantity
        Window duration.
    steptime : quantities.Quantity, optional
        Step size between windows.
        If None, uses timewidth/10.
    func : callable, optional
        Function to apply to each window pair. Should accept
        two array-like arguments and return scalar or 1D array.
        Default: CrossCorr.
    **fkwargs
        Additional keyword arguments passed to func.

    Returns
    -------
    neo.AnalogSignal or tuple
        Result(s) with sampling rate corresponding
        to window step time. If func returns tuple,
        returns tuple of AnalogSignals.
    """
    strides1, time_width, step_time = strides_signal(sig1, timewidth, steptime)
    strides2, _, _ = strides_signal(sig2, timewidth, steptime)

    st = func(strides1, strides2, fs=sig1.sampling_rate, **fkwargs)

    if isinstance(st, tuple):
        ret = []
        for ic, s in enumerate(st):
            ret.append(AnalogSignal(signal=s,
                                    # units=pq.dimensionless,
                                    t_start=sig1.t_start + time_width / 2,
                                    name=f'Corr_{ic}_{sig1.name}-{sig2.name}',
                                    sampling_rate=(1 / step_time).rescale('Hz'),
                                    **sig1.annotations)
                       )
        return tuple(ret)
    else:
        return AnalogSignal(signal=st,
                            # units=pq.dimensionless,
                            t_start=sig1.t_start + time_width / 2,
                            name=f'Corr_{sig1.name}-{sig2.name}',
                            sampling_rate=(1 / step_time).rescale('Hz'),
                            **sig1.annotations)

