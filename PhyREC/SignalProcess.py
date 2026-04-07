#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
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
    derivative_sig = AnalogSignal(
        np.diff(sig.as_quantity(), axis=0) / sig.sampling_period,
        t_start=sig.t_start + sig.sampling_period / 2,
        sampling_period=sig.sampling_period,
        name=sig.name, **sig.annotations)

    return derivative_sig


def DownSampling(sig, Fact, zero_phase=True):
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
    st = np.array(sig)
    st = signal.detrend(st, type=Type, axis=0)
    return sig.duplicate_with_new_data(signal=st * sig.units)


def SetZero(sig, TWind=None):
    if TWind is None:
        TWind = (sig.t_start, sig.t_start + 30 * pq.s)
    st = np.array(sig)
    # offset = np.mean(sig.GetSignal(TWind))
    offset = np.mean(sig.time_slice(TWind[0], TWind[1]), axis=0)
    print(sig.name, offset)
    st_corrected = st - offset.magnitude
    return sig.duplicate_with_new_data(signal=st_corrected * sig.units)


def Gain(sig, Gain):
    return sig * Gain


def Resample(sig, Fs=None, MaxPoints=None):
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
        print(sig.sampling_rate * f, f, uprate, dowrate)
        rs = signal.resample_poly(np.array(sig), uprate, dowrate)
        sig = sig.duplicate_with_new_data(signal=rs * sig.units)
        sig.sampling_rate = sig.sampling_rate * f
        return sig
    else:
        return sig


def Abs(sig):
    st = np.array(sig)
    st = np.abs(st)

    return sig.duplicate_with_new_data(signal=st * sig.units)


def power(sig):  # to solve units
    st = np.array(sig) ** 2
    #    st = st**2
    return sig.duplicate_with_new_data(signal=st * sig.units)


def Filter(sig, Type, Order, Freqs):
    st = np.array(sig)
    Fs = sig.sampling_rate.magnitude
    freqs = Freqs / (0.5 * Fs)

    # b, a = signal.butter(Order, freqs, Type)
    # st = signal.filtfilt(b, a, st, axis=0)

    sos = signal.butter(Order, freqs, Type, output='sos')
    st = signal.sosfiltfilt(sos, st, axis=0)

    # DbgFplt.PlotResponse(a, b, Fs)
    # DbgFplt.PlotResponse(sos, Fs)

    return sig.duplicate_with_new_data(signal=st * sig.units)


def ThresholdTrianGen(sig, RelaxTime=0.4 * pq.s, threshold=None, sign='below'):
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
    s = sig.shape[0] / sFact
    spl = UnivariateSpline(sig.times, sig, s=s)
    return sig.duplicate_with_new_data(signal=spl(sig.times) * sig.units)


def MedianFilt(sig, window_size=None, **kwargs):  # window_size in pq.s
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
    return 20 * np.log10(np.sqrt(np.mean(x ** 2, axis=axis)))


def rms(x, axis=-1):
    return np.sqrt(np.mean(x ** 2, axis=axis))


def power_sliding(x, axis=-1):
    return np.mean(x ** 2, axis=axis)


def strides_signal(sig, timewidth, steptime=None):
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
    pearson_corr = np.ones((x1.shape[0],))
    for ic, (s1, s2) in enumerate(zip(x1, x2)):
        s1 = s1.ravel()
        s2 = s2.ravel()
        res = stats.pearsonr(s1, s2, **fkwargs)
        pearson_corr[ic] = res[0]

    return pearson_corr * pq.dimensionless


def sliding_window_2sigs(sig1, sig2, timewidth, steptime=None, func=CrossCorr, **fkwargs):
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

# def xcorr(x, y, normed=True, detrend=mlab.detrend_none,
#           usevlines=True, maxlags=10, **kwargs):
#     Nx = len(x)
#     if Nx != len(y):
#         raise ValueError('x and y must be equal length')
#
#     x = detrend(np.asarray(x))
#     y = detrend(np.asarray(y))
#
#     correls = np.correlate(x, y, mode="full")
#
#     if normed:
#         correls /= np.sqrt(np.dot(x, x) * np.dot(y, y))
#
#     if maxlags is None:
#         maxlags = Nx - 1
#
#     if maxlags >= Nx or maxlags < 1:
#         raise ValueError('maxlags must be None or strictly '
#                          'positive < %d' % Nx)
#
#     lags = np.arange(-maxlags, maxlags + 1)
#     correls = correls[Nx - 1 - maxlags:Nx + maxlags]
#
#     return lags, correls
#
#
#  def ThresholdInstantRate(sig, RelaxTime=0.1*pq.s, threshold=None,
#                          OutSampling=0.01*pq.s,):
#
#     return elephant.statistics.instantaneous_rate(ThresholdTrianGen(sig,
#                                                                     RelaxTime,
#                                                                     threshold),
#                                                   sampling_period=OutSampling)
#
#
# def HilbertInstantFreq(sig, MaxFreq=20, MinFreq=0):
#     SigH = elephant.signal_processing.hilbert(sig)
#     insfreq = np.diff(np.angle(SigH)[:, 0]) / np.diff(SigH.times)
#
#     return AnalogSignal(signal=np.clip(insfreq.magnitude, 0, 20),
#                         units='Hz',
#                         name=sig.name,
#                         sampling_rate=SigH.sampling_rate,
#                         t_start=SigH.t_start)
#
#
# def HilbertAngle(sig):
#     SigH = elephant.signal_processing.hilbert(sig)
#     return sig.duplicate_with_new_data(signal=np.angle(SigH),
#                                        units=pq.radians)
#
#
# def HilbertAmp(sig):
#     SigH = elephant.signal_processing.hilbert(sig)
#
#     return sig.duplicate_with_new_data(signal=np.array(np.abs(SigH))*sig.units)
