"""
basic_filter_example.py

A short, documented example showing how to generate synthetic analog signals,
apply simple band-pass and band-stop filters using PhyREC.SignalProcess, and
plot results with PhyREC.PlotWaves.

This file is intended for users learning how to:
 - create neo.AnalogSignal objects from numpy data
 - apply filters (bandpass / bandstop) provided by PhyREC
 - display raw and filtered traces with PlotWaves

Notes on units:
 - Signals are generated in volts (V). PlotWaves is instructed to display in mV.
 - Sampling frequency (Fs) is provided in Hz; quantities (pq) are used where
   PhyREC/neo expect a Quantity (sampling_rate=Fs * pq.Hz).

"""

from matplotlib import pyplot as plt
from neo import AnalogSignal
import numpy as np
import PhyREC.PlotWaves as prplt
import PhyREC.SignalProcess as prspro
import quantities as pq
from importlib.resources import files

# ============================================================================
# 1. LOAD PHYREC STYLING
# ============================================================================
# Apply the PhyREC matplotlib style for a consistent, publication-ready
# appearance across example figures. This file is bundled with the package.
plt.style.use(files("PhyREC.style") / "PhyREC.mplstyle")

# ---------------------------------------------------------------------------
# 2. SIMULATION PARAMETERS
# ---------------------------------------------------------------------------
# Fs : sampling frequency in Hz
# T  : total duration in seconds
# f0 : sine wave frequency (Hz)
# A  : sine amplitude in volts (here 250 µV)
# noise_rms : RMS of additive Gaussian noise (in volts)
Fs = 30000  # Hz
T = 300  # seconds
f0 = 10  # Hz
A = 250e-6  # 250 µV → volts
noise_rms = 50e-6  # 50 µV RMS → volts
# time vector (numpy array) sampled at 1/Fs spacing
t = np.arange(0, T, 1 / Fs)

# Number of example channels to create
nSigs = 16
Sigs = []

# ---------------------------------------------------------------------------
# 3. CREATE SYNTHETIC SIGNALS
# ---------------------------------------------------------------------------
# Each channel is a sine wave with a different phase plus Gaussian noise. The
# result is converted to a neo.AnalogSignal so it can be consumed by PhyREC.
for i in range(nSigs):
    # stagger phases so the channels look different when plotted
    phase = i * np.pi / (nSigs - 1) if nSigs > 1 else 0
    # deterministic sinusoid (volts)
    signal = A * np.sin(2 * np.pi * f0 * t + phase)
    # zero-mean Gaussian noise with specified RMS
    noise = np.random.normal(0, noise_rms, size=t.shape)
    # combined trace
    sig = signal + noise

    # Wrap numpy array into neo.AnalogSignal. Neo and PhyREC often accept a
    # Quantity-based sampling_rate; here we provide both 'samplerate' (legacy)
    # and 'sampling_rate' to be robust across versions.
    Sigs.append(AnalogSignal(sig,
                             samplerate=Fs,
                             units='V',
                             sampling_rate=Fs*pq.Hz,
                             name='signal' + str(i)))

# ---------------------------------------------------------------------------
# 4. SETUP PLOTTING SLOTS
# ---------------------------------------------------------------------------
# Create one subplot per signal and collect WaveSlot objects which tell
# PlotWaves how to draw each trace (color, units, label, axis to plot on).
pltFig, pltAxes = plt.subplots(len(Sigs), 1, sharex=True, figsize=(15, 15))

Slots = []
for i, sig in enumerate(Sigs):
    # Raw (noisy) trace in blue, semi-transparent
    Slots.append(prplt.WaveSlot(sig,
                               color='blue',
                               alpha=0.5,
                               Units='mV',  # display units (converted from V)
                               label=sig.name,
                               Ax=pltAxes[i]))

    # -----------------------------------------------------------------------
    # Apply band-pass filter around the sine frequency to isolate oscillation
    # -----------------------------------------------------------------------
    sig_filt = prspro.Filter(sig,
                             Type='bandpass',
                             Order=4,
                             Freqs=(9, 11),  # pass-band in Hz
                             )

    # Filtered trace plotted in red
    Slots.append(prplt.WaveSlot(sig_filt,
                               color='red',
                               Units='mV',
                               label=sig.name + '_filt',
                               Ax=pltAxes[i]))

    # -----------------------------------------------------------------------
    # Apply band-stop (notch) filter at the same frequency range. The result
    # can be interpreted as the 'noise' after removing the oscillatory band.
    # -----------------------------------------------------------------------
    sig_filt = prspro.Filter(sig,
                             Type='bandstop',
                             Order=4,
                             Freqs=(9, 11),
                             )

    # Notch-filtered trace plotted in black
    Slots.append(prplt.WaveSlot(sig_filt,
                                color='black',
                                Units='mV',
                                label=sig.name + '_noise',
                                Ax=pltAxes[i]))

# ---------------------------------------------------------------------------
# 5. COMPOSE AND DISPLAY
# ---------------------------------------------------------------------------
# Pack all slots into a PlotSlots object which manages interactivity and
# convenient plotting helpers.
Splt = prplt.PlotSlots(Slots,
                      Fig=pltFig,
                      LiveControl=True)  # LiveControl enables interactive widgets

# Show a short window of the data (10.0s to 10.5s) so example plots zoom in
Splt.PlotChannels((10*pq.s, 10.5*pq.s))

# Add a legend describing traces
Splt.AddLegend()

# Tidy up spacing so labels and axes don't overlap
pltFig.tight_layout()