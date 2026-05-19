from matplotlib import pyplot as plt
from neo import AnalogSignal
import numpy as np
import PhyREC.PlotWaves as prplt
import PhyREC.SignalProcess as prspro
import quantities as pq

plt.style.use('PhyREC.mplstyle')

# Parameters
Fs = 30000  # Hz
T = 300  # seconds
f0 = 10  # Hz
A = 250e-6  # 250 µV → volts
noise_rms = 50e-6  # 50 µV RMS → volts
t = np.arange(0, T, 1 / Fs)

nSigs = 16
Sigs = []

for i in range(nSigs):
    phase = i * np.pi / (nSigs - 1) if nSigs > 1 else 0
    # Clean sinusoidal signal
    signal = A * np.sin(2 * np.pi * f0 * t + phase)
    # Noise (Gaussian with desired RMS)
    noise = np.random.normal(0, noise_rms, size=t.shape)
    # Final signal
    sig = signal + noise

    Sigs.append(AnalogSignal(sig,
                             samplerate=Fs,
                             units='V',
                             sampling_rate=Fs*pq.Hz,
                             name='signal' + str(i)))


pltFig, pltAxes = plt.subplots(len(Sigs), 1, sharex=True)

Slots = []
for i, sig in enumerate(Sigs):
    Slots.append(prplt.WaveSlot(sig,
                               color='blue',
                               alpha=0.5,
                               Units='mV',
                               label=sig.name,
                               Ax=pltAxes[i]))

    sig_filt = prspro.Filter(sig,
                             Type='bandpass',
                             Order=4,
                             Freqs=(9, 11),
                             )
    Slots.append(prplt.WaveSlot(sig_filt,
                               color='red',
                               Units='mV',
                               label=sig.name + '_filt',
                               Ax=pltAxes[i]))

Splt = prplt.PlotSlots(Slots,
                      Fig=pltFig,
                      LiveControl=True,)


Splt.PlotChannels((10*pq.s, 11*pq.s))

Splt.AddLegend()

