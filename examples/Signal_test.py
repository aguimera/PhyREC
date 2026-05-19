
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import butter, filtfilt

# Paràmetres
fs = 1000  # Hz
t = np.arange(0, 2, 1/fs)

# 1. Soroll pink (aproximació simple)
white = np.random.randn(len(t))
pink = np.cumsum(white)  # integració → espectre ~1/f² (aproximació)

# 2. Oscil·lacions neuronals
theta = 0.5 * np.sin(2*np.pi*6*t)   # 6 Hz
gamma = 0.2 * np.sin(2*np.pi*40*t)  # 40 Hz

# 3. Generar spikes
spikes = np.zeros_like(t)
spike_times = np.random.choice(len(t), size=20, replace=False)

for st in spike_times:
    width = int(0.005 * fs)  # 5 ms
    spike_shape = np.exp(-np.linspace(0, 5, width))
    spikes[st:st+width] += spike_shape[:len(spikes[st:st+width])]

# 4. Combinar
signal = 0.3 * pink + theta + gamma + spikes

# 5. Filtrat passa-baix (simular electròde)
def lowpass(data, cutoff, fs):
    b, a = butter(4, cutoff/(fs/2), btype='low')
    return filtfilt(b, a, data)

signal = lowpass(signal, 100, fs)

# Plot
plt.plot(t, signal)
plt.title("Senyal neuronal sintètica")
plt.xlabel("Temps (s)")
plt.ylabel("Amplitud")
plt.show()
