from matplotlib import pyplot as plt
import PhyREC.PlotWaves as PltW
import quantities as pq
import importlib.resources as pkg_resources
from neo import AnalogSignal
import numpy as np
import quantities as pq

# ============================================================================
# SIGNAL GENERATION BLOCK (Reusable)
# ============================================================================

# ============================================================================
# 1. DEFINE SIGNAL PARAMETERS
# ============================================================================
# Sampling rate (Hz) - typical for electrophysiology recordings
Fs = 30000  # 30 kHz sampling rate

# Signal duration (seconds)
T = 300  # 300 seconds of data

# Frequency of the sinusoidal component we're simulating (Hz)
f0 = 10  # 10 Hz oscillation

# Amplitude of the sinusoidal signal (volts)
A = 250e-6  # 250 µV (micro-volts converted to volts)

# Gaussian noise standard deviation (volts) - simulates background noise
noise_rms = 20e-6  # 20 µV RMS noise level

# Create time vector in seconds
# This defines the time points at which we'll sample our signal
t = np.arange(0, T, 1 / Fs)

# ============================================================================
# 2. GENERATE SYNTHETIC SIGNALS WITH PHASE SHIFTS
# ============================================================================
# Number of channels to simulate
nSigs = 16
Sigs = []

for i in range(nSigs):
    # Define phase shift for each channel (in radians)
    # Each channel has a different phase to simulate multi-electrode recordings
    # with temporal relationships
    phase_shift = (i * 2 * np.pi) / nSigs  # Linear phase distribution across channels

    # Create clean sinusoidal component with phase shift
    # This simulates periodic neural oscillations with phase differences
    # between channels (e.g., from waves propagating across the tissue)
    signal = A * np.sin(2 * np.pi * f0 * t + phase_shift)

    # Generate Gaussian noise to simulate real electrophysiology recording noise
    # Noise is Gaussian-distributed with specified RMS value
    noise = np.random.normal(0, noise_rms, size=t.shape)

    # Combine signal and noise to create realistic synthetic data
    sig = signal + noise

    # Create neo.AnalogSignal object
    # AnalogSignal is the standard container for continuous electrophysiology data
    # It includes both data and metadata (sampling rate, units, channel name)
    Sigs.append(AnalogSignal(
        sig,                           # Signal data array
        units='V',                     # Physical units (volts)
        sampling_rate=Fs*pq.Hz,        # Sampling rate with quantities units
        name='signal_' + str(i),       # Channel identifier
    ))

# ============================================================================
# PLOTTING BLOCK
# ============================================================================
# Note: This block uses Sigs from the signal generation block above

# ============================================================================
# 1. LOAD PHYREC STYLING
# ============================================================================
# Apply the PhyREC matplotlib style for professional-looking plots
# This style includes optimized colors, fonts, and figure layout
with pkg_resources.path("PhyREC.style", "PhyREC.mplstyle") as p:
    plt.style.use(p)


# ============================================================================
# 2. CREATE MATPLOTLIB FIGURE WITH SUBPLOTS
# ============================================================================
# Create figure with one subplot per channel, sharing the same x-axis
# This allows synchronized zooming/panning across all channels
pltFig, pltAxes = plt.subplots(len(Sigs), 1, sharex=True, figsize=(12, 8))

# ============================================================================
# 3. CREATE WAVESLOT OBJECTS FOR EACH CHANNEL
# ============================================================================
# WaveSlot is PhyREC's main visualization class for waveforms
# It handles plot rendering, time navigation, and signal display
Slots = []
for i, sig in enumerate(Sigs):
    # Create a WaveSlot for each signal channel
    Slots.append(PltW.WaveSlot(
        sig,                    # neo.AnalogSignal object containing the data
        color='blue',           # Line color for the waveform
        Units='mV',             # Display units (data will be converted from V to mV)
        label=sig.name,         # Legend label for this channel
        Ax=pltAxes[i]           # Matplotlib axes to plot on
    ))

# ============================================================================
# 4. CREATE PLOTSLOTS CONTAINER WITH INTERACTIVE CONTROLS
# ============================================================================
# PlotSlots manages multiple WaveSlot instances as a group
# It provides unified controls for visualization, formatting, and time navigation
Splt = PltW.PlotSlots(
    Slots,              # List of WaveSlot objects to manage
    Fig=pltFig,         # Matplotlib figure containing all axes
    LiveControl=True,   # Enable interactive controls (time slider, buttons)
)

# ============================================================================
# 5. ADD LEGEND AND PLOT DATA
# ============================================================================
# Add legend to the figure
# Shows the channel names and colors for easy identification
Splt.AddLegend()

# Plot the signals in the time window from 0 to 10 seconds
# Quantities with units (pq.s = seconds) ensure proper time handling
Splt.PlotChannels((10*pq.s, 11*pq.s))

# ============================================================================
# 6. INTERACTIVE FEATURES
# ============================================================================
# With LiveControl=True enabled, the figure now includes:
#   - Time slider for navigating through the signal
#   - Zoom and pan tools via matplotlib toolbar
#   - Synchronized viewing across all channels
#   - Legend display with channel information

plt.tight_layout()
Splt.AddLegend()