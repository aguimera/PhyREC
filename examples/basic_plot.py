"""
Basic Plot Example - PhyREC Library Demonstration

This script demonstrates the fundamental usage of the PhyREC library for 
visualizing electrophysiological signals. It creates synthetic multi-channel 
signal data (16 channels) and displays them using the library's interactive 
plotting capabilities.

The example shows:
    - Creating synthetic signals with phase offsets
    - Converting signals to Neo AnalogSignal format (standard for electrophysiology)
    - Setting up multi-channel plots
    - Using PhyREC's WaveSlot and PlotSlots classes for visualization
    - Interactive navigation with live controls
"""

from matplotlib import pyplot as plt
from neo import AnalogSignal
import numpy as np
import PhyREC.PlotWaves as PltW
import quantities as pq

# Load the PhyREC matplotlib style for consistent appearance
plt.style.use('PhyREC.mplstyle')

# ============================================================================
# Signal Generation Parameters
# ============================================================================

Fs = 30000  # Sampling frequency in Hz (30 kHz)
T = 300  # Total signal duration in seconds (5 minutes)
f0 = 10  # Frequency of the sinusoidal signal in Hz
A = 250e-6  # Amplitude of the signal in volts (250 µV)
noise_rms = 20e-6  # RMS (root mean square) noise level in volts (20 µV)

# Create time vector with sampling interval of 1/Fs seconds
t = np.arange(0, T, 1 / Fs)

# ============================================================================
# Synthetic Signal Creation
# ============================================================================

nSigs = 16  # Number of channels to simulate
Sigs = []  # List to store signal objects

# Generate multi-channel signals with varying phases
for i in range(nSigs):
    # Calculate phase offset for each channel (allows visualization of phase shifts)
    # Creates a phase ramp from 0 to π across all channels
    phase = i * np.pi / (nSigs - 1) if nSigs > 1 else 0

    # Generate clean sinusoidal signal at frequency f0
    signal = A * np.sin(2 * np.pi * f0 * t + phase)

    # Generate Gaussian noise with specified RMS level (simulates real recording noise)
    noise = np.random.normal(0, noise_rms, size=t.shape)

    # Combine clean signal with noise (typical electrophysiology data)
    sig = signal + noise

    # Convert to Neo AnalogSignal format (standard in electrophysiology)
    # This ensures compatibility with analysis pipelines and other Neo-based tools
    Sigs.append(AnalogSignal(sig,
                             units='V',  # Signal in volts
                             sampling_rate=Fs * pq.Hz,  # Specify sampling rate with units
                             name='signal' + str(i)))  # Descriptive name for each channel

# ============================================================================
# Plot Setup and Configuration
# ============================================================================

# Create a figure with one subplot per signal (16 subplots stacked vertically)
# sharex=True links all x-axes together for synchronized zooming/panning
pltFig, pltAxes = plt.subplots(len(Sigs), 1, sharex=True, figsize=(15,15))

# Create WaveSlot objects for each signal
# WaveSlot is a PlotWaves class that handles signal visualization and interactive features
Slots = []
for i, sig in enumerate(Sigs):
    # Each WaveSlot wraps a signal with plotting parameters
    Slots.append(PltW.WaveSlot(sig,
                               color='blue',  # Line color for the waveform
                               Units='mV',  # Display units (converted from volts)
                               label=sig.name,  # Label displayed in legend
                               Ax=pltAxes[i]))  # Assign to corresponding subplot

# ============================================================================
# Interactive Plot Control
# ============================================================================

# PlotSlots is the main container that manages all plot slots and coordinates visualization
Splt = PltW.PlotSlots(Slots,
                      Fig=pltFig,
                      LiveControl=True)  # Enable interactive controls (sliders, buttons, etc.)

# Plot  1 second of data across all channels
# Time is specified using Quantity objects with units from the 'quantities' package
Splt.PlotChannels((10 * pq.s, 11 * pq.s))

# Add a legend showing channel names
Splt.AddLegend()

#
plt.tight_layout()

# The plot is now interactive and can be explored using the live control interface