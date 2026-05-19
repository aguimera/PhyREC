Getting Started
===============

This tutorial introduces basic usage of PhyREC for signal visualization and analysis.

Installation
------------

.. code-block:: bash

    pip install phyrec

Using PhyREC Matplotlib Style
------------------------------

The PhyREC style applies optimized matplotlib settings for electrophysiology data visualization:

- **Minimal visual clutter**: Removes unnecessary axis spines (borders) for cleaner plots
- **Readable labels**: Uses appropriately sized tick labels for multi-channel displays
- **Clean legends**: Removes legend frames for a more professional appearance
- **Consistent formatting**: Ensures uniform styling across all your visualizations

**PhyREC.mplstyle contents:**

.. code-block:: ini

    xtick.labelsize : xx-small
    ytick.labelsize : xx-small
    legend.frameon : False
    legend.loc : lower left

    axes.spines.left   : False
    axes.spines.bottom : False
    axes.spines.top    : False
    axes.spines.right  : False

**To use the style:**

.. code-block:: python

    import matplotlib.pyplot as plt

    # Apply the PhyREC matplotlib style
    plt.style.use('PhyREC.mplstyle')

    # Now create your plots with consistent PhyREC styling
    # ... rest of your plotting code ...

The style can also be applied system-wide by copying ``PhyREC.mplstyle`` to your matplotlib configuration directory.

Key Classes
-----------

- :py:class:`PhyREC.PlotWaves.WaveSlot`: Handles visualization of a single signal channel
- :py:class:`PhyREC.PlotWaves.PlotSlots`: Manages collections of WaveSlot instances with unified controls
- :py:module:`PhyREC.SignalProcess`: Library of signal processing functions (filters, detrending, etc.) that can be applied to signals before plotting

Generating Dummy Signals
------------------------

This reusable code block creates synthetic electrophysiology signals with realistic noise
and phase shifts. This signal generation can be reused across different plotting examples.

.. code-block:: python

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


Basic Plotting Example
----------------------

This example demonstrates visualizing the generated signals using PhyREC's
multi-channel plotting capabilities with interactive controls.

.. code-block:: python

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

.. plot::

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

Basic Signal Filtering Example
-------------------------------------
This example shows how to apply common digital filters to continuous signals using PhyREC.SignalProcess.
It generates synthetic multi‑channel data, then demonstrates a band‑pass filter (to isolate an oscillatory band)
and a band‑stop/notch filter (to remove that band). The example plots raw, band‑passed, and band‑stopped traces
together so you can compare effects visually and verify filter parameters (Freqs, Order). Use LiveControl to navigate
time and adjust the plotted window.

.. code:: python

    import PhyREC.SignalProcess as prspro

    # Note: This code block assumes you have already generated the Sigs list of neo.AnalogSignal objects
    # and set up the plotting figure and axes as shown in the previous example.
    # for i, sig in enumerate(Sigs):
    # Paste the following code inside the loop that iterates over Sigs to apply filters and add new WaveSlots for the filtered signals.

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



.. Plot::

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


