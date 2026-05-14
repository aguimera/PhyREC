Getting Started
===============

This tutorial introduces basic usage of PhyREC for signal visualization.

Installation
------------

.. code-block:: bash

    pip install phyrec

Basic Example
-------------

.. code-block:: python

    import PhyREC as pr

    slot = pr.PlotWaves.WaveSlot(signal)
    slot.PlotSignal((0, 5))

Explanation
-----------

This example creates a waveform visualization slot and plots the signal
between 0 and 5 seconds.

Next Steps
----------

See :py:meth:`PhyREC.PlotWaves.WaveSlot.CalcAvarage` for averaging signals.
``