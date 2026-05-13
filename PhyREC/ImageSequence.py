#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 17 18:42:32 2023

@author: aguimera
"""

import quantities as pq
import numpy as np
from copy import deepcopy

class ImageSequence(pq.Quantity):

    def __new__(cls, signal, sampling_rate,  t_start=0*pq.s, units=None, dtype=None, copy=True,

                name=None, **kwargs):

        signal = cls._rescale(cls, signal=signal, units=units)
        obj = pq.Quantity(signal, units=units,
                          dtype=dtype, copy=copy).view(cls)

        obj.sampling_rate = sampling_rate
        obj.t_start = t_start
        obj.name = name
        obj.annotations = {}
        return obj

    def __array_finalize__(self, obj):

        """

        Finalize the array creation, copying attributes from the original object.



        Parameters

        ----------

        obj : object

            The object from which to copy attributes.

        """

        if obj is None:

            return
        super(ImageSequence, self).__array_finalize__(obj)
        self.sampling_rate = getattr(obj, 'sampling_rate', None)
        self.t_start = getattr(obj, 't_start', None)
        self.name = getattr(obj, 'name', None)
        self.annotations = getattr(obj, 'annotations', {})

    def __repr__(self):
        return("{cls} {frame} frames with {width} pixels of width and {height} pixels of height \n"
               "units {units} \n"
               "sampling_rate {sampling_rate} \n"
               "Duration {t_start} --  {t_stop}\n"
               "datatype {dtype} \n"
               "name {name} ".format(cls=self.__class__.__name__,
                                     frame=self.shape[0],
                                     height=self.shape[1],
                                     width=self.shape[2],
                                     units=self.units,
                                     dtype=self.dtype,
                                     sampling_rate=self.sampling_rate,
                                     t_start=self.t_start,
                                     t_stop=self.t_stop,
                                     name=self.name))

    def __deepcopy__(self, memo):
        """
        Create a deep copy of the data object.

        All attributes and annotations are also deep copied.
        References to parent objects are not kept, they are set to None.

        Parameters
        ----------
        memo : dict
            Objects that have been deep copied already.

        Returns
        -------
        ImageSequence
            Deep copy of the input ImageSequence object.
        """
        cls = self.__class__
        necessary_attrs = {'signal': self,
                           'units': self.units,
                           'sampling_rate': self.sampling_rate,
                           'name': self.name,
                           't_start': self.t_start}
        # Create object using constructor with necessary attributes
        new_obj = cls(**necessary_attrs)
        # Add all attributes

        return new_obj

    def _rescale(self, signal, units=None):
        '''
        Check that units are present, and rescale the signal if necessary.
        This is called whenever a new signal is
        created from the constructor. See :meth:`__new__' in
        :class:`AnalogSignal` and :class:`IrregularlySampledSignal`
        '''
        if units is None:
            if not hasattr(signal, "units"):
                raise ValueError("Units must be specified")
        elif isinstance(signal, pq.Quantity):
            # This test always returns True, i.e. rescaling is always executed if one of the units
            # is a pq.CompoundUnit. This is fine because rescaling is correct anyway.
            if pq.quantity.validate_dimensionality(units) != signal.dimensionality:
                signal = signal.rescale(units)
        return signal

    def rescale(self, units):

        """

        Rescale the signal to new units.



        Parameters

        ----------

        units : str or quantity

            The new units to rescale to.



        Returns

        -------

        ImageSequence

            The rescaled ImageSequence.

        """

        obj = super(ImageSequence, self).rescale(units)

        obj.sampling_rate = self.sampling_rate
        obj.t_start = self.t_start
        obj.name = self.name
        obj.annotations = self.annotations
        return obj.view(ImageSequence)

    def time_index(self, t):
        """
        Return the array index corresponding to the time `t`.

        Parameters
        ----------
        t : quantities.Quantity
            Time value to convert to array index.

        Returns
        -------
        int
            Array index corresponding to the time value.
        """
        i = (t - self.t_start) * self.sampling_rate
        i = int(np.rint(i.simplified.magnitude))
        return i

    def time_slice(self, t_start, t_stop):
        """
        Create a new ImageSequence corresponding to the time slice of the original.

        Creates a time slice between t_start and t_stop. For numerical stability,
        if t_start does not fall exactly on time bins defined by the sampling_period,
        it will be rounded to the nearest sampling bin. The time bin for t_stop will
        be chosen to make the duration of the resultant signal as close as possible to
        t_stop - t_start. This means for a given duration, the slice size is always the same.

        Parameters
        ----------
        t_start : quantities.Quantity
            Start time of the slice. If None, uses the beginning of signal.
        t_stop : quantities.Quantity
            Stop time of the slice. If None, uses the end of signal.

        Returns
        -------
        ImageSequence
            Time-sliced ImageSequence object.

        Raises
        ------
        ValueError
            If t_start or t_stop are outside the signal duration.
        """

        # checking start time and transforming to start index
        if t_start is None:
            i = 0
            t_start = 0 * pq.s
        else:
            i = self.time_index(t_start)

        # checking stop time and transforming to stop index
        if t_stop is None:
            j = len(self)
        else:
            delta = (t_stop - t_start) * self.sampling_rate
            j = i + int(np.rint(delta.simplified.magnitude))

        if (i < 0) or (j > len(self)):
            raise ValueError('t_start, t_stop have to be within the analog \
                              signal duration')

        # Time slicing should create a deep copy of the object
        obj = deepcopy(self[i:j, :, :])

        obj.t_start = self.t_start + i * self.sampling_period

        return obj

    @property
    def duration(self):
        """
        Signal duration.

        Computed as size * sampling_period.

        Returns
        -------
        quantities.Quantity
            Duration of the signal.
        """
        return self.shape[0] / self.sampling_rate

    @property
    def t_stop(self):
        """
        Time when the signal ends.

        Returns
        -------
        quantities.Quantity
            t_start + duration.
        """
        return self.t_start + self.duration

    @property
    def times(self):
        """
        The time points of each sample of the signal.

        Computed as t_start + arange(shape) / sampling_rate.

        Returns
        -------
        quantities.Quantity
            Array of time points for each sample.
        """
        return self.t_start + np.arange(self.shape[0]) / self.sampling_rate

    @property
    def sampling_period(self):
        """
        Interval between two samples.

        Returns
        -------
        quantities.Quantity
            1 / sampling_rate.
        """
        return 1. / self.sampling_rate

    def annotate(self, **annotations):
        """
        Add annotations (non-standardized metadata) to an ImageSequence object.

        Parameters
        ----------
        **annotations
            Arbitrary keyword arguments to store as annotations.

        Returns
        -------
        None

        Examples
        --------
        >>> obj.annotate(key1=value0, key2=value1)
        >>> obj.key2
        value1
        """
        self.annotations.update(annotations)

    def __array_ufunc__(self, ufunc, method, *inputs, out=None, **kwargs):

        """

        Handle universal functions applied to the array.



        This method ensures that operations on ImageSequence preserve the metadata.



        Parameters

        ----------

        ufunc : callable

            The universal function.

        method : str

            The method.

        *inputs

            Input arrays.

        out : array, optional

            Output array.

        **kwargs

            Additional arguments.



        Returns

        -------

        ImageSequence or NotImplemented

            The result or NotImplemented.

        """

        args = []
        in_no = []
        for i, input_ in enumerate(inputs):
            if isinstance(input_, ImageSequence):
                in_no.append(i)
                args.append(input_.view(pq.Quantity))
            else:
                args.append(input_)

        results = super(ImageSequence, self).__array_ufunc__(
            ufunc, method, *args, **kwargs)
        if results is NotImplemented:
            return NotImplemented

        results = results.view(ImageSequence)
        results.sampling_rate = inputs[0].sampling_rate
        results.t_start = inputs[0].t_start
        results.name = inputs[0].name
        results.annotations = inputs[0].annotations

        return results

        # if method == 'at':
        #     if isinstance(inputs[0], A):
        #         inputs[0].info = info
        #     return

        # if ufunc.nout == 1:
        #     results = (results,)

        # results = tuple((np.asarray(result).view(A)
        #                  if output is None else output)
        #                 for result, output in zip(results, outputs))
        # if results and isinstance(results[0], A):
        #     results[0].info = info

        # return results[0] if len(results) == 1 else results
    def duplicate_with_new_data(self, signal, units=None):
        """
        Create a new signal with the same metadata but different data.

        Required attributes of the signal are used. Note that array annotations
        cannot be copied here because the length of data can change.

        Parameters
        ----------
        signal : array-like
            New signal data for the duplicated object.
        units : quantities.Quantity, optional
            Units for the new signal. If None, uses self.units.

        Returns
        -------
        ImageSequence
            New ImageSequence with same metadata but new data.

        Notes
        -----
        Array annotations are not copied here because it is not ensured that the
        same number of samples is used and they would possibly make no sense
        when combined with another signal.
        """
        if units is None:
            units = self.units
        # else:
        #     units = pq.quantity.validate_dimensionality(units)

        # signal is the new signal
        necessary_attrs = {'signal': signal,
                           'units': units,
                           'sampling_rate': self.sampling_rate,
                           'name': self.name,
                           't_start': self.t_start}
        new = self.__class__(**necessary_attrs)
        new.annotations.update(self.annotations)
        # Note: Array annotations are not copied here, because it is not ensured
        # that the same number of signals is used and they would possibly make no sense
        # when combined with another signal
        return new
