# -*- coding: utf-8 -*-
"""This module defines a base class for format specific trajectory classes."""

from os.path import isfile, abspath, split, splitext

import numpy as np

from prody import LOGGER
from prody.utilities import relpath

from .trajbase import TrajBase

__all__ = ['TrajFile']


class TrajFile(TrajBase):

    """A base class for trajectory file classes:

      * :class:`.DCDFile`"""


    def __init__(self, filename, mode='r'):
        """Open *filename* for reading (default, ``mode="r"``), writing
        (``mode="w"``), or appending (``mode="r+"`` or ``mode="a"``).
        Binary mode option will be appended automatically."""

        if not isinstance(filename, str):
            raise TypeError("filename argument must be a string")
        if not isinstance(mode, str):
            TypeError('mode argument must be string')
        if not mode in ('r', 'w', 'a', 'r+'):
            ValueError("mode string must begin with one of 'r', 'w', 'r+', or "
                       "'a'")
        name = splitext(split(filename)[1])[0]
        TrajBase.__init__(self, name)
        self._file = None
        if mode == 'r' and not isfile(filename):
            raise IOError("[Errno 2] No such file or directory: '{0}'"
                          .format(filename))
        self._filename = filename
        if mode in ('a', 'r+'):
            self._file = open(filename, 'r+b')
            self._file.seek(0)
            mode = 'a'
        else:
            if not mode.endswith('b'):
                mode += 'b'
            self._file = open(filename, mode)

        self._mode = mode
        self._bytes_per_frame = None
        self._first_byte = None
        self._dtype = np.float32

        self._timestep = 1
        self._first_ts = 0
        self._framefreq = 1
        self._n_fixed = 0

    def __del__(self):

        if self._file is not None:
            self._file.close()

    def __repr__(self):

        if self._closed:
            return ('<{0}: {1} (closed)>').format(
                        self.__class__.__name__, self._title)

        link = ''
        if self._ag is not None:
            link = 'linked to ' + str(self._ag) + '; '

        if self._mode.startswith('r'):
            next = 'next {0} of {1} frames; '.format(self._nfi,
                                                       self._n_csets)
        else:
            next = '{0} frames written; '.format(self._n_csets)

        if self._indices is None:
            atoms = '{0} atoms'.format(self._n_atoms)
        else:
            atoms = 'selected {0} of {1} atoms'.format(
                                        self.numSelected(), self._n_atoms)

        return '<{0}: {1} ({2}{3}{4})>'.format(
                   self.__class__.__name__, self._title, link, next, atoms)

    def getFilename(self, absolute=False):
        """Return relative path to the current file. For absolute path,
        pass ``absolute=True`` argument."""

        if absolute:
            return abspath(self._filename)
        return relpath(self._filename)

    def getFrame(self, index):
        """Return frame at given *index*."""

        if self._closed:
            raise ValueError('I/O operation on closed file')
        if not isinstance(index, int):
            raise IndexError('index must be an integer')
        if not 0 <= index < self._n_csets:
            raise IndexError('index must be greater or equal to 0 and less '
                             'than number of frames')
        nfi = self._nfi
        if index > nfi:
            self.skip(index - nfi)
        elif index < nfi:
            self.reset()
            if index > 0:
                self.skip(index)
        return next(self)

    getFrame.__doc__ = TrajBase.getFrame.__doc__

    def getCoordsets(self, indices=None):

        if self._closed:
            raise ValueError('I/O operation on closed file')
        if indices is None:
            indices = np.arange(self._n_csets)
        elif isinstance(indices, int):
            indices = np.array([indices])
        elif isinstance(indices, slice):
            indices = np.arange(*indices.indices(self._n_csets))
            indices.sort()
        elif isinstance(indices, (list, np.ndarray)):
            indices = np.unique(indices)
        else:
            raise TypeError('indices must be an integer or a list of integers')

        nfi = self._nfi
        self.reset()

        n_atoms = self.numSelected()
        coords = np.zeros((len(indices), n_atoms, 3), self._dtype)

        prev = 0
        next = self.nextCoordset
        for i, index in enumerate(indices):
            diff = index - prev
            if diff > 1:
                self.skip(diff-1)
            xyz = next()
            if xyz is None:
                LOGGER.warning('Expected {0} frames, but parsed {1}.'
                               .format(len(indices), i))
                self.goto(nfi)
                return coords[:i]
            coords[i] = xyz
            prev = index

        self.goto(nfi)
        return coords

    getCoordsets.__doc__ = TrajBase.getCoordsets.__doc__

    def skip(self, n):

        if self._closed:
            raise ValueError('I/O operation on closed file')
        if not isinstance(n, int):
            raise ValueError('n must be an integer')
        if n > 0:
            left = self._n_csets - self._nfi
            if n > left:
                n = left
            self._file.seek(n * self._bytes_per_frame, 1)
            self._nfi += n

    skip.__doc__ = TrajBase.skip.__doc__

    def goto(self, n):

        if self._closed:
            raise ValueError('I/O operation on closed file')
        if not isinstance(n, int):
            raise ValueError('n must be an integer')
        n_csets = self._n_csets
        if n == 0:
            self.reset()
        else:
            if n < 0:
                n = n_csets + n
            if n < 0:
                n = 0
            elif n > n_csets:
                n = n_csets
            self._file.seek(self._first_byte + n * self._bytes_per_frame)
            self._nfi = n

    goto.__doc__ = TrajBase.goto.__doc__

    def reset(self):

        if self._closed:
            raise ValueError('I/O operation on closed file')
        self._file.seek(self._first_byte)
        self._nfi = 0

    reset.__doc__ = TrajBase.reset.__doc__

    def close(self):

        self._file.close()
        self._nfi = 0
        self._closed = True

    close.__doc__ = TrajBase.close.__doc__

    def getTimestep(self):
        """Return timestep size."""

        return self._timestep

    def getFirstTimestep(self):
        """Return first timestep value."""

        return self._first_ts

    def getFrameFreq(self):
        """Return timesteps between frames."""

        return self._framefreq

    def numFixed(self):
        """Return number of fixed atoms."""

        return self._n_fixed
