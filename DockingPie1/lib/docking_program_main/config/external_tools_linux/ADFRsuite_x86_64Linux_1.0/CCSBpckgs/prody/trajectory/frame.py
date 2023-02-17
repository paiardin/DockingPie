# -*- coding: utf-8 -*-
"""This module defines a class for handling trajectory frames."""

import numpy as np

from prody.measure import getRMSD
from prody.utilities import importLA

__all__ = ['Frame']


class Frame(object):

    """A class for storing trajectory frame coordinates and provide methods
    acting on them."""

    __slots__ = ['_traj', '_index', '_coords', '_unitcell', '_velocs']

    def __init__(self, traj, index, coords, unitcell=None, velocs=None):

        self._traj = traj
        self._index = index

        self._coords = coords
        self._velocs = velocs
        self._unitcell = unitcell

    def __repr__(self):

        sel = ''
        if self._traj._indices is not None:
            sel = 'selected {0} of '.format(self.numSelected())

        return ('<Frame: {0} from {1} ({2}{3} atoms)>'
                ).format(self._index, self._traj.getTitle(), sel,
                         self._traj.numAtoms())

    def __str__(self):

        return 'Frame {0} from {1}'.format(self._index,
                                           self._traj.getTitle())

    def numAtoms(self):
        """Return number of atoms."""

        return self._traj.numAtoms()

    def numSelected(self):
        """Return number of selected atoms."""

        return self._traj.numSelected()

    def getAtoms(self):
        """Return associated/selected atoms."""

        return self._traj.getAtoms()

    def getIndex(self):
        """Return index."""

        return self._index

    def getWeights(self):
        """Return coordinate weights for selected atoms."""

        return self._traj.getWeights()

    def _getWeights(self):

        return self._traj._getWeights()

    def getTrajectory(self):
        """Return the trajectory."""

        return self._traj

    def getCoords(self):
        """Return a copy of coordinates of (selected) atoms."""

        indices = self._traj._indices
        if indices is None:
            return self._getxyz().copy()
        else:
            return self._getxyz()[indices]

    def _getCoords(self):
        """Return coordinates of (selected) atoms."""

        indices = self._traj._indices
        if indices is None:
            return self._getxyz()
        else:
            return self._getxyz()[indices]

    def _getxyz(self):
        """Return coordinates array."""

        ag = self._traj.link()
        if ag is None:
            coords = self._coords
        else:
            coords = ag._getCoords()
        if coords is None:
            raise ValueError('coordinates are not set')
        return coords

    def getVelocities(self):
        """Return a copy of velocities of (selected) atoms."""

        if self._velocs is not None:
            indices = self._traj._indices
            if indices is None:
                return self._velocs.copy()
            else:
                return self._velocs[indices]

    def _getVelocities(self):
        """Return velocities of (selected) atoms."""

        if self._velocs is not None:
            indices = self._traj._indices
            if indices is None:
                return self._velocs
            else:
                return self._velocs[indices]

    def getUnitcell(self):
        """Return a copy of unitcell array."""

        if self._unitcell is not None:
            return self._unitcell.copy()

    def _getUnitcell(self):

        return self._unitcell

    def getDeviations(self):
        """Return deviations from the trajectory reference coordinates."""

        indices = self._traj._indices
        coords = self._getCoords()
        if indices is None:
            return coords - self._traj._coords
        else:
            return coords - self._traj._coords[indices]

    def getRMSD(self):
        """Return RMSD from the trajectory reference coordinates.  If weights
        for the trajectory are set, weighted RMSD will be returned."""

        indices = self._traj._indices
        traj = self._traj
        coords = self._getCoords()
        if indices is None:
            return getRMSD(coords, traj._coords, traj._weights)
        else:
            if traj._weights is None:
                return getRMSD(coords, traj._coords[indices])
            else:
                return getRMSD(coords, traj._coords[indices],
                               traj._weights[indices])

    def superpose(self):
        """Superpose frame onto the trajectory reference coordinates.  Note
        that transformation matrix is calculated based on selected atoms and
        applied to all atoms. If atom weights for the trajectory are set, they
        will be used to calculate the transformation."""

        traj = self._traj
        indices = traj._indices
        mob = mov = self._getxyz()

        weights = traj._weights
        if indices is None:
            tar = traj._coords
            mov = None
        else:
            tar = traj._coords[indices]
            mob = mob[indices]
            if weights is not None:
                weights = weights[indices]

        linalg = importLA()
        if weights is None:
            mob_com = mob.mean(0)
            mob_org = mob - mob_com
            tar_com = tar.mean(0)
            tar_org = tar - tar_com
            matrix = np.dot(tar_org.T, mob_org)
        else:
            weights_sum = weights.sum()
            weights_dot = np.dot(weights.T, weights)
            mob_com = (mob * weights).sum(axis=0) / weights_sum
            mob_org = mob - mob_com
            tar_com = (tar * weights).sum(axis=0) / weights_sum
            tar_org = tar - tar_com
            matrix = np.dot((tar_org * weights).T,
                            (mob_org * weights)) / weights_dot

        U, s, Vh = linalg.svd(matrix)
        Id = np.array([[1, 0, 0],
                       [0, 1, 0],
                       [0, 0, np.sign(linalg.det(matrix))]])
        rotation = np.dot(Vh.T, np.dot(Id, U.T))

        if mov is None:
            np.add(np.dot(mob_org, rotation), tar_com, mob)
        else:
            np.add(np.dot(mov, rotation),
                   (tar_com - np.dot(mob_com, rotation)), mov)
