# -*- coding: utf-8 -*-
"""This module defines customized gamma functions for elastic network model
analysis."""

import numpy as np

from prody.atomic import Atomic

__all__ = ['Gamma', 'GammaStructureBased', 'GammaVariableCutoff']


class Gamma(object):

    """Base class for facilitating use of atom type, residue type, or residue
    property dependent force constants (γ).

    Derived classes:

    * :class:`.GammaStructureBased`
    * :class:`.GammaVariableCutoff`"""

    def __init__(self):
        pass

    def gamma(self, dist2, i, j):
        """Return force constant.

        For efficiency purposes square of the distance between interacting
        atom/residue (node) pairs is passed to this function. In addition,
        node indices are passed."""

        pass


class GammaStructureBased(Gamma):

    """Facilitate setting the spring constant based on the secondary structure
    and connectivity of the residues.

    A recent systematic study [LT10]_ of a large set of NMR-structures analyzed
    using a method based on entropy maximization showed that taking into
    consideration properties such as sequential separation between
    contacting residues and the secondary structure types of the interacting
    residues provides refinement in the ENM description of proteins.

    This class determines pairs of connected residues or pairs of proximal
    residues in a helix or a sheet, and assigns them a larger user defined
    spring constant value.

     DSSP single letter abbreviations are recognized:
       * **H**: α-helix
       * **G**: 3-10-helix
       * **I**: π-helix
       * **E**: extended part of a sheet

    *helix*:
        Applies to residue (or Cα atom) pairs that are in the same helical
        segment, at most 7 Å apart, and separated by at most
        3 (3-10-helix), 4 (α-helix), or 5 (π-helix) residues.

    *sheet*:
        Applies to Cα atom pairs that are in different β-strands and at most
        6 Å apart.

    *connected*:
        Applies to Cα atoms that are at most 4 Å apart.

    Note that this class does not take into account insertion codes.

    .. [LT10] Lezon TR, Bahar I. Using entropy maximization to understand the
       determinants of structural dynamics beyond native contact topology.
       *PLoS Comput Biol* **2010** 6(6):e1000816.

    **Example**:

    Let's parse coordinates and header data from a PDB file, and then
    assign secondary structure to the atoms.

    .. ipython:: python

       from prody import *
       ubi, header = parsePDB('1aar', chain='A', subset='calpha', header=True)
       assignSecstr(header, ubi)

    In the above we parsed only the atoms needed for this calculation, i.e.
    Cα atoms from chain A.

    We build the Hessian matrix using structure based force constants as
    follows;

    .. ipython:: python

       gamma = GammaStructureBased(ubi)
       anm = ANM('')
       anm.buildHessian(ubi, gamma=gamma)

    We can obtain the force constants assigned to residue pairs from the
    Kirchhoff matrix as follows:

    .. ipython:: python

       k = anm.getKirchhoff()
       k[0,1] # a pair of connected residues
       k[0,16] # a pair of residues from a sheet"""

    def __init__(self, atoms, gamma=1.0, helix=6.0, sheet=6.0, connected=10.0):
        """Setup the parameters.

        :arg atoms: A set of atoms with chain identifiers, residue numbers,
            and secondary structure assignments are set.
        :type atoms: :class:`.Atomic`

        :arg gamma: Force constant in arbitrary units. Default is 1.0.
        :type gamma: float

        :arg helix: Force constant factor for residues hydrogen bonded in
            α-helices, 3,10-helices, and π-helices. Default is 6.0, i.e.
            ``6.0`*gamma``.
        :type helix: float

        :arg sheet: Force constant factor for residue pairs forming a hydrogen
            bond in a β-sheet. Default is 6.0, i.e. ``6.0`*gamma``.
        :type sheet: float

        :arg connected: Force constant factor for residue pairs that are
            connected. Default is 10.0, i.e. ``10.0`*gamma``.
        :type connected: float"""

        if not isinstance(atoms, Atomic):
            raise TypeError('atoms must be an Atomic instance')
        n_atoms = atoms.numAtoms()
        if n_atoms < 3:
            raise ValueError('number of atoms must be larger than 2')
        sstr = atoms.getSecstrs()
        assert sstr is not None, 'secondary structure assignments must be set'
        chid = atoms.getChids()
        assert chid is not None, 'chain identifiers must be set'
        rnum = atoms.getResindices()
        assert rnum is not None, 'residue numbers must be set'
        gamma = float(gamma)
        assert gamma > 0, 'gamma must be greater than 0'
        helix = float(helix)
        assert helix > 0, 'helix must be greater than 0'
        sheet = float(sheet)
        assert sheet > 0, 'sheet must be greater than 0'
        connected = float(connected)
        assert connected > 0, 'connected must be greater than 0'

        ssid = np.zeros(n_atoms)
        for i in range(1, n_atoms):
            if (sstr[i-1] == sstr[i] and chid[i-1] == chid[i] and
               rnum[i]-rnum[i-1] == 1):
                ssid[i] = ssid[i-1]
            else:
                ssid[i] = ssid[i-1] + 1
        self._sstr = sstr
        self._chid = chid
        self._rnum = rnum
        self._ssid = ssid
        self._gamma = gamma
        self._helix = gamma * helix
        self._sheet = gamma * sheet
        self._connected = gamma * connected

    def getSecstrs(self):
        """Return a copy of secondary structure assignments."""

        return self._sstr.copy()

    def getChids(self):
        """Return a copy of chain identifiers."""

        return self._chid.socopypy()

    def getResnums(self):
        """Return a copy of residue numbers."""

        return self._rnum.copy()

    def gamma(self, dist2, i, j):
        """Return force constant."""

        if dist2 <= 16:
            return self._connected
        sstr = self._sstr
        ssid = self._ssid
        rnum = self._rnum
        # if residues are in the same secondary structure element
        if ssid[i] == ssid[j]:
            i_j = abs(rnum[j] - rnum[i])
            if ((i_j <= 4 and sstr[i] == 'H') or
                    (i_j <= 3 and sstr[i] == 'G') or
                    (i_j <= 5 and sstr[i] == 'I')) and dist2 <= 49:
                return self._helix
        elif sstr[i] == sstr[j] == 'E' and dist2 <= 36:
            return self._sheet

        return self._gamma


class GammaVariableCutoff(Gamma):

    """Facilitate setting the cutoff distance based on user defined
    atom/residue (node) radii.

    Half of the cutoff distance can be thought of as the radius of a node.
    This class enables setting different radii for different node types.

    **Example**:

    Let's think of a protein-DNA complex for which we want to use different
    radius for different residue types. Let's say, for protein Cα atoms we
    want to set the radius to 7.5 Å, and for nucleic acid phosphate atoms to
    10 Å. We use the HhaI-DNA complex structure :file:`1mht`.

    .. ipython:: python

       hhai = parsePDB('1mht')
       ca_p = hhai.select('(protein and name CA) or (nucleic and name P)')
       ca_p.getNames()

    We set the radii of atoms:

    .. ipython:: python

       varcutoff = GammaVariableCutoff(ca_p.getNames(), gamma=1,
           default_radius=7.5, debug=False, P=10)
       varcutoff.getRadii()

    The above shows that for phosphate atoms radii is set to 10 Å, because
    we passed the ``P=10`` argument.  As for Cα atoms, the default 7.5 Å
    is set as the radius (``default_radius=7.5``).  You can also try this with
    ``debug=True`` argument to print debugging information on the screen.

    We build :class:`.ANM` Hessian matrix as follows:

    .. ipython:: python

       anm = ANM('HhaI-DNA')
       anm.buildHessian(ca_p, gamma=varcutoff, cutoff=20)

    Note that we passed ``cutoff=20.0`` to the :meth:`.ANM.buildHessian`
    method.  This is equal to the largest possible cutoff distance (between
    two phosphate atoms) for this system, and ensures that all of the
    potential interactions are evaluated.

    For pairs of atoms for which the actual distance is larger than the
    effective cutoff, the :meth:`.GammaVariableCutoff.gamma` method returns
    ``0``.  This annuls the interaction between those atom pairs."""

    def __init__(self, identifiers, gamma=1., default_radius=7.5, **kwargs):
        """Set the radii of atoms.

        :arg identifiers: List of atom names or types, or residue names.
        :type identifiers: list or :class:`numpy.ndarray`

        :arg gamma: Uniform force constant value. Default is 1.0.
        :type gamma: float

        :arg default_radius: Default radius for atoms whose radii is not set
            as a keyword argument. Default is 7.5
        :type default_radius: float

        Keywords in keyword arguments must match those in *atom_identifiers*.
        Values of keyword arguments must be :class:`float`."""

        self._identifiers = identifiers
        radii = np.ones(len(identifiers)) * default_radius

        for i, identifier in enumerate(identifiers):
            radii[i] = kwargs.get(identifier, default_radius)
        self._radii = radii
        self._gamma = float(gamma)
        self._debug = bool(kwargs.get('debug', False))

    def getRadii(self):
        """Return a copy of radii array."""

        return self._radii.copy()

    def getGamma(self):
        """Return the uniform force constant value."""

        return self._gamma

    def gamma(self, dist2, i, j):
        """Return force constant."""

        cutoff = (self._radii[i] + self._radii[j])
        cutoff2 = cutoff ** 2

        if dist2 < cutoff2:
            gamma = self._gamma
        else:
            gamma = 0
        if self._debug:
            print(' '.join([self._identifiers[i] + '_' + str(i), '--',
                  self._identifiers[j] + '_' + str(j),
                  'effective cutoff:', str(cutoff), 'distance:',
                  str(dist2**0.5), 'gamma:', str(gamma)]))  # PY3K: OK
        return gamma
