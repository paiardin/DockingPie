# -*- coding: utf-8 -*-
"""This module defines atomic data fields.  You can read this page in
interactive sessions using ``help(fields)``.

Data parsed from PDB and other supported files for these fields are stored in
:class:`.AtomGroup` instances.  Available data fields are listed in the table
below.  :class:`Atomic` classes, such as :class:`.Selection`, offer ``get`` and
``set`` for handling parsed data:

"""

from numpy import array

from prody.utilities import wrapText

from .flags import FIELDS as FLAG_FIELDS

__all__ = ['Field']

READONLY = set()


class Field(object):

    """Atomic data field."""

    __slots__ = ['name', 'dtype',  'doc', 'doc_pl', 'meth', 'meth_pl',
                 'ndim', 'none', 'selstr', 'synonym', 'readonly', 'call',
                 'private', 'depr', 'depr_pl', 'desc', 'flags']

    def __init__(self, name, dtype, **kwargs):

        #: data field name used in atom selections
        self.name = name
        #: data type (primitive Python types)
        self.dtype = dtype
        #: internal variable name used as key for :class:`.AtomGroup` ``_data``
        self.doc = kwargs.get('doc', name)
        #: plural form for documentation
        self.doc_pl = kwargs.get('doc_pl', self.doc + 's')
        #: description of data field, used in documentation
        self.desc = kwargs.get('desc')
        #: expected dimension of the data array
        self.ndim = kwargs.get('ndim', 1)
        #: atomic get/set method name
        self.meth = kwargs.get('meth', name.capitalize())
        #: get/set method name in plural form
        self.meth_pl = kwargs.get('meth_pl', self.meth + 's')
        #: :class:`.AtomGroup` attributes to be set None, when ``setMethod``
        #: is called
        self.none = kwargs.get('none')
        #: list of selection string examples
        self.selstr = kwargs.get('selstr')
        #: deprecated method name
        self.depr = kwargs.get('depr')
        #: deprecated method name in plural form
        self.depr_pl = None
        if self.depr is not None:
            self.depr_pl = kwargs.get('depr_pl', self.depr + 's')
        #: synonym used in atom selections
        self.synonym = kwargs.get('synonym')
        #: read-only attribute without a set method
        self.readonly = kwargs.get('readonly', False)
        #: list of :class:`.AtomGroup` methods to call when ``getMethod`` is
        #: called
        self.call = kwargs.get('call', None)
        #: define only _getMethod for :class:`.AtomGroup` to be used by
        #: :class:`.Select` class
        self.private = kwargs.get('private', False)
        #: **True** when there are flags associated with the data field
        self.flags = kwargs.get('flags', False)

        if self.readonly:
            READONLY.add(self.name)

    def getDocstr(self, meth, plural=True, selex=True):
        """Return documentation string for the field."""

        assert meth in ('set', 'get', '_get'), "meth must be 'set' or 'get'"
        assert isinstance(plural, bool), 'plural must be a boolean'
        assert isinstance(selex, bool), 'selex must be a boolean'

        if meth == 'get':
            if plural:
                docstr = 'Return a copy of {0}.'.format(self.doc_pl)
            else:
                docstr = 'Return {0} of the atom.'.format(self.doc)
        elif meth == 'set':
            if plural:
                docstr = 'Set {0}.'.format(self.doc_pl)
            else:
                docstr = 'Set {0} of the atom.'.format(self.doc)
        else:
            selex = False
            if plural:
                docstr = 'Return {0} array.'.format(self.doc_pl)

        if self.desc:
            docstr += '  ' + self.desc

        selstr = self.selstr
        if selex and selstr:
            if plural:
                doc = self.doc_pl
            else:
                doc = self.doc
            if '(' in doc:
                doc = doc[:doc.index('(')]
            selex = "``, ``".join([repr(s) for s in selstr])
            selex = ("  {0} can be used in atom selections, e.g. "
                     "``{1}``.").format(doc.capitalize(), selex)
            if self.synonym is not None:
                selex = selex + ('  Note that *{0}* is a synonym for '
                                 '*{1}*.').format(self.synonym, self.name)
            return wrapText(docstr + selex)
        else:
            return wrapText(docstr)

HVNONE = ['_hv', 'segindex', 'chindex', 'resindex']

DTYPE = array(['a']).dtype.char  # 'S' for PY2K and 'U' for PY3K

ATOMIC_FIELDS = {
    'name':      Field('name', DTYPE + '6', selstr=('name CA CB',)),
    'type':      Field('type', DTYPE + '6', selstr=('type N.am',)),
    'altloc':    Field('altloc', DTYPE + '1',
                       doc='alternate location indicator',
                       selstr=('altloc A B', 'altloc _'),),
    'anisou':    Field('anisou', float, doc='anisotropic temperature factor',
                       ndim=2),
    'chain':     Field('chain', DTYPE + '2',  doc='chain identifier',
                       meth='Chid', none=HVNONE, synonym='chid',
                       selstr=('chain A', 'chid A B C', 'chain _')),
    'element':   Field('element', DTYPE + '2', doc='element symbol',
                       selstr=('element C O N',)),
    'occupancy': Field('occupancy', float,
                       doc='occupancy value', meth_pl='Occupancies',
                       selstr=('occupancy 1', 'occupancy > 0')),
    'resname':   Field('resname', DTYPE + '6', doc='residue name',
                       selstr=('resname ALA GLY',)),
    'resnum':    Field('resnum', int, doc='residue number', none=HVNONE,
                       selstr=('resnum 1 2 3', 'resnum 120A 120B',
                               'resnum 10 to 20', 'resnum 10:20:2',
                               'resnum < 10'), synonym='resid'),
    'secondary': Field('secondary', DTYPE + '1',
                       doc='secondary structure assignment',
                       meth='Secstr', synonym='secstr',
                       selstr=('secondary H E', 'secstr H E'),
                       ),
    'segment':   Field('segment', DTYPE + '6', doc='segment name',
                       meth='Segname',
                       selstr=('segment PROT', 'segname PROT'),
                       synonym='segname', none=HVNONE),
    'siguij':    Field('siguij', float, doc='standard deviations for '
                       'anisotropic temperature factor',
                       meth='Anistd', ndim=2),
    'serial':    Field('serial', int, doc='serial number (from file)',
                       doc_pl='serial numbers (from file)', none=['_sn2i'],
                       selstr=('serial 1 2 3', 'serial 1 to 10',
                       'serial 1:10:2', 'serial < 10')),
    'beta':      Field('beta', float, doc='β-value (temperature factor)',
                       doc_pl='β-values (or temperature factors)',
                       selstr=('beta 555.55', 'beta 0 to 500', 'beta 0:500',
                       'beta < 500')),
    'icode':     Field('icode', DTYPE + '1', doc='insertion code', none=HVNONE,
                       selstr=('icode A', 'icode _')),
    'type':      Field('type', DTYPE + '6', selstr=('type CT1 CT2 CT3',)),
    'charge':    Field('charge', float, doc='partial charge',
                       selstr=('charge 1', 'abs(charge) == 1', 'charge < 0')),
    'mass':      Field('mass', float, doc_pl='masses',
                       meth_pl='Masses', selstr=('12 <= mass <= 13.5',)),
    'radius':    Field('radius', float, doc='radius',
                       doc_pl='radii', meth_pl='Radii',
                       selstr=('radii < 1.5', 'radii ** 2 < 2.3')),
    'resindex':  Field('resindex', int, doc='residue index',
                       doc_pl='residue indices', meth_pl='Resindices',
                       selstr=('resindex 0',), readonly=True,
                       call=['getHierView'],
                       desc='Residue indices are assigned to subsets of atoms '
                            'with distinct sequences of residue number, '
                            'insertion code, chain identifier, and segment '
                            'name.  Residue indices start from zero, are '
                            'incremented by one, and are assigned in the '
                            'order of appearance in :class:`.AtomGroup` '
                            'instance.'),
    'chindex':   Field('chindex', int, doc='chain index',
                       doc_pl='chain indices', meth_pl='Chindices',
                       selstr=('chindex 0',), readonly=True,
                       call=['getHierView'],
                       desc='Chain indices are assigned to subsets of atoms '
                            'with distinct pairs of chain identifier and '
                            ' segment name.  Chain indices start from zero, '
                            'are incremented by one, and are assigned in the '
                            'order of appearance in :class:`.AtomGroup` '
                            'instance.'),
    'segindex':  Field('segindex', int, doc='segment index',
                       doc_pl='segment indices', meth_pl='Segindices',
                       selstr=['segindex 0'], readonly=True,
                       call=['getHierView'],
                       desc='Segment indices are assigned to subsets of atoms '
                            'with distinct segment names.  Segment indices '
                            'start from zero, are incremented by one, and are '
                            'assigned in the order of appearance in '
                            ':class:`.AtomGroup` instance.'),
    'fragindex': Field('fragindex', int,
                       doc='fragment index', doc_pl='fragment indices',
                       meth_pl='Fragindices',
                       selstr=['fragindex 0', 'fragment 1'],
                       readonly=True, call=['_fragment'], synonym='fragment',
                       desc='Fragment indices are assigned to connected '
                            'subsets of atoms.  Bonds needs to be set using '
                            ':meth:`.AtomGroup.setBonds` method.  Fragment '
                            'indices start from zero, are incremented by '
                            'one, and are assigned in the order of appearance '
                            'in :class:`.AtomGroup` instance.'),
    'numbonds':  Field('numbonds', int, meth_pl='Numbonds',
                       doc='number of bonds',
                       selstr=['numbonds 0', 'numbonds 1'],
                       readonly=True, private=True),
}


__doc__ += """

Many of these data fields can be used to make :ref:`selections`.
Following table lists definitions of fields and selection examples.  Note that
fields noted as *read only* do not have a ``set`` method.

.. glossary::
"""

keys = list(ATOMIC_FIELDS.keys())
keys.sort()


for key in keys:

    field = ATOMIC_FIELDS[key]

    __doc__ += '\n\n   ' + key + '\n'

    if field.synonym:
        __doc__ += '   ' + field.synonym + '\n'

    __doc__ += wrapText(field.doc + (' *(read only)*' if field.readonly
                                     else ''), indent=6)

    if field.selstr:
        sel = "\n      ``'" + "'``,\n      ``'".join(field.selstr) + "'``"
        __doc__ += '\n\n      *E.g.:* ' + sel

for key in FLAG_FIELDS.keys():
    ATOMIC_FIELDS[key].flags = True


def wrapGetMethod(fn):
    def getMethod(self):
        return fn(self)
    return getMethod


def wrapSetMethod(fn):
    def setMethod(self, data):
        return fn(self, data)
    return setMethod
