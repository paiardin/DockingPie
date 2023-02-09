################################################################################
##
## This library is free software; you can redistribute it and/or
## modify it under the terms of the GNU Lesser General Public
## License as published by the Free Software Foundation; either
## version 2.1 of the License, or (at your option) any later version.
## 
## This library is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
## Lesser General Public License for more details.
## 
## You should have received a copy of the GNU Lesser General Public
## License along with this library; if not, write to the Free Software
## Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301 USA
##
## (C) Copyrights Dr. Michel F. Sanner and TSRI 2016
##
################################################################################

#############################################################################
#
# Author: Michel F. SANNER
#
# Copyright: M. Sanner TSRI 2015
#
#
#############################################################################

#
# $Header: /mnt/raid/services/cvs/MolKit2/selection.py,v 1.4.2.1 2017/07/26 22:03:40 annao Exp $
# 
# $Id: selection.py,v 1.4.2.1 2017/07/26 22:03:40 annao Exp $
#
import re, weakref
from time import time

from numpy import concatenate, unique, transpose, concatenate, dot, ones, array, arange

import prody
from prody.atomic.selection import Selection as ProdySelection
from prody import LOGGER
from .tree import TreeObject

from mglutil.util.io import Stream

## the selectionTypes dict is used in PmvGUI to build selection menus
# The key will create a selection submneu
# the value is the (selection argument, menu entry, minimum level, doc string) 
selectionTypes = {
    'Type' : [
        ('protein', 'protein', 2, 'atoms belonging to amino acids'),
        ('stdaa', 'standard amino acids', 2,
         'standard amino acid residues: ALA, ARG, ASN, ASP, CYS, GLN, GLU, GLY, HIS, ILE, LEU, LYS, MET, PHE, PRO, SER, THR, TRP, TYR, and VAL'),
        ('nonstdaa', 'Non standard amino acids', 2,
         'residues: ASX, GLX, CSO, HIP, HSD, HSE, HSP, MSE, SEC, SEP, TPO, PTR, XLE, XAA'),
        ('nucleic', 'nucleic', 2, 'atoms belonging to nucleic acids'),
        'Separator',
        ('hetero', 'hetero', 2, 'atoms other than a protein or a nucleic residue'),
        ('water', 'water', 2, 'atoms in residues named HOH, DOD, WAT, TIP3, H2O, OH2, TIP, TIP2, and TIP4'),
        ('ion', 'ion', 2, 'atoms AL BA CA CD CL CO CS CU CU1 CUA HG IN IOD K MG MN3 NA PB PT RB TB TL WO4 YB ZN CAL CES CLA POT SOD ZN2'),
        ('lipid', 'lipid', 2, 'residues GPE, LPP, OLA, SDS, and STE from PDB, and POPC, LPPC, POPE, DLPE, PCGL, STEA, PALM, OLEO, DMPC from CHARMM force field'),
        ('sugar', 'sugar', 2, 'residues BGC, GLC, and GLO from PDB, and also AGLC from CHARMM'),
        ('heme', 'heme', 2, 'residues 1FH 2FH DDH DHE HAS HDD HDE HDM HEA HEB HEC HEM HEO HES HEV NTE SRM and VER from PDB, and HEMO and HEMR from CHARMM'),
        ],
    'Residue type' : [
        ('acidic', 'acidic', 2, 'residues ASP, GLU, HSP, PTR, SEP, TPO'),
        ('aliphatic', 'aliphatic', 2, 'residues ALA, GLY, ILE, LEU, PRO, VAL, XLE'),
        ('acyclic', 'acyclic', 2, 'residues ALA, ARG, ASN, ASP, ASX, CSO, CYS, GLN, GLU, GLX, GLY, ILE, EU, LYS, MET, MSE, SEC, SEP, SER, THR, TPO, VAL, XLE'),
        ('aromatic', 'aromatic', 2, 'residues HIS, PHE, PTR, TRP, TYR'),
        ('basic', 'basic', 2, 'residues ARG, HIP, HIS, HSD, HSE, LYS'),
        ('charged', 'charged', 2, 'residues ARG, ASP, GLU, HIS, LYS'),
        ('cyclic', 'cyclic', 2, 'residues HIP, HIS, HSD, HSE, HSP, PHE, PRO, PTR, TRP, TYR'),
        ('hydrophobic', 'hydrophobic', 2, 'residues ALA, ILE, LEU, MET, PHE, PRO, TRP, VAL, XLE'),
        'Separator',
        ('nucleobase', 'nucleobase', 2, 'residues ADE (adenine), GUN (guanine), CYT (cytosine), THY (thymine), and URA (uracil)'),
        ('nucleotide', 'nucleotide', 2, 'residues DA DC DG DU A C G T U'),
        ('nucleoside', 'nucleoside', 2, 'residues AMP ADP ATP CDP CTP GMP GDP GTP TMP TTP UMP UDP UTP'),
        ('at', 'ADE A THY T', 2, 'residues ADE A THY T of nucleic acids'),
        ('cg', 'CYT C GUN G', 2, 'residues CYT C GUN G of nucleic acids'),
        ('purine', 'purine', 2, 'residues ADE A GUN G of nucleic acids'),
        ('pyrimidine', 'pyrimidine', 2, 'residues CYT C THY T URA U of nucleic acids'),
        'Separator',
        ('large', 'large', 2, 'residues ARG, GLN, GLU, GLX, HIP, HIS, HSD, HSE, HSP, ILE, LEU, LYS, ET, MSE, PHE, PTR, SEP, TPO, TRP, TYR, XLE'),
        ('medium', 'medium', 2, 'residues ASN, ASP, ASX, CSO, CYS, PRO, SEC, THR, VAL'),
        ('polar', 'polar', 2, 'residues ARG, ASN, ASP, ASX, CSO, CYS, GLN, GLU, GLX, GLY, HIP, HIS, SD, HSE, HSP, LYS, PTR, SEC, SEP, SER, THR, TPO, TYR'),
        ('small', 'small', 2, 'residues ALA, GLY, SER'),
        'Separator',
        ('buried', 'buried', 2, 'residues ALA, CYS, ILE, LEU, MET, MSE, PHE, SEC, TRP, VAL, XLE'),
        ('surface', 'surface', 2, 'residues ARG, ASN, ASP, ASX, CSO, GLN, GLU, GLX, GLY, HIP, HIS, HSD, SE, HSP, LYS, PRO, PTR, SEP, SER, THR, TPO, TYR')
        ],

    'Element' : [
        ('carbon', 'carbon', 0, 'carbon atoms'),
        ('nitrogen', 'nitrogen', 0, 'nitrogen atoms'),
        ('oxygen', 'oxygen', 0, 'oxygen atoms'),
        ('sulfur', 'sulfur', 0, 'sulfur atoms'),
        ('hydrogen', 'hydrogen', 0, 'hydrogen atoms'),
        ('heavy', 'heavy', 0, 'non hydrogen atoms')
        ],
    'Structure' : [
        ('calpha', 'calpha', 2, 'carbon alpha atoms of protein residues'),
        ('backbone', 'backbone', 2, 'non-hydrogen backbone atoms of protein residues'),
        ('backbonefull', 'backbonefull', 2, 'backbone atoms of protein residues (including hydrogens)'),
        ('sidechain', 'sidechain', 2, 'side-chain atoms of protein residues'),
        'Separator',
        ('extended', 'extended', 2, 'extended conformation, same as "secondary E"'),
        ('helix', 'helix', 2, 'alpha-helix conformation, same as "secondary H"'),
        ('helix310', 'helix310', 2, '3_10-helix conformation, same as "secondary G"'),
        ('helixpi', 'helixpi', 2, 'pi-helix conformation, same as "secondary I"'),
        ('turn', 'turn', 2, 'hydrogen bonded turn conformation, same as "secondary T"'),
        ('bridge', 'bridge', 2, 'isolated beta-bridge conformation, same as "secondary B"'),
        ('bend', 'bend', 2, 'bend conformation, same as "secondary S"'),
        ('coil', 'coil', 2, 'not in any particular secondary structure conformation, same as "secondary C"')
        ]
    }

class Selection(ProdySelection):

    def __init__(self, ag, indices, selstr, acsi=None, **kwargs):
        ProdySelection.__init__(self, ag, indices, selstr, acsi=acsi, **kwargs)
        self.updateAgDict() # update _agDict which has atomIndex: None for atoms in selection
        self._molecule = weakref.ref(ag._molecule())

    def getAtomGroup(self):
        return self._molecule()._ag
        
    def updateAgDict(self):
        self.atomsDict = {}.fromkeys(self.getIndices(), None)
        
    def copy(self):
        return self.__class__(
            self.getAtomGroup(),
            self.getIndices(),
            self.getSelstr(),
            acsi = self.getACSIndex())
    
    def getBonds(self):
        # atom indices of atoms with no bonds are returned in bonds[0]
        # atom indices of all single, triple, and aromatic bonds are in bonds[1]
        # atom indices of double bonds are in bonds[2]
        # atom indices of triple bonds are in bonds[3]
        # atom indices of aromatic bonds are in bonds[4]
        ag = self.getAtomGroup()
        # list of: no, single double, triple, aromatic bonds
        if ag._bonds is None: # no bonds in molecule
            if len(self)==0: # no atoms in selection
                return [ [], [], [], [], [] ]
            else:
                return [ self.getIndices().tolist(), [], [], [], [] ]
        elif len(self)==0: # no atoms in selection
            return [ [], [], [], [], [] ]

        bonds = [ [], [], [], [], [] ]
        indices = self.select('not deleted').getIndices()
        inSet = {}.fromkeys(indices)
        bo = self._ag._bondOrder
        t0 = time()
        for i in indices: # loop over atoms indices in this selection
            for bn, j in enumerate(ag._bmap[i]):
                if j==-1:
                    if bn==0:
                       bonds[0].append(i) 
                    break
                if i<j:
                    if inSet.has_key(j): # check if other atom is in this selection
                        if bo is None:
                            bonds[1].append( (i,j) )
                        else:
                            bondOrder = bo['%d %d'%(i,j)]
                            if bondOrder==4 or bondOrder==5: # aromatic and amide bonds will draw as single bonds with added dashed line
                                bonds[1].append( (i,j) )
                                bonds[4].append( (i,j) )
                            elif bondOrder==3: # aromatic will draw as single bonds with added dashed line
                                bonds[1].append( (i,j) )
                                bonds[3].append( (i,j) )
                            else:
                                bonds[bondOrder].append( (i,j) )
                        ## bonds[1].append( (i,j) )
                        ## if bo is not None:
                        ##     bondOrder = bo['%d %d'%(i,j)]
                        ##     if bondOrder > 1:
                        ##         bonds[bondOrder].append( (i,j) )
        #import pdb
        #pdb.set_trace()
        return bonds

    def __or__(self, other):
        """Returns an :class:`.Selection` instance."""

        try:
            ag = other.getAtomGroup()
        except AttributeError:
            raise TypeError('unsupported operand type(s) for +: {0} and '
                            '{1}'.format(repr(type(self).__name__),
                                         repr(type(other).__name__)))

        if ag != self.getAtomGroup():
            raise ValueError('AtomPointer instances must point to the same '
                             'AtomGroup instance')
        acsi = self.getACSIndex()
        if acsi != other.getACSIndex():
            LOGGER.warning('Active coordset indices of atoms are not the same.'
                           ' Result will have ACSI {0}.'.format(acsi))

        indices = unique(concatenate((self._getIndices(),
                                      other._getIndices())))

        return Selection(self.getAtomGroup(), indices, '({0}) + ({1})'
                         .format(self.getSelstr(), other.getSelstr()),
                         acsi, unique=True)

    def __sub__(self, other):
        """Returns an :class:`.Selection` instance."""

        try:
            ag = other.getAtomGroup()
        except AttributeError:
            raise TypeError('unsupported operand type(s) for +: {0} and '
                            '{1}'.format(repr(type(self).__name__),
                                         repr(type(other).__name__)))

        if ag != self.getAtomGroup():
            raise ValueError('AtomPointer instances must point to the same '
                             'AtomGroup instance')
        acsi = self.getACSIndex()
        if acsi != other.getACSIndex():
            LOGGER.warning('Active coordset indices of atoms are not the same.'
                           ' Result will have ACSI {0}.'.format(acsi))

        indices = set(self._getIndices())

        indices = list(indices.difference(other.getIndices()))

        return Selection(self.getAtomGroup(), indices, '({0}) - ({1})'
                         .format(self.getSelstr(), other.getSelstr()),
                         acsi, unique=True)

    def __and__(self, other):
        """Returns an :class:`.Selection` instance."""
        try:
            ag = other.getAtomGroup()
        except AttributeError:
            raise TypeError('unsupported operand type(s) for +: {0} and '
                            '{1}'.format(repr(type(self).__name__),
                                         repr(type(other).__name__)))

        if ag != self.getAtomGroup():
            raise ValueError('AtomPointer instances must point to the same '
                             'AtomGroup instance')
        acsi = self.getACSIndex()
        if acsi != other.getACSIndex():
            LOGGER.warning('Active coordset indices of atoms are not the same.'
                           ' Result will have ACSI {0}.'.format(acsi))

        indices = set(self._getIndices())

        indices = list(indices.intersection(other.getIndices()))
        return Selection(self.getAtomGroup(), indices, '({0}) - ({1})'
                         .format(self.getSelstr(), other.getSelstr()),
                         acsi, unique=True)

    def getPDBRecords(self):
        stream = Stream()
        prody.writePDBStream(stream, self)
        return stream.lines

    def copy(self):
        # override so that it remains a MolKit2 Selection
        return Selection(self.getAtomGroup(), self._indices.copy(),
                         self.getSelstr(), acsi=self.getACSIndex())

    def select(self, what):
        try:
            sel = self.getAtomGroup().select(what)
        except ValueError: # happens for secondary structure when info is not yet assigned
            return None
        if sel:
            sel = sel & self
            if sel:
                return Selection(self.getAtomGroup(), sel.getIndices(), what, acsi=sel.getACSIndex())
            else:
                #return None
                return Selection(self.getAtomGroup(), [], 'empty', acsi=0)
        else:
            #return None
            return Selection(self.getAtomGroup(), [], 'what', acsi=0)

    def  getTransformedCoords(self):
        """ Get transformation matrix of the molecule and apply it to
        coords"""
        # get a handle on the master geometry
        mg = self.getAtomGroup().getMolecule().geomContainer.masterGeom
        matrix = mg.GetMatrix(mg).astype('f')
        # Transpose the transformation matrix.
        trMatrix = transpose(matrix)
        coords = self.getCoords()
        # Transform coords into homogeneous coords.
        hCoords = concatenate((coords, ones( (len(coords), 1), 'f')), 1)
        tCoords = dot(hCoords, trMatrix)
        return tCoords[:,:3]

    def toAtomGroup(self, name=None):
        """creates a new AtomGroup object from selection (self).
        Copies all existing attributes of the selection.""" 
        oldag = self.getAtomGroup()
        if not name:
            name = self.getAtomGroup()._title

        from prody.atomic.atomgroup import AtomGroup
        ag = AtomGroup(name)
        ag._n_atoms = len(self)
        indices = self.getIndices()

        if oldag._n_csets:
            ag._setCoords(oldag._coords[:,indices], overwrite=True)

        from prody.atomic.functions import SAVE_SKIP_POINTER as SKIP
        for key in oldag.getDataLabels():
            if key in SKIP:continue
            ag.setData(key,  oldag.getData(key)[indices] )

        for key in oldag.getFlagLabels():
            if key in SKIP: continue
            ag._setFlags(key, oldag.getFlags(key)[indices])
        
        for key in oldag._bondData:
            ag._bondData[key] = oldag._bondData[key][indices]

        ## if self is not an entire molecule or atoms were deleted
        ## we need to re-index bonds and bond order
        newBo = []
        newBonds = []                
        if len(oldag)==len(self): # whole molecule, no deletes
            # handle bonds and bondOrder
            if oldag._bonds is not None:
                newBonds  = oldag._bonds.copy()
            if oldag._bondOrder is not None:
                newBo = oldag._bondOrder.copy()
        else: # either deleted or sub set ==> re-index
            bonds = self.getBonds()[1]
            if len(bonds)>0:
                size = len(oldag)
                old = indices = arange(0, size)
                a1 = ones(size)*-1
                tokeep = self.getIndices()
                # a1 is an array of old indices with -1 for deleted atoms
                a1[tokeep] = indices[tokeep]
                mapping = ones(size,'i')*-1
                # build the mapping such that mapping[i] wil be index of i
                # in new set
                maxi = 0
                for i in xrange(len(a1)):
                    if a1[i]!=-1:
                        mapping[i] = maxi
                        maxi += 1
                # fix the bonds
                newBonds = []
                newBo = {}
                bo = oldag._bondOrder
                for i,j in bonds:
                    if mapping[i]!=-1 and mapping[j]!=-1:
                        newBonds.append( (mapping[i],mapping[j]) )
                        if bo is not None:
                            v = bo['%d %d'%(i,j)]
                            newBo['%d %d'%(mapping[i],mapping[j])] = v
        if len(newBo):
            ag.setBonds(newBonds, newBo)
        else:
            ag.setBonds(newBonds)
                    
        return ag
        
    def toAtomGroupOLD(self, name=None):
        """creates a new AtomGroup object from selection (self).
        Copies all existing attributes of the selection.""" 
        oldag = self.getAtomGroup()
        if not name:
            name = self.getAtomGroup()._title

        from prody.atomic.atomgroup import AtomGroup
        ag = AtomGroup(name)
        ag._setCoords(self.getCoords())
        # data fields
        for name in ['Resnums', 'Names', 'Chids', 'Altlocs', 'Occupancies',
                     'Elements', 'Betas', 'Icodes', 'Resnames', 'Serials',
                     'Segnames']:
            data = getattr(self, 'get%s'%name)()
            getattr(ag, 'set%s'%name)(data)

        # mgl data fields
        for name in ['atomType', 'hbType', 'atomicNumber']:
            data = self.getData(name)
            if data is None:
                continue
            ag.setData(name, data)
                
        # flags
        for name in ['hetatm', 'pdbter', 'deleted']:
            data = self.getFlags(name)
            if data is not None:
                ag.setFlags(name, data)

        ## if self is not an entire molecule or atoms were deleted
        ## we need to re-index bonds and bond order
        newBo = []
        newBonds = []                
        if len(oldag)==len(self): # whole molecule, no deletes
            # handle bonds and bondOrder
            if oldag._bonds is not None:
                newBonds  = oldag._bonds.copy()
            if oldag._bondOrder is not None:
                newBo = oldag._bondOrder.copy()
        else: # either deleted or sub set ==> re-index
            bonds = self.getBonds()[1]
            if len(bonds)>0:
                size = len(oldag)
                old = indices = arange(0, size)
                a1 = ones(size)*-1
                tokeep = self.getIndices()
                # a1 is an array of old indices with -1 for deleted atoms
                a1[tokeep] = indices[tokeep]
                mapping = ones(size,'i')*-1
                # build the mapping such that mapping[i] wil be index of i
                # in new set
                maxi = 0
                for i in xrange(len(a1)):
                    if a1[i]!=-1:
                        mapping[i] = maxi
                        maxi += 1
                # fix the bonds
                newBonds = []
                newBo = {}
                bo = oldag._bondOrder
                for i,j in bonds:
                    if mapping[i]!=-1 and mapping[j]!=-1:
                        newBonds.append( (mapping[i],mapping[j]) )
                        if bo is not None:
                            v = bo['%d %d'%(i,j)]
                            newBo['%d %d'%(mapping[i],mapping[j])] = v
        if len(newBo):
            ag.setBonds(newBonds, newBo)
        else:
            ag.setBonds(newBonds)
                    
        return ag
        
    def toMolecule(self, name):
        ag = self.toAtomGroup(name=name)
        from MolKit2.molecule import Molecule
        mol = Molecule(name, ag)
        mol._ag.setData('colorsIndices_lines', [0]*len(ag))
        #oldmol = ag.getMolecule()
        #mol._colors = oldmol._colors.copy()
        #mol._multi = oldmol._multi
        #mol._renderingProp = oldmol._renderingProp.copy()
        # this does not copy _bondOrderData, _group .. not sure what to do
        return mol


class SelectionSet(list, TreeObject):

    def nbChildren(self):
        return len(self)
    
    def getChildrenAndNames(self, rootSelection=None):
        children = []
        names = []
        for sel in self:
            mol = sel.getAtomGroup().getMolecule()
            names.append(mol.name)
            children.append( mol )
        return [children, names]

    def __init__(self, itr=[], name='NoName'):
        self._agDict = {} # mapping between AtomGroups representing molecules
                          # and Selection objects in the list. Used to make sure
                          # there is only one Selection object per mol
        self.name = name
        for sel in itr:
            self.append(sel)
        
    def copy(self):
        new = SelectionSet(name=self.name)
        for sel in self:
            new.append(sel.copy())
        return new

    def selectionFor(self, sel):
        # return True is this selectionSet has a selection pointing to the
        # same atom group at sel
        return self._agDict.get(sel.getAtomGroup(), None)

    
    def append(self, obj):
        # insure obj is a Selection object and we only have one Selection
        # object 
        assert isinstance(obj, Selection)
        sel = self._agDict.get(obj.getAtomGroup(), None)
        if sel:
            newsel = sel | obj
            self[self.index(sel)] = newsel
            self._agDict[obj.getAtomGroup()] = newsel
        else:
            list.append(self, obj)
            self._agDict[obj.getAtomGroup()] = obj
            
    def __or__(self, other):
        """Returns an :class:`.SelectionSet` instance."""

        if isinstance(other, SelectionSet):
            newSels = []
            for ag, sel in self._agDict.items():
                otherSel = other._agDict.get(ag, None)
                if otherSel:
                    newSels.append(sel | otherSel)
                else:
                    newSels.append(sel.copy())

            for otherAg, otherSel in other._agDict.items():
                sel = self._agDict.get(otherAg, None)
                if not sel:
                    newSels.append(otherSel.copy())
                
            return SelectionSet(newSels)

        else: # bad other
            raise ValueError('second operand for SelectionSet OR must be SelectionSet. Got %s'%other.__class__)

        return SelectionSet(newSels)
            
    def __and__(self, other):
        """Returns an :class:`.SelectionSet` instance."""

        if isinstance(other, SelectionSet):
            newSels = []
            for ag, sel in self._agDict.items():
                otherSel = other._agDict.get(ag, None)
                if otherSel:
                    sel = sel & otherSel
                    if len(sel):
                        newSels.append(sel)
                else:
                    newSels.append(sel.copy())
            return SelectionSet(newSels)

        else: # bad other
            raise ValueError('second operand for SelectionSet AND must be SelectionSet. Got %s'%other.__class__)

        return SelectionSet(newSels)

    def __sub__(self, other):
        """Returns an :class:`.SelectionSet` instance."""
        
        if isinstance(other, SelectionSet):
            newSels = []
            for ag, sel in self._agDict.items():
                otherSel = other._agDict.get(ag, None)
                if otherSel:
                    sel = sel - otherSel
                    if len(sel):
                        newSels.append(sel)
                else:
                    newSels.append(sel.copy())
            return SelectionSet(newSels)

        else: # bad other
            raise ValueError('second operand for SelectionSet SUB must be SelectionSet. Got %s'%other.__class__)


    def allSelected(self, inSelectionSet):
        # returns True if all atoms in inSelection are in this PmvSelection
        #         False is all atoms in obj are deselected
        #         'partial' else

        if self.nbAtoms()==0: return False # empty selection -> allSelected False

        return self._allSelected(inSelectionSet)
            
    def _allSelected(self, inSelectionSet):
        for inSelection in inSelectionSet:
            # get the selection in self that corresponding to incoming inSelection
            currentSel = self._agDict.get(inSelection.getAtomGroup(), None)

            if currentSel is None:
                # this PmvSelection does not have anything for the molecule inSelection is for
                return False

            # now we have to find out is all indices of inSelection are in currentSel
            if len(inSelection) > len(currentSel): return False

            curIndicesDict = {}.fromkeys(currentSel.getIndices())
            for i in inSelection.getIndices():
                if curIndicesDict.get(i, 1) is 1:
                    return False

        return True

    def isAnySelected(self, inSelectionSet, mode='fast'):
        # returns true is at least one atom in inSelectionSet is in self

        len1 = self.nbAtoms()
        if len1==0: return False

        for selection in self:
             # get the corresponding incoming Selection
            currentSel = inSelectionSet._agDict.get(selection.getAtomGroup(), None)
            if currentSel is None: continue
            if len(currentSel & selection) > 0: return True

        return False

    def update(self):
        # compute self.children and self.atomsDict
        for sel in self:
            sel.updateAgDict()

    def nbAtoms(self):
        # number of atoms in that selection
        size = 0
        for sel in self:
            size += len(sel)
        return size

    def clear(self):
        # remove all selections
        while len(self):
            self.pop()
        self._agDict = {}
        
    def replace(self, newSelections):
        # replace all Selections by the ones in newSelections
        self.clear()
        for sel in newSelections:
            self.append(sel)

    def select(self, what):
        result = SelectionSet()
        for sel in self:
            res = sel.select(what)
            if res:
                result.append(res)
        if len(result): return result
        else: return None

    def numAtoms(self, what=None):
        nb = 0
        for sel in self: nb+= sel.numAtoms(what)
        return nb
    
    def  getTransformedCoords(self):
        """ Get transformation matrix of the molecule and apply it to
        coords"""
        coords = ones( (self.nbAtoms(), 3), 'f')
        off = 0
        for sel in self:
            l = len(sel)
            coords[off:off+l] = sel.getTransformedCoords()
            off += l
        return coords

    def toAtomGroup(self, name=None):
        """creates a new AtomGroup object from selection (self).
        Copies all existing attributes of the selection.""" 
        if not name:
            name = self[0].name#getAtomGroup()._title

        from prody.atomic.atomgroup import AtomGroup
        from prody.atomic.functions import SAVE_SKIP_POINTER as SKIP
        #import pdb; pdb.set_trace()
        ag = AtomGroup(name)
        coords = []
        dataFields = {}
        flags = {}

        bonds = []
        bondOrder = {}
        maxInd = 0
        for i, sel in enumerate(self):
            oldag = sel.getAtomGroup()
            mol = sel.getAtomGroup().getMolecule()
            indices = sel.getIndices()
            if oldag._n_csets:
                if i==0:
                    coords=sel.getCoords()
                else:
                    coords = concatenate((coords, sel.getCoords()))
            # data fields
            for key in oldag.getDataLabels():
                if key in SKIP or key.startswith("colorsIndices"): continue
                if i==0:
                    dataFields[key] =  oldag.getData(key)[indices]
                else:
                    if key == 'Serials':
                        off = dataFields[key][-1]
                        dataFields[key] = concatenate((dataFields[key], oldag.getData(key)[indices] + off))
                    else:
                        dataFields[key] = concatenate((dataFields[key], oldag.getData(key)[indices]))
            # flags
            for key in oldag.getFlagLabels():
                if key in SKIP:continue
                if i == 0:
                    flags[key] = oldag.getFlags(key)[indices] 
                else:
                    flags[key] = concatenate((flags[key], oldag.getFlags(key)[indices]))
        ag._setCoords(coords, overwrite=True)
        for key, data in dataFields.items():
            ag.setData(key, data)
        for key, data in flags.items():
            ag._setFlags(key, data)
            
        for i, sel in enumerate(self):
            oldag = sel.getAtomGroup()
            mol = sel.getAtomGroup().getMolecule()
            ## if sel is not an entire molecule or atoms were deleted
            ## we need to re-index bonds and bond order
            if len(oldag)==len(self) and i == 0: # whole molecule, no deletes
                # handle bonds and bondOrder
                if oldag._bonds is not None:
                    bonds = oldag._bonds.copy()
                if oldag._bondOrder is not None:
                    bondOrder= oldag._bondOrder.copy()
                maxInd = len(sel)
            else: # either deleted or sub set ==> re-index
                _bonds = sel.getBonds()[1]
                if len(_bonds)>0:
                    size = len(oldag)
                    old = indices = arange(0, size)
                    a1 = ones(size)*-1
                    tokeep = sel.getIndices()
                    # a1 is an array of old indices with -1 for deleted atoms
                    a1[tokeep] = indices[tokeep]
                    mapping = ones(size,'i')*-1
                    # build the mapping such that mapping[i] wil be index of i
                    # in new set
                    maxi = maxInd
                    for i in xrange(len(a1)):
                        if a1[i]!=-1:
                            mapping[i] = maxi
                            maxi += 1
                    # fix the bonds
                    newBonds = []
                    newBo = {}
                    bo = oldag._bondOrder
                    for i,j in _bonds:
                        if mapping[i]!=-1 and mapping[j]!=-1:
                            newBonds.append( (mapping[i],mapping[j]) )
                            if bo is not None:
                                v = bo['%d %d'%(i,j)]
                                newBo['%d %d'%(mapping[i],mapping[j])] = v
                    if not len(bonds):
                        bonds = newBonds
                    else:
                        bonds = concatenate((bonds, newBonds))            
                    if len(newBo):
                        bondOrder.update(newBo)
                maxInd = maxInd + len(sel)

        if len(bondOrder):
            ag.setBonds(bonds, bondOrder)
        else:
            ag.setBonds(bonds)
        return ag

    def toMolecule(self, name):
        ag = self.toAtomGroup(name=name)
        from MolKit2.molecule import Molecule
        mol = Molecule(name, ag)
        return mol

    
def getMolsFromName(molName, molSet):
    prog = re.compile(molName)
    mols = {}
    for mol in molSet:
        if prog.match(mol.name):
            mols[mol] = True
    return mols.keys()

def selector(selectionStr, molSet):

    selections = SelectionSet()
    # loop over ';' separated selection strings
    for selStr in selectionStr.split(';'):
        tokens = selStr.split(':')
        molNames = tokens[0].split(',')
        mols = {}
        for molName in molNames:
            mols.update( {}.fromkeys( getMolsFromName(molName, molSet) ))

        if len(tokens)==1: # only molecule names were given
            args = ()
        else:
            args = tokens[1]

        for mol in mols.keys():
            selections.append( mol.select(*args) )
    return selections
