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
# $Header: /mnt/raid/services/cvs/MolKit2/molecule.py,v 1.6.2.1 2017/07/26 22:03:40 annao Exp $
# 
# $Id: molecule.py,v 1.6.2.1 2017/07/26 22:03:40 annao Exp $
#
import prody, os, numpy, re, sys
from time import time, sleep
from math import sqrt

from prody.atomic.atomgroup import AtomGroup
from prody.atomic.chain import Chain as ProDyChain
from prody.atomic.residue import Residue as ProDyResidue
from prody.atomic.atom import Atom as ProDyAtom
from prody import SelectionError

from .selection import Selection, SelectionSet
from .babelElements import babel_elements
from .tree import TreeObject

from bhtree import bhtreelib
from mglutil.math.torsion import torsion
from mglutil.util.packageFilePath import getCacheFolder
import openbabel as ob
import pybel

from MolKit2.PDBresidueNames import AAnames

##########################################################################
##
## functions used to translate MolKit syntax into prody selection strings
##
##########################################################################

# define some regular expression used for selection\
# using the MolKit syntax
is_re = re.compile('.*[\.\^\$\*\+\?\{\}\[\]\\\|\(\)].*')
is_range = re.compile('.*to.*')
splitOpPattern = re.compile('(/[\+\-&]/)')

residueSelKw = {}.fromkeys(['protein', 'aminoacid', 'stdaa', 'nonstdaa', 'acidic', 'acyclic',
                            'SEC,', 'aliphatic', 'aromatic', 'basic', 'buried', 'charged',
                            'cyclic', 'hydrophobic', 'large', 'medium', 'neutral', 'polar',
                            'small', 'surface', 'nucleic', 'nucleobase', 'nucleotide', 'nucleoside',
                            'water', 'lipid', 'sugar', 'heme', 'helix', 'helix310', 'helixpi',
                            'turn', 'bridge', 'bend', 'coil'])

def residueToken(token):
    """token can be:
    - a number (eg 10 for ALA10)
    - #Num for 1 based indices (eg. #1 for first residue)
    - a range of of number 10to27
    - a range of indices #1to#5
    - a regular expression on the residue name
    """
    if token=='': return '','','',''
    resnames = []
    resnums = []
    serials = []
    kw = []
    if is_re.match(token): # regular expression
        resnames.append('"%s"'%token)
    elif is_range.match(token): # range
        if token[0]=='#': # range of 1-based indices
            _from, _to = token.split('to') # make it a range of 0-based serial numbers
            try:
                serials.append('%dto%d'%(int(_from[1:])-1,int(_to[1:])-1))
            except ValueError:
                pass
        else: # keep the range of residues numbers
            resnums.append(token)

    elif AAnames.has_key(token): # AA name such as ALA
        resnames.append(token)

    elif residueSelKw.has_key(token): # key words such as backbone
        kw.append(token)

    elif token[0]=='#': # 1-based index
        try:
            serials.append(str(int(token[1:])-1))
        except ValueError:
            pass
    elif token.isdigit():
        resnums.append(token)

    return resnames, resnums, serials, kw

def getProdyResSelStr(mselStr):
    resnames = []
    resnums = []
    serials = []
    kws = []
    for token in mselStr.split(' '):
        a,b,c,d = residueToken(token)
        resnames.extend(a)
        resnums.extend(b)
        serials.extend(c)
        kws.extend(d)

    selStr = ""
    if len(resnames):
        selStr += 'resname '
        for name in resnames:
            selStr += name + ' '
    if len(resnums):
        if len(selStr):
            selStr += ' or '
        selStr += 'resnum '
        for name in resnums:
            selStr += name + ' '
    if len(serials):
        if len(selStr):
            selStr += ' or '
        selStr += 'resindex '
        for name in serials:
            selStr += name + ' '
    if len(kws):
        if len(selStr):
            selStr += ' and '
        selStr += ' '
        for name in kws:
            selStr += name + ' '
    return selStr

def chainsToken(token, mol):
    """token can be:
    - a letter (eg A for chain A)
    - #Num for 1 based indices (eg. #1 for first chain)
    - a range of indices #1to#5
    """
    # since prody does not support ranges on chains we expand chain IDs
    chainIds = []
    if token=='': return ''
    if token[0]=='#':
        allChainsIds = numpy.unique(mol._ag.getChids())
        match = is_range.match(token)
        if match: # chain range by index
            _from, _to = token.split('to')
            try:
                for chid in allChainsIds[int(_from[1:])-1:int(_to[1:])-1]:
                    chainIds.append(chid)
            except ValueError:
                pass
        else: # keep the range of residues numbers
            try:
                chainIds.append(allChainsIds[int(token[1:])-1])
            except ValueError:
                pass
    else: # single chain by index
        chainIds.append(token)
    return chainIds
    
def getProdyChainSelStr(mselStr, mol):
    chainIds = ""
    for token in mselStr.split(' '):
        chids = chainsToken(mselStr, mol)
        if len(chids):
            for name in chids:
                chainIds += name + ' '
    if len(chainIds):
        return 'chain %s'%chainIds
    else:
        return ''

atomSelKw = {}.fromkeys(['calpha', 'ca', 'backbone', 'bb', 'backbonefull', 'bbfull',
                         'sidechain', 'sc', 'at', 'cg', 'purine', 'pyrimidine', 'hetero',
                         'hetatm', 'ion', 'carbon', 'nitrogen', 'oxygen', 'sulfur', 'hydrogen',
                         'noh'])

def atomToken(token):
    """1 - string e.g. 'CA' ==> 'name CA'
                    '10' ==> 'serial 10'
    2 - atom index #1, ==>'index 1'
    4 - ranges of indices ==> 'index 1to4'
    5 - Regular expressions matched against the atom name
    keywords: see page 58 59 or prody manual
        ::backbone ==> 'backbone' """
    if token=='': return '','','',''
    names = []
    inds = [] # 1 based index in list of atoms
    serials = [] # atom number from PDB file
    kws = []
    if is_re.match(token): # regular expression
        names.append('"%s"'%token)

    elif is_range.match(token): # range
        if token[0]=='#': # range of 1-based indices
            try:
                _from, _to = token.split('to') # make it a range of 0-based serial numbers
                inds.append('%dto%d'%(int(_from[1:])-1,int(_to[1:])-1))
            except ValueError:
                pass
        else: # keep the range of atom indices
            inds.append(token)

    elif atomSelKw.has_key(token): # key words such as backbone
        kws.append(token)

    elif token[0]=='#': # 1-based index
        try:
            inds.append(str(int(token[1:])-1))
        except ValueError:
            pass
    elif token.isalpha():
        names.append(token)
    elif token.isdigit():
        serials.append(token)
    elif token.isalnum():
         names.append(token)
    return names, inds, serials, kws

def getProdyAtomSelStr(mselStr):
    atnames = []
    atinds = []
    atserials = []
    kws = []
    for token in mselStr.split(' '):
        a,b,c,d = atomToken(token)
        atnames.extend(a)
        atinds.extend(b)
        atserials.extend(c)
        kws.extend(d)

    selStr = ""
    if len(atnames):
        selStr += 'name '
        for name in atnames:
            selStr += name + ' '
    if len(atinds):
        if len(selStr):
            selStr += ' or '
        selStr += 'index '
        for name in atinds:
            selStr += name + ' '
    if len(atserials):
        if len(selStr):
            selStr += ' or '
        selStr += 'serial '
        for name in atserials:
            selStr += name + ' '
    if len(kws):
        if len(selStr):
            selStr += ' and '
        selStr += ' '
        for name in kws:
            selStr += name + ' '
    return selStr

##########################################################################
##
## PROSS-based secondary structure assignment
##
##########################################################################
from prody.measure.measure import calcPhi, calcPsi

class AssignSecondaryStructureWtihPross:
    """class to extend GetSecondaryStructure with specific methods to
    get informations on the secondary structure from PROSS This module 
	uses much of the code from the original BIOMOL collection of utilities 
	written by Raj Srinivasani with enhancements by Nick Fitzkee.

    The script was put together by Pat Fleming so that a user would not need
	to have the BIOMOL distribution installed to run PROSS.

    Note that since Raj's time the definitions of mesostates has been superceded
    by the fine grained 30 deg x 30 deg grid for most purposes. Either mesostate
    grid will work for PROSS. Give your choice as an argument (see USAGE below).

    Date: September 2004
    Author: Pat Fleming, pat.fleming@jhu.edu
	http://www.roselab.jhu.edu/utils/pross.html"""

    def __init__(self, default='fgmeso'):
        """ # Setup Fine Grain Mesostate Bins (ala Pat Fleming)"""
        """ see PROSS.py """
        from PROSS import MSDEFS
        self.mode=default
        self.MSDEFS=MSDEFS[default]

    def __call__(self, mol):
        for chain in mol._ag.select('not deleted').getHierView().iterChains():
            proteinAtoms = chain.select("protein and not deleted")
            if proteinAtoms is not None:
                   self.rc_ss(chain)
		
    def computePhiPsi(self, chain):
        phi = []
        psi = []
        for res in chain.iterResidues():
            try:
                phi.append(calcPhi(res))
            except ValueError:
                phi.append(0.)
            try:
                psi.append(calcPsi(res))
            except ValueError:
                psi.append(0.)
        return phi, psi
        
    def rc_ss(self, chain):
        """rc_ss(chain, phi, psi, ome) - calculate secondary structure

        This function calculates the secondary structure using the PROSS method
        with rotamer codes.  Given a chain, and optionally a list of phi,
        psi, and omega, it calculates the backbone secondary structure of
        the chain.  The return value is (phi, psi, ome, sst), where
        phi, psi, and ome are calculated if not specified, and sst is the
        secondary structure codes: H = helix, E = strand, P = PII, C = coil.
        """
        
        ms = self.MSDEFS
        PII    = ms['PII']
        TURNS  = ms['TURNS']
        HELIX  = ms['HELIX']
        STRAND = ms['STRAND']

        nres = chain.numResidues()
        phi, psi = self.computePhiPsi(chain)

        """rc_codes(chain, phi, psi, ome) - return rotamer codes

        Given a protein chain (and optionally phi, psi, omega), this
        function will return a list of mesostate codes that
        applies to the chain, as determined by res_rc.
        """
        codes = [self.res_rc(x, y) for x, y in zip(phi, psi)]

        sst = ['C']*nres

        is_PII = PII.has_key

        for i in xrange(nres-1):
            code = codes[i]
            if is_PII(code):
                sst[i] = 'P'

        is_turn = TURNS.has_key

        for i in xrange(nres-1):
            code = codes[i] + codes[i+1]
            if is_turn(code):		
                sst[i] = sst[i+1] = 'T'

        helices = self._rc_find(codes, HELIX)
        strands = self._rc_find(codes, STRAND)

        for helix in helices:
            i, j = helix
            for k in range(i, j):
                 sst[k] = 'H'
        for strand in strands:
            i, j = strand
            for k in range(i, j):
                if sst[k] in ('C', 'P'): sst[k] = 'E'

        return chain.select("protein and bb and name CA and not deleted").setSecstrs(sst)

    def res_rc(self,r1, r2, r3=180):
        """res_rc(r1, r2, r3) - get mesostate code for a residue

        Given a phi (r1), psi (r2), and omega (r3) torsion angle, calculate
        the mesostate code that describes that residue. 
        A mesostate will be returned in
        all but two cases:  if omega deviates from planarity by more than 90
        degrees, '*' is returned.  Also, if any torsions are greater than
        180.0 (biomol assignes 999.0 degrees to angles that that are
        indeterminate), prosst.INVALID is returned.  Here, r3 defaults to
        180.0.
        """

        ms = self.MSDEFS
        OMEGA   = ms['OMEGA']
        INVALID = ms['INVALID']
        PHI_OFF = ms['PHI_OFF']
        PSI_OFF = ms['PSI_OFF']
        DELTA   = ms['DELTA']
        RC_DICT = ms['RC_DICT']
        if r1 == None or r2 == None : return INVALID 
        if (abs(r3) <= 90.0):
            return OMEGA
        elif r1>180.0 or r2>180.0 or r3>180.0:
            return INVALID

        ir1 = -int(PHI_OFF) + int(round((r1+PHI_OFF)/DELTA )) * int(DELTA)
        ir2 = -int(PSI_OFF) + int(round((r2+PSI_OFF)/DELTA )) * int(DELTA)

        while ir1 <= -180: ir1 = ir1 + 360
        while ir1 >   180: ir1 = ir1 - 360
        while ir2 <= -180: ir2 = ir2 + 360
        while ir2 >   180: ir2 = ir2 - 360

        return RC_DICT[(ir1,ir2)]

    def _rc_find(self, codes, pattern):
        """_rc_find(codes, pat_obj) - find a endpoints of a regexp

        Given a list of mesostate codes, this function identifies a endpoints
        of a match  <pattern>.  <pat_obj> is a compiled regular expression
        pattern whose matches will be returned as pairs indicated start,
        end in <codes>
        """

        CODE_LENGTH = self.MSDEFS['CODE_LENGTH']

        if not type(codes) == type(''):
                codes = ''.join(codes)

        matches = []
        it = pattern.finditer(codes)

        try:
                while 1:
                        mat = it.next()
                        matches.append((mat.start()/CODE_LENGTH, mat.end()/CODE_LENGTH))
        except StopIteration:
                pass

        return matches

##########################################################################
##
## END functions used to translate MolKit syntax into prody selection strings
##
##########################################################################

def getAtomIndicesPerType(atoms):
    d1 = {}
    for i, _type in enumerate(atoms.getTypes()):
        try:
            d1[_type].append(i)
        except KeyError:
            d1[_type] = [i]
    return d1


from mglutil.util.io import Stream

class Molecule(TreeObject):
    """
    The Molecule class is a wrapper of a ProDy AtomGroup object

    mol._ag is the prody atom group

    In addition to the regular ProDy fields molecule have the following
    data fields:
       atomicNumber: (int) atom type periodic table number
       atomType: (string 'S5') open babel atom type (set by .typeAtomsAndBonds())
       hbType:   int 0: None, 1: ACCEPTOR, 2: DONOR, 3: BOTH

    and flags:
       deleted: (bool) turned to true when an atom is deleted

    mol._ag.bondOrder is a dictionary with keys '%d %d' and value bond order
    for the bond between atom i and j
    
    mol.select(str) will select '(str) and not deleted' where str is a prody
       selection string
    """
    vdw_rad = numpy.zeros( 110, 'f' )
    for k,v in babel_elements.items():
        vdw_rad[v['num']] = v['vdw_rad']
        
    bo_rad = numpy.zeros( 110, 'f' )
    for k,v in babel_elements.items():
        bo_rad[v['num']] = v['bond_ord_rad']

    ## def __del__(self):
    ##     del self._ag
    ##     print 'deleting molecule', self.name
    ##     #del self._ag._molecule

    def atomFullName(self, atom):
        name = '%s:%s_%s:%s%d%s:%s'%(
            self.name, atom.getSegname(), atom.getChid(), atom.getResname(),
            atom.getResnum(), atom.getIcode(), atom.getName())
        return name
    
    def __init__(self, name, prodyAtomGroup=None, filename=None):

        self.app = None # will be weakref to app when molecule is loaded into app
        self._indexingThread = None
        assert prodyAtomGroup is None or isinstance(prodyAtomGroup, AtomGroup)
        self.filename = filename
        self.name = name
        if filename is None:
            self._basename = os.path.splitext(os.path.split(name)[1])[0]
        else:
            self._basename = os.path.splitext(os.path.split(filename)[1])[0]
        self.alias = ''
        self._multi = False # could be 'molecules' or 'conformations'
        if prodyAtomGroup is not None:
            self.setAtomgroup(prodyAtomGroup)

        self._bondOrderData = {
            'doubleBonds': None,
            'tripleBonds': None,
            'aromaticRings': None,
            'aromaticArcs': None,
            }
        self._renderingProp = {}
        # create a dictionary of indices for atoms colors. The keys are geoms,
        # values are numeric arrays of indices in to newmol._colors[geom]
        self._colors = {}
        self._automorphisms = None
        
    ## def extend(self, other):
    ##     assert self.__class__ == other.__class__

    ##     import pdb; pdb.set_trace()

    ##     # add the atomSet
    ##     lenBefore = len(self._ag)
    ##     self._ag.extend(other._ag)
    ##     import pdb; pdb.set_trace()
    ##     self.updateBondsByDistance()

    ##     # handle colors
    ##     for name in set(list(self._colors) + list(other._colors)):
    ##         thisColor = self._colors.get(name, None)
    ##         thatColor = other._colors.get(name, None)
    ##         if thisColor is not None or thatColor is not None:
    ##             if thisColor is not None:
    ##                 # offset color indices of other by length of colors of self
    ##                 cols = self._ag._data['colorsIndices_%s'%name]
    ##                 cols[lenBefore:] += len(thisColor)
    ##             self._colors[name] = numpy.concatenate((thisColor, thatColor))
    def getIsomorphisms(self, mol2):
        """
        use openbabel to identify ismorphisms in the molecule (i.e. atoms that are equivalent)
        """
        noH = self.select('not hydrogen')
        #avoid circular imports
        from MolKit2.openBabelInterface import ProdyToOBMol
        obmol = ProdyToOBMol(noH)
        query = pybel.ob.CompileMoleculeQuery(obmol)
        mapper = pybel.ob.OBIsomorphismMapper.GetInstance(query)
        isomorphs = pybel.ob.vvpairUIntUInt()
        mapper.MapAll(obmol, isomorphs)
        ## now renumber atom indices to account for hydrogen atoms
        isomorphisms = []
        for iso in isomorphs:
            isoRenum = [(obmol.obToProdyIndex[i], obmol.obToProdyIndex[j])
                        for i,j in iso]
            isomorphisms.append(isoRenum)
        return isomorphisms

    def getAutomorphisms(self):
        """
        use openbabel to identify automorphisms in the molecule (i.e. atoms that are equivalent)
        """
        if self._automorphisms == None:
            noH = self.select('not hydrogen')
            #avoid circular imports
            from MolKit2.openBabelInterface import ProdyToOBMol
            obmol = ProdyToOBMol(noH)
            gs = pybel.ob.OBGraphSym(obmol)
            symclasses = pybel.ob.vectorUnsignedInt()
            gs.GetSymmetry(symclasses)
            automorphs = pybel.ob.vvpairUIntUInt()
            pybel.ob.FindAutomorphisms(obmol, automorphs, symclasses)
            ## now renumber atom indices to account for hydrogen atoms
            automorphisms = []
            for auto in automorphs:
                autoRenum = [(obmol.obToProdyIndex[i], obmol.obToProdyIndex[j])
                             for i,j in auto]
                automorphisms.append(autoRenum)
            self._automorphisms = automorphisms
        return self._automorphisms

    def _traverse(self, at1, at2, inSubTree, limitTo=None):
        # identified the molecular sub-graph rooted at atom at1
        # while not alloing to move through atom at1.
        # inSubTree is a list intially containing the indices of at1 and at2
        # The method returns a list of atom indices of the sub-graph
        # limitTo is None for no limitation or a dict chid:[resnums]
        for at3 in at2.iterBonded():
            if at3==at1:
                continue
            if limitTo:
                try:
                    resList = limitTo[at3.getChid()]
                    if at3.getResnum() not in resList:
                        continue
                except KeyError:
                    continue # atom is not within limits
            if at3.getIndex() not in inSubTree:
                inSubTree.append(at3.getIndex())
                try:
                    self._traverse(at2, at3, inSubTree, limitTo=limitTo)
                except(RuntimeError, AttributeError):
                    pass
        return inSubTree

    def subTree(self, atom1, atom2, limitTo=None):
        assert isinstance(atom1, (ProDyAtom,prody.atomic.atom.Atom))
        assert isinstance(atom2, (ProDyAtom,prody.atomic.atom.Atom))
        inSubTree = self._traverse(
            atom1, atom2, [atom1.getIndex(), atom2.getIndex()], limitTo)
        return Selection(self._ag, inSubTree, "")

    def setLinewidth(self, width):
        assert isinstance(width, int)
        assert width > 0 
        self._renderingProp['lines']['width'] = width
        
    def setDisplayBondOrder(self, truth):
        assert truth in [True, False, 0, 1]
        assert truth > 0 
        self._renderingProp['lines']['displayBondOrder'] = truth

    def setStipple(self, length, space):
        assert isinstance(length, float)
        assert length > 0.
        self._renderingProp['lines']['stippleLength'] = length
        assert isinstance(space, float)
        assert space > 0.
        self._renderingProp['lines']['stippleSpace'] = space

    def setStippleLength(self, length):
        assert isinstance(length, float)
        assert length > 0.
        self._renderingProp['lines']['stippleLength'] = length

    def setStippleSpace(self, space):
        assert isinstance(space, float)
        assert space > 0.
        self._renderingProp['lines']['stippleSpace'] = space

    def setDoubleBondSep(self, val):
        assert isinstance(val, float)
        assert val > 0.
        self._renderingProp['lines']['doubleBondSep'] = val

    def setTripleBondSep(self, val):
        assert isinstance(val, float)
        self._renderingProp['lines']['TripleBondSep'] = val
        assert val > 0.        

    def setAromaticLineWidth(self, width):
        assert isinstance(width, int)
        assert width > 0 
        self._renderingProp['lines']['aromaticLinewidth'] = width
        
    def numMols(self):
        if not self._multi:
            return 1
        elif isinstance(self, MultiMolecule):
            return len(self.index)
        else:
            return self._ag._coords.shape[0]

    def curMolIndex(self):
        if not self._multi:
            return 0
        elif isinstance(self, MultiMolecule):
            return self.currentMoleculeIndex
        else:
            return self._ag.getACSIndex()

    def addAGBondData(self, name, value, defaultValue):
        ## add a data field to ag._data with default value
        if self._ag._bondData.has_key(name):
            raise ValueError("Molecule.addAGBondData: there already is a field called %s in the AtomGroup._bondData of molecule %s"%(name, self.name))
        self._ag._bondData[name] = value
        #self._ag._customFieldsDefaultValues['bonds'][name] = defaultValue

    def setAtomgroup(self, ag):
        self._ag = ag
        # when adding a new field to the atom groups a default value should
        # added in this dict with to that the addAtoms method can create
        # missing field when adding atoms
        #ag._customFieldsDefaultValues = {
        #    'data':{'secondary':'C'},
        #    'flags':{},
        #    'bonds':{}
        #    }
        self._ag.setMolecule(self)
        # create a list of flags used to inactivate deleted atoms
        self._ag.setFlags('deleted', [False]*len(ag))

        # create atom's data atomicNumber
        atNums = []
        if self._ag.numCoordsets() > 1:
            self._multi = 'conformations'

        #if self._multi=='molecules':
        #    self._basename = '%s_%s'%(os.path.splitext(os.path.split(filename)[1])[0], self.name)
        #else:
        #    self._basename = self.name

        # handle atom.element=='' by using first character of name
        # when it is ''. This version is twice as fast as the one below
        aa = self._ag.getElements()
        inds = numpy.where(aa == '')
        if len(inds[0]):
            aa[inds] = self._ag.getNames().astype('|S2')[inds]
        for a in aa:
            atNums.append(babel_elements.get(a,babel_elements['Xx'])['num'])

        #atNums1 = []
        #t0 = time()
        #for a,b in zip(self._ag.getElements(), self._ag.getNames()):
        #    if a=='':
        #        a = b[0]
        #    atNums1.append(babel_elements[a]['num'])
        #print 'atnums', time()-t0
        #for a,b in zip(atNums, atNums1):
        #    assert a==b
        self._ag.setData('atomicNumber', atNums)
        self.defaultRadii()
        self.selSet = SelectionSet( [self.select()] )
        self._atomSubset = self.select()
        self.children = self._ag.iterChains
        
## TreeObject methods
##
    def nbChildren(self):
        return len(numpy.unique(self.selSet[0]._getChindices()))
    
    def getChildrenAndNames(self, rootSelection=None):
        children = []
        names = []
        ag = self._ag
        cs = ag.getACSIndex()
        for chain in ag.select('not deleted').getHierView().iterChains():
            chid = chain.getChid()
            segment = chain.getSegname()
            if rootSelection:
                obj = Selection(ag, chain.getIndices(), cs, title=chid)
                #obj = self.selSet & rootSelection
                obj = SelectionSet([obj]) & rootSelection
                if obj.nbAtoms():
                    if segment != '':
                        names.append("%s (%s)"%(chid,segment))
                    else:                    
                        names.append(chid)
                    children.append( Chain(obj, chid) )
            else:
                if segment != '':
                    names.append("%s (%s)"%(chid,segment))
                else:
                    names.append(chid)
                children.append( Chain(chain, chid) )
        return [children, names]

## End TreeObject methods
##

    def getAtomGroup(self):
        return self._ag

    def selectMolKit(self, selStr):
        """select atoms in this molecule using the MolKit syntax.
selStr is expected to be ccc:rrr:aaa where

ccc is space separated strings called tokens identifying chains as follows using:
    1 - chain id e.g. A
    2 - chain index e.g. #1
    3 - range using characters e.g. 1-N
    4 - range using indices e.g. #3to#7
    ranges start and ends are optional and default to first and last

rrr is space separated set of tokens for residues that can be:
    1 - numbers ==> e.g. '10' ==> 'resum 10'
    2 - residue index #1, ==>'resindex 1' 
    3 - ranges of indices ==> 'resindex 1to4'
    4 - ranges of residues numbers ==> 'resnum 0to3'
    5 - wildcards 
          all residues starting with A ':A*:' ==> 'resname "A.*"'
    keywords: see page 58 59 or prody manual
        :polar: ==> "polar"

aaa is space separated set of tokens for atoms that can be 
    1 - string e.g. "CA" ==> "CA"
                    "10" ==> "serial 10"
    2 - atom index #1, ==>"index 1" 
    3 - ranges of indices ==> "index 1to4"
    5 - wildcards 
    keywords: see page 58 59 or prody manual
        ::backbone ==> "backbone"

special characters are exscaped using "", e.g. "*" means character *
indices are 1-based
        """

        levels = selStr.split(':')
        psel = ""
        #import pdb
        #pdb.set_trace()
        for i,l in enumerate(levels):
            if i==0:
                cselstr = getProdyChainSelStr(l, self) 
                if len(cselstr):
                    psel += cselstr + ' '
            elif i==1:
                rselstr = getProdyResSelStr(l) 
                if len(rselstr):
                    psel += rselstr + ' '
            elif i==2:
                aselstr = getProdyAtomSelStr(l) 
                if len(aselstr):
                    psel += aselstr + ' '
        #print 'PSEL', psel
        if len(psel):
            return self.select(psel)
        else:
            return self.select()
                    
    def select(self, selstr=None):
        try:
            if selstr == None:
                return Selection(
                    self._ag,
                    self._ag.select('not deleted').getIndices(),
                    'name ".*"',
                    self._ag.getACSIndex())
            else:
                sel = self._ag.select('(%s) and not deleted'%selstr)
                if sel is not None:
                    return Selection(self._ag, sel._indices, selstr,
                                     self._ag.getACSIndex())
        except SelectionError:
            print "Prody Selection Error:", sys.exc_info()[1]

    def emptySelection(self):
        return Selection(self._ag, [], '', None)

    def getBHT(self):
        try:
            return self._bht
        except AttributeError:
            v = self._ag.getCoords().astype('f')
            self._bht = bhtreelib.BHtree( v, None, 10)
            return self._bht
    
    def updateBondsByDistance(self, factor=1.175):
        if self._ag._bonds is None:
            return self.buildBondsByDistance(factor)

        # if there already are bonds we want to keep them
        # build dict of existing bonds
        bonds = self._ag._bonds
        newBonds = []
        bondOrder = []
        _bo = self._ag._bondOrder
        bondDict = {}
        for i, b in enumerate(self._ag._bonds):
            key = '%d %d'%(b[0], b[1])
            bondDict[key] = i
            newBonds.append( (b[0], b[1]) )
            if _bo is not None:
                bondOrder.append(_bo[key])

        # compute bonds by distance
        bo = self.getBondOrdRadius()
        v = self._ag.getCoords().astype('f')
        bht = bhtreelib.BHtree( v, bo, 10)
        pairs = bht.closePointsPairsInTree(1.175)
        del bht

        # find new bonds
        if len(pairs)==len(bonds): # no new bonds
            return []
        else:
            addBonds = []
            for i,j in pairs:
                key = '%d %d'%(i, j)
                ind = bondDict.get(key, None)
                if ind is None: # new bond
                    newBonds.append( (i,j) )
                    addBonds.append( (i,j) )
                    if _bo is not None:
                        bondOrder.append(1)
            if _bo is None:
                self._ag.setBonds(newBonds)
            else:
                self._ag.setBonds(newBonds, bondOrder)
            # make sure geometry for bondOrder gets updated
            #self.setBondorder(self._ag._bondOrder)

        return addBonds
    
    def buildBondsByDistance(self, factor=1.175):
        bo = self.getBondOrdRadius()
        v = self._ag.getCoords().astype('f')
        bht = bhtreelib.BHtree( v, bo, 10)
        pairs = bht.closePointsPairsInTree(1.175)
        del bht
        if len(pairs)>0:
            self._ag.setBonds(pairs)
        return pairs

    def typeAtomsAndBonds(self, obmol=None):
        """Assign atom types, bond orders and HB donor acceptor status
The molecule is expected to be protonated
After running this command, atoms will have:
   - types (atomType),
   - bonds will have bondorders save in mol.bondorder[(i,j)]
     where i,j are 0 based atom indices with i<j
   - HB donor acceptor properties ('hbType')
        """
        bondOrder = {}
        types = numpy.zeros((self._ag.numAtoms(),), dtype='S5')
        hbTypes = numpy.zeros((self._ag.numAtoms(),), dtype='uint8')
        maxLenType = 0
        from MolKit2.openBabelInterface import ProdyToOBMol
        if obmol is None:
            obmolprovided = False
        else:
            obmolprovided = True

        # loop over chains as creating OBMOL for entire proteins tends to crash
        self.initAromaticRings()
        for chId in numpy.unique(self._ag.getChids()):
            if not obmolprovided:
                if chId == ' ':
                    chId = '_'
                chain = self.select('chain %c'%chId)
                obmol = ProdyToOBMol(chain, title='%s: chain %s'%(self.name, chId))

            if self._ag._bondOrder is None:
                obmol.PerceiveBondOrders()
            
            # get atom type 
            for patom, atm in zip(obmol.prodyAtomList, ob.OBMolAtomIter(obmol)):
                pAtmIndex = patom.getIndex()
                types[pAtmIndex] = atm.GetType()
                maxLenType = max(maxLenType, len(atm.GetType()))
                hbTypes[pAtmIndex] = 1 * atm.IsHbondAcceptor() + \
                                     2 * atm.IsHbondDonor()

            # get bond order
            for bond in ob.OBMolBondIter(obmol):
                i = obmol.obToProdyIndex[bond.GetBeginAtomIdx()-1]
                j = obmol.obToProdyIndex[bond.GetEndAtomIdx()-1]
                if bond.IsAromatic():
                    bo = 4
                else:
                    bo = bond.GetBondOrder()
                if i<j:
                    bondOrder['%d %d'%(i,j)] = bo
                else:
                    bondOrder['%d %d'%(j,i)] = bo

            self.addAromaticRings(obmol=obmol)

        self._ag.setData('atomType', types)
        self._ag.setData('hbType', hbTypes)
        self.setBondorder(bondOrder)
        
    def initAromaticRings(self):
        self._aromaticArcs = {}
        self._aromaticAtomsCounter = 0
        
    def setBondorder(self, bondOrder):
        self._ag._bondOrder = bondOrder
        bonds = self.select().getBonds()
        self.addDoubleBonds(bonds[2])
        sep = self._renderingProp['lines']['doubleBondSep']
        v = self.computeDoubleBondsVertices(bonds[2], sep)
        self.geomContainer.geoms['doubleBonds'].vertexSet.vertices.array = v

        self.addTripleBonds(bonds[3])
        sep = self._renderingProp['lines']['tripleBondSep']
        v = self.computeTripleBondsVertices(bonds[3], sep)
        self.geomContainer.geoms['tripleBonds'].vertexSet.vertices.array = v

        if len(bonds[4]):
            from MolKit2.openBabelInterface import ProdyToOBMol
            self.initAromaticRings()
            obmol = ProdyToOBMol(self.select())
            self.addAromaticRings(obmol=obmol)
            vert, norm, vect, angles, rad = self.computeAromaticRingsArcs()
            self.geomContainer.geoms['aromaticBonds'].Set(
                vertices=vert, vnormals=norm, angles=angles,
                vectors=vect, radii=rad, redo=0, tagModified=False,
                transparent='implicit')
    ##
    ## methods for dealing with bond order
    ##
    def addDoubleBonds(self, bonds):
        # created self._bondOrderData['doubleBonds']["1 2"] = 1
        # indicating that we have double bond between atom 1 and 2
        # for each bond we have a key "a1 a2"
        d = {}
        counter = 0
        for i1,i2 in bonds:
            d["%d %d" %(i1, i2)] = counter
            counter +=1
        self._bondOrderData['doubleBonds'] = d

    def addTripleBonds(self, bonds):
        # build self._bondOrderData['tripleBonds']
        # for each bond we have a key "a1 a2"
        tripleBonds = {}
        counter = 0
        for i1,i2 in bonds:
            tripleBonds["%d %d" %(i1, i2)] = counter
            counter +=1
        self._bondOrderData['tripleBonds'] = tripleBonds

    def addAromaticRings(self, obmol=None, append=False):
        # build self._bondOrderData['aromaticRings']
        # for each ring store _center, _normal, _atomIndices
        if obmol is None:
            if atoms is None:
                raise ValueError("ERROR: addAromaticRings() called with no obmol and no atoms")
            else:
                from MolKit2.openBabelInterface import ProdyToOBMol
                obmol = ProdyToOBMol(atoms)
        if append:
            aromatic = self._bondOrderData['aromaticRings']
            if aromatic is None:
                aromatic = []
        else:
            aromatic = []
        for ring in obmol.GetSSSR(): # get the Smallest Set of Smallest Rings
            if not ring.IsAromatic():
                continue
            center = ob.vector3()
            norm1 = ob.vector3()
            norm2 = ob.vector3()
            ring.findCenterAndNormal(center, norm1, norm2)
            _atomIndices = obmol.obToProdyIndex[[x-1 for x in list(ring._path)]]
            _center = numpy.array( (center.GetX(), center.GetY(), center.GetZ()) )
            _normal = numpy.array( (norm1.GetX(), norm1.GetY(), norm1.GetZ()) )
            aromatic.append([_center, _normal, _atomIndices])
        self._bondOrderData['aromaticRings'] = aromatic

    ##
    ## methods to generate vertices for rendering bond order
    ##
    def computeDoubleBond(self, a1, a2, d=0.12):
        # compute 2 points to draw the double bond between a1, a2
        # find a3 that is bonded to a1 and different from a2
        # d is the separation betweenthe 2 parallel lines
        ag = self._ag
        a1 = ag[a1]
        a2 = ag[a2]
        for a3 in a1.iterBonded():
            if a3.getElement()=='H' or a3.getIndex()==a2.getIndex():
                continue
            break
        x1,y1,z1 = a1.getCoords()
        x2,y2,z2 = a2.getCoords()
        x3,y3,z3 = a3.getCoords()
        # compute v(a1-a2) and v(a1-a3)
        vx1, vy1, vz1 = x2-x1, y2-y1, z2-z1
        vx2, vy2, vz2 = x3-x1, y3-y1, z3-z1
        # compute inverse of norms of v1 and v2
        from math import sqrt
        n1 = 1./sqrt((vx1*vx1 + vy1*vy1 + vz1*vz1))
        n2 = 1./sqrt((vx2*vx2 + vy2*vy2 + vz2*vz2))
        # normalize v1 v2
        vx1, vy1, vz1 = vx1*n1, vy1*n1, vz1*n1
        vx2, vy2, vz2 = vx2*n2, vy2*n2, vz2*n2
        # compute cross product: v1xv2
        vx3 = vy1*vz2 - vz1*vy2
        vy3 = vz1*vx2 - vx1*vz2
        vz3 = vx1*vy2 - vy1*vx2
        # check if  v3 vector is 0 or close to 0:
        if (vx3*vx3 + vy3*vy3 + vz3*vz3) < 0.1:
            # find orthogonal vector to (vx1, vy1, vz1):
            if vz1 == 0 and -vx1-vy1 == 0:
                vx3 = -vy1-vz1
                vy3 = vx1
                vz3 = vx1
            else:
                vx3 = vz1
                vy3 = vz1
                vz3 = -vx1-vy1
            n = 1./sqrt((vx3*vx3 + vy3*vy3 + vz3*vz3))
            vx3, vy3, vz3 = vx3*n, vy3*n, vz3*n
        # compute cross product v1xv3
        vx4 = vy1*vz3 - vz1*vy3
        vy4 = vz1*vx3 - vx1*vz3
        vz4 = vx1*vy3 - vy1*vx3

        return [(x1+vx4*d,y1+vy4*d,z1+vz4*d), (x2+vx4*d,y2+vy4*d,z2+vz4*d),
                (x1-vx4*d,y1-vy4*d,z1-vz4*d), (x2-vx4*d,y2-vy4*d,z2-vz4*d)]

    def computeDoubleBondsVertices(self, bonds, separation):
        # compute vertices used to draw all double bonds and store in
        # self._bondOrderData['doubleBondsVertices']
        # a double bond between atom a1 and a2 generates 4 vertices
        #
        #       x3 ---------- x4
        #       a1            a2
        #       x5 ---------- x6
        #
        vertices = []
        for i1,i2 in bonds:
            vertices.extend(self.computeDoubleBond(i1, i2, d=separation))
        return numpy.array(vertices, 'f')

    def computeTripleBond(self, a1, a2, d=0.2):
        # compute vertices for a single tripple bond by computing
        # 2 offset lines like double bounds
        a,b,c,d = self.computeDoubleBond(a1, a2, d)
        return (a, b, c, d)

    def computeTripleBondsVertices(self, bonds, separation):
        # compute vertices used to draw all triple bonds and store in
        # self._bondOrderData['tripleBondsVertices']
        # a triple bond between atom 1 and 2 generates 4 vertices
        # for bonds flanking the central bond drawn as signle bond
        #       x3 ---------- x4
        #       x1 ---------- x2
        #       x5 ---------- x6
        #
        vertices = []
        for i1,i2 in bonds:
            vertices.extend(self.computeTripleBond(i1, i2, separation))
        return numpy.array(vertices, 'f')

    def computeAromaticRingsArcs(self):
        coords = self._ag.getCoords()
        vertices = []
        normals = []
        vectors = []
        angles = []
        radii = []
        aromaticArcs = {}
        # for each atom we self._bondOrderData['aromaticArcs'][ayomIndex] will contain a
        # list of arc parameters, one arc for each aromatic ring the atom belongs to.
        # The parameters or one arc are:
        #     center, normal, start vector, angle, radius
        # This information can be used to create a Fan3D geometry for line drawings
        # or turned into a set of points for S&B rendering of aromatic circle
        #
        for _center, _normal, _atomIndices in self._bondOrderData['aromaticRings']:
            #import pdb; pdb.set_trace()
            #print 'RING', _atomIndices
            # compute vertices, vectors, angles , etc for the ring's arcs
            length = len(_atomIndices)
            angle = 360./length
            rcoords = coords[_atomIndices]
            if length==5:
                radius = .6
            else:
                radius = .8
            for i in range(-1, length-1):
                p1 = rcoords[i]
                p2 = rcoords[i+1]
                # compute arc start vector
                sx1, sy1, sz1 = p1-_center + p2-_center
                # compute inverse of norms of start vector
                n1 = 1./sqrt((sx1*sx1 + sy1*sy1 + sz1*sz1))
                # normalize v1 v2
                sx1, sy1, sz1 = sx1*n1, sy1*n1, sz1*n1
                if aromaticArcs.has_key(_atomIndices[i]):
                    aromaticArcs[_atomIndices[i]].append(
                        (self._aromaticAtomsCounter, _center, _normal, (sx1, sy1, sz1), angle, radius) )
                else:
                    aromaticArcs[_atomIndices[i]] = [
                        (self._aromaticAtomsCounter, _center, _normal, (sx1, sy1, sz1), angle, radius ) ]
                vertices.append(_center)
                normals.append(_normal)
                vectors.append((sx1, sy1, sz1))
                angles.append(angle)
                radii.append(radius)
                self._aromaticAtomsCounter += 1
        self._bondOrderData['aromaticArcs'] = aromaticArcs
        return numpy.array(vertices, 'f'), normals, vectors, angles, radii

    def setDoubleBondSep(self, value):
        assert isinstance(value, float), "ERROR: double bond separation has to be afloating point value, got%g"%value
        bonds = self.select().getBonds()[2]
        self._renderingProp['lines']['doubleBondSep'] = value
        v = self.computeDoubleBondsVertices(bonds, value)
        self.geomContainer.geoms['doubleBonds'].vertexSet.vertices.array = v

    def setTripleBondSep(self, value):
        assert isinstance(value, float), "ERROR: triple bond separation has to be afloating point value, got%g"%value
        bonds = self.select().getBonds()[3]
        self._renderingProp['lines']['tripleBondSep'] = value
        v = self.computeTripleBondsVertices(bonds, value)
        self.geomContainer.geoms['tripleBonds'].vertexSet.vertices.array = v


    ## def addAromaticRings(self, obmol):
    ##     # get aromatic ring info and stores in in self._aromaticRings
    ##     aromatic = []
    ##     coords = self._ag.getCoords()
    ##     vertices = []
    ##     normals = []
    ##     vectors = []
    ##     angles = []
    ##     radii = []
    ##     # for each tom index we dict will contain a list of arc parameters, one arc
    ##     # for each aromatic ring the atom belongs to. The parameters are:
    ##     #     center, normal, start vector, angle, radius
    ##     for ring in obmol.GetSSSR(): # get the Smallest Set of Smallest Rings
    ##         if not ring.IsAromatic():
    ##             continue
    ##         center = ob.vector3()
    ##         norm1 = ob.vector3()
    ##         norm2 = ob.vector3()
    ##         ring.findCenterAndNormal(center, norm1, norm2)
    ##         _atomIndices = obmol.obToProdyIndex[[x-1 for x in list(ring._path)]]
    ##         _center = numpy.array( (center.GetX(), center.GetY(), center.GetZ()) )
    ##         _normal = numpy.array( (norm1.GetX(), norm1.GetY(), norm1.GetZ()) )
    ##         aromatic.append([_center, _normal, _atomIndices])

    ##         # compute vertices, vectors, angles , etc for the ring's arcs
    ##         length = len(_atomIndices)
    ##         angle = 360./length
    ##         rcoords = coords[_atomIndices]
    ##         if length==5:
    ##             radius = .6
    ##         else:
    ##             radius = .8
    ##         for i in range(-1, length-1):
    ##             p1 = rcoords[i]
    ##             p2 = rcoords[i+1]
    ##             # compute arc start vector
    ##             sx1, sy1, sz1 = p1-_center + p2-_center
    ##             # compute inverse of norms of start vector
    ##             n1 = 1./sqrt((sx1*sx1 + sy1*sy1 + sz1*sz1))
    ##             # normalize v1 v2
    ##             sx1, sy1, sz1 = sx1*n1, sy1*n1, sz1*n1
    ##             if self._aromaticArcs.has_key(_atomIndices[i]):
    ##                 self._aromaticArcs[_atomIndices[i]].append(
    ##                     (self._aromaticAtomsCounter, _center, _normal, (sx1, sy1, sz1), angle, radius) )
    ##             else:
    ##                 self._aromaticArcs[_atomIndices[i]] = [
    ##                     (self._aromaticAtomsCounter, _center, _normal, (sx1, sy1, sz1), angle, radius ) ]
    ##             vertices.append(_center)
    ##             normals.append(_normal)
    ##             vectors.append((sx1, sy1, sz1))
    ##             angles.append(angle)
    ##             radii.append(radius)
    ##             self._aromaticAtomsCounter += 1
        
    ##     self._aromaticRings.extend(aromatic)
    ##     return numpy.array(vertices, 'f'), normals, vectors, angles, radii
        

    ## def doubleBondsGeometry(self):
    ##     # compute vertices used to draw all double bonds and store in
    ##     # self._bondOrderData['doubleBondsVertices']
    ##     # a triple bond between atom a1 and a2 generates 4 vertices
    ##     #
    ##     #       x3 ---------- x4
    ##     #       a1            a2
    ##     #       x5 ---------- x6
    ##     #
    ##     bonds = self.select().getBonds()
    ##     vertices = []
    ##     for i1,i2 in bonds[2]:
    ##         vertices.extend(self.computeDoubleBond(
    ##             i1, i2, d=self._bondOrderData['doubleBondLineOffset']))
    ##     self._bondOrderData['doubleBondsVertices'] = vertices

    def getBondOrdRadius(self):
        # return a list of bond or der radii for all atoms
        return self.bo_rad[self._ag.getData('atomicNumber')]

    def defaultRadii(self, atomSet=None):#, united=None, overwrite=0):
        """Assign atom radii to the selected atoms

        radiiLst <- defaultRadii(self, atomSet=None, united=None, overwrite=0)

if atomSet is None, all atoms in the molecule are used
if overwrite is true pre-exiting atom.radius attribute will be overwritten
if united is true large atomic radii are used for atoms which have no hydrogens
"""
        if atomSet==None:
            atomSet = self._ag
        radii = self.vdw_rad[atomSet.getData('atomicNumber')]
        atomSet.setRadii(radii)
        return radii

    def measureCHIs(self, chain, resnum, resname=None, rotlib=None):
        if rotlib is None:
            from MolKit2.rotamerLib import RotamerLib
            rotlib = RotamerLib()
        resAtoms = self.select('chain %c resnum %d'%(chain, resnum))
        if resname is None:
            resname = resAtoms.getResnames()[0]
        if resname.upper() not in rotlib.residueNames:
            return [], resname.upper()
#            raise ValueError, "not rotamer for residue name %s (%d)"%(resname, resnum)
        chiDef = rotlib.getAngleDef(resname)
        angles = []
        for angleDef, movedAtoms in chiDef:
            if len(angleDef):
                try:
                    atoms = resAtoms.select('name '+' '.join(angleDef))
                except KeyError, e:
                    if 'ND1' in angleDef:
                        angleDef1 =  angleDef[:3]+['ND2']
                        atoms = resAtoms.select('name '+' '.join(angleDef1))
                    else:
                        raise
                if len(atoms) == 4:                    
                    angles.append(torsion(*atoms.getCoords()))
                else:
                    angles.append(180.0)
        return angles, resname
    
    ## closest rotamer using CHI angles
    ##
    def closestRotamerCHI(self, angles, resname, rotlib):
        #compute closest rotamer usign CHI angles
        rotAngleList = rotlib.getAngles(resname)
        rotDevList = rotlib.getAnglesDev(resname)
        mini = 360*len(rotAngleList)
        for i, rang in enumerate(rotAngleList):
            dev = []
            devsum = 0.0
            for a1, a2 in zip(angles, rang):
                dev.append(min(abs(a1-a2), abs(a1+360-a2), abs(a1-360-a2)))
                devsum += dev[-1]
            #print devsum, angles, rang, dev
            if devsum < mini:
                mini = devsum
                minidev = dev
                miniind = i
        return miniind, minidev, rotAngleList[miniind], rotDevList[miniind]

    def getPDBRecords(self):
        stream = Stream()
        prody.writePDBStream(stream, self._ag)
        return stream.lines

    ## molecule manipulation methods
    ##
    def deleteAtoms(self, sele):
        flags = self._ag.getFlags('deleted')
        flags[sele.getIndices()] = True
        self._ag.setFlags('deleted', flags)
        deletedIndices = {}.fromkeys(sele.getIndices(), -1)

        # fix the bonds
        newBonds = []
        newBo = []
        bo = self._ag._bondOrder
        for i,j in self._ag._bonds:
            if deletedIndices.get(i, None) is None and \
               deletedIndices.get(j, None) is None: 
                newBonds.append((i,j))
                if bo:
                    k = '%d %d'%(i,j)
                    newBo.append(bo[k])

        if len(newBonds)>0:
            if len(newBo):
                self._ag.setBonds(newBonds, newBo)
            else:
                self._ag.setBonds(newBonds)
                self._ag._bondOrder = None
        else:
            self._ag._bonds = None
            self._ag._bondOrder = None
            self._ag._bmap = None
            
    def undeleteAtoms(self, sele):
        flags = self._ag.getFlags('deleted')
        flags[sele.getIndices()] = False
        self._ag.setFlags('deleted', flags)

    def clone(self, name=None):
        return self.select().toMolecule(name)
    
    def addAtoms(self, sele):
        ## NOT USED ANYMORE
        raise RuntimeError, "Molecule.AddAtoms is deprecated, use _ag.extend(sele) instead"
        ## extend self._ag with atoms from sele
        ## if sele has some data or flags fields missing they are
        ## initialized with self._ag._customFieldsDefaultValues if found
        ## else we raise an exception
        from prody.atomic.functions import SAVE_SKIP_POINTER as SKIP
        
        selag = sele.getAtomGroup()
        ag = self._ag
        length = len(ag)
        ag._n_atoms = length + len(sele)
        
        # FIXME handle multiple coordinate sets
        assert ag.numCoordsets() ==  sele.numCoordsets(), \
               "ERROR: cannto add atoms with %d coordsets to molecule with %d coordsets"%(
            ag.numCoordsets(), sele.numCoordsets())
        nag = prody.AtomGroup(ag._title)
        nag._setCoords(numpy.concatenate((ag.getCoords(), sele.getCoords())))
        #nag._customFieldsDefaultValues = ag._customFieldsDefaultValues
        ##
        ## handle bonds
        ##
        bonds1 = self._ag._bonds
        bonds2 = sele.getBonds()[1]
        if bonds1 is None:
            if len(bonds2)==0:
                pass
            else:
                # the bonds are the bonds added with an offset
                bo = selag._bondOrder
                if bo is not None:
                    newBo = {}
                    for i,j in bonds2:
                        newBo['%d %d'%(i+length,j+length)] = bo['%d %d'%(i,j)]
                else:
                    bo = None
                nag.setBonds( numpy.array(bonds2)+length, bo )

        else: # self has bonds
            if len(bonds2)>0:
                # map bonds indices of sele from 0 to n
                size = len(selag)
                old = indices = numpy.arange(0, size)
                a1 = numpy.ones(size)*-1
                tokeep = sele.getIndices()
                a1[tokeep] = indices[tokeep]
                mapping = numpy.ones(size,'i')*-1
                maxi = 0
                for i in xrange(len(a1)):
                    if a1[i]!=-1:
                        mapping[i] = maxi
                        maxi += 1
                newBonds = mapping[bonds2]+length

                bo = selag._bondOrder
                if bo is not None:
                    newBo = bo.copy()
                    for bn in xrange(len(bonds2)):
                        i,j = bonds[bn]
                        ni,nj = newBonds[bn]
                        newBo['%d %d'%(ni,nj)] = bo['%d %d'%(i,j)]
                else:
                    bo = None

                nag.setBonds( numpy.concatenate((bonds1, newBonds)), bo )
            else:
                bondOrder = []
                bo = self._ag._bondOrder
                if bo is None:
                    nag.setBonds( bonds1)
                else:
                    for b in bonds1:
                        bondOrder.append(bo['%d %d'%(b[0], b[1])])
                    nag.setBonds( bonds1, bondOrder )

        ##
        ## data 
        ##
        for name in ag.getDataLabels():
            if name in SKIP:
                continue
            data1 = ag._getData(name)
            data2 = selag._getData(name)
            if data2 is None:
                # name is missing in atoms
                defaultValue = ag._customFieldsDefaultValues['data'].get(name, None)
                if defaultValue is None:
                    raise RuntimeError, "No default value for field %s"%name
                data2 = numpy.array( [defaultValue]*len(selag), data1.dtype)

            if name=='serials': # serials have to be unique
                data2 = data2 + max(data1)
            elif name=='secondary': # serials have to be unique
                selag.setSecstrs(['C']*len(selag))
            nag.setData(name, numpy.concatenate((data1, data2)) )

        ##
        ## flags
        ##
        for name in ag.getFlagLabels():
            if name in SKIP:
                continue
            data1 = ag._getFlags(name)
            data2 = selag._getFlags(name)
            if data2 is None:
                # name is missing in atoms
                defaultValue = ag._customFieldsDefaultValues['flags'].get(name, None)
                if defaultValue is None:
                    raise RuntimeError, "No default value for field %s"%name
                data2 = numpy.array( [defaultValue]*len(selag), data1.dtype)
            nag._setFlags(name, numpy.concatenate( (data1, data2) ))
                
        ##
        ## bondData 
        ##
        for name in ag._bondData.keys():
            if name in SKIP:
                continue
            data1 = ag._bondData.get(name)
            data2 = selag._bondData.get(name, None)
            if data2 is None:
                # name is missing in atoms
                defaultValue = ag._customFieldsDefaultValues['bonds'].get(name, None)
                if defaultValue is None:
                    raise RuntimeError, "No default value for field %s"%name
                data2 = numpy.array( [defaultValue]*len(selag), data1.dtype)

            nag._bondData[name] = numpy.concatenate((data1, data2))

        self._ag = nag
        nag.setMolecule(self)
        self._atomSubset = self.select()
        self.updateBondsByDistance()

    def isFromPDB(self, filename):
        ext = os.path.splitext(filename)[-1];
        if ext == '.pdb':
            return True
        else:
            return False
        
    def biologicalUnit(self, which='all'):
        """If the molecule was built from a PDB file that contains BIOMT
        records this function will return a new molecule constructed from
        the current one using the BIOMT reecords.
        chains will we rename starting from Chain A
        which specifies the biological unit to be built. By default which is
        set to 'all' and this function returns a list of molecules, one for each
        biological unit.
        """
        if self.filename is None:
            raise ValueError('Molecule has no originating filename information')

        if not self.isFromPDB(self.filename):
            raise ValueError('Molecule is not from a pdb file')
        
        if not hasattr(self, 'pdbHeader'):
            self.pdbHeader = prody.parsePDB(self.filename, model=0, header=True)

        if which is 'all':
            bioMolNames = self.pdbHeader['biomoltrans'].keys()
        else:
            assert str(which) in self.pdbHeader['biomoltrans'].keys(), \
                   "Biomolecule %s not recognized expected one of %s"%(
                str(which), self.pdbHeader['biomoltrans'].keys())
            bioMolNames = [str(which)]

        bioMols = []
        for k in bioMolNames:
            chCount = 65 # 'A' ascii code

            biotrans = self.pdbHeader['biomoltrans'][k]
            # select the atom set to which to apply the transformation
            selStr = "chain "+' '.join(biotrans[0])
            bioUnitAtoms = self.select(selStr)
            nbTrans = (len(biotrans)-1)/3
            nmol = None
            for i in range(nbTrans):
                tmol = bioUnitAtoms.toMolecule('%s_bioUnit_%s'%(self.name, k))
                tmol.buildBondsByDistance()
                tmol.defaultRadii()
                # fix chain ids to be consecutive
                chains = tmol._ag.getChids()
                for chid in biotrans[0]: # loop over original chain ids
                    chains[chains==chid] = chr(chCount)
                    chCount += 1
                tmol._ag._data['chain'][:] = chains

                # transform coordinates
                coords = tmol._ag.getCoords()
                coords1 = numpy.ones( (len(coords), 4), 'f')
                coords1[:, :3] = coords
                mat = numpy.identity(4)
                mat[:3,:] = [ [float(x) for x in biotrans[i*3+1].split()],
                              [float(x) for x in biotrans[i*3+2].split()],
                              [float(x) for x in biotrans[i*3+3].split()]]

                ncoords = numpy.dot(coords1, numpy.transpose(mat))
                tmol._ag._setCoords( ncoords[:, :3] )

                if nmol is None:
                    nmol = tmol
                else:
                    nmol.addAtoms( tmol.select() )
            bioMols.append(nmol)
        return bioMols

    def assignSecondaryStructureWithPross(self):
        """methiod to assign secondary structure information using PROSS
        This method uses much of the code from the original BIOMOL collection of utilities 
	written by Raj Srinivasani with enhancements by Nick Fitzkee.

        The script was put together by Pat Fleming so that a user would not need
	to have the BIOMOL distribution installed to run PROSS.

        Note that since Raj's time the definitions of mesostates has been superceded
        by the fine grained 30 deg x 30 deg grid for most purposes. Either mesostate
        grid will work for PROSS. Give your choice as an argument (see USAGE below).

        Date: September 2004
        Author: Pat Fleming, pat.fleming@jhu.edu
        http://www.roselab.jhu.edu/utils/pross.html"""
        self._ag.setSecstrs(['C']*self._ag.numAtoms())
        SSassigner = AssignSecondaryStructureWtihPross()
        SSassigner(self)

    def getChainCtrlPointsForNA(self, ctrlAtoms, chain, allNAResnums):
        """
        builds a list of control points for a spline

        the ctrlAtoms are a set of consecutive P atoms (i.e. with no gaps) from chain "chain"
        allNAResnums is sthe list of all Nucleotides in chain A

        In addition to the P atoms from ctrlAtoms, the list of ontrol points will include O5' before and O3'
        at the end. We also check for nucleotide before and after the nucleotides corresponding to the
        first and last P in the ctrlAtoms list. If these nucleotides have no P we use the O5' of the nucleotide
        before and the and O3' of the one after. Else we use these atoms in the nucleotice corresponding to
        the first and last P in the ctrlatoms.

        ctrlAtoms returns the list of atoms described above
        ctrl returns the 3D coordinates of ctrlAtoms with 2 points added before and after (used to compute
        cubic splines).
        """
        resnums = ctrlAtoms.getResnums()
        #resnums = chain.select("nucleotide and not deleted").getResnums()
        #if chain.select("nucleotide element P resnum %d and not deleted"%resnums[0]) is None:

        #import pdb; pdb.set_trace()
        indFirstNA = allNAResnums.index(resnums[0])
        # if there is a NA before that has not P check for O5' adn add it
        if indFirstNA > 0 and \
               chain.select("nucleotide name P resnum %d and not deleted"%allNAResnums[indFirstNA-1]) is None:
           
            o5 = chain.select("nucleotide name O5' resnum %d and not deleted"%allNAResnums[indFirstNA-1])
            if o5:
                o5c = o5.getCoords()[0]
                ctrlAtoms._indices = numpy.array( o5.getIndices().tolist() + ctrlAtoms._indices.tolist())

        else: # check for O5' in resnums[0] and add it
            o5 = chain.select("nucleotide name O5' resnum %d and not deleted"%resnums[0])
            if o5:
                o5c = o5.getCoords()[0]
                ctrlAtoms._indices = numpy.array( o5.getIndices().tolist() + ctrlAtoms._indices.tolist())
            
        indLastNA = allNAResnums.index(resnums[-1])
        # if there is a NA after and it has not P check for O3'
        if indLastNA < len(allNAResnums) and \
               chain.select("nucleotide name P resnum %d and not deleted"%allNAResnums[indLastNA]) is None:

            o3 = chain.select("nucleotide name O3' resnum %d and not deleted"%allNAResnums[indLastNA])
            if o3:
                o3c = o3.getCoords()[0]
                ctrlAtoms._indices = numpy.array( ctrlAtoms._indices.tolist() + o3.getIndices().tolist())

        else: # check for O3' in resnums[-1] and add it
            o3 = chain.select("nucleotide name O3' resnum %d and not deleted"%resnums[-1])
            if o3:
                o3c = o3.getCoords()[0]
                ctrlAtoms._indices = numpy.array( ctrlAtoms._indices.tolist() + o3.getIndices().tolist())

        ## add a ctrl point before and after
        ctrl = numpy.zeros( (len(ctrlAtoms)+2, 3), 'd')
        ctrl[1:-1] = ctrlAtoms.getCoords()
            
        # add ctrl points before O5' and after O3'
        vec = [ctrl[2][0]-ctrl[1][0],
               ctrl[2][1]-ctrl[1][1],
               ctrl[2][2]-ctrl[1][2]]
        ctrl[0] = [ctrl[1][0]-vec[0], ctrl[1][1]-vec[1], ctrl[1][2]-vec[2]]

        n = len(ctrl)
        vec = [ctrl[n-2][0]-ctrl[n-3][0],
               ctrl[n-2][1]-ctrl[n-3][1],
               ctrl[n-2][2]-ctrl[n-3][2]]
        ctrl[-1] = [ctrl[n-2][0]+vec[0], ctrl[n-2][1]+vec[1], ctrl[n-2][2]+vec[2]]
        return ctrl, ctrlAtoms
        
import mmap, weakref, threading
MOLECULE_ATOM_LABEL = 1
MOLECULE_RESIDUE_LABEL = 2

class MultiMolecule(Molecule):
    """
    The MultiMolecules class is a class for representing multi molecule files as a Molecule

    It add method for indicating the number of molecules and jumping to a molecule
    """

    # rendering bits, set when the entire molecule is rendered as:
    _LINES = 1    # lines
    _SB = 2       # sticks and balls
    _CPK = 4      # CPK
    _CARTOON = 8  # cartoon
    _SURFACE = 16 # molecular surface
    _COARSE_SURFACE = 32 # coarse molecular surface
    _SECONDARY_STRUCTURE = 64  # ???
    
    # labeling bits
    _NO_LABELS = 0
    _ATOM_LABEL = MOLECULE_ATOM_LABEL
    _RESIDUE_LABEL = MOLECULE_RESIDUE_LABEL
    _CHAIN_LABEL = 4
    _MOLECULE_LABEL = 8
    
    def __del__(self):
        if self._fp is not None:
            self._fp.close()

    def __init__(self, filename, fileFormat='auto', group=None, indexingMode='foreground'):

        self._renderingBits = self._LINES
        self._labelingBits = self._NO_LABELS
        
        self._fp = None
        if fileFormat is 'auto':
            self.fileFormat = os.path.splitext(filename)[1][1:]
        else:
            self.fileFormat = fileFormat
        assert self.fileFormat in ['mol2','sdf'], "ERROR: bad fileFormat parameter, got%s expected 'mol2' or 'sdf'"%self.fileFormat
        self._fp = open(filename)
        self.data = mmap.mmap(self._fp.fileno(), 0, access=mmap.ACCESS_READ)
        self.numberOfMolecules = -1
        self.index = []
        self._hasIndex = False

        self.group = group

        ag = self.getFirstAtomgroup()
        Molecule.__init__(self, ag.getTitle(), ag, filename=filename)
        self.currentMoleculeIndex = 0
        self.filename = filename
        #self.name will be molecule name in multi molecule file
        #self._basename will be the name of the file without path or extension
        self.alias = '%s 1/?'%(self.name,)
        #bi = BuildIndex(self)
        #thread.start_new_thread( bi.run, ())
        if indexingMode=='background':
            t1 = threading.Thread(target=self.buildIndex)
            t1.start()
            self._indexingThread = (t1, "indexing file %s"%filename)
        else:
            self.buildIndex()
        
    def getFirstAtomgroup(self):
        if self.fileFormat=='mol2':
            # get the first molecule 'by hand', i.e. without using index
            end = self.data.find("@<TRIPOS>MOLECULE\n", len('@<TRIPOS>MOLECULE\n'))-1
            source = self.data[0:end]
            return self.agFromOBmolSource(source) # prody molecule for current index
        elif self.fileFormat=='sdf':
            # molecule ends with "$$$$\n"
            end = self.data.find("$$$$")
            source = self.data[0:end]
            return self.agFromOBmolSource(source)

    def agFromOBmolSource(self, source):
        self.obmol = obmol = ob.OBMol()
        obconversion = ob.OBConversion()
        if self.fileFormat=='mol2':
            formatok = obconversion.SetInFormat('mol2')
        elif self.fileFormat=='sdf':
            formatok = obconversion.SetInFormat('sdf')
        success = obconversion.ReadString(obmol, source)
        if not success:
            raise IOError("Failed to convert '%s' to format %s" % (
                source, self.fileFormat))

        from MolKit2.openBabelInterface import OBMolToProdyAG
        return OBMolToProdyAG(obmol)

    def gotoMolecule(self, idx):
        idx = min(idx, len(self.index)-1)
        idx = max(0, idx)
        if idx==len(self.index)-1:
            end = -1
        else:
            end = self.index[idx+1]-1
        source = self.data[self.index[idx]:end]
        ag = self.agFromOBmolSource(source) # prody molecule for current index
        self.setAtomgroup(ag)
        self.currentMoleculeIndex = idx
        
    def buildIndex(self):
        #FIXME change to use file format
        t0 = time()
        #print 'FUGU', resourceFolder, self.filename
        cacheFolder = getCacheFolder()
        indexFilename = os.path.join(cacheFolder, os.path.splitext(os.path.basename(self.filename))[0]+'.idx')
        buildIndex= False
        if os.path.exists(indexFilename):
            f = open(indexFilename)
            modeTimeFromCache = float(f.readline())
            if abs(os.path.getmtime(self.filename)-modeTimeFromCache)>0.01:
                #print 'modeTimeFromCache', modeTimeFromCache, os.path.getmtime(self.filename), modeTimeFromCache-os.path.getmtime(self.filename)
                print 'WARNING: cache file %s found but modification time does not match. Re-indexing'%self.filename
                buildIndex = True
            else:
                self.index = [int(x) for x in f.readline().split()]
                self._hasIndex = True
        else:
            buildIndex = True
            
        if buildIndex: # build cache
            start = 0
            if self.fileFormat=='mol2':
                # each molecule starts with "@<TRIPOS>MOLECULE\n"
                ct = 0
                while True:
                    start = self.data.find("@<TRIPOS>MOLECULE\n", start)
                    if start == -1: break
                    self.index.append(start)
                    start += 18
                    ct += 1
                    if ct%1000==0:
                        sleep(0.001)
            elif self.fileFormat=='sdf':
                # each molecule ends with "$$$$\n"
                ct = 0
                self.index = [0]
                while True:
                    end = self.data.find("$$$$\n", start)
                    if end==-1: break
                    self.index.append(end+5)
                    start = end+5
                    ct += 1
                    if ct%1000==0:
                        sleep(0.001)                
            else:
                return
            self._hasIndex = True
            print 'MSG: time to index', time()-t0

            # save index in file cache
            try:
                f = open(indexFilename, 'w')
                print 'MSG: writing index file:', indexFilename, os.path.getmtime(self.filename)
            except IOError:
                print 'WARNING: Failed to save index cache %', indexFilename
                return # failed to write cache
            modTime = os.path.getmtime(self.filename)
            f.write('%f\n'%modTime) # write the file last modification time for sanity check
            dum = [f.write('%i '%x) for x in self.index]
            f.close()
        else:
            self._hasIndex = True

class MoleculeSet(list):

    def __init__(self, name, itr=None):
        if itr:
            list.__init__(self, [x for x in itr])
            
        self.name = name

    def getMolecules(self, selStr):
        """selStr can be a comma separated list ot tokens where tokens can be
    - molecule names
    - 1-based indices in the molecule list (#1)
    - 1-based indices range in the molecule list (#4to#6)
    - regular expressions on molecule names
    if selStr is '' all molecules are selected
        """
        tokens = selStr.split(' ')
        molDict = {}
        for token in tokens:
            if is_re.match(token):
                pat = re.compile(token)
                for mol in self:
                    if pat.match(mol.name):
                        molDict[mol] = 1
            elif token[0]=='#':
                match = is_range.match(token)
                if match: # range
                    _from, _to = token.split('to') # make it a range of 0-based serial numbers
                    for mol in self[int(_from[1:])-1:int(_to[1:])]:
                        molDict[mol] = 1
                else:
                    molDict[self[int(token[1:])-1]] = 1
            else: # keep the range of residues numbers
                for mol in self:
                    if mol.name==token:
                        molDict[mol] = 1
                        break
        return MoleculeSet('selection', molDict.keys())

    def selectMolKit(self, selstr=''):
        """select molecular fragments. Returns a SelectionSet
        """
        finalResult = SelectionSet(name=selstr)
        # split using operators /+/ /-/ /&/ and /^/
        i = -1
        words = splitOpPattern.split(selstr)
        #import pdb
        #pdb.set_trace()
        for selector in words[::2]: # skip operators
            if i==-1:
                operator = '+'
            else:
                operator = words[i][1]
            selector = selector.strip()
            levels = selector.split(':')
            result = SelectionSet(name=selstr)

            # first get the molecules specidfied at the top level
            if levels[0]=='':
                mols = self
            else:
                mols = self.getMolecules(levels[0])

            for mol in mols:
                mselStr = ''
                if len(levels)>1:
                    for l in levels[1:]:
                        mselStr+= l+':'
                    mselStr = mselStr[:-1]
                sel = mol.selectMolKit( mselStr )
                if sel:
                    result.append( sel) 
            if operator=='+': 
                finalResult = finalResult.__or__(result)
            elif operator=='-': 
                finalResult = finalResult.__sub__(result)
            elif operator=='&': 
                finalResult = finalResult.__and__(result)
            #elif operator=='^': 
            #    finalResult = finalResult ^ result
            i += 2
        return finalResult

    def select(self, selStr):
        result = []
        selections = selStr.split(",")
        for item in selections:
            ss = item.split(":")
            if len(ss) == 2:
                molname, sel = ss
            else:
                molname = None
                sel = ss[0]
            for mol in self:
                molsel = []
                if molname is not None:
                    if molname == mol.name:
                        molsel = mol.select(sel)
                else:
                    molsel = mol.select(sel)
                if len(molsel):
                    result.append(molsel)        
        return SelectionSet(result, name=selStr)

        
class Chain(TreeObject):

    def __init__(self, pchain, name):
        self.name = name
        self.alias = ''
        if isinstance(pchain, SelectionSet):
            self.selSet = pchain
            self._atomSubset = pchain[0].getHierView().iterChains().next()
        elif isinstance(pchain, ProDyChain):
            self._atomSubset = pchain
            self.selSet = SelectionSet( [self.select()] )
        else:
            raise ValueError("Chain needs a SelectionSet or prody.chain.Chain instance")

    def select(self, selstr=None):
        if selstr:
            sel = self._atomSubset.select(selstr)
            return Selection(self._atomSubset.getAtomGroup(),
                             sel._indices, selstr,
                             self._atomSubset.getACSIndex())
        else:
            return Selection(self._atomSubset.getAtomGroup(),
                             self._atomSubset.getIndices(),
                             self._atomSubset.getSelstr(),
                             self._atomSubset.getACSIndex())

    def getAtomGroup(self):
        return self.mol.getAtomGroup()
    
## TreeObject methods
##
    def nbChildren(self):
        return len(numpy.unique(self.selSet[0]._getResindices()))

    def getChildrenAndNames(self, rootSelection=None):
        children = []
        names = []
        ag = self._atomSubset
        cs = ag.getACSIndex()
        for res in ag.select('not deleted').getHierView().iterResidues():
            resname = res.getResname()
            resnum = res.getResnum()
            name = '%s%d'%(resname, resnum)
            if rootSelection:
                obj = Selection(ag.getAtomGroup(), res.getIndices(), cs, title=resname)
                obj = SelectionSet([obj]) & rootSelection
                if obj.nbAtoms():
                    names.append(name)
                    children.append( Residue(obj, name) )
            else:
                names.append(name)
                children.append( Residue(res, name) )
        return children, names

class ChainSet(list):
    pass

class Residue(TreeObject):

    def __init__(self, presidue, name):
        self.name = name
        self.alias = ''
        if isinstance(presidue, SelectionSet):
            self.selSet = presidue
            self._atomSubset = presidue[0].getHierView().iterResidues().next()
        elif isinstance(presidue, ProDyResidue):
            self._atomSubset = presidue
            self.selSet = SelectionSet( [self.select()] )
        else:
            raise ValueError("Residue needs a SelectionSet or prody.residue.Resdiue instance")

    def getAtomGroup(self):
        return self.mol.getAtomGroup()
    
    def select(self, selstr=None):
        if selstr:
            sel = self._atomSubset.select(selstr)
            return Selection(self._atomSubset.getAtomGroup(),
                             sel._indices, selstr,
                             self._atomSubset.getACSIndex())
        else:
            return Selection(self._atomSubset.getAtomGroup(),
                             self._atomSubset.getIndices(),
                             self._atomSubset.getSelstr(),
                             self._atomSubset.getACSIndex())

    def nbChildren(self):
        return len(self.selSet[0])

    def getChildrenAndNames(self, rootSelection=None):
        children = []
        names = []
        ag = self._atomSubset
        cs = ag.getACSIndex()
        for atom in ag.select('not deleted'):#.getHierView().iterAtoms():
            name = atom.getName()
            if rootSelection:
                obj = Selection(ag.getAtomGroup(), atom.getIndices(), cs, title=name)
                obj = SelectionSet([obj]) & rootSelection
                if obj.nbAtoms():
                    names.append(name)
                    children.append( Atom(obj) )
            else:
                names.append(name)
                children.append( Atom(atom) )
        return children, names

class ResidueSet(list):
    pass
 
class Atom(TreeObject):

    def __init__(self, patom):
        if isinstance(patom, SelectionSet):
            self._atomSubset = ProDyAtom(patom[0].getAtomGroup(),
                                         patom[0].getIndices()[0],
                                         patom[0].getACSIndex())
            self.name = patom[0].getNames()[0]
            self.alias = ''
            self.selSet = SelectionSet( [self.select()] )
        elif isinstance(patom, ProDyAtom):
            self._atomSubset = patom
            self.name = patom.getName()
            self.selSet = SelectionSet( [self.select()] )
        else:
            raise ValueError("Atom needs a SelectionSet or prody.atom.Atom instance")

    def select(self):
        ind = self._atomSubset.getIndex()
        return Selection(self._atomSubset.getAtomGroup(),
                         [ind], 'index (%d)'%ind,
                         self._atomSubset.getACSIndex())
                         
    def getAtomGroup(self):
        return self.mol.getAtomGroup()
                         
    def nbChildren(self):
        return 0

    def getChildrenAndNames(self, rootSelection=None):
        return [], []

#MoleculeGroup = MoleculeSet
