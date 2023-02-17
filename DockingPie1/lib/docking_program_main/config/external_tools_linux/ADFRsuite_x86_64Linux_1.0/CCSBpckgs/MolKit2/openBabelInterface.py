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
# $Header: /mnt/raid/services/cvs/MolKit2/openBabelInterface.py,v 1.2.4.1 2017/07/26 22:03:40 annao Exp $
# 
# $Id: openBabelInterface.py,v 1.2.4.1 2017/07/26 22:03:40 annao Exp $
#
from .molecule import Molecule

import prody
from prody.atomic.atomgroup import AtomGroup
from prody.atomic import ATOMIC_FIELDS
import numpy as np
import openbabel as ob
from MolKit2.selection import Selection

def ProdyToOBMol(sel, title = "ProdyToOBMOL"):
    """sele is an isntance of MolKit.sleection.Selection"""
    #pmol = molecule._ag
    assert isinstance(sel, Selection), "ERROR: bad argument"

    obmol = ob.OBMol()
    obmol.BeginModify()
    eTable = ob.OBElementTable()
    atoms = []
    #hv = pmol.getHierView()
    hv = sel.getHierView()

    # bidirectional atom indices lookup between OB and ProDy mol
    prodyToObIndex = np.zeros( len(sel.getAtomGroup()), 'int' )
    obToProdyIndex = np.zeros( len(sel), 'int' )
    nbOBAtoms = 0
    for chain in hv.iterChains():
        cId = chain.getChid()
        for res in chain.iterResidues():
            obres = obmol.NewResidue()
            obres.SetName(res.getResname() )
            obres.SetNum( int(res.getResnum()) )
            obres.SetChain(cId[0])
            #print "RESNAME", obres.GetName(), res.getResname()
            #print "RESNUM", obres.GetNum(), res.getResnum()
            for atom in res.iterAtoms():
                obatom = obmol.NewAtom()
                #pind = ob.OBPairData()
                #pind.SetAttribute('prodyIndex')
                #pind.SetValue(atom.getIndex())
                #obatom.SetData(pind)
                #obatom.SetIdx(atom.getIndex()+1)
                obatom.SetAtomicNum( int(atom.getData('atomicNumber') ))
                prodyToObIndex[atom.getIndex()] = nbOBAtoms
                obToProdyIndex[nbOBAtoms] = atom.getIndex()
                nbOBAtoms += 1
                obres.AddAtom(obatom)
                obatom.SetResidue(obres)
                obres.SetAtomID(obatom, atom.getName())
                atoms.append(atom)
                #obatom.SetAtomicNum( int(atom.getData('atomicNumber') ))
                coords = atom.getCoords()
                vec = ob.vector3()
                vec.SetX(coords[0])
                vec.SetY(coords[1])
                vec.SetZ(coords[2])
                obatom.SetVector(vec)
                obres.SetHetAtom(obatom, bool(atom.getFlag('hetatm')))
                #obres.SetSerialNum(obatom, int(atom.getSerial()))
    
    ag = sel.getAtomGroup()
    if ag._bonds is not None:
        bonds = sel.getBonds()
        for bo in [1,2,3]:
            for i,j in bonds[bo]:
                if ag._bondOrder:
                    bondOrder  = ag._bondOrder['%d %d'%(i,j)]
                else:
                    bondOrder = 1
                #print i, j, int(prodyToObIndex[i]+1), int(prodyToObIndex[j]+1), bondOrder
                obmol.AddBond(int(prodyToObIndex[i]+1),
                              int(prodyToObIndex[j]+1), bondOrder)
            #print i, j,bondOrder
            #obmol.AddBond(i+1, j+1, bondOrder)
            #obmol.AddBond(int(i+1), int(j+1), bondOrder)

    ## if pmol._bonds is not None:
    ##     for n in range(len(pmol._bonds)):
    ##         i,j = pmol._bonds[n]
    ##         if pmol._bondOrder:
    ##             bondOrder  = pmol._bondOrder['%d %d'%(i,j)]
    ##         else:
    ##             bondOrder = 1
    ##         obmol.AddBond(int(i+1), int(j+1), bondOrder)

    #print 'FAGA', obmol.NumResidues(), obmol.NumBonds()

    resdat = ob.OBResidueData()
    resdat.AssignBonds(obmol,ob.OBBitVec())
    obmol.EndModify(True)
    # MS: adds a bond in chain A of 1jff
    #obmol.ConnectTheDots() # build bonds based on covalent radii
    obmol.PerceiveBondOrders()

    obmol.prodyToObIndex = prodyToObIndex
    obmol.obToProdyIndex = obToProdyIndex
    obmol.prodyAtomList = atoms

    #typer = ob.OBAtomTyper()
    #typer.Init()
    #typer.AssignHyb(obmol)

    #for res in ob.OBResidueIter(obmol):
    #    print res.GetName(), type(res.GetName()), len(res.GetName())

    #obmolNew = ob.OBMol(obmol)
    #obmolNew.prodyToObIndex = prodyToObIndex
    #obmolNew.obToProdyIndex = obToProdyIndex
    #import pybel
    #pybel.Molecule(obmol).write('pdb', 'debug.pdb', overwrite=1)
    #return obmolNew
    return obmol


def OBMolToProdyAG(obmol):
    # return Prody AtomGroup-Based MolKit.molecule.Moelcule object
    assert isinstance(obmol, ob.OBMol), "ERROR: bad argument"
    name = obmol.GetTitle()
    atomgroup = AtomGroup(name)

    asize = obmol.NumAtoms()
    prodyToObIndex = np.zeros( asize, 'int' )
    obToProdyIndex = np.zeros( asize, 'int' )

    # create numpy arrays for
    coordinates = np.zeros((asize, 3), dtype=float)
    atomnames = np.zeros(asize, dtype=ATOMIC_FIELDS['name'].dtype)
    atomtypes = np.zeros(asize, dtype=ATOMIC_FIELDS['type'].dtype)
    resnames = np.zeros(asize, dtype=ATOMIC_FIELDS['resname'].dtype)
    resnums = np.zeros(asize, dtype=ATOMIC_FIELDS['resnum'].dtype)
    chainids = np.zeros(asize, dtype=ATOMIC_FIELDS['chain'].dtype)
    hetero = np.zeros(asize, dtype=bool)
    serials = np.zeros(asize, dtype=ATOMIC_FIELDS['serial'].dtype)
    segnames = np.zeros(asize, dtype=ATOMIC_FIELDS['segment'].dtype)
    elements = np.zeros(asize, dtype=ATOMIC_FIELDS['element'].dtype)

    #termini = np.zeros(asize, dtype=bool)
    #altlocs = np.zeros(asize, dtype=ATOMIC_FIELDS['altloc'].dtype)

    charges = np.zeros(asize, dtype=ATOMIC_FIELDS['charge'].dtype)
    radii = np.zeros(asize, dtype=ATOMIC_FIELDS['radius'].dtype)

    #bfactors = np.zeros(asize, dtype=ATOMIC_FIELDS['beta'].dtype)
    #occupancies = np.zeros(asize, dtype=ATOMIC_FIELDS['occupancy'].dtype)

    #print obmol.NumAtoms()
    
    #maxAtypLen = 0
    elementTable = ob.OBElementTable()
    atnum = 0

    for atom in ob.OBMolAtomIter(obmol):
        res = atom.GetResidue()
        resName = res.GetName()
        resNum = res.GetNum()
        aIndex = atom.GetIdx()-1
        #print 'ATMO', atnum, atom.GetIdx(), resName, resNum, residue.GetAtomID(atom)
        prodyToObIndex[atnum] = aIndex
        obToProdyIndex[aIndex] = atnum

        coordinates[aIndex] = atom.GetX(),atom.GetY(),atom.GetZ()
        atomnames[aIndex] = res.GetAtomID(atom)
        resnames[aIndex] = resName[:3]
        resnums[aIndex] = resNum
        chainids[aIndex] = 'A'
        hetero[aIndex] = atom.IsHeteroatom()
        #termini[aIndex] = False
        serials[aIndex] = aIndex
        #segnames[aIndex] = 'noSegm'
        atype = atom.GetType()
        #if len(atype)> maxAtypLen:
        #    maxAtypLen = len(atype)
        atomtypes[aIndex] = atype
        #print atype # different from file :( "." is missing
        elements[aIndex] = elementTable.GetSymbol(atom.GetAtomicNum())

        charges[aIndex] = atom.GetPartialCharge()
        radii[aIndex] = 1.0

        #altlocs[aIndex] = ' '
        #bfactors[aIndex] = 0.0
        #occupancies[aIndex] = 0.0

        #atId = atom.GetIdx()
        #name = res.GetAtomID(atom)
        #coords = atom.GetCoordinate()
        #atomicNumber = atom.GetAtomicNum()
        #valence = atom.GetValence()
        #hyb = atom.GetHyb()
        #atype = atom.GetType()
        #pcharge = atom.GetPartialCharge()
        #print resName, resNum, name, atId, coords, atomicNumber, valence, hyb, atype, pcharge
        atnum += 1

    atomgroup._setCoords(coordinates)
    atomgroup.setNames(atomnames)
    atomgroup.setTypes(atomtypes)
    atomgroup.setData('atomType', atomtypes)
    atomgroup.setResnames(resnames)
    atomgroup.setResnums(resnums)
    atomgroup.setChids(chainids)
    atomgroup.setFlags('hetatm', hetero)
    #atomgroup.setFlags('pdbter', termini)
    #atomgroup.setAltlocs(altlocs)
    #atomgroup.setIcodes(np.char.strip(icodes))
    atomgroup.setSerials(serials)
    #atomgroup.setBetas(bfactors)
    #atomgroup.setOccupancies(occupancies)
    #atomgroup.setSegnames(np.char.strip(segnames))
    atomgroup.setElements(np.char.strip(elements))
    atomgroup.setCharges(charges)
    atomgroup.setRadii(radii)

    # add bonds and bondOrders
    pairs = []
    bo = []
    #print 'AAAA', obmol.NumAtoms(), obmol.NumBonds()
    #print 'FAGA', obmol.NumAtoms(), obmol.NumBonds(), id(obmol)
    for bond in ob.OBMolBondIter(obmol):
        #i = bond.GetBeginAtomIdx()-1
        #j = bond.GetEndAtomIdx()-1
        i = bond.GetBeginAtomIdx()-1
        j = bond.GetEndAtomIdx()-1
        #print 'BO', i, j, bond.GetBondOrder()
        if i>j:
            a=i
            i=j
            j=a
        pairs.append( (i,j) )
        if bond.IsAromatic():
            bondOrder = 4
        elif bond.IsAmide():
            bondOrder = 5
        else:
            bondOrder = bond.GetBondOrder()
        bo.append(bondOrder)
        #print 'BOND', bond.GetBeginAtomIdx()-1, bond.GetEndAtomIdx()-1, bondOrder
        #atom1 = obmol.GetAtom(bond.GetBeginAtomIdx())
        #atomname1 = atom1.GetResidue().GetAtomID(atom1)
        #atom2 = obmol.GetAtom(bond.GetEndAtomIdx())
        #atomname2 = atom2.GetResidue().GetAtomID(atom2)
        #print atomname1, atomname2, bond.GetBondOrder()
    #print 'FAGA END'
    #print pairs
    atomgroup.setBonds(pairs, bo)
    obmol.prodyToObIndex = prodyToObIndex
    obmol.obToProdyIndex = obToProdyIndex
    obmol.prodyAtomList = atomgroup
    return atomgroup

def OBMolToPrody(obmol):
    return Molecule(obmol.GetTitle(), OBMolToProdyAG(obmol))
