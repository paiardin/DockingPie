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
# $Header: /mnt/raid/services/cvs/MolKit2/__init__.py,v 1.2.2.2 2017/11/10 18:49:08 sanner Exp $
# 
# $Id: __init__.py,v 1.2.2.2 2017/11/10 18:49:08 sanner Exp $
#
from __future__ import absolute_import
import openbabel
import mmtf as MMTF
from .babelElements import babel_elements as elements

for k,v in elements.items():
   if len(k)==2:
       elements[k.upper()] = v

import os
import prody
# put error as default level in ProDyDist/prody/__init__.py
prody.confProDy(verbosity='error', auto_secondary=True, silent=True)
from .molecule import Molecule, MultiMolecule

_KNOW_FORMATS = ["pdb", "pdbqt", "ent", "pdb.gz",
                 "mmtf", "mol2", "sdf"]#, "pqr", "cif"]

def Read(filename, fileFormat='auto', **kw):
   if fileFormat == 'auto':
      fileFormat = os.path.splitext(filename)[1][1:]

   indexingMode = kw.get('indexingMode', 'foreground')
   assert fileFormat in _KNOW_FORMATS, 'ERROR: bad fileFormat parameter, got ".%s" expected one of %s'%(
      fileFormat,_KNOW_FORMATS)
   if fileFormat in ["pdb", "pdbqt", "ent", "pdb.gz"]:
      header = kw.get('header', False)
      model = kw.get('model',None)
      mol = readWithPrody(filename, header=header, model=model)
   elif fileFormat == "mol2":
      mol = readWithOB(filename, 'mol2', indexingMode=indexingMode)
   elif fileFormat == "sdf":
      mol = readWithOB(filename, 'sdf', indexingMode=indexingMode)
   elif fileFormat == "mmtf":
      decoder = kw.get('decoder', None)
      mol = readMMTF(filename, decoder)
   ## elif fileExt == "pqr":
   ##    mols = readPQR(filename)
   ## elif fileExt == "cif":
   ##    mols = readMMCIF(filename)
   return mol

def readMMTF(filename, decoder=None):
    if decoder is None:
        mmtfDecoder = MMTF.parse(filename)
    else:
        mmtfDecoder = decoder
    from .mglmmtf import MMTFtoPrody
    name = os.path.splitext(os.path.basename(filename))[0]
    atGroup = MMTFtoPrody(mmtfDecoder, name=name)
    mol = Molecule(name, atGroup, filename)
    mol.name = name
    return mol

def readWithPrody(filename, header=False, model=None):
    #prody.confProDy(verbosity='error')
    if not os.path.exists(filename):
        raise AssertionError , "%s doesn't exist" %filename
    ext = os.path.splitext(filename)[1]
    name = os.path.splitext(os.path.split(filename)[1])[0]
    if ext.lower()=='.pdb':
        ag = prody.parsePDB(filename, model=model)
        mol = Molecule(name, ag, filename=filename)
        mol.buildBondsByDistance()
        mol.defaultRadii()
        if header is True:
            try:
                mol.pdbHeader = prody.parsePDB(filename, model=0, header=True)
            except ValueError:
                mol.pdbHeader = {}
    elif ext.lower()=='.pdbqt':
        ag = prody.parsePDB(filename, format='PDBQT', model=model)
        mol = Molecule(name, ag, filename=filename)
        mol.buildBondsByDistance()
        mol.defaultRadii()
        if header is True:
            try:
               mol.pdbHeader = prody.parsePDB(filename, model=0, header=True)
            except ValueError:
               mol.pdbHeader = {}
        findHbAcceptorsDonors(mol)
    return mol

def readWithOB(filename, _format, group=None, indexingMode='foreground'):
   assert _format in ["mol2", "sdf"], "ERROR: bad file format, expected mol2 or sdf got %s"%_format
   obconv = openbabel.OBConversion()
   obconv.SetInFormat(_format)
   obmol = openbabel.OBMol() 
   thereIsMore = obconv.ReadFile(obmol, filename)
   # check if there are more molecules in the file
   obmolDum = openbabel.OBMol()
   thereIsMore = obconv.Read(obmolDum)
   if thereIsMore: # MultiMolecule
      mol = MultiMolecule(filename, indexingMode=indexingMode)
      mol._multi = 'molecules'

   else: # single molecule in Mol2 file
      from MolKit2.openBabelInterface import OBMolToPrody
      mol = OBMolToPrody(obmol)
      mol._multi = 'False'
      mol.name = obmol.GetTitle()
      mol._basename = os.path.splitext(os.path.basename(filename))[0]
   return mol

## def readPQR(filename, group=None):
##    from MolKit2.pdbParser import PQRParser
##    newparser = PQRParser(filename, modelsAs=modelsAs)
##    mols = newparser.parse()
##    if mols is None :
##       del newparser
##       return
##    newmol = []

##    for m in mols:
##       mol = self.addMolecule(m, group=group, filename=filename)
##       if mol is None:
##          del newparser
##          return mols.__class__([])
##       newmol.append(mol)
##    return mols.__class__(newmol)


## def readMMCIF(filename, group=None):
##    from MolKit2.mmcifParser import MMCIFParser
##    newparser = MMCIFParser(filename)
##    mols = newparser.parse()
##    if mols is None: return
##    newmol = []
##    for m in mols:
##       mol = self.addMolecule(m, group=group, filename=filename)
##       if mol is None:
##          del newparser
##          return mols.__class__([])
##       newmol.append(mol)
##    return mols.__class__(newmol)

from MolKit2.PDBresidueNames import AAnames
def getSequence(atoms):
   """return sequence string with 1 character amino acids names"""
   try:
       return "".join([AAnames[x] for x in atoms.select('ca').getResnames()]) 
   except KeyError:
       return None

## HBOND functions for PDBQT files
from .selection import Selection
import numpy
def findHbAcceptorsDonors(mol):
   
    acceptor_types = {}.fromkeys(['OA', 'NA', 'SA'])
    donor_types = {}.fromkeys(['N', 'O', 'OA', 'NA', 'SA'])
    donors = [False]*len(mol._ag)
    acceptors = [False]*len(mol._ag)
    for n, atom in enumerate(mol._ag):
        atype = atom.getType().strip()
        index = atom.getIndex()
        if acceptor_types.has_key(atype):
            acceptors[n] = True
        if donor_types.has_key(atype):
           hasH = False
           for neighbor in atom.iterBonded():
               if neighbor.getType()=='HD':
                   donors[n] = True

    mol._ag.setFlags('donor', donors)
    mol._ag.setFlags('acceptor', acceptors)

def hbGeom(atoms1, atoms2, cut=3.21, foff=0):
    cut2 = cut*cut
    hbCoords = []
    faces = []
    n = 0
    if len(atoms2) > len(atoms1):
       tmp = atoms1
       atoms1 = atoms2
       atoms2 = tmp
    ag1 = atoms1.getAtomGroup()
    ag2 = atoms2.getAtomGroup()
    d1 = atoms1.select('donor')
    if d1:
        #rcoords = d1.getCoords() this give the coord set fro when atom was created
        rcoords = ag1.getCoords()[d1.getIndices()]
        a2 = atoms2.select('acceptor')
        if a2:
            for lc in ag2.getCoords()[a2.getIndices()]:
                #lc = a.getCoords()
                delta = lc-rcoords
                dist2 = numpy.sum(delta*delta, 1)
                for i in numpy.where(dist2 < cut2)[0]:
                    hbCoords.append( lc )
                    hbCoords.append( rcoords[i] )
                    faces.append( (n+foff, n+1+foff) )
                    n += 2

    a1 = atoms1.select('acceptor')
    if a1:
        #rcoords = a1.getCoords()
        rcoords = ag1.getCoords()[a1.getIndices()]
        d2 = atoms2.select('donor')
        if d2:
            #for a in d2:
                #lc = a.getCoords()
            for lc in ag2.getCoords()[d2.getIndices()]:
                delta = lc-rcoords
                dist2 = numpy.sum(delta*delta, 1)
                for i in numpy.where(dist2 < cut2)[0]:
                    hbCoords.append( lc )
                    hbCoords.append( rcoords[i] )
                    faces.append( (n+foff, n+1+foff) )
                    n += 2
    return hbCoords, faces
  
