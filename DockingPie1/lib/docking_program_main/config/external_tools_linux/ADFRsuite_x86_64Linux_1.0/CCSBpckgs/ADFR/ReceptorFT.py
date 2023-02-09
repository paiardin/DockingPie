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
# Copyright: M. Sanner and TSRI 2015
#
#########################################################################
#
# $Header: /mnt/raid/services/cvs/ADFR/ReceptorFT.py,v 1.10 2017/04/15 02:04:38 sanner Exp $
#
# $Id: ReceptorFT.py,v 1.10 2017/04/15 02:04:38 sanner Exp $
#

from ADFR.AARotamers import RotamerLib

from ADFRcc.adfr import AtomSet, AtomSetStatic, FTDeclareAtomSet, FTAtomsNode

from mglutil.math.torsion import torsion
from ADFR.LigandFT import MotionBase, PyFTTorsion
from random import random, uniform, gauss
from MolKit2.selection import Selection

class FTSoftRotamer(MotionBase):

    def __init__(self, name, angles, deviations):
        self._deltaAmplitude = []
        self.deltaAmplitude = []
        self.name = name
        self.angles = angles
        self.deviations = deviations
        MotionBase.__init__(self, name)
        self.torsions = [] # will hold PyFTTorsion objects
        self.nodes = []
        
    def __repr__(self):
        return '%s %s'%(self.__class__, self.name)
    
    def getNumVariables(self):
        return len(self.angles[0])

    def getDeltaAmplitudeBase(self):
        deltas = []
        for i, tors in enumerate(self.torsions):
            deltas.append(tors.getDeltaAmplitudeBase())
        return deltas
        
    def setDeltaAmplitude(self, deltas):
        for i, tors in enumerate(self.torsions):
            tors.setDeltaAmplitude(deltas[i])
        self.deltaAmplitude = deltas[:]
        self._deltaAmplitude = deltas[:]

    def resetDeltaAmplitude(self):
        for tors in self.torsions:
            tors.resetDeltaAmplitude()
        self.deltaAmplitude = self._deltaAmplitude

    def getGenesPy(self):
        values = []
        for tors in self.torsions:
            values.append(tors.getGenesPy()[0])
        return values
    
    def getIdentityGenesPy(self):
        values = []
        for tors in self.torsions:
            values.append(tors.getIdentityGenesPy()[0])
        return values
    
    def getVariablesPy(self):
        values = []
        for tors in self.torsions:
            values.append(tors.getVariablesPy()[0])
        return values
    
    ## def randomValues(self):
    ##     values = []
    ##     for tors in self.torsions:
    ##         tors.randomize()
    ##         #print 'R', m.__class__.__name__, m.getVariablesPy()
    ##         values.extend( tors.getGenes() )
    ##     return values
    
    def initialize(self):
        for tors in self.torsions:
            tors.setGenes(tors.getIdentityGenesPy())

    ##
    ## for Solis Wets
    ##
    ## def genDeviate(self, searchRate):
    ##     """
    ##     return a list of deltas for each motion
    ##     """
    ##     deltas = []
    ##     ct = 0
    ##     for tors in self.torsions:
    ##         ctm, devs = tors.genDeviate(searchRate)
    ##         ct += ctm
    ##         deltas.extend( devs )
    ##     return ct, deltas
        
    ## def initialBias(self):
    ##     bias = []
    ##     for tors in self.torsions:
    ##         bias.extend( tors.initialBias())
    ##     return bias
        
    ## def biasedDevs(self, dev, bias):
    ##     """
    ##     let each motion object apply the bias to the deviation
    ##     """
    ##     off = 0
    ##     biasedDev = []
    ##     for tors in self.torsions:
    ##         biasedDev.extend( tors.biasedDevs(dev[off:off+1], bias[off:off+1]))
    ##         off += 1
    ##     return biasedDev
        
    ## def applyDelta(self, values, delta, direction='+'):
    ##     """
    ##     return a list of deltas for each motion
    ##     """
    ##     off = 0
    ##     nv = []
    ##     for tors in self.torsions:
    ##         nv.extend(tors.applyDelta(values[off:off+1], delta[off:off+1],
    ##                                direction))
    ##         off += 1
    ##     return nv

    ## def scaleBias(self, bias, factor):
    ##     off = 0
    ##     newBias = []
    ##     for tors in self.torsions:
    ##         newBias.extend(tors.scaleBias(bias[off:off+1], factor))
    ##         off += 1
    ##     return newBias
        
    ## def resetDeltaAmplitude(self):
    ##     for tors in self.torsions:
    ##         tors.resetDeltaAmplitude()
        
    ## def newBias(self, coef1, bias, coef2, dev):
    ##     off = 0
    ##     newBias = []
    ##     for tors in self.torsions:
    ##         newBias.extend(tors.newBias(coef1, bias[off:off+1],
    ##                                  coef2, dev[off:off+1]))
    ##         off += 1
    ##     return newBias

    ## def scaleUpAmplitude(self, factor):
    ##     import pdb; pdb.set_trace()
    ##     for tors in self.torsions:
    ##         tors.scaleUpAmplitude(factor)
            
    ## def scaleDownAmplitude(self, factor):
    ##     import pdb; pdb.set_trace()
    ##     for tors in self.torsions:
    ##         terminate = tors.scaleDownAmplitude(factor)
    ##         if terminate is True:
    ##             return True
    ##     return False

    ## def getVariablesPy(self):
    ##     import pdb; pdb.set_trace()
    ##     values = []
    ##     for tors in self.torsions:
    ##         values.extend(tors.getVariablesPy())
    ##     return values

    ## def setValuesFromGenes(self, genes):
    ##     import pdb; pdb.set_trace()
    ##     for gene, tors in zip(genes, self.torsions):
    ##         tors.setVariables([tors.mini[0] + gene*tors.length[0]])
    ##         tors.setGenes([gene])
    
    ## def getGenesForValues(self, values):
    ##     import pdb; pdb.set_trace()
    ##     genes = []
    ##     for val, tors in zip(values, self.torsions):
    ##         genes.append( tors.getGenesForValuesPy([val])[0] )
    ##     return genes

    ## def setVariables(self, angles):
    ##     for a, m in zip(angles, self.torsions):
    ##         tors.setVariables([a])

    ## def initialize(self):
    ##     for tors in self.torsions:
    ##         tors.initialize()
    ##         oAngle = tors.getOrigAngle()
    ##         if oAngle < 0:
    ##             oAngle += 360.            
    ##         tors.setGenes([oAngle/360.])
                
    def chi1Rotamers(self, originalChi1Angle):
        # find rotamers with Xtal Chi 1
        self.XtalChi1Rotamers = []
        xtalChi1 = originalChi1Angle
        if originalChi1Angle < 0.: originalChi1Angle += 360.

        off = 0.0
        while len(self.XtalChi1Rotamers) == 0:
            for i, angles in enumerate(self.angles):
                rotChi1 = angles[0]
                if rotChi1 < 0.: rotChi1 += 360.
                #if self.atoms[0].parent.name == 'VAL32':
                #    import pdb
                #    pdb.set_trace()
                #if xtalChi1 < 0.: xtalChi1 += 360.
                if abs(rotChi1-originalChi1Angle) < 50.+off:
                    self.XtalChi1Rotamers.append(i)
                self.confNBXtalChi1 = len(self.XtalChi1Rotamers)
            off += 1
        #print 'chi1Rotamers', self.name, originalChi1Angle, self.XtalChi1Rotamers
        self.confNBXtalChi = len(self.XtalChi1Rotamers)
        
    def randomize(self):
        #import pdb; pdb.set_trace()
        rvalue = random()
        if rvalue < 0:#0.33:# 0.66: # restore Xtal conformation for this side chain
            for i, gene in enumerate(genes):
                #angle = gauss(self.rotamer.originalAngles[i], 0.2)
                angle = self.rotamer.originalAngles[i]
                #if angle < 0: angle += 360.
                #elif angle > 360: angle -= 360.
                if angle != gene._value:
                    gene._value = angle
                    #print 'ROTARAND SETTING XTAL ROT for', self.rotamer.atoms[0].parent.name
                else:
                    #print 'ROTARAND WAS ALREADY XTAL'
                    return -1
            return 0
        else:
            if rvalue < 0:#0.5:#0.85: # use rotamer hi1Prob:
                n = int(uniform(0., 0.999999)*self.rotamer.confNBXtalChi1)
                index = self.rotamer.XtalChi1Rotamers[n]
                #print 'ROTARAND KEEPING XTAL CHI1 ROT for', self.rotamer.atoms[0].parent.name
                returnValue = 1
            else:
                index = int(uniform(0., 0.999999)*len(self.angles))
                returnValue = 2
            dev = self.deviations[index]
            angles = self.angles[index]
            for i, tors in enumerate(self.torsions):
                #gene._value = 0.0
                angle = gauss(angles[i], dev[i])
                if angle < 0: angle += 360.
                elif angle > 360: angle -= 360.
                #tors.setVariables([angle])
                tors.setGenes([angle/360.])
            #print 'ROTARAND CHANGING ROTAMER for', self.getName(), index, self.getGenesPy()
            return returnValue

    def mutate(self, genes, mutationRate, dev=0.2):
        mutated = 0
        status = -1
        #if random() < mutation_rate:#*len(genes):
        status = self.randomize()
        if status > -1:
            mutated += len(genes)

        return mutated, status

from ADFR import getDonorTypes

class MkReceptorFT:

    def __init__(self):
        self.rotlib = RotamerLib()
        self.motions = [] # will hold a list of FTSoftRotamer. One per felxible side chain
        self.allMotions = [] # will hold a list of all FTTorsion objects across all FTSoftRotamer
        self.ft = None    # will be a FT for all receptor FSCs
        
    ## def getDonorTypes(self, atoms):
    ##     """
    ##     for all atoms with HD neighbors we need to assign a new type

    ##     N-HD becomes NX-HD
    ##     NA-HD becomes N2-HD
    ##     OA-HD becomes OX-HD
    ##     SA-HD becomes SX-HD
    ##     """
    ##     def hasHD(atom, atypes):
    ##         for neighbor in atom.iterBonded():
    ##             if atypes[neighbor.getIndex()]=='HD':
    ##                 return True
    ##         return False

    ##     mol = atoms.getAtomGroup().getMolecule()
    ##     adTypes = mol._ag.getData("AD_element")
    ##     subsetADTypes = []
    ##     for a in atoms:
    ##         i = a.getIndex()
    ##         # skip N from backbone else it becomes NX
    ##         # and the AtomSet complains that NX has no H
    ##         # attached since HN is not in the FRatoms
    ##         if adTypes[i]=='N' and a.getName()!='N':
    ##             if hasHD(a, adTypes):
    ##                 subsetADTypes.append('NX' )
    ##             else:
    ##                 subsetADTypes.append(adTypes[i])
    ##         elif adTypes[i]=='NA':
    ##             if hasHD(a, adTypes):
    ##                 subsetADTypes.append('N2' )
    ##             else:
    ##                 subsetADTypes.append(adTypes[i])
    ##         elif adTypes[i]=='OA':
    ##             if hasHD(a, adTypes):
    ##                 subsetADTypes.append('OX')
    ##             else:
    ##                 subsetADTypes.append(adTypes[i])
    ##         elif adTypes[i]=='SA':
    ##             if hasHD(a, adTypes):
    ##                 subsetADTypes.append('SX')
    ##             else:
    ##                 subsetADTypes.append(adTypes[i])
    ##         else:
    ##             subsetADTypes.append(adTypes[i])
    ##     return subsetADTypes

    def mkFT(self, mol, flexRes):
        """
        build FT for receptor with flexible side chains
        flexres is a list of residue names and numbers eg.
            [ ['A', ('GLU', '167'), ('TYR', '628'), ('ARG', '87')],
              ['B', ('GLU', '167')] ]
        """
        #build an atom set for flexible receptor atoms
        from ADFR.utils.MakeGrids import splitFlexRes
        receptorAtoms, sideChainAtoms = splitFlexRes(mol, flexRes, exclude='C O')
        # build a look up allowing to get a relative index in FRatoms from 
        # absolute atom index in receptor molecule. Also build list of atoms 
        # not scorable ingrid scorer and pairwise scorers
        indexLookUp = {}
        notGridScorable = []
        notPairwiseScorable = []

        # sort sideChainAtoms
        indices = sideChainAtoms.getIndices()
        indices.sort()
        sideChainAtoms._indices[:] = indices

        print 'Flexible Receptor Atoms:'
        for i, a in enumerate(sideChainAtoms):
            print "%c:%s%d:%s"%(a.getChid(),a.getResname(), a.getResnum(),
                                a.getName())
            
        #import pdb; pdb.set_trace()
        for i, a in enumerate(sideChainAtoms):
            indexLookUp[a.getIndex()] = i
            name = a.getName()
            if name in ['N', 'CA', 'CB']:
                print 'NOT GRID SCORABLE', a, a.getName(), a.getResname(),
                a.getResnum(), a.getChid()
                notGridScorable.append(i)
                if name == 'N':
                    print 'NOT PAIREWISE SCORABLE', a, a.getName(), \
                          a.getResname(), a.getResnum(), a.getChid()
                    notPairwiseScorable.append(i)

        atomSetStatic = AtomSetStatic(len(sideChainAtoms),  'FRatoms')
        for n in notGridScorable:
            atomSetStatic.setScorableGrid(n, False)
        for n in notPairwiseScorable:
            atomSetStatic.setScorablePairwise(n, False)
            
        adtypes = getDonorTypes(sideChainAtoms)
        self._adtypes = adtypes
	print 'RECEPTOR ATOM TYPES', adtypes
        atomSetStatic.setAtomTypes(adtypes)
        atomSetStatic.setCharges(sideChainAtoms.getCharges())
        atomSetStatic.setOrigCoords(sideChainAtoms.getCoords())
        bonds = sideChainAtoms.getBonds()
        bl = []
        for bds in bonds[1:]:
            for bond in bds:
                bl.extend((indexLookUp[bond[0]], indexLookUp[bond[1]]))
        bl.append(-1)
        atomSetStatic.setCovalentBonds(bl)

        atomSet = AtomSet(atomSetStatic)

        # remove non bonded pairs within rigid body
        #from ADFR.torTreeFromPDBQT import weedBonds
        # FIXME .. bond weeding for side chains ? TRP
        #pairs = weedBonds(mol, torTree)

        #for i,j in pairs:
        #    self.atomSetStatic.setPairScorable(i, j, False)

        # create FT root
        ftRRoot = FTDeclareAtomSet()
        ftRRoot.setAtomSet(atomSet)

        adTypes = mol._ag.getData("AD_element")
        for chid, residues in flexRes:
            for resname, resnum in residues:
                print resname, resnum
                angleDef, angleList, angleDev = self.rotlib.get(resname)

                # create dict of rotatble bonds used to add H atoms
                rotBonds = {}
                for adef in angleDef:
                    rotBonds['%s-%s'%(adef[0][1], adef[0][2])] = True

                # create FTSoftRotamer motion object to be placed in Genome
                softRot = FTSoftRotamer('%s_%s'%(resname,resnum), angleList, angleDev)
                self.motions.append(softRot)
                ftParent = ftRRoot
                for i, adef in enumerate(angleDef):

                    if chid == ' ':
                        tatoms = mol.select('chid "%s" resname %s resnum %s name %s %s %s %s'%(
                            (chid, resname, resnum)+tuple(adef[0])))
                    else:
                       tatoms = mol.select('chid %s resname %s resnum %s name %s %s %s %s'%(
                            (chid, resname, resnum)+tuple(adef[0])))
                       if len(tatoms)<4:
                           msg = 'ERROR: Chi angle defining atoms not found in receptor\n'
                           msg += '    looking for %s:%s%d: %s %s %s %s but only got'%(
                               chid, resname, resnum,
                               adef[0][0], adef[0][1], adef[0][2], adef[0][3])
                           for name in tatoms.getNames():
                               msg += ' '+name
                           raise ValueError, msg
                    a, b, c, d = [indexLookUp[x] for x in tatoms.getIndices()]
                    na, nb, nc, nd = adef[0]

                    ftTorsion = PyFTTorsion('Torsion [%s]-%s--%s-[%s]'%(na,nb,nc,nd))
                    ftTorsion.setDihedralAtomIndexes(a,b,c,d, atomSet)
                    c1, c2, c3, c4 = tatoms.getCoords()  
                    #ftTorsion.originalValue = torsion(c1, c2, c3, c4)
                    print 'TOR_rec', na, nb, nc, nd, ftTorsion, ftTorsion.getOrigAngle(), ftTorsion.getIdentityGenesPy()
                    if i==0:
                        softRot.chi1Rotamers(ftTorsion.getOrigAngle())
                    ftParent.addChild(ftTorsion)
                    softRot.torsions.append(ftTorsion)

                    child = FTAtomsNode()
                    softRot.nodes.append(child)
                    ftTorsion.addChild(child)
                    if chid==' ':
                        selStr = 'chid "%s" resname %s resnum %s name '%(chid, resname, resnum)
                    else:
                        selStr = 'chid %s resname %s resnum %s name '%(chid, resname, resnum)
                        
                    for atName in adef[1]:
                        selStr += '%s '%atName
                    matoms = mol.select(selStr)
                    #print '  moving', matoms.getNames(),
                    atIndices = [indexLookUp[x] for x in matoms.getIndices()]
                    ## add hydrogen atoms attached the second atom in rotatable bond
                    ## and attached to moving atoms
                    # H bound second atom in rotatble bond always move with that torsion
                    atom3 = mol.select('resname %s resnum %s name %s'%(resname, resnum, nc))
                    for atom in atom3:
                        for neighbor in atom.iterBonded():
                            if adTypes[neighbor.getIndex()]=='HD':
                                atIndices.append( indexLookUp[neighbor.getIndex()] )
                                #print neighbor.getName(),

                    # H atoms attached to moving atoms all move if bond c-d in a-b-c-d torsion
                    # is NOT rotatable
                    if not rotBonds.has_key('%s-%s'%(nc, nd)):
                        for a in matoms:
                            for neighbor in a.iterBonded():
                                if adTypes[neighbor.getIndex()]=='HD':
                                    atIndices.append( indexLookUp[neighbor.getIndex()] )
                                    #print neighbor.getName(),
                    #print
                    #print 'for', atIndices,
                    #print 'relative to', [torTree.serialToIndex[x] for x in ttchild.parent().atoms]
                    child.setAtomIndexes(atIndices)
                    ftParent = child
                softRot.setDeltaAmplitude([0.01]*len(angleDef))
                #print 'FUGU', softRot.identityGenes()
                self.allMotions.extend(softRot.torsions)
                
        #ftRRoot.update() # force evaluation of torsion so that
                        # ftTorsion.getOrigAngle() returns the value
        self.ftRoot = ftRRoot
        self.indexLookUp = indexLookUp # key: prody full receptor index, value: self.atoms index
        self.atoms = sideChainAtoms # prody atom set for FR
        self.atomSet = atomSet # ADFRcc atomset for FR
        
if __name__=='__main__':
    from MolKit2 import Read
    rec = Read('Astex/receptorsuniq-85/1n1m_rec.pdbqt')
    RFTGen = MkReceptorFT()
    flexRes = [('GLU', '167'), ('TYR', '628'), ('ARG', '87')]
    RFTGen.mkFT(rec, flexRes)

    #for m in RFTGen.motions:
        #for tors in m.torsions:
        #    tors.idGenes = tors.getGenesForValues([tors.getOrigAngle()])

    coords = RFTGen.atomSet.getCoordsPy()

    #for m in RFTGen.motions:
    #    for tors in m.torsions:
    #        tors.setVariables([180.])

    for m in RFTGen.motions:
        m.randomize()

    RFTGen.ft.update()
    RFTGen.atoms.setCoords( RFTGen.atomSet.getCoordsPy() )
    import prody
    prody.writePDB('1n1m_frec.pdb', RFTGen.atoms)
