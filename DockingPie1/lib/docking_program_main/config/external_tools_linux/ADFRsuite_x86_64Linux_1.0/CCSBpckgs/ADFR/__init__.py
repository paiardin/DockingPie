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
# $Header: /mnt/raid/services/cvs/ADFR/__init__.py,v 1.32.2.1 2017/08/01 00:35:32 annao Exp $
#
# $Id: __init__.py,v 1.32.2.1 2017/08/01 00:35:32 annao Exp $
#

import os, sys
from math import sqrt

_intelContest_ = False
__version__ = 1
__subversion__ = 0
__revision__ = ""

class PDBQTWriter:
    """
    Class for writing a PDBQT file from a torsion tree as obtained from
    ADFR.torTreeFromPDBQT
    """
    # we use columns 68:69 for 4 letter segment 
    PDBQT_FMT = "%-6s%5d %4s%c%3s %c%4d%c   %8.3f%8.3f%8.3f%6.2f%6.2f%4s%5.2f  %-2s"
    def getLines(self, root):
        ag = root.mol._ag
        self.names = ag.getNames()
        self.resnames = ag.getResnames()
        self.resnums = ag.getResnums()
        self.chids = ag.getChids()
        self.altoc = ag.getAltlocs()
        self.icodes = ag.getIcodes()
        self.coords = ag.getCoords()
        self.occup = ag.getOccupancies()
        self.tfact = ag.getBetas()
        self.charges = ag.getCharges()
        self.segnames = ag.getSegnames()
        self.elems = ag.getData('AD_element')
        self.lines = ['ROOT']
        self.getNodeLines(root, root.atoms)
        self.lines.append('ENDROOT')
        #import pdb; pdb.set_trace()
        for child in root.children:
            self._processChild(root, child)
        self.lines.append('TORSDOF %d'%root.torsdof)
        return self.lines

    def _processChild(self, root, node):
        self.lines.append('BRANCH %d %d'%tuple(node.bond))
        self.getNodeLines(root, node.atoms)
        for child in node.children:
            self._processChild(root, child)
        self.lines.append('ENDBRANCH %d %d'%tuple(node.bond))
    
    def getNodeLines(self, root, atIndices):
        for ind in atIndices:
            serial = root.serialToIndex[ind]
            name = self.names[serial]
            if len(name)==4 and name[-1].isdigit():
                name = name[3]+name[:3]
            elif len(self.elems[serial])==2:
                name = '%-4s'%name
            else:
                name = ' %-3s'%name
            icode = self.icodes[serial]
            if icode=='':
                icode = ' '
            self.lines.append(
                self.PDBQT_FMT%(
                    'ATOM', ind, name, self.altoc[serial],
                    self.resnames[serial], self.chids[serial],
                    self.resnums[serial], icode, self.coords[serial][0],
                    self.coords[serial][1], self.coords[serial][2],
                    self.occup[serial], self.tfact[serial],
                    self.segnames[serial], self.charges[serial], self.elems[serial]))

def getDonorTypes(atoms):
    """
    for all atoms with HD neighbors we need to assign a new type

    N-HD becomes NX-HD
    NA-HD becomes N2-HD
    OA-HD becomes OX-HD
    SA-HD becomes SX-HD
    """
    def hasHD(atom, atypes):
        for neighbor in atom.iterBonded():
            if atypes[neighbor.getIndex()]=='HD':
                return True
        return False

    mol = atoms.getAtomGroup().getMolecule()
    adTypes = mol._ag.getData("AD_element")
    subsetADTypes = []
    for a in atoms:
        i = a.getIndex()
        # skip N from backbone else it becomes NX
        # and the AtomSet complains that NX has no H
        # attached since HN is not in the FRatoms
        if adTypes[i]=='N' and a.getName()!='N':
            if hasHD(a, adTypes):
                subsetADTypes.append('NX' )
            else:
                subsetADTypes.append(adTypes[i])
        elif adTypes[i]=='NA':
            if hasHD(a, adTypes):
                subsetADTypes.append('N2' )
            else:
                subsetADTypes.append(adTypes[i])
        elif adTypes[i]=='OA':
            if hasHD(a, adTypes):
                subsetADTypes.append('OX')
            else:
                subsetADTypes.append(adTypes[i])
        elif adTypes[i]=='SA':
            if hasHD(a, adTypes):
                subsetADTypes.append('SX')
            else:
                subsetADTypes.append(adTypes[i])
        else:
            subsetADTypes.append(adTypes[i])
    return subsetADTypes

def getTORSDOF(filename):
    if os.path.splitext(filename)[1]!='.pdbqt':
        return -1
    f = open(filename)
    lines = f.readlines()
    f.close()
    return int(lines[-1].split()[1])
    
def checkLigandFile(filename):
    # check if there is a torsion tree in the file by looking a ROOT record
    if os.path.splitext(filename)[1]!='.pdbqt':
        return False
    f = open(filename)
    lines = f.readlines()
    f.close()
    for line in lines:
        if line.startswith('TORSDOF'):
            return 1
    return False

def getLigandFromFile(filename):
    if not checkLigandFile(filename):
        return None, 'ERROR: %s does not seem to contain a ligand prepared for ADFR'%filename
    from MolKit2 import Read
    mol = Read(filename)
    return mol, None
    
def ligandFTfromFile(filename, mol=None, covalentLigandAnchorAtoms=None,
                     fixedRoot=False, neighborSearchCutoff=-1.):
    """Build a flexibility tree for a prody molecule"""

    from ADFR.torTreeFromPDBQT import TorTreeFromPDBQT     
    mkFt = TorTreeFromPDBQT()
    root = mkFt(filename, mol=mol)

    from ADFR.LigandFT import LigandFT
    return LigandFT(root.mol, root,
                    covalentLigandAnchorAtoms=covalentLigandAnchorAtoms,
                    fixedRoot=fixedRoot, neighborSearchCutoff=neighborSearchCutoff)


import numpy, os, prody
from glob import glob
from MolKit2 import Read
from time import time

class ADFR:

    def __init__(self, ligand, mapsFolder, receptor=None,
                 flexibleResidues=None, mapFilesRoot=None,
                 tpointsFilename=None, seedValue=-1, logFile=None,
                 covalentIndices=None, FTRecsrc=None, fixedRoot=False,
                 neighborSearchCutoff=-1.):
        """Class allowing to dock a ligand into a rigid receptor
        ligand has to be a ProDy molecule with a .filename attribute and come
        from a valid pdbqt file of a ligand
        Affinity mapFiles are to be provided using either the mapsFolder or
        the mapFiles argument. The mapsFolder argument is a path to a folder
        containing AutoDock affinity maps. If mapFiles is None, we will try
        to read files matching the '*.$ADTYPE$.map'. Else mapFiles has to
        provide a file name for each atom type as a dict with key atom type
        and value a filename.
        if tpointsFilename is None we will look for a unique file matching
        *anchor*.npy in mapsFolder. Else if tpointsFilename starts with a path
        separator we will read this file, else read the file in mapsFolders.
        If receptorFilename is provided, it means we have a flexible receptor and
        flexibleResidues has to provide a list of flexible side chains e.g.
        [('GLU', '167'), ('TYR', '628'), ('ARG', '87')]
        FTRecsrc is None ot a path to a .py file which when executed defines 'rft'
        which has to be an RFT instance.
        """

        self.folder = '.'

        self.ligand = ligand
        self.receptor = receptor
        # handle covalent indices
        if covalentIndices is not None:
            l1, l2, l3, r1, r2, r3 = covalentIndices
            al1 = self.ligand.select('serial %d'%l1)
            al2 = self.ligand.select('serial %d'%l2)
            al3 = self.ligand.select('serial %d'%l3)
            ar1 = self.receptor.select('serial %d'%r1)
            ar2 = self.receptor.select('serial %d'%r2)
            ar3 = self.receptor.select('serial %d'%r3)
            al1c = al1.getCoords()[0]
            al2c = al2.getCoords()[0]
            al3c = al3.getCoords()[0]
            ar1c = ar1.getCoords()[0]
            ar2c = ar2.getCoords()[0]
            ar3c = ar3.getCoords()[0]
            print 'covalent ligand fitting Lig (%d) %c:%s%d:%s to Rec (%d) %c:%s%d:%s and Lig (%d) %c:%s%d:%s to Rec (%d) %c:%s%d:%s'%(
                al2.getIndices()[0], al2.getChids()[0], al2.getResnames()[0], al2.getResnums()[0], al2.getNames()[0],
                ar2.getIndices()[0], ar2.getChids()[0], ar2.getResnames()[0], ar2.getResnums()[0], ar2.getNames()[0],
                al3.getIndices()[0], al3.getChids()[0], al3.getResnames()[0], al3.getResnums()[0], al3.getNames()[0],
                ar3.getIndices()[0], ar3.getChids()[0], ar3.getResnames()[0], ar3.getResnums()[0], ar3.getNames()[0])

            # get indices for atoms al1, al2, al3 in the prody structure
            covLigAnchorAt = [self.ligand._ag[al1.getIndices()[0]],
                              self.ligand._ag[al2.getIndices()[0]],
                              self.ligand._ag[al3.getIndices()[0]]]
            
            # compute 4th point using vector normal the al1-al2-al3 plan
            v1 = numpy.cross( al1c-al2c, al3c-al2c)
            v1 = v1/sqrt(numpy.sum(v1*v1))
            l4c = al2c+v1
            # compute 4th point using vector normal the ar1-ar2-ar3 plan
            v2 = numpy.cross( ar1c-ar2c, ar3c-ar2c)
            v2 = v2/sqrt(numpy.sum(v2*v2))
            r4c = ar2c+v2

            # align l1,l2,l3,l4 and r1,r2,r3,r4
            from mglutil.math.rigidFit import RigidfitBodyAligner
            aligner = RigidfitBodyAligner( numpy.array( [ar1c, ar2c, ar3c, r4c] ) )
            aligner.rigidFit( numpy.array( [al1c, al2c, al3c, l4c] ) )
            #rmsd = aligner.rmsdAfterSuperimposition(numpy.array( [al1c, al2c, al3c, l4c] ))
            tcoords = aligner.transformCoords(self.ligand._ag.getCoords())[:, :3]
            
            # set ligand co0rdinate to have atoms l1 l2 l3 be superimposed to
            # receptor atoms r1 r2 r3
            self.ligand._ag.setCoords(tcoords)
        else:
            covLigAnchorAt = None

        # build FT for ligand
        from ADFR import ligandFTfromFile
        ligFT = ligandFTfromFile(self.ligand.filename, mol=self.ligand,
                                 covalentLigandAnchorAtoms=covLigAnchorAt,
                                 fixedRoot=fixedRoot, neighborSearchCutoff=neighborSearchCutoff)
        self.ligandFT = ligFT
        if FTRecsrc is not None:
            d = {}
            execfile(FTRecsrc, d, d)
            self.RFTGen = d['rft']
            self.receptorFT = self.RFTGen.ftRoot
           
            # create a single FT containing both the ligand and receptor FTs
            from ADFRcc.adfr import FTBase
            self.FT = FTBase('root')
            self.FT.nbSoftRotamers = len(self.RFTGen.motions)-1
            self.FT.addChild(self.ligandFT.ftRoot)
            self.FT.addChild(self.receptorFT)

            self.FT.motions = self.ligandFT.motions + self.RFTGen.motions
            # self.FT.allMotions contains all FTTorsion motion objects (i.e. flattened for C++)
            self.FT.allMotions = self.ligandFT.motions + self.RFTGen.allMotions
            self.FT.ftRoot = self.FT
            self.FT.receptorAtomSet = self.RFTGen.atomSet
            
        elif flexibleResidues is not None:
            assert self.receptor is not None, "ERROR: flexible residues require providing the receptor pdbqt file"
            # FIXME here we should check that the specified flexbile residues
            # do exist
            self.flexibleResidues = flexibleResidues
            from .ReceptorFT import MkReceptorFT
            self.RFTGen = MkReceptorFT()
            try:
                self.RFTGen.mkFT(self.receptor, flexibleResidues)
            except ValueError, e:
                sys.stderr.write(e.message)
                sys.exit(1)
            self.receptorFT = self.RFTGen.ftRoot
            # create a single FT containing both the ligand and receptor FTs
            from ADFRcc.adfr import FTBase
            self.FT = FTBase('root')
            self.FT.nbSoftRotamers = len(self.RFTGen.motions)
            self.FT.addChild(self.ligandFT.ftRoot)
            self.FT.addChild(self.receptorFT)

            # self.RFTGen.motions contains SoftRotamer which each contains torsions
            self.FT.motions = self.ligandFT.motions + self.RFTGen.motions
            # self.FT.allMotions contains all FTTorsion motion objects (i.e. flattened for C++)
            self.FT.allMotions = self.ligandFT.motions + self.RFTGen.allMotions
            self.FT.ftRoot = self.FT
            self.FT.receptorAtomSet = self.RFTGen.atomSet
        else:
            self.FT = self.ligandFT
            self.FT.allMotions = self.ligandFT.motions
            self.FT.nbSoftRotamers = 0
            self.receptor = None
            self.flexibleResidues = None
            self.RFTGen = None
            self.receptorFT = None
            self.FT.receptorAtomSet = None

        if len(self.FT.allMotions)==0:
            sys.stderr.write("ERROR: no motion object in Flexibility Tree. There are no variables to be optimized\n")
            sys.exit(1)
        
        self.FT.ligandAtomSet = self.ligandFT.atomSet
        
        if logFile is not None:
            sys.stdout = open(logFile, 'w')
            self.folder = os.path.split(logFile)[0]

        # build scorer
        #if self.receptorFT is None:
        if self.RFTGen is None:
            from ADFRcc.adfr import RigidReceptorScorer, GridMap
            scorer = RigidReceptorScorer()
            scorer.setLigandAtomSet(ligFT.atomSet)
            #print 'HACKKKKKKKK' 
            #scorer.setLlCoefficient(0.1)
        else:
            from ADFRcc.adfr import FlexibleReceptorScorer, GridMap
            scorer = FlexibleReceptorScorer()
            scorer.setLigandAtomSet(ligFT.atomSet)
            scorer.setFrAtomSet(self.RFTGen.atomSet)
            #nbfr = numpy.sum([len(x[1]) for x in self.flexibleResidues])
            nbfr = self.FT.nbSoftRotamers
            coef = 1.0/nbfr
            scorer.setFrrrCoefficient(coef)
            scorer.setFrfrCoefficient(coef)
        self.scorer = scorer

        self.maps = []
        if mapFilesRoot is None:
            mapFiles = {}
	    atypes = self.ligand._adtypes#.tolist()
	    if self.receptorFT:
		atypes += self.RFTGen._adtypes
	    atypes = numpy.unique(atypes).tolist()
            for atype in atypes + ['ELECTROSTATIC', 'DESOLVATION']:
                fatype = atype
                if atype=='NX': fatype='N'
                elif atype=='N2': fatype='NA'
                elif atype=='OX': fatype='OA'
                elif atype=='SX': fatype='SA'
                elif atype=='ELECTROSTATIC': fatype='e'
                elif atype=='DESOLVATION': fatype='d'
                mapfilename = glob(os.path.join(mapsFolder,"*.%s.map"%fatype))
                assert len(mapfilename)==1
                mapFiles[atype] = os.path.split(mapfilename[0])[1]
        else:
            mapFiles = {}
            atypes = self.ligand._adtypes
            if self.receptorFT:
                atypes += self.RFTGen._adtypes#atoms.getData('AD_element')
            atypes = numpy.unique(atypes).tolist()

            for atype in atypes + ['ELECTROSTATIC', 'DESOLVATION']:
                fatype = atype
                if atype=='NX': fatype='N'
                elif atype=='N2': fatype='NA'
                elif atype=='OX': fatype='OA'
                elif atype=='SX': fatype='SA'
                elif atype=='ELECTROSTATIC': fatype='e'
                elif atype=='DESOLVATION': fatype='d'
                mapfilename = os.path.join(mapsFolder,"%s.%s.map"%(mapFilesRoot, fatype))
                mapFiles[atype] = os.path.split(mapfilename)[1]
            
        for mapType, mapFilename in mapFiles.items():
          _map = GridMap()
          _map.loadFromMapFile(mapType, mapsFolder, mapFilename)
          self.maps.append(_map)
          scorer.addGridMap(_map)

        #if tpointsFilename is None:
        #    tpointsFilename = glob(os.path.join(mapsFolder,
        #                                        "*anchorPoints*.npy"))
        #    assert len(tpointsFilename)==1
        #   tpointsFilename = tpointsFilename[0]
        from ADFRcc.adfr import FTDiscreteTranslation
        if not covLigAnchorAt and tpointsFilename is not None and fixedRoot is False and neighborSearchCutoff < 0:
            goodPoints = numpy.load(tpointsFilename)
            print 'Using %d translation points from %s'%(len(goodPoints),
                                                     tpointsFilename)
            ligFT.motions[0].setPreferredPoints(goodPoints)

        elif isinstance(ligFT.motions[0],FTDiscreteTranslation) and  neighborSearchCutoff < 0:
            ligFT.motions[0].setPreferredPoints([[0,0,0]])
        
        if not covLigAnchorAt and not fixedRoot:
            cx, cy, cz = self.maps[0].getCenterPy()
            spacing = self.maps[0].getDistBetweenGridPoints()
            sx, sy, sz = self.maps[0].getNumGridPointsPy()
            mini = [ cx-((sx-1)/2)*spacing,
                     cy-((sy-1)/2)*spacing,
                     cz-((sz-1)/2)*spacing]
            maxi = [ cx+((sx-1)/2)*spacing,
                     cy+((sy-1)/2)*spacing,
                     cz+((sz-1)/2)*spacing]
            ligFT.setBox(mini, maxi)

        # handle random seed
        import random
        if seedValue == -1:
            seedValue = random.randint(1,999999)
        random.seed(seedValue)
        self.initialSeed = seedValue
        from ADFRcc.adfr import Randomize
        Randomize.setRandSeed(seedValue)
        
    def createPopulation(self, size, neighborSearchCutoff=-1. ): 
        # create Genome an population
        #print "Builing Genome and Population..."
        from ADFR.GA import GenomePy, IndividualPy, Population, GA
        from .ReceptorFT import FTSoftRotamer
        
        if self.RFTGen:
            scaleRE= 1.0/len(self.RFTGen.motions)
        else:
            scaleRE=None

        genome = GenomePy(self.FT, self.scorer, scaleRE=scaleRE)
        
        #ind = IndividualPy(genome)
        pop = Population()
        from time import time
        t0 = time()
        for i in range(size):
            ind = IndividualPy(genome)
            ind.randomize()
            genes = []
            for motion in self.FT.motions:
                if isinstance(motion, FTSoftRotamer):
                    motion.initialize()
                genes.extend( motion.getGenesPy())
            ind.setGenes(genes)
            ind.score()
            pop.append( ind )
            #print 'ind %i, %f'%(i, ind._score)
        print 'Build pop of %d individuals in %.2f second'%(len(pop), time()-t0)
        return pop
        
    def savePopulation(self, pop, torTree, filename=None):
        if filename is None:
            filename = 'sol_lig'
        for i, ind in enumerate(pop):
            fname = filename+'%04d'%i
            self.writeSolution(ind, torTree, fname, i)
            
    def writeSolution(self, solution, torTree, filename, solNum,
                      recFilename=None, rmsdCalculators=None,
                      status='Unknown'):

        # hack to replace the coordinates in the pdbqt file for ligand
        mol = self.ligand
        score = solution.score() # update tree and coordinates
        # save coordinates of input atomm set
        origCoords = mol._ag.getCoords()
        # set the kligand atom coordinates to current solution
        mol._ag.setCoords(solution.phenotype[1])
        # get lines for PDBQT file
        writer = PDBQTWriter()
        lines = writer.getLines(torTree)
        # restore atom set coordinates
        mol._ag.setCoords(origCoords)
        torsdof = torTree.torsdof

        # add the receptor atoms with segment 'REC '
        class stringBuffer:
            def __init__(self):
                self.lines = []
            def write(self, line):
                if line.startswith('ATOM'):
                    line = line[:67] + 'REC' + line[70:][:-1] # remove new line
                self.lines.append(line)

        # add receptor atom records and move TORSDOF to the end
        if self.RFTGen:
            recAtoms = self.RFTGen.atoms
            recAtoms.setCoords(solution.phenotype[0])
            buf = stringBuffer()
            prody.writePDBStream(buf, self.RFTGen.atoms)
            lastLine = lines[-1]
            lines = lines[:-1]
            lines.extend(buf.lines[1:])
            lines.append(lastLine)

        if not filename.endswith('.pdbqt'):
            filename += '.pdbqt'
        f = open(filename, 'w')

        f.write("USER: ADFR SOLUTION from run %d\n"%solNum)
        import ADFRcc
        f.write("USER: SCORE %f LL: %7.3f LR: %7.3f RR: %7.3f FEB: %7.3f\n"%(
            score, solution.energies['LL'],
            solution.energies['RRL']+solution.energies['FRL'],
            solution.energies['RRFR']+solution.energies['FRFR'],
            solution.energies['RRL'] + solution.energies['FRL'] + \
            torsdof*ADFRcc._parameters.feCoeffTors))

        if rmsdCalculators and len(rmsdCalculators):
            f.write("USER: RMSD %f\n"%rmsdCalculators[0].computeRMSD(solution.phenotype[1]))
        else:
            f.write("USER: RMSD -1\n")

        f.write("USER: STATUS %s\n"%status)
            
        nbMotions = solution.genomePy.getNumMotionNodes()
        f.write("USER: GENES %d\n"%nbMotions)
        for i in range(nbMotions):
            motion = solution.genomePy.getMotionNodeAtIndex(i)
            f.write("USER: genes %d %s |==| "%(i, motion.getName()))
            for gene in motion.getGenesPy():
                f.write(" %.5f"%gene)
            f.write("\n")
        [f.write(l+'\n') for l in lines]
        f.close()

    def createNeighborPopulation(self, size, neighborSearchCutoff):
        from ADFR.GA import GenomePy, IndividualPy, Population, GA
        #from .ReceptorFT import FTSoftRotamer
        from .LigandFT import PyFTTorsion, PyFTDiscreteTranslation, \
             PyFTRotationAboutPointQuat
        from mglutil.math.rotax import rotax, mat_to_quat
        from quat import axisAngleToQuat, qmult
        from math import pi
        from random import random, gauss, uniform
        if self.RFTGen:
            scaleRE= 1.0/len(self.RFTGen.motions)
        else:
            scaleRE=None

        genome = GenomePy(self.FT, self.scorer, scaleRE = scaleRE , neighborSearchCutoff = neighborSearchCutoff)
        ident = genome.getIdentityGenesPy()
        pop = Population()
        neighborPop = 0 # set to True when ind within RMSD cutoff found
        from mglutil.math.rigidFit import RigidfitBodyAligner
        refCoords = self.ligand._ag.getCoords()
        aligner = RigidfitBodyAligner()
        from time import time
        t0 = time()
        ct = 0
        fttrans = ftrot = None
        off = 0
        fttrans = ftrot = None
        for motion in self.FT.motions:
            if isinstance(motion, PyFTDiscreteTranslation):
                tb = off # gene offset for translation
                fttrans = motion
                delta = 1./(numpy.array(fttrans.maxi)-numpy.array(fttrans.mini))		
            elif isinstance(motion, PyFTRotationAboutPointQuat):
                rb = off # gene offset for rotation
                ftrot = motion
            off += motion.getNumVariables()
        assert fttrans is not None, "Error: NeighborhoodPopulation: not Transaltion node found"
        mat4 = numpy.identity(4)
        to = None
        ind = IndividualPy(genome)
        ind.setGenes(ident)
        ind.score()
        #pop.append( ind )
        while len(pop) < size:
            ind = IndividualPy(genome)
            ind.setGenes(ident)
            ind.score() # for all motion object in genome to be configured
            #import pdb;pdb.set_trace()
            if to is None:
                # get translation vector for identity genes
                to = fttrans.getVariablesPy()
            genes = []
            off = 0
            for motion in self.FT.motions:
                nbg = motion.getNumVariables()
                if isinstance(motion, PyFTTorsion):
                    if neighborPop < size * 0.1 or random()< 0.2:
                        motion.randomize()
                    genes.extend( motion.getGenesPy())
                elif isinstance(motion, PyFTDiscreteTranslation):
                    
                    #fttrans.mutate(fttrans.getGenesPy(), 0.3)
                    dx = numpy.array(([uniform(-1,1), uniform(-1,1), uniform(-1,1)]) * delta * neighborSearchCutoff / 4)
                    genes.extend(ident[off:off+nbg]+dx)
                else:
                    genes.extend(ident[off:off+nbg])
                off += nbg
            ind.setGenes(genes)
            ind.score()
            #print ind.genomePy.scorer.getRMSD()
            if ind._neighborRMSD < neighborSearchCutoff:
                neighborPop += 1
                pop.append(ind)
            elif len(pop) < size * 0.9 or neighborPop > size * 0.1:
                pop.append(ind)
            #print len(neighborPop),len(pop)
            #    print ct, len(pop), rmsd, aligner.translationMatrix
        print "built neigborhood population of %d in %.2f seconds"%(size, time()-t0)
        return pop


    def createLocalPopulation(self, size, radius=2.0): 
        # create Genome an population
        #print "Builing Genome and Population..."
        from ADFR.GA import GenomePy, IndividualPy, Population, GA
        genome = GenomePy(self.FT, self.scorer)
        
        ind = IndividualPy(genome)
        pop = Population()
        from random import random, gauss, uniform
        from time import time
        t0 = time()
        idGenes = genome.identityGenes()
        transg = idGenes[:3]
        rotg = idGenes[3:7]
        torsg = idGenes[7:]
        radius2 = 0.5*radius
        while len(pop) < size:
            ind = IndividualPy(genome)#, values=idGenes)            
            trans = genome.motions[0]
            rot = genome.motions[1]
            ct = 0
            while ct < 1:
                delta = []
                for i in range(3):
                    if random()< 0.3:
                        delta.append( uniform(0.0, 0.15) )
                        ct += 1
                    else:
                        delta.append( 0.0 )
                ngenes = trans.applyDelta(transg, delta)

                delta = []
                for i in range(4):
                    if random()< 0.3:
                        delta.append( gauss(0.0, 0.01) )
                        ct += 1
                    else:
                        delta.append( 0.0 )
                ngenes.extend(rot.applyDelta(rotg, delta ))

                for n, tors in enumerate(genome.motions[2:]):
                    if random()< 0.3:
                        ngenes.extend(
                            tors.applyDelta([torsg[n]], [uniform(0.0, 1.)]))
                        ct += 1
                    else:
                        ngenes.append( torsg[n] )#m.getGenesForValues([m.getOrigAngle()]) )
            #import pdb; pdb.set_trace()
            ind.values[:] = ngenes
            ind.score()
            pop.append( ind )
        print 'Build pop of %d individuals in %.2f second'%(len(pop), time()-t0)
        return pop
    
    def createGA(self, pop, referenceLigFile=None, quickmini=None, hardmini=None, seedValue=-1, RMSDMatching='hungarian', neighborSearchCutoff=-1.0):
        from ADFR.GA import GA
        if referenceLigFile:
            refLig = Read(referenceLigFile)
        else:
            refLig = None
        self.referenceLigand = refLig
        
        self.ga = ga = GA(pop, self.ligand, refLig, folder=self.folder,
                          RMSDMatching=RMSDMatching, neighborSearchCutoff=neighborSearchCutoff)
        if quickmini:
            ga.quickmini = quickmini
        if hardmini:
            ga.hardmini = hardmini

        # handle random seed
        import random
        if seedValue == -1:
            seedValue = self.initialSeed
        random.seed(seedValue)

        ga.seed = seedValue
        self.dockingSeed = seedValue
        from ADFRcc.adfr import Randomize
        Randomize.setRandSeed(seedValue)

## def writeSolution(solution, filename, rmsd=None):
##     # hack to replace the coordinates in the pdbqt file
##     mol = solution.genome.ft.mol
##     atoms = mol._ag
##     score = solution.score() # update tree and coordinates
##     f = open(mol.filename)
##     lines = f.readlines()
##     f.close()

##     atomByName = {}
##     coords = solution.phenotype[1]
##     for i, atom in enumerate(atoms):
##         name = atom.getName()
##         assert not atomByName.has_key(name), "ERROR atoms with identical names %s"%name
##         atomByName[atom.getName()] = coords[i]

##     ct = 0
##     for i in range(len(lines)):
##         if lines[i].startswith('HETATM') or lines[i].startswith('ATOM'):
##             atName = lines[i][12:16].strip()
##             x,y,z = atomByName[atName]
##             lines[i] = lines[i][:30] + '%8.3f%8.3f%8.3f'%(x,y,z) + lines[i][54:]

##     f = open(filename, 'w')
##     f.write("USER: ADFR SOLUTION 0\n")
##     f.write("USER: ADFR SCORE %f\n"%(0.-score))
##     if rmsd is not None:
##         f.write("USER: ADFR RMSD %f\n"%rmsd)
##     f.write("USER: ADFR genes %s\n"%str(solution.values))
##     for l in lines:
##         if l.startswith('USER: '): continue
##         f.write(l) 
         
##     f.close()

class ADFRRigidReceptor(ADFR):
    pass
