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

########################################################################
#
# Date: 2016 Authors: Michel Sanner
#
#    sanner@scripps.edu
#
#       The Scripps Research Institute (TSRI)
#       Molecular Graphics Lab
#       La Jolla, CA 92037, USA
#
# Copyright: Michel Sanner and TSRI
#
#########################################################################
#
# $Header: /mnt/raid/services/cvs/ADFR/dro.py,v 1.12.2.1 2017/10/26 23:23:22 annao Exp $
#
# $Id: dro.py,v 1.12.2.1 2017/10/26 23:23:22 annao Exp $
#

import tempfile, tarfile, os, shutil, numpy, pickle
from ADFR.utils.maps import MapsFile
from ADFR.utils.maps import flexResStr2flexRes
from MolKit2 import Read
from MolKit2.selection import Selection

from ADFRcc.adfr import RigidReceptorScorer, Parameters
parameters = Parameters.getParameters()
feCoeffVdw = parameters.feCoeffVdw
feCoeffHbond = parameters.feCoeffHbond
feCoeffEstat = parameters.feCoeffEstat
feCoeffDesolv = parameters.feCoeffDesolv
feCoeffTors = parameters.feCoeffTors

class DockingResultsObject:
    """
    proxy for a docking result file produced by ADFR
    """

    def __init__(self):
        self.filename = None
        self._tmpFolder = None # temp folder into which dro is extracted
        self._topFolder = None # top ofrlder of .dro fiel
        self._resultsFolder = None # folder holding solutions and dlgs of GAs
        self._inputFolder = None # folder containing input files
        self._dockingData = None
        self._gridData = None
        self._mapTypes = []
        self._mapOrigin = None
        self._mapSize = None
        self._mapSpacing = None
        self._mapCenter = None
        self._rigidDocking = False
        self._lig = None
        self._rec = None
        self._solutions = None

    def load(self, filename):
        """
        open and extract the .dro file
        """
        # FIXME we should use MapsFile object here
        self.filename = filename
        
        tar = tarfile.open(filename)
        self._tmpFolder = tempfile.mktemp()
        tar.extractall(path=self._tmpFolder)
        folder = self._topFolder = os.path.join(self._tmpFolder,
                                                os.listdir(self._tmpFolder)[0])
        # get docking data
        with open(os.path.join(folder, 'data.pkl')) as f:
            self._dockingData = pickle.load(f)
        self._resultsFolder = os.path.join(folder, 'dockingDetails')
        self._inputFolder = os.path.join(folder, 'input')

        self._lig = Read(os.path.join(self._inputFolder,
                                      self._dockingData['ligandFilename']))
        self._rec = Read(os.path.join(self._inputFolder,
                                      self._dockingData['fullReceptorFilename']))
        self._rec._ag.setSegnames(['REC ']*len(self._rec._ag))

        solutionsFile = os.path.join(self._topFolder,"%s_out.pdbqt"%
                                     self._dockingData['basename'])
        self._solutions = Read(solutionsFile)
        self._solGenes = self.getSolutionGenes(solutionsFile)
        # get grid data

        mf = MapsFile(os.path.join(self._inputFolder, self._dockingData['receptorMapsFile']))
        mf.unzipMaps()
        self.mf = mf
        self._gridData = mf.getData()
        self._mapSize = self._gridData['boxSize']
        self._mapLength = self._gridData['boxLengths']
        self._mapCenter = self._gridData['boxCenter']
        self._mapOrigin = numpy.array(self._mapCenter) - 0.5*(numpy.array(self._mapLength))
        self._mapSpacing = self._gridData['spacing']
        self._flexRes = self._gridData.get('flexResStr', None)
        if self._flexRes is not None:
            self._flexRes = flexResStr2flexRes(self._flexRes)

        if self._gridData['covalentLigandAtomIndices']:
            # tag receptor atoms forming covalent ligand starting from second
            # index which is the the second atom in the covalent bond
            indices = self._gridData['covalentLigandAtomIndices'][1:]
            self._rec._ag._data['segment'][indices] = 'CLIG'

        # add conformations to the receptor
        if self._flexRes:
            FRAtInd = self._dockingData['FRAtomsIndices']
            coords = self._rec._ag.getCoords()
            solRecAtInd = self._solutions.select('segment REC').getIndices()
            for i in range(self._solutions.numMols()):
                if i < self._solutions.numMols()-1:
                    self._rec._ag.addCoordset(coords,
                                              label='solution %d')
                self._rec._ag._coords[i][FRAtInd] = \
                      self._solutions._ag._coords[i][solRecAtInd]

            self._rec._multi = 'conformations'
            self._solutions._ag._flags['deleted'][solRecAtInd] = True
            # make atom set for FR Atoms
            atoms = Selection(self._rec._ag, FRAtInd, 'FlexRec')
            #self._FRAtoms = atoms.select('not backbone and not name CB')
            self._FRAtoms = atoms.select('sidechain or ca')

    def getSolutionGenes(self, filename):
        if self._solutions is None: return
        solGenes = []
        with open(filename) as f: 
            for line in f:
                if line.startswith('MODEL '):
                    genes = []
                elif line.startswith('ENDMDL'):
                    solGenes.append(genes)
                elif line.startswith('USER: genes '):
                    index = line.find("|==|")+4
                    genes.extend( [float(x) for x in line[index:].split()] )
        return solGenes
    
    def makeScorer(self):
        from ADFR.utils.scorer import ADFRscorer
        folder = self._inputFolder
        receptor = None
        if self._flexRes is not None:
            receptor = self._rec

        if self._gridData.has_key('covalentBond') and self._gridData['covalentBond']:
            # build l1 l2 l3 r1 r2 r3
            covalentIndices = self._dockingData['covalentLigand']
            l3 = int(self._gridData['covalentBondTorsionAtom'].split()[1][1:-1])
            covalentIndices.append(l3)
            covalentIndices.extend(self._gridData['covalentBond'])
            receptor = self._rec
        else:
            covalentIndices=None
        adfr, ind, _score = ADFRscorer(
            self._lig, self.mf.getMapsFolder(), mapFilesRoot='rigidReceptor',
            receptor=receptor, flexRes=self._flexRes,
            covalentIndices=covalentIndices)

        self._adfr = adfr
        self._ind = ind
        self._pop = adfr.createPopulation(1)
        score = self.scoreSolution(0)

        ## ADelem = adfr.ligandFT.mol._ag._data['AD_element']
        ## RRL = adfr.scorer.getLrrGridScorer().getScoreArray()
        ## for i, a1 in enumerate(adfr.ligandFT.mol.select()):
        ##     x, y ,z = a1.getCoords()
        ##     print "%4s %2s %9.3f %9.3f %9.3f %9.3f"%(
        ##         a1.getName(), ADelem[a1.getIndex()], RRL[i],x, y, z)
        #print 'SCORE', score
        
    def scoreSolution(self, num):
        if self._solutions is None: return
        if num<0 or num> self._solutions.numMols():
            return
        self._ind.setGenes(self._solGenes[num])
        return self._ind.score()

    def perAtomInteractionEnergy(self):
        return self._ind.genomePy.scorer.getLrrGridScorer().getScoreArray()

    def getLigandComposition(self):
        allTypes = self._lig._ag.getData("AD_element")
        types = numpy.unique(allTypes)
        d = {}.fromkeys(types, 0)
        for a in allTypes:
            d[a] += 1
        return d

    def dictToString(self, d):
        s = ''
        for k,v in d.items():
            s += '%s:%s '%(k,str(v))
        return s

    def printInfo(self):
        print 'docking result file'
        print '  date       : %s'%self._dockingData['date']
        print '  node       : %s'%self._dockingData['node']
        print '  platform   : %s'%self._dockingData['platform']
        print '  ncores     : %s'%self._dockingData['ncores']
        print
        print '  receptor   : %s'%self._dockingData['receptorMapsFile']
        self.mf.printInfo(indent='   ')
        print
        print '  ligand     : %s'%self._dockingData['ligandFilename']
        print '     nbAtoms : %d'%len(self._lig._ag)
        print '     types   : %s'%self.dictToString(self.getLigandComposition())
        print '  lig ref    : %s'%self._dockingData['referenceLigName']
        print '  summary    : %s'%self._dockingData['summaryFile']
        print '  GA params  : '
        print '    maxEvals : %s'%self._dockingData['maxEvals']
        print '    nbRuns   : %s'%self._dockingData['nbRuns']
        print '    jobName  : %s'%self._dockingData['jobName']
        print '    seed     : %s'%self._dockingData['seed']
        print '    covLig   : %s'%self._dockingData['covalentLigand']
         
if __name__=="__main__":
    name = '../../tests/ADFR/rigid_and_flexible/4EK4_random_rigid.dro'
    dro = DockingResultsObject()
    dro.load(name)
    print dro._data
    print dro._mapTypes
    print dro._mapSize
    print dro._mapOrigin
    print dro._mapCenter
    print dro._mapSpacing
    print dro._flexRes
    print dro._solutions.numMols()
    print dro._solGenes
    dro.makeScorer()
    print dro.scoreSolution(0)
    print dro.scoreSolution(1)
    #genes = [0.52521, 0.51207, 0.49969, 0.75501, 0.58870, 0.09668, 0.37984,
    #         0.55289, 0.26428]
    #dro._ind.setGenes(genes)
    #_score = dro._ind.score()
    RRL =  dro._ind.genomePy.scorer.getLrrGridScorer().getScoreArray()
    ADelem = dro._adfr.ligandFT.mol._ag._data['AD_element']
    for i, a1 in enumerate(dro._adfr.ligandFT.mol.select()):
        x, y ,z = a1.getCoords()
        print "%4s %2s %9.3f %9.3f %9.3f %9.3f"%(
            a1.getName(), ADelem[a1.getIndex()], RRL[i],x, y, z)
        
