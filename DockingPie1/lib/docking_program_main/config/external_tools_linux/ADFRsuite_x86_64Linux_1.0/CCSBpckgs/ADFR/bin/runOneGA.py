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
# $Header: /mnt/raid/services/cvs/ADFR/bin/runOneGA.py,v 1.20.2.1 2017/08/21 22:23:11 annao Exp $
#
# $Id: runOneGA.py,v 1.20.2.1 2017/08/21 22:23:11 annao Exp $
#

import sys, traceback

from ADFR.utils.optParser import ArgParser
parser = ArgParser('OneGA')
args = parser.parse_args()

# use 4.2 by default
#from ADFRcc import setForceFieldVersion
#parameters = setForceFieldVersion('4.2')

if len(args)<1:
    parser.print_help()
    sys.exit(1)
    
from ADFR import ADFR
import os

seedValue = args['seedValue']

from time import time
from ADFR import getLigandFromFile
ligand, error = getLigandFromFile(args['ligandFile'])
if error:
    print error
    sys.exit(1)

mapsFolder = None
reference = args['reference']

neighborSearchCutoff = args['neighborSearchCutoff']
neighborSearchGroup = args['neighborSearchGroup']
if neighborSearchGroup is None:
    neighborSearchGroup = "all"
    
flexRes = args['flexRes']
fixedRoot = args['fixedRoot']
if flexRes is not None:
    from ADFR.utils.maps import flexResStr2flexRes
    flexRes = flexResStr2flexRes(flexRes)
    
covalentLigand = args.get('covalentLigand', None)
if args['receptorMapsFile'] is not None:
    from ADFR.utils.maps import MapsFile, flexResStr2flexRes
    mf = MapsFile(args['receptorMapsFile'])
    mf.unzipMaps()
    mapsFolder = mf.getMapsFolder()
    receptorFilename = os.path.join(mf.getMapsFolder(),
                                    mf.getReceptorFilename())
    flexResStr = mf.getFlexResStr()
    covalentRec = mf.getCovalentBond()
    #mapsFolder, receptorFilename, flexRes, msg, amap = unzipMaps(
    #    args['receptorMapsFile'])
    #flexRes = flexResStr2flexRes(flexRes)
    #f = open(os.path.join(mapsFolder, 'data.pkl'))
    #gridData = pickle.load(f)
    #f.close()
    #covalentRec = d.get('covalentBond', None)
    if covalentRec is not None:
        covalentRec += [d.get('covalentBondTorsionAtoms')[-1]]
    else:
        covalentRec = None    
    if mapsFolder is None:
        print msg
        sys.exit(1)
    tPtsFile = os.path.join(mapsFolder, "translationPoints.npy")
    if os.path.exists(tPtsFile):
        tpointsFilename = os.path.join(mapsFolder, "translationPoints.npy")
    else:
        tpointsFilename = None
    mapFilesRoot = 'rigidReceptor'

else:
    mapsFolder = args['mapFilesFolder']
    receptorFilename = args['receptorFile']
    mapFilesRoot = args['mapFilesRoot']
    tpointsFilename = args.get('tpointsFilename', None)
    cr = args.get('covalentReceptor', None)
    covalentRec = args.get('covalentRec', None)

if covalentLigand is not None:
    if covalentRec is None:
        raise RuntimeError("ERROR: covalentLigand specified but no covalentRec")
elif covalentRec is not None:
    if covalentLigand is None:
        raise RuntimeError("ERROR: covalentRec specified but covalentLigand")

if covalentLigand is not None and covalentRec is not None:
    covalentIndices = covalentLigand + covalentRec
else:
    covalentIndices = None

if receptorFilename is not None:
    from MolKit2 import Read
    receptor = Read(receptorFilename)
else:
    receptor = None

t00= time()
try:
    adfr = ADFR(ligand, mapsFolder, logFile=args['logfile'],
            receptor=receptor,
            flexibleResidues=flexRes,
            mapFilesRoot=mapFilesRoot,
            tpointsFilename=tpointsFilename,
            seedValue=seedValue,
            covalentIndices=covalentIndices,
            FTRecsrc=args.get('FTRecsrc', None),
            fixedRoot=fixedRoot,
            neighborSearchCutoff=neighborSearchCutoff)
except Exception, e:
    print traceback.format_exc()
    sys.exit(1)
#           torsionInfoFilename=args['torsInfo)

print 'receptor', receptorFilename
print 'flexRes', flexRes
print 'mapFilesRoot', mapFilesRoot
print 'tpointsFilename', tpointsFilename
print 'covalentIndices', covalentIndices
print 'RMSDMatching', args['RMSDMatching']
print 'MAXGENS', args['maxGens']
print 'MAXEVALS', args['maxEvals']
print 'NOIMPROVESTOP', args['noImproveStop']
if args['popSize'] is 'auto':
    nLigGenes = 7+adfr.ligandFT.torTree.torsdof
    popSize = 50 + nLigGenes*10
else:
    popSize = args['popSize']
if neighborSearchCutoff>0:
    adfr.scorer.setNeighborRMSDcutoff(neighborSearchCutoff,
                         ligand._ag.select(neighborSearchGroup).getIndices())
    pop = adfr.createNeighborPopulation(popSize, neighborSearchCutoff)
else:
    pop = adfr.createPopulation(popSize, neighborSearchCutoff)
#adfr.savePopulation(pop, 'test')
#import pdb; pdb.set_trace()

#ind = pop[0]
#if args['swBiasCoeff1:
#    sw.biasCoef1 = args['swBiasCoeff1
#if args['swBiasCoeff2:
#    sw.biasCoef2 = args['swBiasCoeff2

adfr.createGA(pop, reference, RMSDMatching=args['RMSDMatching'], neighborSearchCutoff=neighborSearchCutoff)
adfr.ga.setDebug(args['debug'])

from ADFRcc import getFFParameters
parameters = getFFParameters()
print "ForceField cutoffs: vdw %4.2f estat %4.2f hbond %4.2f desolv %4.2f\n"%(
    parameters.maxDistVdw, parameters.maxDistEstat,
    parameters.maxDistHbond, parameters.maxDistSolv)

stat = adfr.ga.evolve(maxGens=args['maxGens'],
                      maxEvals=args['maxEvals'],
                      noImproveStop=args['noImproveStop'])

outputpop=[]

if neighborSearchCutoff>0:
    adfr.scorer.setNeighborRMSDcutoff(-1.0,ligand._ag.select(neighborSearchGroup).getIndices())
    for finalpop in pop:
        finalpop.genomePy.neighborSearchCutoff=-1.0
        #sscore=finalpop._score
        finalpop.score()
        #finalpop.genomePy.neighborSearchCutoff=neighborSearchCutoff
        if finalpop._neighborRMSD<neighborSearchCutoff:
            outputpop.append(finalpop)
    if len(outputpop)>0:
        #pop=outputpop
        print "Rescore populations with rmsd larger than RMSD cutoff left with %3d populations"%(len(outputpop))
    else:
        print "Nothing left"
    pop.sort()

# build Hungarian matching RMSD calculator regarless of matching method used
# in GA
if adfr.referenceLigand:
    atoms = adfr.referenceLigand.select()
    from MolKit2.molecule import getAtomIndicesPerType
    from mglutil.math.rmsd import HungarianMatchingRMSD_prody
    atoms = adfr.referenceLigand.select()
    d1 = getAtomIndicesPerType(atoms)
    d2 = getAtomIndicesPerType(adfr.ligand.select())
    rmsdCalc = HungarianMatchingRMSD_prody(atoms.getCoords(), d1, d2)
    rmsd = rmsdCalc.computeRMSD(pop[0].phenotype[1])
else:
    rmsd = -1
    
print 'evolution stopped with status %s after %s seconds'%(stat, time()-t00)
print 'END bestScore: %.3f RMSD: %.2f nbgen: %3d nbEvals %8d time %.3f'%(
    -pop[0]._score, rmsd, adfr.ga.gen, adfr.scorer.getNumEvals(), time()-t00)
topSol = pop[0]
print pop[0].energies

if args.get('logfile', None):
    import os
    folder = os.path.split(args['logfile'])[0]
else:
    folder = '.'

if args.get('jobNumber',None) is None:
    jobNum=adfr.initialSeed
else:
    jobNum = args['jobNumber']

torTree = adfr.ligandFT.torTree
if adfr.referenceLigand:
    adfr.writeSolution(topSol, torTree, os.path.join(
        folder, '%s%04d'%(args['jobName'],jobNum)),
                       jobNum, rmsdCalculators=[rmsdCalc], status=stat)
else:
    adfr.writeSolution(topSol,  torTree, os.path.join(
        folder, '%s%04d'%(args['jobName'],jobNum)), jobNum, status=stat)

print 'GA %d terminated successfully'%jobNum

