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
# $Header: /mnt/raid/services/cvs/ADFR/bin/ade.py,v 1.7 2016/12/07 00:38:32 sanner Exp $
#
# $Id: ade.py,v 1.7 2016/12/07 00:38:32 sanner Exp $
#
import os, sys
from time import time

from optparse import OptionParser
parser = OptionParser(usage="usage: %prog [options] filename",
                      version="%prog 0.1")
parser.add_option("-l", "--ligand",
                  action="store", # optional because action defaults to "store"
                  dest="ligand",
                  help="docked ligand filename",)
parser.add_option("--flexibleReceptor",
                  action="store", # optional because action defaults to "store"
                  dest="flexibleReceptor",
                  help="pdb file of flexible receptor",)
parser.add_option("-t", "--target",
                  action="store", # optional because action defaults to "store"
                  dest="receptorMapsFile",
                  help=".zip file with receptor and map",)
parser.add_option("--mapsFolder",
                  action="store", # optional because action defaults to "store"
                  dest="mapsFolder",
                  help="folder containing receptor affinity maps",)
parser.add_option("--receptor",
                  action="store", # optional because action defaults to "store"
                  dest="receptorFile",
                  help="receptor PDBQT file",)
parser.add_option("--flexRes",
                  action="store", # optional because action defaults to "store"
                  dest="flexRes",
                  help="list of flexbile receptor side chains e.g. [('GLU', '167'), ('TYR', '628')]",)
parser.add_option("--mapFilesRoot",
                  action="store", # optional because action defaults to "store"
                  dest="mapFilesRoot",
                  help="list of flexbile receptor side chains e.g. [('GLU', '167'), ('TYR', '628')]",)
parser.add_option("-c", "--covalentLigand",
                  action="store", 
                  dest="covalentLigand",
                  type="int",
                  nargs=3,
                  help="specify l1 l2 l3 indices of ligand and receptor atoms. Ligand atoms l2 and l3 will be overlapped with receptor atoms r2 and r3. ligand atom l1 has to correspond to receptor atom r1 and is use to define first ligand torsion as l1-l2-l3-neighbor(l2).")
(options, args) = parser.parse_args()


if options.receptorMapsFile is not None:
    #from ADFR.utils.maps import unzipMaps
    #mapsFolder, receptorFilename, flexRes, msg, amap = unzipMaps(options.receptorMapsFile)
    from ADFR.utils.maps import MapsFile
    mf = MapsFile(options.receptorMapsFile)
    mf.unzipMaps()
    mapsFolder  = mf.getMapsFolder()
    receptorFilename = receptorFilename = os.path.join(mf.getMapsFolder(),
                                                       mf.getReceptorFilename())
    flexResStr = mf.getFlexResStr()
    if flexResStr:
        from ADFR.utils.maps import flexResStr2flexRes
        flexRes = flexResStr2flexRes(flexResStr)
    else:
        flexRes=None
        
    if mapsFolder is None:
        print msg
        sys.exit(1)
        
    tpointsFilename = os.path.join(mapsFolder, "translationPoints.npy")
    mapFilesRoot = 'rigidReceptor'
else:
    if options.mapsFolder is not None:
        mapsFolder = options.mapsFolder

    if options.receptorFile is not None:
        receptorFilename = options.receptorFile
    else:
        receptorFilename = None

    if options.tpointsFilename is not None:
        tpointsFilename = options.tpointsFilename
    else:
        from glob import glob
        tpointsFilename = glob(os.path.join(mapsFolder, "*anchorPoints*.npy"))
        if len(tpointsFilename)==1:
            tpointsFilename = tpointsFilename[0]
        elif options.covalentLigand is None:
            print "ERROR: no valid translation points"
            sys.exit(1)
    if options.mapFilesRoot is not None:
        mapFilesRoot = options.mapFilesRoot
    else:
        mapFilesRoot = None

    if options.flexRes is not None:
        from ADFR.utils.maps import flexResStr2flexRes
        flexRes = flexResStr2flexRes(options.flexRes)
    else:
        flexRes = None

from MolKit2 import Read
ADRecElem = []
recAtomsIndices = []
if receptorFilename is not None:
    receptor = Read(receptorFilename)
    if options.flexibleReceptor is not None:
        conf = Read(options.flexibleReceptor)
        hv = conf._ag.getHierView()
        for res in hv.iterResidues():
            chid = res.getChid()
            resname = res.getResname()
            resnum = res.getResnum()
            for atom in res.iterAtoms():
                selStr = 'chid "%c" resname %s resnum %d name %s'%(
                    chid, resname, resnum, atom.getName())
                atomInRec = receptor.select(selStr)
                print selStr, atom.getCoords(), atomInRec.getCoords()
                atomInRec.setCoords([atom.getCoords()])
                ind = atomInRec.getIndices()[0]
                ADRecElem.append(receptor._ag._data['AD_element'][ind])
                recAtomsIndices.append( ind )

#from ADFRcc.adfr import Parameters
#parameters = Parameters.getParameters()
#parameters.setUseTables(False)

# use 4.2 by default
#from ADFRcc import setForceFieldVersion
#parameters = setForceFieldVersion('4.2')

covalentIndices = options.covalentLigand
if covalentIndices is not None:
    assert options.receptorFile is not None or options.receptorMapsFile is not None, "ERROR: covalent docking requires providing the receptor molecule"
    covalentIndices = list(covalentIndices)
    covalentIndices.append(int(mf.getCovalentBondTorsionAtom().split()[1][1:-1]))
    covalentIndices.extend(mf.getCovalentBond())
    
from ADFR import getLigandFromFile
ligand, error = getLigandFromFile(options.ligand)
if error:
    print error
    sys.exit(1)
    
from ADFR import ADFR
t00= time()
adfr = ADFR(ligand, mapsFolder, 
            receptor=receptor,
            flexibleResidues=flexRes,
            mapFilesRoot=mapFilesRoot,
            tpointsFilename=tpointsFilename,
            covalentIndices=covalentIndices)

#adfr = ADFRRigidReceptor(options.ligand, options.mapsFolder)
#ind = adfr.createPopulation(1)[0]
from ADFR.GA import GenomePy, IndividualPy, Population
if adfr.RFTGen:
    scaleRE= 1.0/len(adfr.RFTGen.motions)
else:
    scaleRE=None
    
genome = GenomePy(adfr.FT, adfr.scorer, scaleRE=scaleRE)

ind = IndividualPy(genome)
ind.setGenes(ind.genomePy.getIdentityGenesPy())
_score = ind.score()

import numpy
diffCoords = numpy.sum(adfr.ligand._ag.getCoords() - ind.phenotype[1])
assert diffCoords < 0.0001, "ValueError: phenotype coords do not match input ligand coordinates"
#adfr.writeSolution(ind, 'adeLig', 'adeRec')

#print os.path.split(options.ligand)[1][:4],'SCORE', -ind.score(),ind.energies
#print ind.energies

from ADFRcc.adfr import RigidReceptorScorer, GridMap, Parameters
parameters = Parameters.getParameters()
feCoeffVdw = parameters.feCoeffVdw
feCoeffHbond = parameters.feCoeffHbond
feCoeffEstat = parameters.feCoeffEstat
feCoeffDesolv = parameters.feCoeffDesolv
feCoeffTors = parameters.feCoeffTors

def LLpairs(ind, ligAtoms):
    #_score = ind.score()
    scorer = genome.scorer.getLlPairwiseScorer()
    energiesLL = scorer.getTotalScore()
    distance = scorer.getDistanceMatrix()
    eArray = scorer.getEstatEnergyMatrix()
    vdwArray = scorer.getVdwEnergyMatrix()
    hArray = scorer.getHbEnergyMatrix()
    dsArray = scorer.getSolvEnergyMatrix()
    atomSet = scorer.getAtomSet1().atomSetStatic
    
    print "	       	   Ligand Intramolecular Energy Analysis"
    print "		   ====================================="
    print
    print "Non-bond      Atom1-Atom2      Distance   Total        Elec     vdW        Hb     Desolv"
    print "________  ___________________  ________   ______  _________  ________ ________  ________"

    ct = 1
    esum = 0.
    vsum = 0.
    hsum = 0.
    dsum = 0.
    #import pdb; pdb.set_trace()
    for i, a1 in enumerate(ligAtoms):
        n1 = a1.getName()
        for j, a2 in enumerate(ligAtoms):
            if j<=i: continue
            n2 = a2.getName()
            if atomSet.getPairScorable(i,j):
                e  = eArray[i][j]*feCoeffEstat
                v = vdwArray[i][j]*feCoeffVdw
                h = hArray[i][j]*feCoeffHbond
                d = dsArray[i][j]*feCoeffDesolv
                print "  %5d %4d-%-4d %4s %4s     %7.3f   %7.4f   %7.4f   %7.4f   %7.4f   %7.4f"%(
                    ct, i+1, j+1, n1, n2, distance[i][j], e+v+h+d, e, v, h, d)
                ct += 1
                esum += e
                vsum += v
                hsum +=h
                dsum += d
    print "--------------------------------------------------------------------------------------------------------"
    print "          %25s      %7.4f   %7.4f   %7.4f   %7.4f   %7.4f\n\n\n"%("Sum", esum+vsum+hsum+dsum, esum, vsum, hsum, dsum)

LLpairs(ind, adfr.ligandFT.mol.select())

for atom in adfr.ligandFT.mol._ag:
    x, y ,z = atom.getCoords()
    print "%4s %9.3f %9.3f %9.3f"%(atom.getName(), x, y, z)

print "RR-L ---------------------------------------------------------------"
ADelem = adfr.ligandFT.mol._ag._data['AD_element']
RRL = genome.scorer.getLrrGridScorer().getScoreArray()
for i, a1 in enumerate(adfr.ligandFT.mol.select()):
    x, y ,z = a1.getCoords()
    print "%4s %2s %9.3f %9.3f %9.3f %9.3f"%(
        a1.getName(), ADelem[a1.getIndex()], RRL[i],x, y, z)
print "RR-L END------------------------------------------------------------"

#import pdb; pdb.set_trace()
if flexRes:          
    print "RR-FR ---------------------------------------------------------------"
    FRRR = genome.scorer.getFrrrGridScorer().getScoreArray()
    coords = conf._ag.getCoords()
    for i, atype in enumerate(ADRecElem):
        x, y ,z = coords[i]
        print "%4s %2s %9.3f %9.3f %9.3f %9.3f"%(
            receptor._ag.getNames()[recAtomsIndices[i]], atype, FRRR[i],x, y, z)
    print "RR-L END------------------------------------------------------------"

print ind.energies
print "SCORE:", -ind._score

