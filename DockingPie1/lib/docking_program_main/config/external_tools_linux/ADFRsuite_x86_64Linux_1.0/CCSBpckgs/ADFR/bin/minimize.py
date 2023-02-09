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
# $Header: /mnt/raid/services/cvs/ADFR/bin/minimize.py,v 1.2 2016/12/07 00:38:32 sanner Exp $
#
# $Id: minimize.py,v 1.2 2016/12/07 00:38:32 sanner Exp $
#
import os
from time import time

from optparse import OptionParser
parser = OptionParser(usage="usage: %prog [options] filename",
                      version="%prog 0.1")
parser.add_option("-l", "--ligand",
                  action="store", # optional because action defaults to "store"
                  dest="ligand",
                  help="docked ligand filename",)
parser.add_option("-r", "--reference",
                  action="store", # optional because action defaults to "store"
                  dest="reference",
                  help="reference ligand filename",)
parser.add_option("-m", "--mapsFolder",
                  action="store", # optional because action defaults to "store"
                  dest="mapsFolder",
                  help="folder containing receptor affinity maps",)
parser.add_option("--swBiasCoeff1",
                  action="store", # optional because action defaults to "store"
                  dest="swBiasCoeff1",
                  type="float",
                  default=0.2,
                  help="SW bias coef 1",)
parser.add_option("--swBiasCoeff2",
                  action="store", # optional because action defaults to "store"
                  dest="swBiasCoeff2",
                  type="float",
                  default=0.4,
                  help="SW bias coef 2",)
(options, args) = parser.parse_args()

from ADFR import ADFRRigidReceptor

syst = os.path.split(options.ligand)[1][:4]
if options.mapsFolder is None:
    mapsFolder ='Astex/fixMapsAstex85varSC/%s/fixMaps/'%syst 
else:
    mapsFolder = args[1]

if options.reference is None:
    reference = 'Astex/ligands/%s_lig.pdbqt'%syst
else:
    reference = options.reference
t00= time()
adfr = ADFRRigidReceptor(options.ligand, mapsFolder)
print 'seed', adfr.initialSeed
pop = adfr.createPopulation(1)
ind = pop[0]
ind.sw.biasCoef1 = options.swBiasCoeff1
ind.sw.biasCoef2 = options.swBiasCoeff2

ind.values[:] = ind.genome.identityGenes()
s0 = ind.score()
rl0 = ind.energies['RRL']
ll0 = ind.energies['LL']

from MolKit2.molecule import getAtomIndicesPerType
from mglutil.math.rmsd import HungarianMatchingRMSD_prody
if reference:
    from MolKit2 import Read
    refLig = Read(reference)   
    atomsRef = refLig.select()
    d1 = getAtomIndicesPerType(atomsRef)
    rmsdCalcRef = HungarianMatchingRMSD_prody(atomsRef.getCoords(), d1, d1)
else:
    rmsdCalcRef = None
    
atomsIn = adfr.ligand.select()
d2 = getAtomIndicesPerType(atomsIn)
rmsdCalcIn = HungarianMatchingRMSD_prody(atomsIn.getCoords(), d2, d2)

from ADFR.GA import GA
ga = GA(pop)
t0 = time()
#sc, nb = ga.minimize(ind, nbSteps=500, noImproveStop=50,
#                     max_steps=10000, MAX_FAIL=100, searchRate=.1)
sc, nb = ga.minimize(ind, nbSteps=100, noImproveStop=20,
                     max_steps=3000, MAX_FAIL=100, searchRate=.1)
rmsdi = rmsdCalcIn.computeRMSD(ind.phenotype[1])
if rmsdCalcRef:
    rmsdx = rmsdCalcRef.computeRMSD(ind.phenotype[1])
else:
    rmsdx = -1
    
print os.path.split(options.ligand)[1][:4],'SCORE', sc-s0, -s0, rl0, ll0, nb, -ind.score(), ind.energies['RRL'], ind.energies['LL'], rmsdi, rmsdx, time()-t0 

from ADFR import writeSolution
writeSolution(ind, 'minimized.pdbqt')
    
## # get angle genes
## agenes = []
## for m in ind.genome.motions[2:]:
##     agenes.extend(m.getGenesForValues([m.getOrigAngle()]))
    
## v0 = [0.06392117353891653, 0.059391144389635234, 0.005844099857817785, 0.40730072356123265, 0.4537822533982281, 0.48742041161038324, 0.5202989077752548, 0.18932262273006206, 0.2058284489880567, 0.06294398525324646, 0.6565562613033274, 0.04688135031824053, 0.04001205895042704, 0.4195974429534351, 0.4072710547225021]
## vr = [0.5, 0., 0., 0.06392117353891653, 0.059391144389635234, 0.005844099857817785, 0.40730072356123265] + agenes
## vradfr = [0.06392117353891653, 0.059391144389635234, 0.005844099857817785, 0.40730072356123265, 0.5, 0., 0.] + [0.]*len(agenes)

## from math import pi, sin, cos, sqrt
## a,b,c,d = [sin(pi*.25), 0, 0, cos(pi*.25)]
## n1 = 1/sqrt(a*a + b*b + c*c + d*d)
## qx90 = [a/n1, b/n1, c/n1, d/n1]

## vr = [0.5, 0., 0.] + [0.06392117353891653, 0.059391144389635234, 0.005844099857817785, 0.40730072356123265] + agenes
## #ind.genome.motions[0].getGenesForValues((24.390,  16.395,  63.165))


## ind.values[:] = vr[:]
## print ind.score()

## writeSolution(ind, 'abcd.pdbqt')

## sc0 = score
## hardmini = {'eCut':0.01, 'maxSteps':100, 'SW_nbSteps':1, 'SW_MAX_FAIL':1}
## nbMaxSteps = 0
## t0 = time()
## while True:
##     sc, nb = ind.minimize(**hardmini)
##     #print sc, nb
##     if nb==hardmini['maxSteps']:
##         nbMaxSteps += 1
##         if nbMaxSteps == 5:
##             break
##     else:
##         nbMaxSteps = 0
##     sc0 = sc

## print -ind._score, time()-t0
