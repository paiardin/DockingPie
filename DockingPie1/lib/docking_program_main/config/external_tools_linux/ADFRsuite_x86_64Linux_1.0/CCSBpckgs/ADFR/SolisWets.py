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
# $Header: /mnt/raid/services/cvs/ADFR/SolisWets.py,v 1.2 2016/12/07 00:38:32 sanner Exp $
#
# $Id: SolisWets.py,v 1.2 2016/12/07 00:38:32 sanner Exp $
#
## NOTE: does not seem to work :(
## the fastest way to minimize is calling search many times with 1 step
## i.e. the bias is never used but it is very fast !
##

# DEBUG
## from mglutil.math.rmsd import RMSDCalculator
#end debug

class SolisWets:

    ## allRmsdList = []
    ## successRmsdList = []
    ## failRmsdList = []
    ## devCnts = None
    ## devSum = None

    def search(self, ind, max_steps=100,
               MAX_SUCCESS=4,
               MAX_FAIL=10,
               FACTOR_EXPANSION=2,
               FACTOR_CONTRACTION=0.5,
               searchRate=0.3):
        success = 0
        fail = 0 
        steps = 0 
        genome = ind.genome
        genome.resetDeltaAmplitude()
        bias = genome.initialBias()
        #score = genome.score( ind.values )
        score = ind._score # we assuem this individual has a score
        best = ind.values[:]
        terminate = False
        # DEBUG
        ## if self.devCnts==None:
        ##     self.devCnts = [0.]* genome.nbVar
        ##     self.devSum = [0.]* genome.nbVar
        ##     self.contract = 0
        ##     self.expand = 0
        ##     self.nbCalls = 0
        ## self.nbCalls += 1
        ## refCoords = genome.getTransformedCoords()
        ## rmsdCalc = RMSDCalculator(refCoords)        
        # END DEBUG
        while (steps < max_steps and not terminate):
            ct = 0
            while ct==0:
                ct, dev = genome.genDeviate(searchRate)
            # DEBUG
            ## for i, v in enumerate(dev):
            ##     self.devSum[i] += v
            ##     if v !=0.0:
            ##         self.devCnts[i]+=1
            # END DEBUG
            biasedDev = genome.biasedDevs(dev, bias)
            nv = genome.applyDelta(best, biasedDev, '+')
            nscore = genome.score( nv )
            # DEBUG
            ## coords = genome.getTransformedCoords()
            ## rmsd = rmsdCalc.computeRMSD(coords)
            ## self.allRmsdList.append(rmsd)
            # END DEBUG
            if nscore > score:
                ## self.successRmsdList.append(rmsd)
                best[:] = nv[:]
                score = nscore
                success += 1 
                fail = 0 
                #bias = genome.newBias(0.4, biasedDev, 0.2, bias)
                bias = genome.newBias(0.2, bias, 0.4, biasedDev)
                
            else:
                ## self.failRmsdList.append(rmsd)
                nv = genome.applyDelta(best, biasedDev, '-')
                nscore = genome.score(nv)
                # DEBUG
                ## coords = genome.getTransformedCoords()
                ## rmsd = rmsdCalc.computeRMSD(coords)
                ## self.allRmsdList.append(rmsd)
                # END DEBUG
                if nscore > score:
                    ## self.successRmsdList.append(rmsd)
                    best[:] = nv[:]
                    score = nscore
                    success += 1 
                    fail = 0 
                    bias = genome.newBias(0.2, bias, -0.4, biasedDev)
                
                else: #  Score still isn't favorable = fail
                    ## self.failRmsdList.append(rmsd)
                    bias = genome.scaleBias(bias, 0.5)
                    fail += 1
                    success = 0 
                
            if(success >= MAX_SUCCESS):
                genome.scaleUpAmplitude(FACTOR_EXPANSION)
                success = 0 
                #self.expand += 1

            # If you have made a X steps in a row that are unfavorable,
            #take a smaller step
            elif(fail >= MAX_FAIL):
                terminate = genome.scaleDownAmplitude(FACTOR_CONTRACTION)
                #self.contract += 1
                if terminate:
                    break 

            steps+=1 

        return best, score, steps
