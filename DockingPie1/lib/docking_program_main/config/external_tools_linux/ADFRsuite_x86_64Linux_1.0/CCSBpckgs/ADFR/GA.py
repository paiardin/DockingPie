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

############################################################################
#
# Author: Michel F. SANNER
#
# Copyright: M. Sanner and TSRI 2015
#
#########################################################################
#
# $Header: /mnt/raid/services/cvs/ADFR/GA.py,v 1.14 2017/04/14 17:50:25 sanner Exp $
#
# $Id: GA.py,v 1.14 2017/04/14 17:50:25 sanner Exp $
#
import os, sys, numpy
from random import random, gauss

from LigandFT import MotionBase
from ADFR.ReceptorFT import FTSoftRotamer

## class Gene:
##     """object representing 1 or more variables"""

##     def __init__(self, motionObject):
##         assert isinstance(motion, MotionBase)
##         self.motion = motionObject
##         self.nbVar = motionObject.getNumVariables()
##         self.values = [0]*self.nbVar # from 0 to 1
        
##     def setValues(self, values):
##         self.values[:] = values

from ADFRcc.adfr import Genome, Individual

class GenomePy(Genome):
    """list of motion objects"""

    def __init__(self, ft, scorer, scaleRE=None, neighborSearchCutoff=-1.0):
        super(GenomePy, self).__init__(ft.allMotions, ft.ftRoot, scorer)
        self.ft = ft
        self.scaleRE = scaleRE

        #self.ftRoot = ft.ftRoot
        self.scorer = scorer
        #self.motions = ft.motions
        self.nbVar = 0
        self.softRotMotions = []

        # negative values disable neighborsearch, a positive value is the RMSD cutoff
        self.neighborSearchCutoff=neighborSearchCutoff
        
        cutPoints = []
        offset = 0 # offset into genes values for each motion object in the genome
        for m in ft.motions:
            m._offset = offset
            n = m.getNumVariables()
            self.nbVar += n
            offset += n
            cutPoints.append(self.nbVar)
            if isinstance(m, FTSoftRotamer):
                self.softRotMotions.append(m)
                
        if len(cutPoints) > 2:
            cutPoints = cutPoints[:-1]
        self.cutPoints = cutPoints
        print 'CUTPOINTS', self.cutPoints

    def getIdentityGenesPy(self):
        values = []
        for m in self.ft.motions:
            values.extend(m.getIdentityGenesPy())
        return values
    ## def identityGenes(self):
    ##     values = []
    ##     for m in self.motions:
    ##         values.extend(m.idGenes)
    ##     return values
    
    ## def randomValues(self):
    ##     values = []
    ##     for m in self.motions:
    ##         m.randomize()
    ##         #print 'R', m.__class__.__name__, m.getVariablesPy()
    ##         values.extend( m.getGenesForValues(m.getVariablesPy()) )
    ##     return values
    
    ## def initialValues(self):
    ##     values = []
    ##     for m in self.motions:
    ##         m.initialize()
    ##         #print 'I', m.__class__.__name__, m.getVariablesPy()
    ##         values.extend( m.getGenesForValues(m.getVariablesPy()))
    ##     return values

    ##
    ## for Solis Wets
    ##
    ## def genDeviate(self, searchRate):
    ##     """
    ##     return a list of deltas for each motion
    ##     """
    ##     deltas = []
    ##     ct = 0
    ##     for m in self.motions:
    ##         ctm, devs = m.genDeviate(searchRate)
    ##         ct += ctm
    ##         deltas.extend( devs )
    ##     return ct, deltas
        
    ## def initialBias(self):
    ##     bias = []
    ##     for m in self.motions:
    ##         bias.extend( m.initialBias())
    ##     return bias
        
    ## def biasedDevs(self, dev, bias):
    ##     """
    ##     let each motion object apply the bias to the deviation
    ##     """
    ##     off = 0
    ##     biasedDev = []
    ##     for m in self.motions:
    ##         nb = m.getNumVariables()
    ##         biasedDev.extend( m.biasedDevs(dev[off:off+nb], bias[off:off+nb]))
    ##         off += nb
    ##     return biasedDev
        
    ## def applyDelta(self, values, delta, direction='+'):
    ##     """
    ##     return a list of deltas for each motion
    ##     """
    ##     off = 0
    ##     nv = []
    ##     ## backwards compatible with ADFR
    ##     #m = self.motions[1]
    ##     #nvr = m.applyDelta(values[3:7], delta[0:4], direction)
    ##     #nv = m.applyDelta(values[0:3], delta[4:7], direction)
    ##     #nv.extend(nvr)
    ##     #off = 7
    ##     #for m in self.motions[2:]:
    ##     for m in self.motions:
    ##         nb = m.getNumVariables()
    ##         nv.extend(m.applyDelta(values[off:off+nb], delta[off:off+nb],
    ##                                direction))
    ##         off += nb
    ##     return nv

    ## def scaleDeltaAmplitudeBase(self, factor):
    ##     print 'SCALEDELTA',
    ##     for m in self.motions:
    ##         m.scaleDeltaAmplitudeBase(factor)
    ##         print m._deltaAmplitude,
    ##     print

    ## def scaleBias(self, bias, factor):
    ##     off = 0
    ##     newBias = []
    ##     for m in self.motions:
    ##         nb = m.getNumVariables()
    ##         newBias.extend(m.scaleBias(bias[off:off+nb], factor))
    ##         off += nb
    ##     return newBias
        
    ## def resetDeltaAmplitude(self):
    ##     off = 0
    ##     for m in self.motions:
    ##         m.resetDeltaAmplitude()
        
    ## def newBias(self, coef1, bias, coef2, dev):
    ##     off = 0
    ##     newBias = []
    ##     for m in self.motions:
    ##         nb = m.getNumVariables()
    ##         newBias.extend(m.newBias(coef1, bias[off:off+nb],
    ##                                  coef2, dev[off:off+nb]))
    ##         off += nb
    ##     return newBias

    ## def scaleUpAmplitude(self, factor):
    ##     for m in self.motions:
    ##         m.scaleUpAmplitude(factor)
            
    ## def scaleDownAmplitude(self, factor):
    ##     for m in self.motions:
    ##         terminate = m.scaleDownAmplitude(factor)
    ##         if terminate is True:
    ##             return True
    ##     return False
    
    ##
    ## END for Solis Wets
    ##
    def getLigandTransformedCoords(self):
        return self.ft.ligandAtomSet.getCoordsPy()

    def getReceptorTransformedCoords(self):
        return self.ft.receptorAtomSet.getCoordsPy()
    
    ## def score(self, values):
    ##     off = 0
    ##     for m in self.motions:
    ##         nb = m.getNumVariables()
    ##         m.setValuesFromGenes(values[off:off+nb])
    ##         off += nb
    ##     self.ftRoot.update()
    ##     return -self.scorer.calculateScores()

    ## def boundGenes(self, values):
    ##     off = 0
    ##     bvalues = []
    ##     for m in self.motions:
    ##         nbg = m.getNumVariables()
    ##         for i, value in enumerate(values[off:off+nbg]):
    ##             isCyclic = m.isCyclic[i]
    ##             if value > 1.0:
    ##                 if isCyclic:
    ##                     value =  (value - 1.)%1.
    ##                 else:
    ##                     value = 1.
    ##             elif value < 0.:
    ##                 if isCyclic:
    ##                     value = 1. - (0. - value)%1.
    ##                 else:
    ##                     value = 0.
    ##             bvalues.append(value)
    ##         off += nbg
    ##     return bvalues
    
#from ADFR.SolisWets import SolisWets
from time import time
class IndividualPy(Individual):

    #sw = SolisWets()
    
    def __init__(self, genomePy, values=None):

        super(IndividualPy, self).__init__(genomePy, values)
        assert isinstance(genomePy, GenomePy)
        self.genomePy = genomePy

        ## if values is not None:
        ##     assert len(values) == self.genome.nbVar
        ##     self.values = values
        ## else:
        ##     self.values = self.genome.initialValues()
        self.phenotype = [None, None]
        self.energies = {'LL':0., 'RRL':0., 'FRL':0., 'FRFR':0., 'RRFR':0., 'wRR':0.}
        self._score = None # -total energy
        self._fitness = None
        self._neighborRMSD = None
        
    ## def setGenes(self, values):
    ##     assert len(values)==len(self.values)
    ##     self.values[:] = self.genome.boundGenes(values)
            
    def clone(self):
        ind = self.__class__(self.genomePy)
        #import pdb; pdb.set_trace()
        ind.setGenes(self.getGenesPy())
        ind.energies = self.energies.copy()
        ind._score = self._score
        ind._neighborRMSD = self._neighborRMSD
        ind._fitness = self._fitness
        ind.phenotype = [None, None]
        if self.phenotype[0] is not None:
            ind.phenotype[0] = self.phenotype[0].copy()
        if self.phenotype[1] is not None:
            ind.phenotype[1] = self.phenotype[1].copy()
        return ind
    
    def score(self):
        genome = self.genomePy
        #genome.setGenes(self.getGenesPy())
        self._score = self.calculateScores() # positive number when good (i.e -energy)
        # we do not use this as this energy cannot scale down receptor internal energy yet
        if genome.neighborSearchCutoff>0.:
            self._neighborRMSD = genome.scorer.getNeighborRMSD()
        self.phenotype[1] = self.genomePy.getLigandTransformedCoords()
        internalLigandE = self.energies['LL'] = genome.scorer.getLlPairwiseScorer().getTotalScore()
        interactionE = self.energies['RRL'] = genome.scorer.getLrrGridScorer().getTotalScore()
        if genome.ft.receptorAtomSet:
            self.phenotype[0] = self.genomePy.getReceptorTransformedCoords()
            a = self.energies['FRFR'] = genome.scorer.getFrfrPairwiseScorer().getTotalScore()
            b = self.energies['RRFR'] = genome.scorer.getFrrrGridScorer().getTotalScore()
            self.energies['FRL'] = genome.scorer.getLfrPairwiseScorer().getTotalScore()
            self.energies['wRR'] = a*genome.scorer.getFrfrCoefficient() + b*genome.scorer.getFrrrCoefficient()
            #self._score = -(internalLigandE + interactionE + self.energies['FRL'] + self.energies['wRR'])
        #else:
        #    self._score = -(internalLigandE + interactionE)
        #print 'ENERGIES', self.energies
        return self._score

    ## def randomize(self):
    ##     self.values[:] = self.genome.randomValues()

    ## def initialize(self):
    ##     self.values[:] = self.genome.initalValues()

    ## def NEWBADminimize(self, eCut=0.1, maxSteps=100, SW_nbSteps=1, SW_MAX_FAIL=3):
    ##     t0 = time()
    ##     totalSteps = 0
    ##     i = 0
    ##     end = False
    ##     oldScore = origScore = self._score
    ##     while not end and i<maxSteps:
    ##         #best, score, steps = self.sw.search(self, nbSteps=100, MAX_FAIL=3)
    ##         best, score, steps = self.sw.search(self, nbSteps=SW_nbSteps, MAX_FAIL=SW_MAX_FAIL)
    ##         totalSteps += steps
    ##         if score > oldScore:
    ##             self.values[:] = best[:]
    ##             if score-oldScore < eCut:
    ##                 end = True
    ##             oldScore = score
    ##         i += 1
    ##     self.score()
    ##     #if self._score > -100 and abs(self._score+(self.energies['RRL']+self.energies['LL'])) > 0.1:
    ##     #    import pdb; pdb.set_trace()
    ##     #print 'MINI %f -> %f in %d %d steps %f(s)'%(origScore, score, i, totalSteps, time()-t0)
    ##     return score, i
        
def score_maximize(x,y): 
	"""Maximization comparator for raw score."""
	#return cmp(y.score(),x.score())
	#return cmp(y.evaluate(),x.evaluate())
	return cmp(y._score, x._score)

from ADFR.scaling import sigma_truncation_scaling
from ADFR.selection import srs_selector

class Population(list):

    def __init__(self, values=[]):
        self.extend(values)
        self.scaler = sigma_truncation_scaling()
        self.selector = srs_selector()
        self.scaled = False
        self.select_ready = False
        self.sorted = False

    def touch(self):
        self.scaled = False
        self.sorted = False
        self.select_ready = False
        
    def sort(self):
        list.sort(self, score_maximize)
        self.sorted = True

    def scores(self):
        # return list of scores
        return [x._score for x in self]

    def fitnesses(self):
        # return list of fitness values
        return [x._fitness for x in self]

    def scale(self, force=0):
        if not self.scaled or force:
            self.scaler.scale(self)			
        self.scaled = True

    def select(self, cnt=1):
        """Calls the selector and returns *cnt* individuals.
        Arguments:
        cnt -- The number of individuals to return.
        """
        if not self.select_ready:
            self.selector.update(self)
            self.select_ready = True
        if len(self.selector.choices):
            return self.selector.select(self,cnt)
        else:
            return None, None

    def clone(self):
        cind = [x.clone() for x in self]
        return self.__class__(values=cind)

from random import randint, random
from MolKit2.molecule import getAtomIndicesPerType
from mglutil.math.rmsd import HungarianMatchingRMSD_prody, RMSDCalculator
from .ga_util import my_mean, my_std, flip_coin

def mutate(ind, mutationRate):
    ## self.mutation_rate is GA_mutation  parameter
    mutated = 0
    values = ind.getGenesPy()
    offset = 0

    #mutatedSD = []
    mutation_rate = 0.05
    nb = max(0, int(gauss(2, 0.5)))
    genome = ind.genomePy
    indices = [int(random()*len(genome.softRotMotions)) for i in range(nb)]
    status = None
    for i, motion in enumerate(genome.softRotMotions):
        nbg = motion.getNumVariables()
        if not i in indices:
            continue
        offset = motion._offset
        # randomly pick a side chain and mutate it
        lmutated = motion.mutate(values[offset:offset+nbg], mutation_rate)
        lmutated, status = lmutated
        mutated += lmutated
        offset += nbg
        #if status > -1:
        #    self.softRotXtal[i] = status
            #mutatedSD.append( (i,status) )

    for motion in genome.ft.motions:
        if isinstance(motion, FTSoftRotamer):
            continue
        nbg = motion.getNumVariables()
        offset = motion._offset
        lmutated = motion.mutate(values[offset:offset+nbg], mutationRate)
        mutated += lmutated
        offset += nbg
    #print 'mutating INDIVIDUAL ====', mutated, mutatedSD, self.softRotXtal, id(self)
    if mutated > 0:
        ind._score = None
        ind._fitness = None
    return mutated, status

from ADFRcc.adfr import SolisWets
from ADFR import _intelContest_

class GA:

    def setDebug(self, level):
        self.debug = level
        
    def __init__(self, pop, ligand, reflig=None, folder='.',
                 RMSDMatching='hungarian', neighborSearchCutoff=-1.0):
        self.debbug = False
        self.gen = 0
        self.pop = pop
        self.sw = SolisWets.getSingleton()

        self.folder = folder
        # declared here, set by ga.evolve()
        self.MaxGens = None
        self.p_stdev = None
        self.p_crossover = None
        self.p_mutation = None
        self.maxEvals = None
        self.useClustering = True
        self.clusteringCutoff = 2.0
        self.neighborSearchCutoff = neighborSearchCutoff
        self.bestNeighborScore = None
        self.bestNeighborInd = None
        self.mutate = mutate

        #quickmini = {'eCut':0.1, 'maxSteps':100, 'SW_nbSteps':5} # good
        #self.quickmini = {'eCut':0.1, 'maxSteps':20, 'SW_nbSteps':3}
        #self.hardmini = {'eCut':0.01, 'maxSteps':200, 'SW_nbSteps':10}
        self.quickmini = {'nbSteps':1, 'noImproveStop':1, 'max_steps':100, 'MAX_FAIL':3, 'searchRate':.3}
        #self.quickmini = {}
        self.hardmini = {'nbSteps':10, 'noImproveStop':3, 'max_steps':300, 'MAX_FAIL':4, 'searchRate':0.05}
        
        # rmsd calculator used by GA to cluster population
        atoms = ligand.select()
        if RMSDMatching=='hungarian':
            d1 = getAtomIndicesPerType(atoms)
            self.rmsdCalc = HungarianMatchingRMSD_prody(atoms.getCoords(), d1, d1)
        else:
            self.rmsdCalc = RMSDCalculator(atoms.getCoords())

        # other RMSD calcualtors
        self.rmsdCalculators = []
        if reflig:
            atoms = reflig.select()
            if RMSDMatching=='hungarian':
                d1 = getAtomIndicesPerType(atoms)
                rmsdCalc = HungarianMatchingRMSD_prody(atoms.getCoords(), d1, d1)
            else:
                rmsdCalc = RMSDCalculator(atoms.getCoords())                
            self.rmsdCalculators.append(rmsdCalc) # list of RMSD calculators

        self.rmsdCalculatorsRec = [] # list of RMSD calculators

        # clustering Energy cutoff
        self.clusterEnergyCut = 2.0

        from ADFRcc.adfr import FlexibleReceptorScorer
        #import pdb; pdb.set_trace()
        if isinstance(pop[0].genomePy.scorer, FlexibleReceptorScorer):
            self.rigidReceptorDocking = False
            self.genStatsOutput = self.genStatsOutputFR
            self.crossover = self.rigidAndflex_crossover
        else:
            self.rigidReceptorDocking = True
            self.genStatsOutput = self.genStatsOutputRR
            self.crossover = self.singlepoint_crossover

    def SWstats(self):
        print 'SolisWets: calls: %d with gain: %d without gain:%d eGain: %f eGain/call: %f'%(
            self.nbSWcalls, self.nbSWcallsWithGain, self.nbSWcallsWithoutGain,
            self.totalEgain, self.totalEgain/ self.nbSWcallsWithGain)
        egain = numpy.sum(self.MinimizeEGains[:self.gen])
        print 'Minimizer: calls: %d eGain: %f eGain/call: %f'%(
            self.nbMinimize, egain, egain/self.nbMinimize)
                                       
    def minimize(self, individual, nbSteps=5, noImproveStop=2, max_steps=20, MAX_FAIL=3, searchRate=0.3):
        if _intelContest_ and self._nbWritten< len(self.pop):
            [self._f.write("%.3f "%x) for x in individual.getGenesPy()]
            self._f.write("\n")
            self._nbWritten += 1
        
        self.sw.minimize(individual, nbSteps, noImproveStop, max_steps, MAX_FAIL, searchRate)
        return individual.score(), self.sw.getSteps()
    
        ## #minimize_param = self.configure_minimize
        ## # MS sept 3 2014 ,, I think the _fitness_score is already known
        ## ##last_score = individual.score()
        ## s0 = -individual._score
        ## last_score = individual._score
        ## noImprovement = 0
        ## totalSteps = 0
        ## best = individual.values[:]
        ## for i in range(nbSteps):
        ##     #neighbor, nbSteps = solisWets.search(individual, **kw)
        ##     #totalSteps += nbSteps
        ##     miniVals, score, steps = individual.sw.search(individual, **kw)
        ##     totalSteps += steps

        ##     #self.nbSWcalls += 1
        ##     #self.nbSWsteps += steps
        ##     #self.SWNSteps[self.nbSWcallsWithGain] = steps
        ##     if score > last_score:
        ##         #print 'round %d: %f -> %f %4d %f'%(i, -individual._score, -neighbor._score, nbSteps, kw['MIN_VAR'])
        ##         #eGain = score-last_score
        ##         #self.SWEGains[self.nbSWcallsWithGain] = eGain
        ##         #self.totalEgain += score-last_score
        ##         #self.nbSWcallsWithGain +=1

        ##         last_score = score
        ##         best[:] = miniVals[:]
        ##         #individual.setValues(best)
        ##     else:
        ##         #self.nbSWcallsWithoutGain +=1

        ##         noImprovement += 1
        ##         if noImprovement > noImproveStop:
        ##             break
        ## individual.values[:] = best[:]
        ## individual.score()
        ## individual._totalStepInLastMinimize = totalSteps
        ## #self.MinimizeEGains[self.nbMinimize] = score-last_score
        ## #self.nbMinimize += 1
        ## #print 'mini', i, totalSteps, s0, -last_score, kw
        ## return individual._score, totalSteps

    def minimizePop(self, pop, **kw):
        t0 = time()
        mini = -999999999
        for individual in pop:
            self.minimize(individual, **kw)
            if individual._score>mini:
                mini = individual._score
        print 'Minimized pop of %d individuals in %.2f second with miniE: %f'%(len(pop), time()-t0, -mini)

    def pop_deviation(self):
        #compute the coefficient of variation (CV): STDV/mean
        scores = self.pop.scores()
        denom = my_mean(scores)
        if denom == 0.: denom = .0001  # what should I do here?
        return abs(my_std(scores) / denom)

    def iteration_output(self):
        # Return the (-) of the score.  Indicating more negative, now more favorable
        score = -self.pop[0]._score

        output = ('\ngen: ' + `self.gen` + ' ' 
                  + 'min: ' + str(score)  + ' ' 
                  + 'pop: %d'%len(self.pop) + ' ' 
                  + 'dev: %.4e'%self.p_dev + ' ' 
                  + 'seed: %d '%self.seed + ' ' 
                  + '#evals: %7d '%(self.pop[0].genomePy.scorer.getNumEvals()) + ' ' 
                  + 'Total Time:  %7.2f'%(time()-self.beginTime))
        print( output )

    def genStatsOutputRR(self):
        """
        None <- GenStatsOutput()

        This function will take a genome (collection of individual
        genes that make up the population) and determine the rmsd and 
        score of each gene. Used to track the minimum RMSD of the 
        population after each step of the GA.
        """
        # Handle to GA instance

        # loop over all RMSD calculators
        if len(self.rmsdCalculators) == 0 or self.debug==0: # simple output, no RMSD
            bestInd = self.pop[0]
            bie = bestInd.energies
            print " _Gen%04d Score: %7.3f LL: %7.3f LR: %7.3f evals: %7d   %2d"%(
                self.gen, -bestInd._score, bie['LL'], bie['RRL'],
                self.pop[0].genomePy.scorer.getNumEvals(),
                self.dStarcnt)
            return # output of iteration has the info
        
        print "           RMSD   score    LL      RRL     #  RMSD   score     LL     RRL    evals  ndStar"
        for rmsdcl in self.rmsdCalculators:
            rmsd_lst = []
            for i, gene in enumerate(self.pop):
                FR_coords, L_coords = gene.phenotype
                rmsd = rmsdcl.computeRMSD(L_coords)
                rmsd_lst.append(rmsd)
            minRMSD = min(rmsd_lst)
            minInd = rmsd_lst.index(minRMSD)
            maxRMSD = max(rmsd_lst)
            maxInd = rmsd_lst.index(maxRMSD)
            
            indMin = self.pop[minInd] # min rmsd individual
            bestInd = self.pop[0] # best score individual
            mie = indMin.energies
            bie = bestInd.energies
            print " _Gen%04d %5.2f %7.3f %7.3f %7.3f |  %d %5.2f %7.3f %7.3f %7.3f %7d   %2d"%(
                self.gen, rmsd_lst[0], -bestInd._score, bie['LL'], bie['RRL'],
                minInd, minRMSD, -indMin._score, mie['LL'], mie['RRL'],
                self.pop[0].genomePy.scorer.getNumEvals(),
                self.dStarcnt)
        
        return rmsd_lst

    def genStatsOutputFR(self):
        """
        None <- GenStatsOutput()

        This function will take a genome (collection of individual
        genes that make up the population) and determine the rmsd and 
        score of each gene. Used to track the minimum RMSD of the 
        population after each step of the GA.
        """
        # Handle to GA instance
        if len(self.rmsdCalculators) == 0 or self.debug==0: # simple output, no RMSD
            bestInd = self.pop[0]
            bie = bestInd.energies
            print " _Gen%04d Score: %7.3f LL: %7.3f LR: %7.3f RR: %7.3f evals: %7d   %2d"%(
                self.gen, -bestInd._score, bie['LL'], bie['RRL']+ bie['FRL'],
                bie['FRFR']+bie['RRFR'],
                self.pop[0].genomePy.scorer.getNumEvals(),
                self.dStarcnt)
            return # output of iteration has the info

        # loop over all RMSD calculators
        print "           RMSD   score    LL      LR      RR       #  RMSD   score     LL     RL     RR    evals  ndStar"

        for rmsdcl in self.rmsdCalculators:
            rmsd_lst = []
            for i, gene in enumerate(self.pop):
                FR_coords, L_coords = gene.phenotype
                rmsd = rmsdcl.computeRMSD(L_coords)
                rmsd_lst.append(rmsd)
            minRMSD = min(rmsd_lst)
            minInd = rmsd_lst.index(minRMSD)
            maxRMSD = max(rmsd_lst)
            maxInd = rmsd_lst.index(maxRMSD)
            
            indMin = self.pop[minInd] # min rmsd individual
            bestInd = self.pop[0] # best score individual
            mie = indMin.energies
            bie = bestInd.energies
            print " _Gen%04d %5.2f %7.3f %7.3f %7.3f %7.3f |  %d %5.2f %7.3f %7.3f %7.3f %7.3f %7d   %2d"%(
                self.gen, rmsd_lst[0], -bestInd._score, bie['LL'], bie['RRL']+ bie['FRL'],bie['FRFR']+bie['RRFR'],
                minInd, minRMSD, -indMin._score, mie['LL'], mie['RRL']+mie['FRL'], mie['FRFR']+mie['RRFR'],
                self.pop[0].genomePy.scorer.getNumEvals(),
                self.dStarcnt)

        #import pdb; pdb.set_trace()
            
        return rmsd_lst
        ## if len(self.rmsdCalculators):
        ##     rmsdcl = self.rmsdCalculators[0]
        ##     rmsd_lst = []
        ##     for i, gene in enumerate(self.pop):
        ##         FR_coords, L_coords = gene.phenotype
        ##         ##print "\n",FR_coords,"\n","\n"
        ##         rmsd = rmsdcl.computeRMSD(L_coords)
        ##         rmsd_lst.append(rmsd)
        ##     minRMSD = min(rmsd_lst)
        ##     minInd = rmsd_lst.index(minRMSD)
        ##     maxRMSD = max(rmsd_lst)
        ##     maxInd = rmsd_lst.index(maxRMSD)
        ## else:
        ##     rmsdcl = None
        ##     minRMSD = -1
        ##     minInd = -1
        ##     maxRMSD = -1
        ##     maxInd = -1
        ##     rmsd_lst = [-1]
            
        ## if len(self.rmsdCalculatorsRec):
        ##     rmsdcr = self.rmsdCalculatorsRec[0]
        ##     rmsd_lstr = []
        ##     for i, gene in enumerate(self.pop):
        ##         FR_coords, L_coords = gene.phenotype
        ##         rmsd = rmsdcr.computeRMSD(FR_coords)
        ##         rmsd_lstr.append(rmsd)
        ##     minRMSDR = min(rmsd_lstr)
        ##     minIndR = rmsd_lstr.index(minRMSDR)
        ##     maxRMSDR = max(rmsd_lstr)
        ##     maxIndR = rmsd_lstr.index(maxRMSDR)
        ## else:
        ##     rmsdcr = None
        ##     minRMSDR = -1
        ##     minIndR = -1
        ##     maxRMSDR = -1
        ##     maxIndR = -1
        ##     rmsd_lstr = [-1]
            
        ## indMin = self.pop[minInd]
        ## bestInd = self.pop[0]
        ## mie = indMin.energies
        ## bie = bestInd.energies
        ## print " _Gen%03d  %5.2f %5.2f   %3d   ( %9.3f %9.3f %9.3f %9.3f ) |  %5.2f %5.2f  |      %5.2f %5.2f ( %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f  %9.3f %9.3f ) | %d  | %d"%(
        ##     self.gen, minRMSD, minRMSDR, minInd,
        ##     mie['RRL']+mie['FRL'], mie['FRFR']+mie['RRFR'], mie['LL'], mie['RRL']+mie['FRL'],
        ##     maxRMSD, maxRMSDR, rmsd_lst[0], rmsd_lstr[0],
        ##     bie['RRL'], bie['FRL'], bie['FRFR'], bie['RRFR'], bie['wRR'], bie['LL'],
        ##     -bestInd._score, bie['RRL']+bie['FRL'],self.pop[0].genome.scorer.getNumEvals(),
        ##     self.dStarcnt)
        ## #print "BestInd Rec-Rec Unbound : ",bie['RRUb, "weighted Rec-Rec Energy : ",bie['wRR

        ## for rmsdc in self.rmsdCalculators[1:]:
        ##     pass
        ## for rmsdc in self.rmsdCalculatorsRec[1:]:
        ##     pass
        
        ## return rmsd_lst

    def evolve(self, maxGens=1000, p_crossover=0.8, p_mutation=0.05,
               maxEvals=1500000, popStdev=1e-4, noImproveStop=5):

        self.MaxGens = maxGens
        self.p_crossover = p_crossover
        self.p_mutation = p_mutation
        self.maxEvals = maxEvals
        self.popStdev = popStdev
        self.stopNoImprove = noImproveStop
        self.nbScaleAmpBase = 0
        print 'quick', self.quickmini
        print 'hard', self.hardmini
        print 'stopNoImprov', self.stopNoImprove
        print 'deltas',
        for m in self.pop[0].genomePy.ft.motions:
            print m.getDeltaAmplitudeBase(),
        print
        ##
        ## Search Stats
        self.nbSWcalls = 0
        self.nbSWcallsWithGain = 0
        self.nbSWcallsWithoutGain = 0
        self.nbSWsteps = 0 # total number of steps performed by SW
        self.totalEgain = 0 # total number of steps performed by SW
        self.nbMinimize = 0
        self.SWEGains = numpy.zeros( (maxEvals,) , 'f')
        self.SWNSteps = numpy.zeros( (maxEvals,) , 'f')
        self.MinimizeEGains = numpy.zeros( (maxEvals,) , 'f')
        self.bestScore = numpy.zeros( (maxEvals,) , 'f')

        self.beginTime = time()
        self.nbClusters = -1
        self.nbGenWithCstClusters = 0
        self.clustersBestE = [] # list of best energy in cluster in previous
                                # generation
        self.nbNoClusterEimprovement = 0
        self.dStarcnt = 0
        self.starMinData = []

        self.p_dev = self.pop_deviation()
        status = 'searching'
        self.gen = 1
        self.pop.sort()
        self.pop.scale()
        self.iteration_output()

        if _intelContest_:
            self._f = open('states.txt', 'w')
            self._nbWritten = 0
            self._f.write('POPSIZE %d\n'%len(self.pop))
            for ind in self.pop:
                [self._f.write("%.3f "%x) for x in ind.getGenesPy()]
                self._f.write("\n")
                
        while (self.gen <= self.MaxGens and self.pop[0].genomePy.scorer.getNumEvals()<=self.maxEvals) \
                  and status=='searching':

            self.lastGenBest = self.pop[0]._score
            self.bestScore[self.gen-1] =self.lastGenBest

            if self.neighborSearchCutoff>0:                    
                for ind in self.pop:
                    if ind._neighborRMSD < self.neighborSearchCutoff:
                        if self.bestNeighborScore is None or ind._score > self.bestNeighborScore:
                            self.bestNeighborScore=ind._score
                            self.bestNeighborInd=ind.clone()
            
            # GA step where crossover, mutation, replacement, local search can occur
            status = self.step()

            # Population Coefficient of variation (CV): SD/mean
            self.p_dev = self.pop_deviation()

            # Write out the generation information
            self.iteration_output()

            # print out stats for this Gen
            rmsd_lst = self.genStatsOutput()

            self.gen += 1
            
            if self.p_dev<self.p_stdev:
                status='population deviation is %f'%p_stdev
                break
            sys.stdout.flush()
        # end of evolution
        if _intelContest_:
            self._f.close()

        #atoms.setCoords(topSol.phenotype[1])
        #import prody
        #prody.writePDB('sol0.pdb', atoms) 
        return status
            
    def cluster(self, remainder):
        clusters = []
        seedInd = remainder[0]
        pop = self.pop
        RMSDcalc = self.rmsdCalc
        while len(remainder):
            #print '%d left to cluster seed=%d'%(len(remainder), seedInd), remainder
            ref = pop[seedInd]
            dum2, L_coordsRef = ref.phenotype
            RMSDcalc.setRefCoords(L_coordsRef)
            cluster = [seedInd]
            bestScore = min([pop[i]._score for i in remainder])
            notSelected = []
            seed = seedInd # seed will not change in the for loop but seedInd will
            for i in remainder:
                if i==seed:
                    continue
                FR_coords, L_coords = pop[i].phenotype
                # RMSD calc
                rmsd = RMSDcalc.computeRMSD(L_coords)
                if rmsd<=self.clusteringCutoff:
                    cluster.append(i)
                else:
                    notSelected.append(i)
                    if pop[i]._score >= bestScore:
                        seedInd = i
                        bestScore = pop[i]._score
            #print 'found cluster', cluster, seedInd
            clusters.append( cluster )
            remainder = notSelected
        return clusters

    def singlepoint_crossover(self, mom, dad):
        if len(mom.genomePy.cutPoints) > 1:
            crosspoint = mom.genomePy.cutPoints[randint(0,len(mom.genomePy.cutPoints)-1)]
        else: 
            crosspoint = mom.genomePy.cutPoints[0]

        klass = mom.__class__
        dadGenes = dad.getGenesPy().tolist()
        momGenes = mom.getGenesPy().tolist()
        brother = klass(mom.genomePy)
        brother.setGenes(momGenes[:crosspoint] + dadGenes[crosspoint:])
        sister = klass(mom.genomePy)
        sister.setGenes(dadGenes[:crosspoint] + momGenes[crosspoint:])
        return brother, sister

    def rigidAndflex_crossover(self, mom, dad):
        """
        For rigid docking a single point crossover is done, and for a flexible docking a
        double point crossover is done (one crosspoint picked randomly from the cutpoints of
        flexible sidechains (genes), and the other crosspoint is picked randomly from the cutpoints
        of ligand genes)
        """
        nbSoftRot = mom.genomePy.ft.nbSoftRotamers
        cutPoints = mom.genomePy.cutPoints
        if nbSoftRot==0: # rigid docking
            crosspoint = cutPoints[randint(0,len(cutPoints)-1)]
        else:
            crosspoint1 = cutPoints[randint(0,len(cutPoints)-nbSoftRot-1)]
            crosspoint2 = cutPoints[randint(len(cutPoints)-nbSoftRot,len(cutPoints)-1)]

	#print nbSoftRot, crosspoint1,crosspoint2
        klass = mom.__class__
        dadGenes = dad.getGenesPy().tolist()
        momGenes = mom.getGenesPy().tolist()
        if nbSoftRot == 0: # single point cross over
            brother = klass(mom.genomePy)
            brother.setGenes(momGenes[:crosspoint] + dadGenes[crosspoint:])
            sister = klass(mom.genomePy)
            sister.setGenes(dadGenes[:crosspoint] + momGenes[crosspoint:])
        else:
            brother = klass(mom.genomePy)
            brother.setGenes(momGenes[:crosspoint1] + dadGenes[crosspoint1:crosspoint2] + momGenes[crosspoint2:])
            sister = klass(mom.genomePy)
            sister.setGenes(dadGenes[:crosspoint1] + momGenes[crosspoint1:crosspoint2] + dadGenes[crosspoint2:])
        return brother, sister

    ## def writeSolution(self, solution, filename):
    ##     # hack to replace the coordinates in the pdbqt file
    ##     mol = solution.genomePy.ft.mol
    ##     atoms = mol._ag
    ##     score = solution.score() # update tree and coordinates
    ##     f = open(mol.filename)
    ##     lines = f.readlines()
    ##     f.close()

    ##     atomByName = {}
    ##     coords = solution.phenotype[1]
    ##     for i, atom in enumerate(atoms):
    ##         name = atom.getName()
    ##         assert not atomByName.has_key(name), "ERROR atoms with identical names %"%name
    ##         atomByName[atom.getName()] = coords[i]

    ##     ct = 0
    ##     for i in range(len(lines)):
    ##         if lines[i].startswith('HETATM') or lines[i].startswith('ATOM'):
    ##             atName = lines[i][12:16].strip()
    ##             x,y,z = atomByName[atName]
    ##             lines[i] = lines[i][:30] + '%8.3f%8.3f%8.3f'%(x,y,z) + lines[i][54:]

    ##     f = open(filename, 'w')
    ##     f.write("USER: ADFR SOLUTION 0\n")
    ##     f.write("USER: SCORE %f\n"%(0.-score))
    ##     #if len(self.rmsdCalculators):
    ##     #    f.write("USER: RMSD %f\n"%self.rmsdCalculators[0].computeRMSD(solution.phenotype[1]))
    ##     f.write("USER: genes %s\n"%str(solution.getGenesPy()))
    ##     [f.write(l) for l in lines]
    ##     f.close()

    def loadPopFromADFRGenes(self, pop, filename):
        d = {}
        execfile(filename, d, d)
        genes = d['pop']
        assert len(genes) == len(pop)
        assert len(genes[0][0]) == pop[0].genomePy.nbVar
        motions = pop[0].genomePy.motions
        gv = [0.]*pop[0].genomePy.nbVar
        for ind, g in zip(pop, genes):
            adfrg, rmsd, ADFRscore = g
            gv[3:7] = adfrg[0:4] # rotation
            gv[0:3] = adfrg[4:7] # translation
            for n in range(7, pop[0].genomePy.nbVar):
                m = motions[n-5]
                gv[n] = m.idGenes[0]+adfrg[n]
            ind.setGenes(gv)
            ns = ind.score()
            # overwrite score with ADFR score
            ind._score = -ADFRscore
            #self.writeSolution(ind, 'adfrSol0.pdbqt')
            #import pdb; pdb.set_trace()
            #assert ns == ADFRscore
        
    def step(self):
        # count similar individuals (within clusterEcut of the top scored ind)

        #self.loadPopFromADFRGenes(self.pop, '../ADFR_adfr/ADFR/1l7f/1l7f_gen2.py')
        #import pdb; pdb.set_trace()

        if _intelContest_:
            self._nbWritten = 0

        sanityCheck = False
        self._nbMut = 0
        ref = self.pop[0]
        refScore = ref._score # lowest score in pop al all times
        remainder = []
        for i,x in enumerate(self.pop):
            if refScore-x._score > self.clusterEnergyCut:
                break
            remainder.append(i)

        gen = self.gen
        popSize = len(self.pop)

        selectionPopulation = []
        selectionPopulationd = {}
        if sanityCheck:
            for ind in self.pop:
                if ind._neighborRMSD < self.neighborSearchCutoff:
                    print "bestNeighborScore",ind._score,ind._neighborRMSD
                    break
                
        incluster = {} # keep track who is in a cluster
        tokeep = [] # clones of best in cluster

        if len(remainder)>2 and self.useClustering:
            print '\n  CLUSTERING %d individual(s) with %.2f'%(len(remainder), self.clusterEnergyCut)
            clusters = self.cluster(remainder)

            print "   CNUM  len best  Rmsd        Score         FEB        <Score>  stdev cluster"
            bestEinClust = []
            ediffCluster = 0.0
            ediffsCluster = [0.]
            for cnum, c in enumerate(clusters):
                maxi = self.pop[c[0]]._score
                rep = c[0]
                scores = [-self.pop[c[0]]._score]
                incluster[c[0]] = True
                for ind in c[1:]:
                    scores.append(-self.pop[ind]._score)
                    incluster[ind] = True
                    if self.pop[ind]._score > maxi:
                        maxi = self.pop[ind]._score
                        rep = ind
                bestEinClust.append(self.pop[rep]._score)
                if cnum < len(self.clustersBestE):
                    diff = self.pop[rep]._score - self.clustersBestE[cnum]
                    ediffCluster += diff
                    ediffsCluster.append(diff)
                #self.pop[rep].minimize(**self.hardmini)
                tokeep.append(self.pop[rep].clone())

                # every best in cluster get best score to guarantee offsprings
                self.pop[rep]._real_score = self.pop[rep]._score
                self.pop[rep]._score = refScore

                # add best in cluster to selection pop
                selectionPopulation.append(self.pop[rep])
                selectionPopulationd[rep] = True
                
                # compute RMSD
                print "   %3d  %3d  %3d"%(cnum, len(c), rep),
                if len(self.rmsdCalculators) and self.debug>0:
                    b, L_coords = self.pop[rep].phenotype
                    rmsdList = []
                    for rmsdc in self.rmsdCalculators:
                        rmsd = rmsdc.computeRMSD(L_coords)
                        rmsdList.append(rmsd)
                        print "%6.2f"%rmsd,
                else:
                    print "%6.2f"%-1,

                cstr = repr(c)
                if len(cstr)> 20:
                    cstr = cstr[:8]+'....'+cstr[-8:]
                print " %12.3f %12.3f %12.3f %6.3f %s"%(-self.pop[rep]._real_score,\
                      -self.pop[rep]._real_score-self.pop[rep].energies['LL'], numpy.mean(scores), numpy.std(scores), cstr),
                #print " %12.3f %12.3f %12.3f %6.3f %s"%(-self.pop[rep]._score,\
                #      -self.pop[rep]._score-self.pop[rep].energies['LL'], numpy.mean(scores), numpy.std(scores), cstr),
                print
                #if rmsdList[0] < 2.0:
                #    print '&CL1'
                #    self.writeSolution(
                #        self.pop[rep], 
                #        os.path.join(self.folder, 'isol_gen%d_%d.pdbqt'%(self.gen,rep)))
                #else:
                #    print
            #print out genes for best in cluster
            #if len(clusters):
            #    for cnum, c in enumerate(clusters):
            #        rep = selectionPopulation[cnum]
                    #print "#best in cluster %d %f"%(cnum, rep._real_score)
                    #print rep.getGenesPy()
                
            self.clustersBestE = bestEinClust

            if len(clusters) and ediffCluster==0.0:
                self.nbNoClusterEimprovement += 1
                if self.nbNoClusterEimprovement==self.stopNoImprove:
                    return "no improvement in best in %d clusters for %d generations"%(
                        len(clusters),self.stopNoImprove)
            else:
                self.nbNoClusterEimprovement=0
                    
            if len(clusters)==self.nbClusters and ediffCluster==0.0:
                self.nbGenWithCstClusters += 1
                
                if self.nbGenWithCstClusters > self.stopNoImprove:
                    if self.clusterEnergyCut >= 1.0:
                        self.clusterEnergyCut = self.clusterEnergyCut-0.1
            else:
                self.nbClusters = len(clusters)
                self.nbGenWithCstClusters = 0
        else: # nothing to cluster -> keep best
            tokeep.append(self.pop[0].clone())
           
        if len(remainder)>2 and self.useClustering:
            print "   CEDiff %f maxEDIFF %.3f noImprov: %3d CCst: %5.2f CEcut: %5.2f"%(
                ediffCluster, max(ediffsCluster), self.nbNoClusterEimprovement,
                self.nbGenWithCstClusters, self.clusterEnergyCut)

        # fill selectionPopulation using current population and skipping ones in cluster
        for num, ind in enumerate(self.pop):
            if incluster.get(num, None):
                continue # skip individuals in clusters
            #if incluster[num] != True: continue# skip individuals in clusters
            selectionPopulation.append(ind)
            selectionPopulationd[num] = True
            if len(selectionPopulation)==popSize: break
        #print selectionPopulation[0]._score
        print "Population Size", popSize, "Selection Population", len(selectionPopulation)
        if len(selectionPopulation) == len(tokeep): # only best in clusters are in selection population
            if len(selectionPopulation) <= 2:
                return '1 or 2 clusters left'
                        
        #selectionPopulation.sort(sc_maximize) # so that printing our choices makes sense
        #self.pop[:] = selectionPopulation[:]
        selPop = Population(selectionPopulation)
        # set a flag used to prevent individual from being minimized hard more tha once
        for ind in selPop:
            ind._hasBeenMinimizedInThisRound = False
        # prepare the population for selection
        selPop.touch()
        selPop.scale()

        # restore real score for best in cluster
        for ind in selPop:
            if hasattr(ind, '_real_score'):
                ind._score = ind._real_score
                delattr(ind, '_real_score')

        print ' ',
        
        # create crossover and mutate
        gaind = []
        gaindSize = popSize - len(tokeep)

        #if None in selPop.fitnesses():
        #import pdb; pdb.set_trace()
            
        while len(gaind) < gaindSize:
            newInds = []
            mom, dad = selPop.select(2)
            if mom is None: # ran out of choices
                break
            momScore = mom._score
            dadScore = dad._score
            # indices of selected individuals
            #p1Index, p2Index = self.pop.selector._selected
            p1Index, p2Index = selPop.selector._selected
            #print ' DEBUG: selected', p1Index, p2Index, id(mom), id(dad)
            minMstar = minDstar = ''
            #if self.nbNoClusterEimprovement==4 and self.nbScaleAmpBase<1:
            #    mom.genomePy.scaleDeltaAmplitudeBase(0.5)
            #    self.nbScaleAmpBase += 1

            if self.nbNoClusterEimprovement==3:
                if random() < 1.:
                    if not mom._hasBeenMinimizedInThisRound:
                        momc = mom.clone()
                        momcScore, nbStepsM = self.minimize(momc, **self.hardmini)
                        mom._hasBeenMinimizedInThisRound = True
                        if momcScore>refScore:
                            self.starMinData.append( ('*MIN MOM', gen, p1Index, p2Index, -momScore, -momcScore, 
                                                      -(momcScore-momScore)) )
                            #newInds.append(mom) 
                            mom = momc
                            minMstar = '***'
                            refScore = momcScore
                            newInds.append(momc)
                            #print ' DEBUG1: added mom1', id(momc)
                        elif incluster.get(p1Index, None) is None:
                            newInds.append(mom) 
                            incluster[p1Index] = True
                            #print ' DEBUG1: added mom2', id(mom)
                    elif incluster.get(p1Index, None) is None:
                        newInds.append(mom)
                        incluster[p1Index] = True
                        #print ' DEBUG1: added mom3', id(mom)

                    if not dad._hasBeenMinimizedInThisRound:
                        dadc = dad.clone()
                        dadcScore, nbStepsD = self.minimize(dadc, **self.hardmini)
                        dad._hasBeenMinimizedInThisRound = True
                        if dadcScore>refScore:
                            self.starMinData.append( ('*MIN DAD', gen, p1Index, p2Index, -dadScore, -dadcScore, 
                                                      -(dadcScore-dadScore)) )
                            dad = dadc
                            minDstar = '***'
                            refScore = dadScore
                            newInds.append(dadc)
                            #print ' DEBUG2: added dad1', id(dadc)
                        elif incluster.get(p2Index, None) is None:
                            newInds.append(dad) 
                            incluster[p2Index] = True
                            #print ' DEBUG2: added dad2', id(dad)
                    elif incluster.get(p2Index, None) is None:
                        newInds.append(dad)
                        incluster[p2Index] = True
                        #print ' DEBUG2: added dad3', id(dad)

                    if minMstar:
                        print '\n  (%3d)      MMI: %12.3f -> %12.3f%1s (%4d)'%(
                            p1Index, -momScore, -momcScore, minMstar, nbStepsM)
                    if minDstar:
                        print '\n  (%3d)       DMI: %12.3f -> %12.3f%1s (%4d)'%(
                            p2Index, -dadScore, -dadcScore, minDstar, nbStepsD)
            else: # no improvement is les than 3
                if incluster.get(p1Index, None) is None:
                    newInds.append(mom)
                    incluster[p1Index] = True
                if incluster.get(p2Index, None) is None:
                    newInds.append(dad)
                    incluster[p2Index] = True

            #if mom._score > dad._score:
            #    best2 = [mom, dad] # keep track of best individuals for this parents
            #else:
            #    best2 = [dad, mom] # keep track of best individuals for this parents
            ##
            ## version where we minimize cross and mutate
            ## in this version the lowest score found after cross
            ## might get lost if we mutate later and it does not get better
            ## yet it seems to work the best
            ##
            o1star = o2star = o1star1 = o2star1 = o1starmin = o2starmin = o1star1min = o2star1min = ''
            if p1Index != p2Index and flip_coin(self.p_crossover):
                bro, sis = self.crossover(mom,dad)
                broScore = bro.score()
                bronScore, nbSteps1 = self.minimize(bro, **self.quickmini)
                if bronScore>refScore:
                    o1star = '*'
                    refScore = bronScore
                    minbroScore, nbSteps1min = self.minimize(bro, **self.hardmini)
                    if minbroScore>refScore:
                        self.starMinData.append( ('*MIN CRO O1', gen, p1Index, p2Index, -bronScore, -minbroScore, 
                                                  -(minbroScore-bronScore)) )
                        o1starmin = '**'
                        self.dStarcnt +=1
                        refScore = minbroScore
                newInds.append(bro)
 
                sisScore = sis.score()
                sisnScore, nbSteps2 = self.minimize(sis, **self.quickmini)
                if sisnScore>refScore:
                    o2star = "*"
                    refScore = sisnScore
                    minsisScore, nbSteps2min = self.minimize(sis, **self.hardmini)
                    if minsisScore>refScore:
                        self.starMinData.append( ('*MIN CRO O2', gen, p1Index, p2Index, -sisnScore, -minsisScore, 
                                                  -(minsisScore-sisnScore)) )
                        o2starmin = '**'
                        self.dStarcnt +=1
                        refScore = minsisScore
                newInds.append(sis)

                if o1star or o2star:
                    print '\n  (%3d, %3d) CRO O1: %12.3f -> %12.3f%1s (%4d) O2: %12.3f -> %12.3f%1s (%4d)'%(
                        p1Index, p2Index,
                        -broScore, -bronScore, o1star, nbSteps1,
                        -sisScore, -sisnScore, o2star, nbSteps2),
                #if o1starmin or o2starmin:
                if o1starmin:
                    print '\n               SMI O1: %12.3f -> %12.3f%1s (%4d)'%(
                        -bronScore, -minbroScore, o1starmin, nbSteps1min),
                if o2starmin:
                    print '\n               SMI O2: %12.3f -> %12.3f%1s (%4d)'%(
                        -sisnScore, -minsisScore, o2starmin, nbSteps2min),           
                
            else:
                bro = dad
                sis = mom

            a1 = bro._score
            broc = bro.clone()
            mutated1, status = self.mutate(broc, self.p_mutation)
            self._nbMut += mutated1
            if mutated1:
                c = broc.score()
                broScore, nbStepsmut = self.minimize(broc, **self.quickmini)
                if broScore>refScore:
                    o1star1 = '*'
                    refScore = broScore
                    minbroScore, nbStepsmutmin = self.minimize(broc, **self.hardmini)
                    if minbroScore>refScore:
                        self.starMinData.append( ('*MIN MUT O1', gen, p1Index, p2Index, -broScore, -minbroScore, 
                                                  -(minbroScore-broScore)) )
                        o1star1min = '**'
                        self.dStarcnt +=1
                        refScore = minbroScore
                newInds.append(broc)

                if o1star1:
                    if not (o1star or o2star):
                        print '\n  (%3d, %3d)'%(p1Index, p2Index),
                    print '\n MUT O1: %12.3f -> %12.3f -> %12.3f%1s (%4d)'%(
                        -a1, -c, -minbroScore, o1star1, nbStepsmut),# self._totalStepInLastMinimize),
                if o1star1min:
                    #if not (o1star or o2star):
                    #    print '\n  (%3d, %3d)'%(p1Index, p2Index),
                    print '\n SMI O1: %12.3f -> %12.3f -> %12.3f%1s (%4d)'%(
                        -c, -broScore,-minbroScore, o1star1min, nbStepsmutmin),# self._totalStepInLastMinimize),
             
            a = sis._score
            sisc = sis.clone()
            mutated, status = self.mutate(sisc, self.p_mutation)
            self._nbMut += mutated
            if mutated:
                c = sisc.score()
                sisScore, nbStepsmut = self.minimize(sisc, **self.quickmini)
                if sisScore>refScore:
                    o2star1 = '*'
                    refScore = sisScore
                    minsisScore, nbStepsmutmin = self.minimize(sisc,**self.hardmini)
                    if minsisScore>refScore:
                        self.starMinData.append( ('*MIN MUT O2', gen, p1Index, p2Index, -sisScore, -minsisScore, 
                                                  -(minsisScore-sisScore)) )
                        o2star1min = '**'
                        self.dStarcnt +=1
                        refScore = minsisScore
                newInds.append(sisc)

                if o2star1:
                    if not (o1star or o2star or o1star1):
                        print '\n  (%3d, %3d)'%(p1Index, p2Index),
                    print '\n MUT O2: %12.3f -> %12.3f -> %12.3f%1s (%4d)'%(
                        -a, -c, -sisScore, o2star1, nbStepsmut),# self._totalStepInLastMinimize),
                if o2star1min:
                    #if not (o1star or o2star or o1star1):
                    #    print '\n  (%3d, %3d)'%(p1Index, p2Index),
                    print '\n SMI O2: %12.3f -> %12.3f -> %12.3f%1s (%4d)'%(
                        -c, -sisScore,-minsisScore, o2star1min, nbStepsmutmin),# self._totalStepInLastMinimize),
               
            ## if not (o1star or o2star or o1star1 or o2star1):
            ##    print '(%3d, %3d)'%(p1Index, p2Index),
            ## else:
            ##    print '\n ',
            if (o1star or o2star or o1star1 or o2star1):
                print '\n ',
                
            if len(newInds):
                newInds.sort(score_maximize)
                #if id(newInds[0]) in {}.fromkeys([id(x) for x in gaind]):
                #if self.neighborSearchCutoff>0.:
                    #if newInds[0]._score>self.bestNeighborScore:
                    #    self.bestNeighborInd=newInds[0].clone()
                    #    self.bestNeighborScore=newInds[0]._score
                gaind.append(newInds[0])
                if len(newInds)>1:
                    #if id(newInds[1]) in {}.fromkeys([id(x) for x in gaind]):
                    #    import pdb; pdb.set_trace()
                    gaind.append(newInds[1])
                    #if newInds[1]._score>self.bestNeighborScore:
                    #    self.bestNeighborInd=newInds[1].clone()
                    #    self.bestNeighborScore=newInds[1]._score
                ## nbAdded = 0
                ## for ind in newInds:
                ##     gaind.append(ind)
                ##     incluster[ind] = True
                ##     nbAdded += 1
                ##     if nbAdded==2:
                ##         break
        print

        # FIXME probably not needed
        #self.pop.sort()#score_maximize)
        #gaind.sort(score_maximize)
        pop = tokeep + gaind 
        pop.sort(score_maximize)
        if sanityCheck:
            assert self.lastGenBest <= pop[0]._score, "ERROR: lost best score %f -> %f"%(
                self.lastGenBest, pop[0]._score)
        if len(pop) < popSize:
            return 'ran out of choices'

        if self.neighborSearchCutoff>0:
            newBestNeighbor=False                 
            for ind in pop:
                if ind._neighborRMSD < self.neighborSearchCutoff:
                    if self.bestNeighborScore is None or ind._score >= self.bestNeighborScore:
                        self.bestNeighborScore=ind._score
                        self.bestNeighborInd=ind.clone()
                        newBestNeighbor=True
            if pop[popSize-1]._score>=self.bestNeighborScore or newBestNeighbor is False:
                pop[popSize-1]=self.bestNeighborInd
                
        self.pop[:] = pop[:(popSize)]

        if sanityCheck:
            NnewBestNeighbor=False 
            for ind in self.pop:
                if ind._neighborRMSD < self.neighborSearchCutoff:
                    print "bestNeighborScore",ind._score,ind._neighborRMSD
                    NnewBestNeighbor=True
                    break
            if NnewBestNeighbor is False:
                raise
            #    import pdb; pdb.set_trace()
        #import pdb; pdb.set_trace()

        # sanity check: check for duplicates in population
        #xdict = {}.fromkeys([id(x) for x in self.pop])
        #assert len(xdict) == popSize

        #for ind in self.pop:
        #    if ind._score > -100 and abs(ind._score+(ind.energies['RRL']+ind.energies['LL'])) > 0.1:
        #        import pdb; pdb.set_trace()
        print 'Mutation Rate %.4f'%(float(self._nbMut) / (popSize*len(self.pop[0].getGenesPy())))

        #import pdb; pdb.set_trace()

        return 'searching'
