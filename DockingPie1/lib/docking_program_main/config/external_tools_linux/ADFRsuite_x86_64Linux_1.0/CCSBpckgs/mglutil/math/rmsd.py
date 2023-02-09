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
# Author: Sophie I. COON, William LINDSTROM, Michel F. SANNER
#
# Copyright: M. Sanner TSRI 2000
#
#############################################################################
#
# Last modified on Thu Aug  9 20:00:31 PDT 2001 by lindy
#
# $Header: /mnt/raid/services/cvs/python/packages/share1.5/mglutil/math/rmsd.py,v 1.10.4.1 2017/07/26 22:35:42 annao Exp $
#

import math

from mglutil.math.munkres import Munkres
import numpy

def getAtomIndicesPerType(atoms, typename='autodock_element'):
    d1 = {}
    for i, a in enumerate(atoms):
        try:
            d1[getattr(a, typename)].append(i)
        except KeyError:
            d1[getattr(a, typename)] = [i]
    return d1

    
class HungarianMatchingRMSD:
    """
    class to compute RMSD between 2 poses of the same molecule with pairing
    calculated using the Hungarian matching algorithm.

    typeIndicesRef are dictionary of where the key is an atom type and the value
    is a 0-based list of indices for atoms of that type in the list of atoms provided
    to the constructor (i.e. the reference atoms).

    the
    """

    def __init__(self, atoms, typeIndicesRef, typeIndicesMoving,
                 ignoreTypes=['HD']):
        # create rmsd calculator
        self.sortedRefAts = atoms
        self.typeIndicesRef = typeIndicesRef
        self.typeIndicesMoving = typeIndicesMoving
        self.atypes = typeIndicesRef.keys()
        for typeName in ignoreTypes:
            if typeName in self.atypes:
                self.atypes.remove(typeName)
        self.matching = None # will hold a list of computed pairs after matching
        

    def setRefCoords(self, coords):
        """
        set the reference atoms
        """
        self.sortedRefAts.updateCoords(coords)
        

    def computeRMSD(self, coords):
        """
        compute RMSD with reference atoms. coords are assumed to be in the same order
        as self.sortedRefAts
        """

        # use the Hungarian matching algorithm to pair up atoms
        # of the same time while minimizing the sum of the distances squared
        matching = []
        total = 0 # sum up square of distances

        # loop over atoms types
        for atype in self.atypes:
            inds1 = self.typeIndicesRef[atype]
            inds2 = self.typeIndicesMoving.get(atype, None)
            if inds2 is None:
                continue

            if len(inds1)==1 and len(inds2)==1: # only one atom of this type, matching is obvious
                matching.append( (inds1[0], inds2[0]) )
                x1, y1, z1 = self.sortedRefAts[inds1[0]].coords
                x2, y2, z2 = coords[inds2[0]]
                total += (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) + (z1-z2)*(z1-z2)
               
            elif len(inds1)==2 and len(inds2)==2: # only two atoms of this type, matching is obvious
                x1i,y1i,z1i = self.sortedRefAts[inds1[0]].coords
                x1j,y1j,z1j = self.sortedRefAts[inds1[1]].coords
                x2i,y2i,z2i = coords[inds2[0]]
                x2j,y2j,z2j = coords[inds2[1]]
                # compute dist(i,i) + dist(j,j)
                sum1 = ((x1i-x2i)*(x1i-x2i) + (y1i-y2i)*(y1i-y2i) + (z1i-z2i)*(z1i-z2i) +
                        (x1j-x2j)*(x1j-x2j) + (y1j-y2j)*(y1j-y2j) + (z1j-z2j)*(z1j-z2j))
                sum2 = ((x1i-x2j)*(x1i-x2j) + (y1i-y2j)*(y1i-y2j) + (z1i-z2j)*(z1i-z2j) +
                        (x1j-x2i)*(x1j-x2i) + (y1j-y2i)*(y1j-y2i) + (z1j-z2i)*(z1j-z2i))
                if sum1 < sum2:
                    matching.append( (inds1[0],inds1[0]) )
                    matching.append( (inds1[1],inds1[1]) )
                    total += sum1
                else:
                    matching.append( (inds1[0],inds1[1]) )
                    matching.append( (inds1[1],inds1[0]) )
                    total += sum2
                   
            else: # use Hungarian matching algorithm to find best assignment
                l1 = len(inds1)
                l2 = len(inds2)
                #print atype, inds
                matrix = numpy.zeros( (l1,l2), 'f')
                for i, n1 in enumerate(inds1):
                    x, y, z = self.sortedRefAts[n1].coords
                    for j, n2 in enumerate(inds2):
                        x1, y1, z1 = coords[n2]
                        matrix[i][j] = (x-x1)*(x-x1) + (y-y1)*(y-y1) + (z-z1)*(z-z1)

                # compute best assignment
                m = Munkres()
                indexes = m.compute(matrix.tolist())
                ltotal = 0.
                for row, column in indexes:
                    value = matrix[row][column]
                    ltotal += value
                    matching.append( (inds1[row], inds2[column]) )
                #print matrix
                #print 'total cost: %s %f %d' % (atype, ltotal, len(inds))
                total += ltotal

        self.matching = matching

        from math import sqrt
        rmsd = sqrt(total/len(matching))

        # sanity check. Recompute RMSD
        ## coords1 = []
        ## coords2 = []
        ## for pair in matching:
        ##     i, j = pair
        ##     #print i, j, atoms1[i].name, atoms2[j].name
        ##     coords1.append(self.sortedRefAts[i].coords)
        ##     coords2.append(coords[j])

        ## from mglutil.math.rmsd import RMSDCalculator
        ## rmsdCalc = RMSDCalculator(coords1)
        ## rmsd1 = rmsdCalc.computeRMSD(coords2)
        ## print "matched RMSD %.2f %.2f"%(rmsd, rmsd1)
        
        return rmsd


    
class HungarianMatchingRMSD_prody:
    """
    same as HungarianMatchingRMSD but uses sortedRefAtsCoords instead of
    sortedRefAts to be independent of MolKit and usable with Prody
    """

    def __init__(self, atomCoords, typeIndicesRef, typeIndicesMoving,
                 ignoreTypes=['HD']):
        # create rmsd calculator
        self.sortedRefAtsCoords = atomCoords
        self.typeIndicesRef = typeIndicesRef
        self.typeIndicesMoving = typeIndicesMoving
        self.atypes = typeIndicesRef.keys()
        for typeName in ignoreTypes:
            if typeName in self.atypes:
                self.atypes.remove(typeName)
        self.matching = None # will hold a list of computed pairs after matching
        

    def setRefCoords(self, coords):
        """
        set the reference atoms
        """
        self.sortedRefAtsCoords = coords
        

    def computeRMSD(self, coords):
        """
        compute RMSD with reference atoms. coords are assumed to be in the same order
        as self.sortedRefAtsCoords
        """

        # use the Hungarian matching algorithm to pair up atoms
        # of the same time while minimizing the sum of the distances squared
        matching = []
        total = 0 # sum up square of distances

        # loop over atoms types
        for atype in self.atypes:
            inds1 = self.typeIndicesRef[atype]
            inds2 = self.typeIndicesMoving.get(atype, None)
            if inds2 is None:
                continue

            if len(inds1)==1 and len(inds2)==1: # only one atom of this type, matching is obvious
                matching.append( (inds1[0], inds2[0]) )
                x1, y1, z1 = self.sortedRefAtsCoords[inds1[0]]
                x2, y2, z2 = coords[inds2[0]]
                total += (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) + (z1-z2)*(z1-z2)
               
            elif len(inds1)==2 and len(inds2)==2: # only two atoms of this type, matching is obvious
                x1i,y1i,z1i = self.sortedRefAtsCoords[inds1[0]]
                x1j,y1j,z1j = self.sortedRefAtsCoords[inds1[1]]
                x2i,y2i,z2i = coords[inds2[0]]
                x2j,y2j,z2j = coords[inds2[1]]
                # compute dist(i,i) + dist(j,j)
                sum1 = ((x1i-x2i)*(x1i-x2i) + (y1i-y2i)*(y1i-y2i) + (z1i-z2i)*(z1i-z2i) +
                        (x1j-x2j)*(x1j-x2j) + (y1j-y2j)*(y1j-y2j) + (z1j-z2j)*(z1j-z2j))
                sum2 = ((x1i-x2j)*(x1i-x2j) + (y1i-y2j)*(y1i-y2j) + (z1i-z2j)*(z1i-z2j) +
                        (x1j-x2i)*(x1j-x2i) + (y1j-y2i)*(y1j-y2i) + (z1j-z2i)*(z1j-z2i))
                if sum1 < sum2:
                    matching.append( (inds1[0],inds1[0]) )
                    matching.append( (inds1[1],inds1[1]) )
                    total += sum1
                else:
                    matching.append( (inds1[0],inds1[1]) )
                    matching.append( (inds1[1],inds1[0]) )
                    total += sum2
                   
            else: # use Hungarian matching algorithm to find best assignment
                l1 = len(inds1)
                l2 = len(inds2)
                #print atype, inds
                matrix = numpy.zeros( (l1,l2), 'f')
                for i, n1 in enumerate(inds1):
                    x, y, z = self.sortedRefAtsCoords[n1]
                    for j, n2 in enumerate(inds2):
                        x1, y1, z1 = coords[n2]
                        matrix[i][j] = (x-x1)*(x-x1) + (y-y1)*(y-y1) + (z-z1)*(z-z1)

                # compute best assignment
                m = Munkres()
                indexes = m.compute(matrix.tolist())
                ltotal = 0.
                for row, column in indexes:
                    value = matrix[row][column]
                    ltotal += value
                    matching.append( (inds1[row], inds2[column]) )
                #print matrix
                #print 'total cost: %s %f %d' % (atype, ltotal, len(inds))
                total += ltotal

        self.matching = matching

        from math import sqrt
        rmsd = sqrt(total/len(matching))

        # sanity check. Recompute RMSD
        ## coords1 = []
        ## coords2 = []
        ## for pair in matching:
        ##     i, j = pair
        ##     #print i, j, atoms1[i].name, atoms2[j].name
        ##     coords1.append(self.sortedRefAts[i].coords)
        ##     coords2.append(coords[j])

        ## from mglutil.math.rmsd import RMSDCalculator
        ## rmsdCalc = RMSDCalculator(coords1)
        ## rmsd1 = rmsdCalc.computeRMSD(coords2)
        ## print "matched RMSD %.2f %.2f"%(rmsd, rmsd1)
        
        return rmsd


class RMSDCalculator:
    """
    This class implements method to compute RMSD and distance vector
    between two given lists of coordinates.
    """
    def __init__(self, refCoords = None):
        self.refCoords = refCoords

    def setRefCoords(self, refCoords):
        self.refCoords = numpy.array(refCoords)
        
    def computeRMSD(self, listCoords):
        """rmsd <- computRMSD(listCoords)
        rmsd returns the overall root mean square distance (rmsd) and
        also sets self.distVect as the vector of distances between each
        pair of points.
        """
        if self.refCoords is None:
            raise ValueError("no reference coordinates set")
        if len(self.refCoords) != len(listCoords):
            raise ValueError("input vector length mismatch")

        deltaVect = self.refCoords - numpy.array(listCoords)
        distSquaredVect = numpy.sum(numpy.transpose(deltaVect*deltaVect))
        self.distVect = numpy.sqrt(distSquaredVect)
        self.rmsd = math.sqrt(numpy.sum(distSquaredVect)/len(self.refCoords))
        return self.rmsd

    def computeRMSDfast(self, coordsArray):
         delta = self.refCoords - coordsArray
         d2 = numpy.sum(delta*delta, 0)
         self.rmsd = math.sqrt( numpy.sum(d2)/len(self.refCoords) )
         return self.rmsd

class MorpicRMSD:
    """
    Calculate RMSD as the minimum RMSD over a set of morphisms
    Each mophism is provided as a list of of 0-based atom indices paires of
    length matching the length of self.refCoords
    """
    def __init__(self, morphisms, refCoords = None):
        self.morphisms = None
        self.refCoords = refCoords
        self.setMorphisms(morphisms)
        
    def setMorphisms(self, morphisms):
        #if self.refCoords:
        #    for iso in morphisms:
        #        assert len(self.refCoords)==len(iso)
        self.morphisms = numpy.array(morphisms)

    def setRefCoords(self, refCoords):
        #if self.morphisms:
        #    assert len(refCoords)==len(self.morphisms[0])
        self.refCoords = numpy.array(refCoords)
        
    def computeRMSD(self, listCoords):
        """rmsd <- computRMSD(listCoords)
        rmsd returns the overall root mean square distance (rmsd) and
        also sets self.distVect as the vector of distances between each
        pair of points.
        """
        if self.refCoords is None:
            raise ValueError("no reference coordinates set")
        if len(self.refCoords) != len(listCoords):
            raise ValueError("input vector length mismatch")
        coords = numpy.array(listCoords)

        rmsds = []
        #import pdb; pdb.set_trace()
        for morph in self.morphisms:
            #import pdb; pdb.set_trace()
            rCoords = self.refCoords[morph[:, 0]]
            mCoords = coords[morph[:, 1]]
            deltaVect = rCoords - numpy.array(mCoords)
            distSquaredVect = numpy.sum(numpy.transpose(deltaVect*deltaVect))
            self.distVect = numpy.sqrt(distSquaredVect)
            rmsds.append(math.sqrt(numpy.sum(distSquaredVect)/len(rCoords)))

        self._rmsds = rmsds
        #print 'RMSD auto', rmsds
        #import pdb; pdb.set_trace()
        return min(self._rmsds)
