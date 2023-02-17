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
# Author: Pradeep Anand Ravindranath
#
# Copyright: Pradeep Anand Ravindranath and TSRI 2016
#
#########################################################################
#
# $Header: /opt/cvs/AutoSite/fillBuriedness.py,v 1.3 2016/12/07 00:38:33 sanner Exp $
#
# $Id: fillBuriedness.py,v 1.3 2016/12/07 00:38:33 sanner Exp $
#
import numpy

class Buriedness:

    def __init__(self, centersL, atypesL, centersR,radiiR,origin):
        self.centersR = centersR
        self.radiiR = radiiR
        self.centersL = centersL
        self.atypesL = atypesL
        from bhtree import bhtreelib
        bht = bhtreelib.BHtree( self.centersR, None, 10)
        results = numpy.zeros(len(self.centersR), 'i')
        dist2 = numpy.zeros(len(self.centersR), 'f')
        radiiL = []
        for atype in atypesL:
            if atype=='C': radiiL.append(1.85)
            if atype=='O': radiiL.append(1.4)
            if atype=='H': radiiL.append(1.2)
        acResCenter = []
        acResRadii = []
        for cen in centersL:
            nb = bht.closePointsDist2(tuple(cen), 5.0, results, dist2)
            if nb:
                for ind, d2 in zip(results[:nb],dist2[:nb]):
                    if self.centersR[ind] not in acResCenter:
                        acResCenter.append(self.centersR[ind])
                        acResRadii.append(self.radiiR[ind])
        maxr = max(max(radiiL), max(radiiR))
        mini = numpy.min( (numpy.min(centersL, 0), numpy.min(centersR, 0)), 0 ) - maxr
        if origin is None:
            self.origin = mini
        else:
            self.origin = origin
        #import pdb;pdb.set_trace()
        self.centersR = acResCenter
        self.radiiR = acResRadii
        self.radiiL = radiiL
        
        

    def spheres2Grid(self,centers, radii, origin, spacing):
        """Compute a list of (i,j,k) grid point indices covered by spheres
        centered at 'centers' and with radii 'radii' indices on a grid with
        a given orgin and spacing

        [(i,j,k)], [(x,y,z)] = spheres2Grid(centers, radii, origin, spacing)
        """
        from math import ceil
        gpts = [] # list of (i,j,k) grid points covered by molecule
        pts = [] # list of (x,y,z) grid points covered by molecule
        ox, oy, oz = origin
        spacing1 = 1.0/spacing
        used = {}
        for n in xrange(len(centers)):
            x,y,z = centers[n]
            r = radii[n]
            x0 = int((x-ox)*spacing1)
            y0 = int((y-oy)*spacing1)
            z0 = int((z-oz)*spacing1)
            r0 = int(ceil(r*spacing1))
            r2 = r*r
            for i in range(x0-r0, x0+1+r0):
                gx = ox + i*spacing
                for j in range(y0-r0, y0+1+r0):
                    gy = oy + j*spacing
                    for k in range(z0-r0, z0+1+r0):
                        gz = oz + k*spacing
                        d2 = (gx-x)*(gx-x) + (gy-y)*(gy-y) + (gz-z)*(gz-z)
                        if d2<r2 and not used.has_key((i,j,k)):
                            used[(i,j,k)] = True
                            gpts.append((i,j,k))
                            pts.append((gx,gy,gz))

        return gpts, pts

    def inflatePocketandBuriedness(self,padding=1.0, spacing=1.0, origin=None):
        """Compute a numerical approximation of the 'burriedness' of
        spheres (centersL, raddiL).

        The 'burrideness' is a number ranging from 0.0 to 1.0 indicating
        what percentage of the surface of the first set of spheres is buried
        by the second set of spheres. The approximation is calculated by
        stamping onto a grid the first step of spheres with radii increased by padding,
        then unstamp the first set of spheres yeliding a shell around the first set of
        spheres comprising N grid points. Next the second set of spheres is unstamped
        leaving M set of points from the shell. The burridness is M/N. A value of 1.0
        means the first set is entirely buried by the second set of spheres.
        the spacing parameters control the accuracy of the result.
        If origin is provided it is used as grid lower left corner, else it is derived
        from the sets of spheres in order to make sure we encompass the spheres
        """

        centersR = self.centersR
        centersL = self.centersL
        radiiR = self.radiiR
        radiiL = self.radiiL
        origin = origin
        radiiL = [spacing+0.1] * len(centersL)
        # compute grid box dimensions as bounding box of the union of set of spheres
        if origin is None:
            maxr = max(max(radiiL), max(radiiR))
            mini = numpy.min( (numpy.min(centersL, 0), numpy.min(centersR, 0)), 0 ) - maxr
        else:
            mini = numpy.array(origin)

        # get (i,j,k) for enlarged first set of sphere
        set1LargeInd, set1LargeCoords = self.spheres2Grid(centersL, numpy.array(radiiL)+padding, mini, spacing)

        # remove indices for first set of sphere with original radii
        set1Ind, set1Coords = self.spheres2Grid(centersL, radiiL, mini, spacing)
        shellDict = {}
        for ind, coords in zip(set1LargeInd, set1LargeCoords):
            shellDict[ind] = coords 

        for inds in set1Ind:
            if shellDict.has_key(inds):
                del shellDict[inds]

        # remember the number of point in the shell

        # remove  indices for second set of sphere with original radii
        # we increase the raddi by 1.0 to get of points in the layer between close spheres
        #import pdb;pdb.set_trace()
        set2Ind, set2Coords = self.spheres2Grid(centersR, numpy.array(radiiR)+1, mini, spacing)
        for inds in set2Ind:
            if shellDict.has_key(inds):
                del shellDict[inds]
        
        #numpy.save('buriedNum/shellKept%d.xyzr'%num, shellDict.values())
        #numpy.save('buriedNum/shellRemoved%d.xyzr'%num, removedCoords)

        #Recalculate for inflated pocket
        for ind,coord in shellDict.iteritems():
            centersL = numpy.vstack((centersL, coord))
        radiiL = [spacing+0.1] * len(centersL)
        # get (i,j,k) for enlarged first set of sphere
        set1LargeInd, set1LargeCoords = self.spheres2Grid(centersL, numpy.array(radiiL)+padding, mini, spacing)

        # remove indices for first set of sphere with original radii
        set1Ind, set1Coords = self.spheres2Grid(centersL, radiiL, mini, spacing)
        shellNewDict = {}
        for ind, coords in zip(set1LargeInd, set1LargeCoords):
            shellNewDict[ind] = coords 

        for inds in set1Ind:
            if shellNewDict.has_key(inds):
                del shellNewDict[inds]

        # remember the number of point in the shell
        M = len(shellNewDict)

        # remove  indices for second set of sphere with original radii
        # we increase the raddi by 1.0 to get of points in the layer between close spheres
        #import pdb;pdb.set_trace()
        for inds in set2Ind:
            if shellNewDict.has_key(inds):
                del shellNewDict[inds]


        return 1.0 - (len(shellDict)/float(M)), shellDict


    def NumericalBurriedness(self,padding=1.0, spacing=0.5, origin=None):
        """Compute a numerical approximation of the 'burriedness' of
        spheres (centersL, raddiL).

        The 'burrideness' is a number ranging from 0.0 to 1.0 indicating
        what percentage of the surface of the first set of spheres is buried
        by the second set of spheres. The approximation is calculated by
        stamping onto a grid the first step of spheres with radii increased by padding,
        then unstamp the first set of spheres yeliding a shell around the first set of
        spheres comprising N grid points. Next the second set of spheres is unstamped
        leaving M set of points from the shell. The burridness is M/N. A value of 1.0
        means the first set is entirely buried by the second set of spheres.
        the spacing parameters control the accuracy of the result.
        If origin is provided it is used as grid lower left corner, else it is derived
        from the sets of spheres in order to make sure we encompass the spheres
        """

        centersR = self.centersR
        centersL = self.centersL
        radiiR = self.radiiR
        radiiL = self.radiiL
        origin = self.origin
        # compute grid box dimensions as bounding box of the union of set of spheres
        if origin is None:
            maxr = max(max(radiiL), max(radiiR))
            mini = numpy.min( (numpy.min(centersL, 0), numpy.min(centersR, 0)), 0 ) - maxr
        else:
            mini = numpy.array(origin)

        # get (i,j,k) for enlarged first set of sphere
        set1LargeInd, set1LargeCoords = self.spheres2Grid(centersL, numpy.array(radiiL)+padding, mini-2*padding, spacing)

        # remove indices for first set of sphere with original radii
        set1Ind, set1Coords = self.spheres2Grid(centersL, radiiL, mini-2*padding, spacing)
        shellDict = {}
        for ind, coords in zip(set1LargeInd, set1LargeCoords):
            shellDict[ind] = coords 

        for inds in set1Ind:
            if shellDict.has_key(inds):
                del shellDict[inds]

        # remember the number of point in the shell
        M = len(shellDict)

        # remove  indices for second set of sphere with original radii
        # we increase the raddi by 1.0 to get of points in the layer between close spheres
        #import pdb;pdb.set_trace()
        set2Ind, set2Coords = self.spheres2Grid(centersR, numpy.array(radiiR)+1, mini-2*padding, spacing)
        removedCoords = []
        for inds in set2Ind:
            if shellDict.has_key(inds):
                removedCoords.append(shellDict[inds])
                del shellDict[inds]
        #numpy.save('buriedNum/shellKept%d.xyzr'%num, shellDict.values())
        #numpy.save('buriedNum/shellRemoved%d.xyzr'%num, removedCoords)

        

        return 1.0 - (len(shellDict)/float(M))

        
