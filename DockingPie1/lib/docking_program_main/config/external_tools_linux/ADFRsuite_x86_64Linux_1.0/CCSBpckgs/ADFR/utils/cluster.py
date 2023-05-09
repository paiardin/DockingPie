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
# Copyright: M. Sanner and TSRI 2016
#
#########################################################################
#
# $Header: /mnt/raid/services/cvs/ADFR/utils/cluster.py,v 1.3 2016/12/07 00:38:33 sanner Exp $
#
# $Id: cluster.py,v 1.3 2016/12/07 00:38:33 sanner Exp $
#

def oneCluster(coords, seed, rmsdCalc, cutOff): 
    # compute cluster seeded at 'seed' for a given RMSD cutoff
    # coords is an array of (n, natoms,3) floats
    cluster = [seed]
    bestScore = -1000000.0
    notSelected = []
    for i, L_coords in enumerate(coords):
        if i==seed:
            continue
        # RMSD calc
        rmsd = rmsdCalc.computeRMSD(L_coords)
        if rmsd<=cutOff:
            cluster.append(i)
    return cluster

def clusterPoses(coords, order, rmsdCalc, cutOff=2.0):
    # cluster a set of solutions in the AutoDock fashion, i.e. use the best
    # solution as a seed and add all solution within cutOff RMSD to this cluster
    # then re-seed the algorithm with the best energy solution not yet clustered
    # until all solution given in "order" are clustered
    #
    # the coordinates are in coords (nsol, natoms, 3)
    # order id a list of indices into coords indicating the subset of
    # solutions to be clustered. These indices are expected to point to
    # solutions sorted by decreasing GA scores
    #
    remainder = order[:]
    clusters = []
    while len(remainder):
        seed = remainder[0]
        #print '%d left to cluster seed=%d'%(len(remainder), seed), remainder
        rmsdCalc.setRefCoords(coords[seed])
        cluster = [seed]
        notSelected = []
        for i in remainder:
            if i==seed:
                continue
            rmsd = rmsdCalc.computeRMSD(coords[i])
            #print '   %d %f'%(i, rmsd)
            if rmsd<=cutOff:
                cluster.append(i)
            else:
                notSelected.append(i)
        #print 'found cluster', cluster, seed
        clusters.append( cluster )
        remainder = notSelected
    return clusters

#clSize = len(oneCluster(coords, seed, rmsdCalc, 2.0))
