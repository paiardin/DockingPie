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
# $Header: /opt/cvs/AutoSite/utils/clusterTPoints.py,v 1.4 2017/06/09 22:58:45 zhangyq Exp $
#
# $Id: clusterTPoints.py,v 1.4 2017/06/09 22:58:45 zhangyq Exp $
#
import numpy, random
from math import sqrt

class Kmeans:
    def __init__(self, pts, k, cutoff=0.0001, returnCentroids=False):
        self.pts = pts
        self.nclusters = k
        self.tolerance = cutoff
        self.centroids = returnCentroids
        self.kmeansCluster()
        
    def kmeansCluster(self):
        """Apply kmeans clustering to a set of 3D points.
        k is the desired number of clusters
        cutoff is the shift threshold for centroids
        Returns a list of clusters (as 3D points) or cluster centroids depending
        on the values of returnCentroids
        """
        pts = self.pts
        # Randomly sample k Points from the points list, build Clusters around them
        initialpts = random.sample(pts, self.nclusters)

        # Enter the program loop
        centroids = []
        for p in initialpts:
            centroids.append((p))

        while True:
            clusters = []
            for p in initialpts: clusters.append([])
            # For each Point:
            for p in pts:
                # Figure out which Cluster's centroid is the nearest
                ccx, ccy, ccz = centroids[0]

                smallest_distance = sqrt( (p[0]-ccx)*(p[0]-ccx) +
                                          (p[1]-ccy)*(p[1]-ccy) + 
                                          (p[2]-ccz)*(p[2]-ccz) )
                index = 0
                for i in range(len(clusters[1:])):
                    cx, cy, cz = centroid = centroids[i+1]
                    distance = sqrt( (p[0]-cx)*(p[0]-cx) +
                                     (p[1]-cy)*(p[1]-cy) +
                                     (p[2]-cz)*(p[2]-cz) )
                    if distance < smallest_distance:
                        smallest_distance = distance
                        index = i+1
                # Add this Point to that Cluster's corresponding list
                clusters[index].append(p)
            # Update each Cluster with the corresponding list
            # Record the biggest centroid shift for any Cluster
            biggest_shift = 0.0
            for i in range(len(clusters)):
                ncx, ncy, ncz = numpy.mean(clusters[i],axis=0)
                shift = sqrt( (ncx-centroids[i][0])*(ncx-centroids[i][0]) +
                              (ncy-centroids[i][1])*(ncy-centroids[i][1]) +
                              (ncz-centroids[i][2])*(ncz-centroids[i][2]) )
                centroids[i] = [ncx, ncy, ncz]
                biggest_shift = max(biggest_shift, shift)
            # If the biggest centroid shift is less than the cutoff, stop
            if biggest_shift < self.tolerance: break

        if self.centroids:
            self.cCentroids = centroids
        else:
            self.cClusters = clusters

    
class DensityClustering:

    def __init__(self, spacing, neighborPts=14):
        self._indices = None
        self.neighborpts = neighborPts
        self.spacing = spacing
        self._clusters = None
        self._clen = []
        
    def findClustersD(self, ptsIndices, ptsTypes=None, cVolcut = 50):
        """returns list of contiguous clusters of filter points with every point in a
        cluster having the specified number of neighbors. """
	    #print "wathch me", self.spacing
        # get number of neighbors needed for density clustering
        nneighbor = self.neighborpts
        sx,sy,sz = self.spacing
        si = int(round(1.0/sx))
        sj = int(round(1.0/sy))
        sk = int(round(1.0/sz))        
        clusters = [] # list of clusters identified
        cluster = [] # current cluster to build

        # keycache
        keycache = {}
        indcache = {}
        typesCOcache = {}
        n = 0
        for i,j,k in ptsIndices:
            #if 
            key = '%d,%d,%d'%(i,j,k)
            keycache[(i,j,k)] = key
            indcache[(i,j,k)] = n
            #if ptsTypes[n]!='H':
            #    typesCOcache[(i,j,k)] = ptsTypes[n]
            n += 1
            #if i%si == 0 and j%sj==0 and k%sk==0:
            #    newptsIndices.append((i,j,k))
            
        
        indices = {} # indices of points not clustered yet
        #n = 0
        for ii in range(len(ptsIndices)):
            i,j,k = ptsIndices[ii]
            if i%si == 0 and j%sj==0 and k%sk==0:
                indices[ii] = ii
                #n += 1
        #import pdb
        #pdb.set_trace()
        ci = 0 # index of point in cluster
        index = indices.values()[0]
        i, j, k = ptsIndices[index] # first point
        cluster.append( index )
        indices.pop(index)
        #indices.pop(indices.keys()[0])
        incluster = {index:True}
        while True:
            neighbors = [] # list of neighbors of current point

            # build list of grid points neighboring index and in ptsIndices
            if keycache.has_key( (i-si,j,  k  ) ): neighbors.append( indcache[(i-si,j,  k  )] )
            if keycache.has_key( (i+si,j,  k  ) ): neighbors.append( indcache[(i+si,j,  k  )] )
            if keycache.has_key( (i,  j-sj,k  ) ): neighbors.append( indcache[(i,  j-sj,k  )] )
            if keycache.has_key( (i,  j+sj,k  ) ): neighbors.append( indcache[(i,  j+sj,k  )] )
            if keycache.has_key( (i,  j,  k-sk) ): neighbors.append( indcache[(i,  j  ,k-sk)] )
            if keycache.has_key( (i,  j,  k+sk) ): neighbors.append( indcache[(i,  j  ,k+sk)] )

            if keycache.has_key( (i-si,j+sj,  k  ) ): neighbors.append( indcache[(i-si,j+sj,  k  )] )
            if keycache.has_key( (i+si,j+sj,  k  ) ): neighbors.append( indcache[(i+si,j+sj,  k  )] )
            if keycache.has_key( (i-si,j-sj,  k  ) ): neighbors.append( indcache[(i-si,j-sj,  k  )] )
            if keycache.has_key( (i+si,j-sj,  k  ) ): neighbors.append( indcache[(i+si,j-sj,  k  )] )
            if keycache.has_key( (i-si,  j ,k+sk ) ): neighbors.append( indcache[(i-si,  j ,k+sk )] )
            if keycache.has_key( (i+si,  j ,k+sk ) ): neighbors.append( indcache[(i+si,  j ,k+sk )] )
            if keycache.has_key( (i-si,  j ,k-sk ) ): neighbors.append( indcache[(i-si,  j ,k-sk )] )
            if keycache.has_key( (i+si,  j ,k-sk ) ): neighbors.append( indcache[(i+si,  j ,k-sk )] )
            if keycache.has_key( (i   ,j-sj,k+sk ) ): neighbors.append( indcache[(i   ,j-sj,k+sk )] )
            if keycache.has_key( (i   ,j+sj,k+sk ) ): neighbors.append( indcache[(i   ,j+sj,k+sk )] )
            if keycache.has_key( (i   ,j-sj,k-sk ) ): neighbors.append( indcache[(i   ,j-sj,k-sk )] )
            if keycache.has_key( (i   ,j+sj,k-sk ) ): neighbors.append( indcache[(i   ,j+sj,k-sk )] )

            if keycache.has_key( (i+si,j+sj,k+sk ) ): neighbors.append( indcache[(i+si,j+sj,k+sk )] )
            if keycache.has_key( (i-si,j+sj,k+sk ) ): neighbors.append( indcache[(i-si,j+sj,k+sk )] )
            if keycache.has_key( (i+si,j-sj,k+sk ) ): neighbors.append( indcache[(i+si,j-sj,k+sk )] )
            if keycache.has_key( (i+si,j+sj,k-sk ) ): neighbors.append( indcache[(i+si,j+sj,k-sk )] )
            if keycache.has_key( (i-si,j-sj,k+sk ) ): neighbors.append( indcache[(i-si,j-sj,k+sk )] )
            if keycache.has_key( (i-si,j+sj,k-sk ) ): neighbors.append( indcache[(i-si,j+sj,k-sk )] )
            if keycache.has_key( (i+si,j-sj,k-sk ) ): neighbors.append( indcache[(i+si,j-sj,k-sk )] )
            if keycache.has_key( (i-si,j-sj,k-sk ) ): neighbors.append( indcache[(i-si,j-sj,k-sk )] )
            
            #neighCO = [x for x in neighbors if not typesCOcache.has_key(x)]
    
            if len(neighbors)>=nneighbor: # if point index (i,j,k) has enough neighbors
            #if len(neighCO)>=nneighbor:
                neigh = [x for x in neighbors if not incluster.has_key(x)]
                cluster.extend(neigh)
                for n in neigh:
                    incluster[n] = True
                    indices.pop(n)             

            ci += 1
            if ci<len(cluster): # more points in this cluster can be checked for neighbors
                i,j,k = ptsIndices[cluster[ci]]
            else: # all points in cluster have been cheked for neighbors
                #print 'CLUSTER', len(cluster)
                clusters.append(cluster)
                if len(indices)==0:
                     break
                cluster = []
                ci = 0
                index = indices.values()[0]
                i,j,k = ptsIndices[index]
                cluster.append( index )
                indices.pop(index)
                #indices.pop(indices.keys()[0])
                incluster[index] = True
                
        clen = [len(x) for x in clusters]  
        cl = []
        for i in numpy.argsort(clen)[::-1]:
            cl.append(clusters[i])#numpy.array(clusters[i], 'i'))
        clv = [x for x in cl if len(x) >= cVolcut]
        cl = clv

        #import pdb;pdb.set_trace()
        if (si,sj,sk) != (1.0,1.0,1.0):
            nclusters = []
            for c in clv:
                #c = c.tolist()
                c2 = []
                for index in c:
                    i, j, k = ptsIndices[index]
                    for ind in range(1,si):
                        if keycache.has_key( (i-ind,j,  k  ) ) and indcache[(i-ind, j,  k  )] not in c and indcache[(i-ind, j,  k  )] not in c2: c2.append( indcache[(i-ind, j,  k  )] )
                        if keycache.has_key( (i+ind,j,  k  ) ) and indcache[(i+ind, j,  k  )] not in c and indcache[(i+ind, j,  k  )] not in c2: c2.append( indcache[(i+ind, j,  k  )] )
                    for ind in range(1,sj):
                        if keycache.has_key( (i,  j-ind,k  ) ) and indcache[(i, j-ind,  k  )] not in c and indcache[(i, j-ind,  k  )] not in c2: c2.append( indcache[(i,  j-ind,k  )] )
                        if keycache.has_key( (i,  j+ind,k  ) ) and indcache[(i, j+ind,  k  )] not in c and indcache[(i, j+ind,  k  )] not in c2: c2.append( indcache[(i,  j+ind,k  )] )
                    for ind in range(1,sk):
                        if keycache.has_key( (i,  j,  k-ind) ) and indcache[(i, j,  k-ind  )] not in c and indcache[(i, j,  k-ind  )] not in c2: c2.append( indcache[(i,  j  ,k-ind)] )
                        if keycache.has_key( (i,  j,  k+ind) ) and indcache[(i, j,  k+ind  )] not in c and indcache[(i, j,  k+ind  )] not in c2: c2.append( indcache[(i,  j  ,k+ind)] )
                    for indi in range(1,si):
                        for indj in range(1,sj):
                            if keycache.has_key( (i-indi,j+indj,  k  ) ) and indcache[(i-indi,j+indj,  k  )] not in c and indcache[(i-indi,j+indj,  k  )] not in c2: c2.append( indcache[(i-indi,j+indj,  k  )] )
                            if keycache.has_key( (i+indi,j+indj,  k  ) ) and indcache[(i+indi,j+indj,  k  )] not in c and indcache[(i+indi,j+indj,  k  )] not in c2: c2.append( indcache[(i+indi,j+indj,  k  )] )
                            if keycache.has_key( (i-indi,j-indj,  k  ) ) and indcache[(i-indi,j-indj,  k  )] not in c and indcache[(i-indi,j-indj,  k  )] not in c2: c2.append( indcache[(i-indi,j-indj,  k  )] )
                            if keycache.has_key( (i+indi,j-indj,  k  ) ) and indcache[(i+indi,j-indj,  k  )] not in c and indcache[(i+indi,j-indj,  k  )] not in c2: c2.append( indcache[(i+indi,j-indj,  k  )] )
                    for indi in range(1,si):
                        for indk in range(1,sk):
                            if keycache.has_key( (i-indi,j,  k+indk  ) ) and indcache[(i-indi,j,  k+indk  )] not in c and indcache[(i-indi,j,  k+indk  )] not in c2: c2.append( indcache[(i-indi,j,  k+indk  )] )
                            if keycache.has_key( (i+indi,j,  k+indk  ) ) and indcache[(i+indi,j,  k+indk  )] not in c and indcache[(i+indi,j,  k+indk  )] not in c2: c2.append( indcache[(i+indi,j,  k+indk  )] )
                            if keycache.has_key( (i-indi,j,  k-indk  ) ) and indcache[(i-indi,j,  k-indk  )] not in c and indcache[(i-indi,j,  k-indk  )] not in c2: c2.append( indcache[(i-indi,j,  k-indk  )] )
                            if keycache.has_key( (i+indi,j,  k-indk  ) ) and indcache[(i+indi,j,  k-indk  )] not in c and indcache[(i+indi,j,  k-indk  )] not in c2: c2.append( indcache[(i+indi,j,  k-indk  )] )
                    for indj in range(1,sj):
                        for indk in range(1,sk):
                            if keycache.has_key( (i,j-indj,  k+indk  ) ) and indcache[(i,j-indj,  k+indk  )] not in c and indcache[(i,j-indj,  k+indk  )] not in c2: c2.append( indcache[(i,j-indj,  k+indk  )] )
                            if keycache.has_key( (i,j+indj,  k+indk  ) ) and indcache[(i,j+indj,  k+indk  )] not in c and indcache[(i,j+indj,  k+indk  )] not in c2: c2.append( indcache[(i,j+indj,  k+indk  )] )
                            if keycache.has_key( (i,j-indj,  k-indk  ) ) and indcache[(i,j-indj,  k-indk  )] not in c and indcache[(i,j-indj,  k-indk  )] not in c2: c2.append( indcache[(i,j-indj,  k-indk  )] )
                            if keycache.has_key( (i,j+indj,  k-indk  ) ) and indcache[(i,j+indj,  k-indk  )] not in c and indcache[(i,j+indj,  k-indk  )] not in c2: c2.append( indcache[(i,j+indj,  k-indk  )] )
                    
                    for indi in range(1,si):
                        for indj in range(1,sj):
                            for indk in range(1,sk):
                                if keycache.has_key( (i+indi,j+indj,k+indk ) ) and indcache[(i+indi,j+indj,k+indk )] not in c and indcache[(i+indi,j+indj,k+indk )] not in c2: c2.append( indcache[(i+indi,j+indj,k+indk )] )
                                if keycache.has_key( (i-indi,j+indj,k+indk ) ) and indcache[(i-indi,j+indj,k+indk )] not in c and indcache[(i-indi,j+indj,k+indk )] not in c2: c2.append( indcache[(i-indi,j+indj,k+indk )] )
                                if keycache.has_key( (i+indi,j-indj,k+indk ) ) and indcache[(i+indi,j-indj,k+indk )] not in c and indcache[(i+indi,j-indj,k+indk )] not in c2: c2.append( indcache[(i+indi,j-indj,k+indk )] )
                                if keycache.has_key( (i+indi,j+indj,k-indk ) ) and indcache[(i+indi,j+indj,k-indk )] not in c and indcache[(i+indi,j+indj,k-indk )] not in c2: c2.append( indcache[(i+indi,j+indj,k-indk )] )
                                if keycache.has_key( (i-indi,j-indj,k+indk ) ) and indcache[(i-indi,j-indj,k+indk )] not in c and indcache[(i-indi,j-indj,k+indk )] not in c2: c2.append( indcache[(i-indi,j-indj,k+indk )] )
                                if keycache.has_key( (i-indi,j+indj,k-indk ) ) and indcache[(i-indi,j+indj,k-indk )] not in c and indcache[(i-indi,j+indj,k-indk )] not in c2: c2.append( indcache[(i-indi,j+indj,k-indk )] )
                                if keycache.has_key( (i+indi,j-indj,k-indk ) ) and indcache[(i+indi,j-indj,k-indk )] not in c and indcache[(i+indi,j-indj,k-indk )] not in c2: c2.append( indcache[(i+indi,j-indj,k-indk )] )
                                if keycache.has_key( (i-indi,j-indj,k-indk ) ) and indcache[(i-indi,j-indj,k-indk )] not in c and indcache[(i-indi,j-indj,k-indk )] not in c2: c2.append( indcache[(i-indi,j-indj,k-indk )] )
                #c.extend(c2)
                nc2 = [x for x in c2 if not incluster.has_key(x)]
                c.extend(nc2)
                for idx in nc2:
                    incluster[idx] = True
                #nclusters.append(c)
                nclusters.append(numpy.array(c, 'i'))
            cl = nclusters
        
        self._clusters = cl # clusters sorted by size
        self._clen = [len(x) for x in cl]   # sorted clusters lengths
