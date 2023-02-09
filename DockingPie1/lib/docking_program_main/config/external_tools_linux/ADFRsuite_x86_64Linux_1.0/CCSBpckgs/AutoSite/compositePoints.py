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
# Author: Pradeep Anand Ravindranath
#
# Copyright: Pradeep Anand Ravindranath and TSRI 2016
#
#########################################################################
#
# $Header: /opt/cvs/AutoSite/compositePoints.py,v 1.4 2017/06/10 01:02:48 zhangyq Exp $
#
# $Id: compositePoints.py,v 1.4 2017/06/10 01:02:48 zhangyq Exp $
#

from ADFR.utils.MakeGrids import CalculateAD4Grids
import platform, os, subprocess, tempfile, numpy, shutil
from ADFR.utils.MakeGpf import ADGPF
from ADFRcc.adfr import GridMap
from MolKit2.selection import Selection
from time import time
from MolKit2 import Read
from utils.clusterTPoints import DensityClustering
from scoreClusters import scoreClusters
from clusterNode import clusterNode


class CompositePoints(CalculateAD4Grids):
    
    def getASPoints(self, carbon_cutoff=-0.30, oxygen_cutoff=-0.66, hydrogen_cutoff=-0.5):
        cmapFile = self.mapFiles['C']
        cmap = GridMap()
        cmap.loadFromMapFile('C', '', cmapFile)
        omapFile = self.mapFiles['OA']
        omap = GridMap()
        omap.loadFromMapFile('OA', '', omapFile)
        hmapFile = self.mapFiles['HD']
        hmap = GridMap()
        hmap.loadFromMapFile('HD', '', hmapFile)

        cdata = cmap.getGridDataPy()#[::3, ::3, ::3]
        ox, oy, oz = cmap.getOriginPy().round(3)
        sx = sy = sz = cmap.getDistBetweenGridPoints()
        odata = omap.getGridDataPy()
        hdata = hmap.getGridDataPy()
        cmap = numpy.select( [cdata<carbon_cutoff,cdata>=carbon_cutoff],[cdata,[999]])
        omap = numpy.select( [odata<oxygen_cutoff,odata>=oxygen_cutoff],[odata,[999]])
        hmap = numpy.select( [hdata<hydrogen_cutoff,hdata>=hydrogen_cutoff],[hdata,[999]])
        tmp = numpy.zeros( cdata.shape, 'f')
        mini = numpy.zeros( cdata.shape, 'f')
        numpy.minimum(cmap,omap,tmp)
        numpy.minimum(tmp,hmap,mini)
        i,j,k = numpy.where(mini!=999)
        amin = numpy.argmin([cmap,omap,hmap],axis=0)
        coords = numpy.zeros((len(i),3),'f')
        coords[:,0] = ox + i*sx
        coords[:,1] = oy + j*sy
        coords[:,2] = oz + k*sz
        ati = amin[[i,j,k]]
        atypes = numpy.select( [ati==0, ati==1, ati==2], ['C', 'O', 'H'])
        indices = numpy.zeros((len(i),3),'i')
        indices[:,0] =  i
        indices[:,1] =  j
        indices[:,2] =  k
        self._indices = indices   # list of (i,j,k) in grid
        self._coords = coords.round(3)    # list of (x,y,z) of selected grid points
        self._potential = mini[[i,j,k]] # list of potential values at selected grid points
        self._atype = atypes     # list of atom types at selected grid points
        
        #import pdb;pdb.set_trace()

    #Find the best cutoff to build the translational map.
    def bestCutoffClustering(self,receptor,spacing,carbon_cutoff=-0.30,oxygen_cutoff=-0.66,hydrogen_cutoff=-0.5,nbSteps=10,nNeighbor=14):
        rangeC =  carbon_cutoff *0.5
        rangeO =  oxygen_cutoff *0.5
        rangeH =  hydrogen_cutoff *0.5
        stepC= rangeC/nbSteps
        stepO= rangeO/nbSteps
        stepH= rangeH/nbSteps
        max_eff = 0.0
        max_effi = None
        cluster_dict={}
        clusters_to_keep=[]
        propclusters_to_keep=[]
        clusterid=0
        scanStopflag=False


        for i in range(nbSteps):
            ccut = carbon_cutoff - stepC*i
            ocut = oxygen_cutoff - stepO*i
            hcut = hydrogen_cutoff - stepH*i
            print "Scanning at:",ccut,ocut,hcut            
            self.getASPoints(carbon_cutoff=ccut, oxygen_cutoff=ocut,hydrogen_cutoff=hcut)
            dcl = DensityClustering([spacing,spacing,spacing], neighborPts=nNeighbor)
            if len(self._indices)>0:
                dcl.findClustersD(self._indices, self._atype)
            else:
                continue
            clusters, clProp = scoreClusters(receptor, dcl, self, inflate=False)
            clustersorted = sorted(clusters,key=lambda x:x[4],reverse=True)
            if len(clusters)>0:
                clPropsorted = sorted(clProp,key=lambda x:x[5],reverse=True)
            else:
                print "%d #: %d"%(i, len(clusters))
                continue
            print 'clust.| Energy| # of |Rad. of | energy |   bns    |score   |'
            print 'number|       |points|gyration|per vol.|buriedness|v*b^2/rg|'
            print '------+-------+------+--------+--------+----------+--------|'
            clN = 1
            clustNum = 1
            for cl, clp in zip(clustersorted, clPropsorted):
                if clp[1]>3000:
                    scanStopflag=True
                    print "Flooded with more than 3000 points in one pocket, stop after this scan"
                if clp[1]>5000:
                    print "Flooded with more than 5000 points in one pocket, ignore this scan"
                    break
                clusterid += 1
                #build dictionary of all the coords in the current cluster
                newClustDict = {}.fromkeys([str(x) for x in cl[0]])
                contains = []
                registercluster = False
                clusters_to_rm=[]
                newTreeNode=clusterNode(id=clusterid,size=clp[1],score=clp[5],rg=clp[3],totalE=clp[0],buriedness=clp[4],gen=i)                          
                if i==0:                    
                    cluster_dict[str(cl[0][10])]=newTreeNode
                    print '%5d %8.2f %5d %7.2f %8.2f     %.3f  %8.2f '%(clN,clp[0],clp[1],clp[3],clp[2],clp[4],clp[5])
                    clustNum += 1
                    clN=clN+1
                    clusters_to_keep.append(cl)
                    propclusters_to_keep.append(clp)
                    continue

                for key, value in cluster_dict.iteritems():
                    if newClustDict.has_key(key):
                        newTreeNode.add_child(value)
                        clusters_to_rm.append(key)
                        registercluster = True

                for delkey in clusters_to_rm:
                    del cluster_dict[delkey]

                if registercluster or ccut<-0.299 or len(cluster_dict)<10:
                    print '%5d %8.2f %5d %7.2f %8.2f     %.3f  %8.2f '%(clN,clp[0],clp[1],clp[3],clp[2],clp[4],clp[5])
                    newTreeNode.id=len(clusters_to_keep)+1
                    clusters_to_keep.append(cl)
                    propclusters_to_keep.append(clp)
                    cluster_dict[str(cl[0][10])]=newTreeNode

                clN=clN+1

            if scanStopflag:
                break
        
        headNode=clusterNode(id=9999,size=999999,score=0)
        for key, value in cluster_dict.iteritems():
            headNode.add_child(value)
        #import pdb; pdb.set_trace()
        #import json
        #with open('pockettree.txt', 'w') as outfile:
        #    json.dump(headNode.writeJSON(), outfile)

        dcl = DensityClustering([spacing,spacing,spacing], neighborPts=nNeighbor)
        clls=[]
        cllens=[]
        self._coords = self._coords.round(3)
        gc_dict = {}
        for ii in range(len(self._indices)):
            gc_dict[str(self._coords[ii])] = ii


        #import pdb; pdb.set_trace()
        for keptcluster in clusters_to_keep:
            cll=[]
            for clcoor in keptcluster[1]:
                cll.append(gc_dict[str(clcoor.round(3))])
            clls.append(cll)
            cllens.append(len(cll))

        dcl._clusters=clls
        dcl._clen=cllens
           
        return dcl,headNode

        
if __name__=='__main__':
    from MolKit2 import Read
    rec = Read('CDK2/receptors/4EK3_rec.pdbqt')  
    gc = CalculateAD4Grids(rec, (25.834, 27.584, 27.528), (70,70,70), ['C'],
                           flexibleResidues=[
                               ('A', [('ILE','10'), ('VAL', '18'), ('LYS', '33'), ('VAL', '64'),
                                      ('PHE', '80'), ('PHE', '82'), ('GLN', '85'), ('ASP', '86'),
                                      ('LYS', '89'), ('ASN', '132'), ('LEU', '134'),
                                      ('ASP', '145')]
                                )])
    
    #rec = Read('CDK2/receptors/2CCH_rec.pdbqt')
    #gc = CalculateAD4Grids(rec, (25.834, 27.584, 27.528), (70,70,70), ['A', 'C', 'HD', 'N', 'NA', 'OA', 'P'],
    #                       flexRes=[('ILE','10'), ('VAL', '18'), ('LYS', '33'), ('VAL', '64'),
    #                                ('PHE', '80'), ('PHE', '82'), ('GLN', '85'), ('ASP', '86'),
    #                                ('LYS', '89'), ('ASN', '132'), ('LEU', '134'), ('ASP', '145')])
    #gc = CalculateAD4Grids(rec, (25.834, 27.584, 27.528), (70,70,70), ['C'],
    #                       flexRes=[('ILE','10'), ('VAL', '18'), ('LYS', '33'), ('VAL', '64'),
    #                                ('PHE', '80'), ('PHE', '82'), ('GLN', '85'), ('ASP', '86'),
    #                                ('LYS', '89'), ('ASN', '132'), ('LEU', '134'), ('ASP', '145')])
    gc.getTPoints()
    print len(gc._coords)
    
