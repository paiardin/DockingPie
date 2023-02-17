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
# $Header: /mnt/raid/services/cvs/AutoSite/ASfeaturepoints.py,v 1.3 2016/12/07 00:38:33 sanner Exp $
#
# $Id: ASfeaturepoints.py,v 1.3 2016/12/07 00:38:33 sanner Exp $
#
import numpy

class featurePts:       
    def __init__(self, receptor, ccoords, cpot, ocoords, opot, hcoords, hpot,scale=1):
        self.rec = receptor
        self.cc = ccoords
        self.cp = cpot
        self.oc = ocoords
        self.op = opot
        self.hc = hcoords
        self.hp = hpot
        self.scale = scale
        self.getFeaturePts()


    def getFeaturePts(self):
        from AutoSite.utils.clusterTPoints import Kmeans
        nbc = len(self.cc)
        #import pdb
        #pdb.set_trace()
        self.cfp, self.cfpp,self.cBres = self.recPharmacophores(self.cc, self.cp, ['NA', 'OA', 'SA','HD','C','A','N','O','S','Br','BR','F','Cl','CL'], 4.0)
        if nbc > 14*self.scale:
            ###cfpK = Kmeans(self.cc, int(round(nbc/14.3))*self.scale, returnCentroids=True)
            cfpK = Kmeans(self.cc, int(round(nbc/(14.3*self.scale))), returnCentroids=True)
            self.cfp = cfpK.cCentroids
            #import pdb;pdb.set_trace()
            self.cfpp = [-0.43]*len(self.cfp)

        else:
            self.cfp = self.cc
            self.cfpp = self.cp


        self.ofp,self.ofpp,self.oBres= self.recPharmacophores(self.oc, self.op, ['HD'], 2.5)

        self.hfp, self.hfpp, self.hBres= self.recPharmacophores(self.hc, self.hp, ['NA', 'OA', 'SA'], 2.5)

        self.bsres = set(self.cBres+self.hBres+self.oBres)
        #return cfp, ofp, hfp, cfpp, ofpp, hfpp,set(cBres+hBres+oBres)
        
    def recPharmacophores(self, oc, op, recat, dcutoff=2.5):
        """
        Returns the coordinate of the fill point with minimum potential for every unique
        receptor atom picked (of provided atom types) at the binding site.
        """
        from bhtree import bhtreelib
        recAtoms = self.rec._ag
        recCoords = recAtoms.getCoords().tolist()
        rbht = bhtreelib.BHtree( recCoords, None, 10)
        pCcenter = []
        pCres = []
        recPdict = {}
        for i,(pC,pP) in enumerate(zip(oc,op)):
            results = numpy.zeros(5000, 'i')
            dist2 = numpy.zeros(5000, 'f')
            nb = rbht.closePointsDist2(tuple(pC), dcutoff, results, dist2)
            if nb:                
                for ind, d2 in zip(results[:nb],dist2[:nb]):
                    if recAtoms.getData('AD_element')[ind] in recat:
                        resname = '%s:%s%d'%(recAtoms.getChids()[ind],recAtoms.getResnames()[ind], recAtoms.getResnums()[ind])
                        recPdict.setdefault(ind, [])
                        recPdict[ind].append((pC,pP))
                        if resname not in pCres:
                            pCres.append(resname)
        pC = []
        pP = []
        for k,v in recPdict.iteritems():
            sortrPclusters = sorted(v, key=lambda x:x[1])
            pC.append(sortrPclusters[0][0])
            pP.append(sortrPclusters[0][1])
        return pC, pP, pCres
        

 
