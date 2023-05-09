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
# $Header: /opt/cvs/AutoSite/scoreClusters.py,v 1.3 2017/04/13 00:52:43 annao Exp $
#
# $Id: scoreClusters.py,v 1.3 2017/04/13 00:52:43 annao Exp $
#
import numpy

from math import sqrt
from ADFRcc.adfr import GridMap
from AutoSite.fillBuriedness import Buriedness

def scoreClusters(rec, dcl, gc, inflate=False, pepScore=False):
    coords = rec._ag.getCoords()
    radii = rec._ag.getRadii()
    clusters = []
    clProp = []
    if inflate==True:
        cmapFile = gc.mapFiles['C']
        cmap = GridMap()
        cmap.loadFromMapFile('C', '', cmapFile)
        omapFile = gc.mapFiles['OA']
        omap = GridMap()
        omap.loadFromMapFile('OA', '', omapFile)
        hmapFile = gc.mapFiles['HD']
        hmap = GridMap()
        hmap.loadFromMapFile('HD', '', hmapFile)
        cdata = cmap.getGridDataPy()#[::3, ::3, ::3]
        odata = omap.getGridDataPy()
        hdata = hmap.getGridDataPy()
    gridSpacing = 0
    ox = oy = oz = 0
    for clinds in dcl._clusters:
        cli = gc._indices[clinds]
        clc = gc._coords[clinds]
        clp = gc._potential[clinds]
        cla = gc._atype[clinds]
        if gridSpacing==0:
            gridSpacing = sqrt((clc[0][0]-clc[1][0])*(clc[0][0]-clc[1][0])+(clc[0][1]-clc[1][1])*(clc[0][1]-clc[1][1])+(clc[0][2]-clc[1][2])*(clc[0][2]-clc[1][2]))/sqrt((cli[0][0]-cli[1][0])*(cli[0][0]-cli[1][0])+(cli[0][1]-cli[1][1])*(cli[0][1]-cli[1][1])+(cli[0][2]-cli[1][2])*(cli[0][2]-cli[1][2]))
            ox = clc[0][0]-cli[0][0]*gridSpacing
            oy = clc[0][1]-cli[0][1]*gridSpacing
            oz = clc[0][2]-cli[0][2]*gridSpacing

        if len(cli)==1:
            metric = Rg = fburied = 0.
        else:
            #import pdb;pdb.set_trace()
            dmetric = Buriedness(clc, cla, coords.tolist(), radii,origin=[ox,oy,oz] )
            if inflate==True:
                fburied,inflateIndices = dmetric.inflatePocketandBuriedness(origin=[ox,oy,oz],spacing=1.0) 
                for ind,coord in inflateIndices.iteritems():
                    #import pdb;pdb.set_trace()
                    if ind[0]>=odata.shape[0] or ind[1]>=odata.shape[1] or ind[2]>=odata.shape[2]:
                        continue
                    clc = numpy.vstack((clc, coord))
                    cli = numpy.vstack((cli, ind))
                    potential = 0.0
                    potential = min(hdata[ind[0]][ind[1]][ind[2]],cdata[ind[0]][ind[1]][ind[2]],odata[ind[0]][ind[1]][ind[2]],0.0)
                    clp = numpy.append(clp, potential)
                    cla = numpy.append(cla, 'C')
            else:
                fburied = dmetric.NumericalBurriedness(origin=[ox,oy,oz],spacing=1.0)
            cx,cy,cz = numpy.mean(clc,axis=0)
            dst2 = 0.0
            for px,py,pz in clc:
                dst2 += (px-cx)*(px-cx)+(py-cy)*(py-cy)+(pz-cz)*(pz-cz)
            complen = len(clc)
            Rg2 = (dst2/complen)
            Rg = sqrt(Rg2)

            #import pdb;pdb.set_trace()
            if not pepScore:
                metric = (len(clc)*fburied*fburied)/Rg
            else:
                metric = (len(clc)*fburied*sqrt(fburied))

        clusters.append( [cli, clc, clp, cla, metric] )
        clProp.append([numpy.sum(clp), len(clc), numpy.sum(clp)/len(clc), Rg, fburied, metric])
    return clusters, clProp

