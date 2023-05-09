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
# $Header: /mnt/raid/services/cvs/ADFR/utils/analyze.py,v 1.2 2016/12/07 00:38:33 sanner Exp $
#
# $Id: analyze.py,v 1.2 2016/12/07 00:38:33 sanner Exp $
#
from ADFRcc import getFFParameters
parameters = getFFParameters()
feCoeffHbond = parameters.feCoeffHbond

def getHBPairs(ind, ligAtoms, cutOffEne=0.0):
    scorer = ind.genomePy.scorer.getLlPairwiseScorer()
    hArray = scorer.getHbEnergyMatrix()
    atomSet = scorer.getAtomSet1().atomSetStatic
    hbPairs = []
    hbEne = []
    for i, a1 in enumerate(ligAtoms):
        n1 = a1.getName()
        for j, a2 in enumerate(ligAtoms):
            if j<=i: continue
            n2 = a2.getName()
            if atomSet.getPairScorable(i,j):
                h = hArray[i][j]*feCoeffHbond
                if h < cutOffEne:
                    hbPairs.append( (i, j) )
                    hbEne.append( h)
    return hbPairs, hbEne

def addHBlines(viewer, hbPairs, hbEne, coords):
    from DejaVu2.IndexedPolylines import IndexedPolylines
    hbLines = IndexedPolylines('hbonds', visible=1,
                               vertices=coords,
                               faces=hbPairs,
                               materials=[[1,1,0]], inheritMaterial=False,
                               stippleLines=True, inheritStippleLines=False)
    viewer.AddObject(hbLines)
    
    from DejaVu2.glfLabels import GlfLabels
    pts = []
    labels = []
    for n in range(len(hbPairs)):
        i,j = hbPairs[n]
        labels.append('%.5f'%hbEne[n])
        pts.append(0.5*(coords[i] + coords[j]))
    import numpy
    g = GlfLabels('hbEneLablel', fontStyle='solid3d',
                  inheritMaterial=1,
                  materials=[([.8,.8,0])],
                  vertices=pts, labels=labels,
                  fontTranslation=(0,0,.1),
                  fontScales=(0.1, 0.1, 0.05),
                  pickable=0,
                  labelTranslation=numpy.zeros( (len(pts), 3), 'f')
                  )
    viewer.AddObject(g)
