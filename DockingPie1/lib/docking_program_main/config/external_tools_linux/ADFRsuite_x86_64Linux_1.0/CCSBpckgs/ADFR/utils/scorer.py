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
# $Header: /mnt/raid/services/cvs/ADFR/utils/scorer.py,v 1.4 2016/12/07 00:38:33 sanner Exp $
#
# $Id: scorer.py,v 1.4 2016/12/07 00:38:33 sanner Exp $
#
from ADFR import ADFR

def ADFRscorer(ligand, mapsFolder, mapFilesRoot=None,
               receptor=None, flexRes=None, tpointsFilename=None,
               covalentIndices=None):

    adfr = ADFR(ligand, mapsFolder, receptor=receptor,
                flexibleResidues=flexRes, mapFilesRoot=mapFilesRoot,
                covalentIndices=covalentIndices)

    from ADFR.GA import GenomePy, IndividualPy, Population
    if adfr.RFTGen:
        scaleRE= 1.0/len(adfr.RFTGen.motions)
    else:
        scaleRE=None
    
    genome = GenomePy(adfr.FT, adfr.scorer, scaleRE=scaleRE)
    ind = IndividualPy(genome)
    ind.setGenes(ind.genomePy.getIdentityGenesPy())
    _score = ind.score()
    return adfr, ind, _score

