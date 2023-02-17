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
# Copyright: M. Sanner TSRI 2015
#
#
#############################################################################

#
# $Header: /mnt/raid/services/cvs/MolKit2/tree.py,v 1.1.1.1.4.1 2017/07/26 22:03:40 annao Exp $
# 
# $Id: tree.py,v 1.1.1.1.4.1 2017/07/26 22:03:40 annao Exp $
#

class TreeObject:
    ## base class for obejct thata can be added to a tree widget
    
    ## def __init__(self):
    ##     self.selSet = None

    def getChildrenAndNames(self):
        #Return a list of children which all are DashboardTreeObject instances
        return [], []

    def nbChildren(self):
        return len(self.selSet)
