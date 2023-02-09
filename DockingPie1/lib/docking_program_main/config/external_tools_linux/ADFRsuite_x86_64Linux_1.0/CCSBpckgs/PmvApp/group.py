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
# Copyright: M. Sanner TSRI 2016
#
#########################################################################
#
# $Header: /mnt/raid/services/cvs/PmvApp/group.py,v 1.1.4.1 2017/07/13 20:52:52 annao Exp $
#
# $Id: group.py,v 1.1.4.1 2017/07/13 20:52:52 annao Exp $
#
class Group:

    def __init__(self, name, locked=False, parent=None, validTypes=None):
        self._basename = name
        self.name = name
        self._children = []
        self._locked = locked # set to True to prevent adding or removing 
                              # from/to this group
        self._validTypes = validTypes # None or a set of classes of classes
                                      # acceptable for this group
        self._group = parent # called _group to be consistent wit mol._group
        self._multi = None
        self._currentIndex = 0
        
    def numMols(self):
        return max(0, max(
            [x.numMols() for x in self._children if hasattr(x, '_multi')]))

    def curMolIndex(self):
        return self._currentIndex
    
    def validChild(self, child):
        if self._validTypes is None:
            return True
        else:
            return isinstance(child, tuple(self._validTypes))

    def nbChildren(self):
        return len(self._children)
    
    def getChildrenAndNames(self):
        #Return a list of children which all are DashboardTreeObject instances
        return self._children, [x.name for x in self._children]

    def add(self, child, atIndex=None):
        if not self.validChild(child):
            return False
        if atIndex is None:
            self._children.append(child)
        else:
            self._children.insert(atIndex, child)
        child._group = self
