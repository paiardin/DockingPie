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
# Copyright: M. Sanner and TSRI 2015
#
#########################################################################
#
# $Header: /mnt/raid/services/cvs/ADFR/bin/about.py,v 1.4 2016/12/07 00:38:32 sanner Exp $
#
# $Id: about.py,v 1.4 2016/12/07 00:38:32 sanner Exp $
#
#
import os, sys
from ADFR.utils.maps import MapsFile
from ADFR.dro import DockingResultsObject

if len(sys.argv)<2:
    print "usage: about filename"
    sys.exit(1)

filename = sys.argv[1]
name, ext = os.path.splitext(filename)

if ext=='.trg':
    mf = MapsFile(filename)
    mf.printInfo()

elif ext=='.dro':
    d = DockingResultsObject()
    d.load(filename)
    d.printInfo()

else:
    print 'Nothing known about %s files'%ext
    sys.exit(1)

