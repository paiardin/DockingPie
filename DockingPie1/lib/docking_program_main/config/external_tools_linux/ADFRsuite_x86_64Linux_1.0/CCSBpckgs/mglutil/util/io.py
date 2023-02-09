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

########################################################################
#
# Date: 2015 Authors: Michel Sanner
#
#    sanner@scripps.edu
#
#       The Scripps Research Institute (TSRI)
#       Molecular Graphics Lab
#       La Jolla, CA 92037, USA
#
# Copyright: Michel Sanner and TSRI 2015
#
#########################################################################
#
# $Header: /mnt/raid/services/cvs/python/packages/share1.5/mglutil/util/io.py,v 1.1.4.1 2017/07/26 22:31:38 annao Exp $
#
# $Id: io.py,v 1.1.4.1 2017/07/26 22:31:38 annao Exp $
#
class Stream:
    def __init__(self):
        self.lines = []
    def write(self, line):
        self.lines.append(line)

# helper class to make stdout set of lines look like a file that ProDy can parse
class BufferAsFile:
    def __init__(self, lines):
        self.lines = lines
    def readlines(self):
        return self.lines
