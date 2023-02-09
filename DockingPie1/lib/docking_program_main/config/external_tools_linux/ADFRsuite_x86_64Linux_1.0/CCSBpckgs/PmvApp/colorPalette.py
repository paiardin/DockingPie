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
# Date: April 2006 Authors: Guillaume Vareille, Michel Sanner
#
#    vareille@scripps.edu
#    sanner@scripps.edu
#
#       The Scripps Research Institute (TSRI)
#       Molecular Graphics Lab
#       La Jolla, CA 92037, USA
#
# Copyright: Guillaume Vareille, Michel Sanner and TSRI
#
#########################################################################
#
# $Header: /mnt/raid/services/cvs/PmvApp/colorPalette.py,v 1.2.4.1 2017/07/13 20:55:28 annao Exp $
#
# $Id: colorPalette.py,v 1.2.4.1 2017/07/13 20:55:28 annao Exp $
#


from DejaVu2.colorMap import ColorMap

class ColorPalette(ColorMap):
    FLOAT = 0

    def __init__(self, name, colorDict={}, readonly=0, colortype=None,
                 info='', sortedkeys=None, lookupMember=None):

        if len(colorDict) > 0:
            if sortedkeys is None:
                labels = colorDict.keys()
                values = colorDict.values()
            else:
                labels = sortedkeys
                values = []
                for label in labels:
                    values.append(colorDict[label])
        else:
            labels = None
            values = None

        ColorMap.__init__(self, name=name, ramp=values, labels=labels)

        self.readonly = readonly
        self.info = info
        #self.viewer = None
        self.sortedkeys = sortedkeys
        if colortype is None:
            self.colortype = self.FLOAT
        self.lookupMember = lookupMember

    def _lookup(self, name):
        #print "_lookup", name, type(name)
        try:
            col = ColorMap._lookup(self, name)
            if len(col) == 4:
                return col[:3]
            else:
                return col
        except:
            return (0., 1., 0.)


    def lookup(self, objects):
        # Maybe should try that first in case all the objects don't have the
        # lookup member
        names = objects.getAll(self.lookupMember)
        return map( self._lookup, names)


class ColorPaletteFunction(ColorPalette):

    def __init__(self, name, colorDict={}, readonly=0, colortype=None,
                 info='', sortedkeys=None, lookupFunction = None):
        """ lookupFunction : needs to be function or a lambda function"""
        ColorPalette.__init__(self, name, colorDict, readonly,colortype,
                               info, sortedkeys)
        from types import FunctionType
        if not type(lookupFunction) is FunctionType:
            self.lookupFunction = None

        self.lookupFunction = lookupFunction
                     
          
    def lookup(self, objects):
        # maybe should do that in a try to catch the exception in case it
        # doesnt work
        names = map(self.lookupFunction, objects)
        return map(self._lookup, names)


