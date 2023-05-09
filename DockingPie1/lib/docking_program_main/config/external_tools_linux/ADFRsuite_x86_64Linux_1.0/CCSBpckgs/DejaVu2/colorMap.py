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
# Copyright: M. Sanner TSRI 2000
#
# Revision: Guillaume Vareille
#
#############################################################################
#
# $Header: /mnt/raid/services/cvs/DejaVu2/colorMap.py,v 1.1.1.1.4.1 2017/07/13 22:28:33 annao Exp $
#
# $Id: colorMap.py,v 1.1.1.1.4.1 2017/07/13 22:28:33 annao Exp $
#

import numpy
import os
import warnings
import types

from DejaVu2 import viewerConst
from DejaVu2.colorTool import RGBRamp, ToHSV, TkColor

class ColorMap:
    """A Colormap is an object designed to provide a mapping from a 
property value to an entry in a 'ramp', which is a list of 
lists. Each entry defines a color for that value.  The property 
range is set by the attributes 'mini' to 'maxi'.  The corresponding 
colors are obtained from the Colormap's ramp, which is n by 4 
matrix of floats: rgba or red, green, blue, alpha values.
If it is initialized with only rgb, the alpha values are added.

A Colormap has always a name and a ramp, if these values are not provided
the creator will try to provide both if they are missing

A Colormap may be created by reading a file and may saved in a file.
   
A Colormap can return its ramp as hsva or as TkColors via its 
methods 'asHSV' and 'asTKCOL'

A Colormap can be associated to a ColormapGUI object.
"""
    def __init__(self, name=None, ramp=None, labels=None, 
                 mini='not passed', maxi='not passed',
                 cmap=None, filename=None, **kw):

        self.name = None
        self.ramp = None
        self.labels = None
        self.mini = None
        self.maxi = None
        self.lastMini = None
        self.lastMaxi = None

        if cmap is not None:
            assert filename is None, \
                "color map creation with both cmap and filename doesn't make sense"
            assert ramp is None, \
                "color map creation with both cmap and ramp doesn't make sense"
            if name is not None:
               warnings.warn("color map creation with cmap and name: name won't be used")
            if mini is not None:
               warnings.warn("color map creation with cmap and mini: mini won't be used")
            if maxi is not None:
               warnings.warn("color map creation with cmap and maxi: maxi won't be used")
            name = cmap.name
            ramp = cmap.ramp
            labels = cmap.labels
            mini = cmap.mini
            maxi = cmap.maxi
            ColorMap.configure(self, name=name, ramp=ramp, labels=labels, mini=mini, maxi=maxi)
        elif filename is not None:
            assert ramp is None, \
                "color map creation with both ramp and filename doesn't make sense"
            if ColorMap.read(self, filename) is False:
                # bad file, but now the colormap is instanced, 
                # we configure it with good values
                #ColorMap.__init__(self, name=name, ramp=ramp, labels=labels, mini=mini, maxi=maxi)
                ColorMap.configure(self, name=name, ramp=ramp, labels=labels, mini=mini, maxi=maxi)
                #but we still raise an error
                raise RuntimeError('invalid Colormap file')
                #return
            else:   
                if name is not None:
                    self.name = name
                elif self.name is None: 
                    # filename without path and without extension
                    self.name = filename.split(os.sep)[-1].split('.')[0]
        else:
            if ramp is None:    
                ramp = RGBRamp(size=16)
                if name is None:
                    name = 'RGBRamp_Colormap'
            if name is None:
                name = 'Untitled_Colormap'
            ColorMap.configure(self, name=name, ramp=ramp, labels=labels, mini=mini, maxi=maxi)


    def configure(self, name=None, ramp=None, labels=None, mini='not passed', maxi='not passed', **kw):
        """Configure the colormap with the given values.
"""
        if name is None:
            assert self.name is not None, "colormap.configure() needs a name"
        else:
            self.name = name
            assert name != ''

        if ramp is None:
            assert self.ramp is not None, "colormap.configure() needs a ramp"
        else:
            assert len(ramp) > 0, \
                "colormap.configure() needs a ramp with at least one element"
            self.ramp = self.checkRamp(ramp)[:]

        if hasattr(self.ramp, 'tolist'):
            self.ramp = self.ramp.tolist()

        assert type(self.ramp) is types.ListType, type(self.ramp)

        lenRamp = len(self.ramp)
        if labels is not None:
            if len(labels) == 0:
                # passing an empty list is a shortcut for a numeric list
                self.labels = range(lenRamp)
            else:
                assert len(labels) == lenRamp
                self.labels = labels
        elif self.labels is not None:
            lenLabels = len(self.labels)
            if lenLabels == 0:
                self.labels = range(lenRamp)
            elif lenLabels < lenRamp:
                lMissingLength = lenRamp - lenLabels
                lMissingLabels = range(lenLabels, lenRamp)
                #print "lMissingLabels", lMissingLabels
                self.labels += lMissingLabels
            else:
                while len(self.labels) > lenRamp:
                    self.labels.pop()
            assert len(self.labels) == lenRamp

        if self.labels is not None:
            for i in range( len(self.labels) ):
                self.labels[i] = str(self.labels[i])

        if mini != 'not passed':
            if mini is None:
                self.mini = None
            else:
                if maxi is not None and (maxi != 'not passed'):
                    if mini >= maxi: 
                        mini = None
                elif self.maxi is not None:
                    if mini >= self.maxi:
                        mini = None
                self.mini = mini
            self.lastMini = self.mini

        if maxi != 'not passed':
            if maxi is None:
                self.maxi = None
            else:
                if mini is not None and (mini != 'not passed'):
                    if mini >= maxi:
                        maxi = None
                elif self.mini is not None:
                    if self.mini >= maxi:
                        maxi = None    
                self.maxi = maxi
            self.lastMaxi = self.maxi


    def set(self, value):
        """set the colormap ramp with the given values.
"""
        if isinstance(value, ColorMap):
            self.configure(name=value.name, ramp=value.ramp, labels=value.labels, 
                           mini=value.mini, maxi=value.maxi)
        else:
            self.configure(ramp=value)


    def get(self):
        """return the colormap ramp
"""
        return self.ramp


    def getDescr(self):
        cfg = {}
        cfg['name'] = self.name            
        cfg['ramp'] = self.ramp
        cfg['labels'] = self.labels
        cfg['mini'] = self.mini
        cfg['maxi'] = self.maxi
        return cfg
        

    def checkRamp(self, ramp):
        """Method to check the given ramp. If only rgb values given then 1
is added for the alpha values.
"""
        #if no alpha values, add 1's
        if len(ramp[0]) == 4:
            return ramp
        if len(ramp[0]) == 3:
            lenRgb = len(ramp)
            _ones = numpy.ones(lenRgb, 'f')
            _ones.shape = (lenRgb, 1)
            ramp = numpy.concatenate((numpy.array(ramp), _ones),1)
            return ramp.tolist()

        
    def read(self, fileName):
        """Reinitialize colormap with data from file
"""
        l = {}
        g = {}
        try:
            execfile(fileName, g, l)
        except:
            return False
        
        cm = None
        for name, object in l.items():
            if isinstance(object, ColorMap):
                cm = object
                break
        else: # we didn't break
            return False
        
        if cm.maxi > cm.mini:
            cfg = {'name':cm.name, 'ramp':cm.ramp, 'mini':cm.mini, 'maxi':cm.maxi}
        else:            
            cfg = {'name':cm.name, 'ramp':cm.ramp}

        if hasattr(cm,'labels'):
            cfg['labels'] = cm.labels

        apply( ColorMap.configure, (self,), cfg)

        if self.ramp is None:
            return False
        else:
            return True


    def write(self, fileName):
        """Write the colormap's source code to a file
"""
        code = self.sourceCode()
        f = open(fileName, 'w')
        f.write("%s"%code)
        f.close()
        
                        
    def sourceCode(self):
        """returns python code that will recreate this object
"""
        cm = self
        code = """from %s import %s\n"""%(ColorMap.__module__, 'ColorMap')
        code = code + """cm = %s(name='%s')\n"""%('ColorMap', cm.name)
        cfg = cm.getDescr()
        code = code + "cfg = "+str(cfg)+"""\n"""
        code = code + """apply( cm.configure, (), cfg)\n"""
        return code


    def Map(self, values, mini='not passed', maxi='not passed'):
        """Get colors corresponding to values in a colormap.
if mini or maxi are provided, self.mini and self.maxi are not
used and stay inchanged.
if mini or maxi are not provided or set to None, 
self.mini and self.maxi are used instead.
if mini or maxi are set to None, self.mini and self.maxi are
ignored
"""
        #print "Map", mini, maxi, self.mini, self.maxi, values
        values = numpy.array(values)

        if len(values.shape)==2 and values.shape[1]==1:
            values.shape = ( values.shape[0], )
        elif len(values.shape) > 1:
            raise ValueError('ERROR! values array has bad shape')

        ramp = numpy.array(self.ramp)
        if len(ramp.shape) !=2 or ramp.shape[1] not in (3,4):
            raise ValueError('ERROR! colormap array has bad shape')

        if mini == 'not passed':
             mini = self.mini
        if maxi == 'not passed':
             maxi = self.maxi

        if mini is None:
            mini = min(values)
        else:
            # All the values < mini will be set to mini
            values = numpy.maximum(values, mini)
            
        if maxi is None:
            maxi = max(values)
        else:
            values = numpy.minimum(values, maxi)       

        # mini and maxi are now set
        if mini >= maxi:
            txt = 'mini:%f must be < maxi:%f'%(mini, maxi)
            warnings.warn( txt )

        valrange = maxi - mini
        if valrange < 0.0001:
            ind = numpy.ones( values.shape )
        else:
            colrange = ramp.shape[0]
            ind = (values - mini) * (colrange/float(valrange))
            ind = ind.astype(viewerConst.IPRECISION)
            ind = numpy.minimum(ind, colrange - 1)

        col = numpy.take(self.ramp, ind, axis=0)
        
        self.lastMini = mini
        self.lastMaxi = maxi

        return col


    #alternative representations of the ramp
    def asHSV(self, redo=0):
        self.hsv = map(lambda x, conv=ToHSV: conv(x), self.ramp)
        #if redo or not hasattr(self, 'hsv'):
        #    self.hsv = map(lambda x, conv=ToHSV: conv(x), self.ramp)
        return self.hsv


    def asHSL(self, redo=0):
        from mglutil.util.colorUtil import RGBA2HSLA_list
        return map(lambda x, conv=RGBA2HSLA_list: conv(x), self.ramp)


    def asTKCOL(self, redo=0):
        if redo or not hasattr(self, 'tkcol'):
            self.tkcol = map(lambda x, conv=TkColor: conv(x), self.ramp)
        return self.tkcol


    def _lookup(self, name):
        #name = str(name)
        if self.labels is None:
            for i in range( len(self.ramp) ):
                self.labels[i] = str(i)

        if type(name) != types.StringTypes:
            name = str(name)

        return self.ramp[self.labels.index(name)]


    def lookup(self, names):
        return map( self._lookup, names)
