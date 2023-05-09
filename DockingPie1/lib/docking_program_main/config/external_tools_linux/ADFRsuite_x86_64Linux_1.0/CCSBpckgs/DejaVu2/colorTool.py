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
#############################################################################
#
# $Header: /mnt/raid/services/cvs/DejaVu2/colorTool.py,v 1.1.1.1.4.1 2017/07/13 22:28:33 annao Exp $
#
# $Id: colorTool.py,v 1.1.1.1.4.1 2017/07/13 22:28:33 annao Exp $
#


import numpy
import viewerConst
from opengltk.extent import _gllib as gllib

from math import fabs
from opengltk.OpenGL import GL
# properties are GL_AMBIENT, GL_DIFFUSE, GL_SPECULAR, GL_EMISSION, GL_SHININESS
materialMemory = {
    GL.GL_FRONT: [[-1.,-1.,-1.,-1], [-1.,-1.,-1.,-1], [-1.,-1.,-1.,-1],
                  [-1.,-1.,-1.,-1], -1],
    GL.GL_BACK: [[-1.,-1.,-1.,-1], [-1.,-1.,-1.,-1], [-1.,-1.,-1.,-1],
                 [-1.,-1.,-1.,-1], -1]
    }

from math import fabs

def resetMaterialMemory():
    global materialMemory
    #print 'RESET'
    for i in range(4):
        for j in range(4):
            materialMemory[GL.GL_FRONT][i][j] = -1
            materialMemory[GL.GL_BACK][i][j] = -1
        materialMemory[GL.GL_FRONT][4] = -1
        materialMemory[GL.GL_BACK][4] = -1


def glAllMaterialWithCheck(face, mat, num, eps=0.001, check=True, debug=0):
    p = mat.prop
    f = face
    if mat.binding[0] == viewerConst.OVERALL:
        glMaterialWithCheck(f, GL.GL_AMBIENT, p[0][0], eps=eps, check=check)
    else:
        if debug:
            print 'ambi', p[0][num]
        glMaterialWithCheck(f, GL.GL_AMBIENT, p[0][num], eps=eps, check=check)

    if mat.binding[1] == viewerConst.OVERALL:
        glMaterialWithCheck(f, GL.GL_DIFFUSE, p[1][0], eps=eps, check=check)
    else:
        if debug:
            print 'diff', p[1][num]
        glMaterialWithCheck(f, GL.GL_DIFFUSE, p[1][num], eps=eps, check=check)

    if mat.binding[2] == viewerConst.OVERALL:
        glMaterialWithCheck(f, GL.GL_SPECULAR, p[2][0], eps=eps, check=check)
    else:
        if debug:
            print 'spec', p[2][num]
        glMaterialWithCheck(f, GL.GL_SPECULAR, p[2][num], eps=eps, check=check)

    if mat.binding[3] == viewerConst.OVERALL:
        glMaterialWithCheck(f, GL.GL_EMISSION, p[3][0], eps=eps, check=check)
    else:
        if debug:
            print 'emis', p[3][num]
        glMaterialWithCheck(f, GL.GL_EMISSION, p[3][num], eps=eps, check=check)

    if mat.binding[4] == viewerConst.OVERALL:
        glMaterialWithCheck(f, GL.GL_SHININESS, p[4][0], eps=eps, check=check)
    else:
        if debug:
            print 'shini', p[4][num]
        glMaterialWithCheck(f, GL.GL_SHININESS, p[4][num], eps=eps, check=check)

    
def glMaterialWithCheck(face, property, material, eps=0.001, check=True):
    """
    Only calls glMaterial if the material is different from he current value
    face can be GL_FRONT or GL_BACK
    propety can be GL_AMBIENT, GL_DIFFUSE, GL_SPECULAR or GL_EMISSION
    material is a 3-sequence RGBA, Alpha values are only test for ambient and
    diffuse properties
    """
    global materialMemory

    #print 'glMaterialWithCheck', property, material
    
    if face==GL.GL_FRONT_AND_BACK:
        face=GL.GL_FRONT
    if property==GL.GL_SHININESS:
        matMem = materialMemory[face]
        if not check or fabs(matMem[4]-material) > eps:
            matMem[4] = material
            #GL.glMaterialfv( face, property, matMem[4] )
            GL.glMaterialf( face, property, float(matMem[4]) )
    else:
        if not material.flags.contiguous:
            material = numpy.array(material,copy=1)
        propNum = viewerConst.propNum[property]
        matMem = materialMemory[face][propNum]
##         print 'DIFFERENCE'
##         print id(matMem), id(materialMemory[face][propNum])
##         print matMem
##         print material
##         print fabs(matMem[0]-material[0]) > eps
##         print fabs(matMem[1]-material[1]) > eps
##         print fabs(matMem[2]-material[2]) > eps
##         print (fabs(matMem[3]-material[3]) > eps and propNum in (0,1))
        if check:
            newCol = fabs(matMem[0]-material[0]) > eps or \
                     fabs(matMem[1]-material[1]) > eps or \
                     fabs(matMem[2]-material[2]) > eps or \
                     (fabs(matMem[3]-material[3]) > eps and propNum in (0,1))
        else:
            newCol = True
            
        if newCol:
            #GL.glMaterialfv( face, property, material.tolist() )
            #print 'SETTING'
            gllib.glMaterialfv( face, property, material)
            matMem[0] = material[0]
            matMem[1] = material[1]
            matMem[2] = material[2]
            matMem[3] = material[3]
            #print matMem
##         else:
##             curcol = GL.glGetMaterialfv(face, property)
##             print 'SKIPPING'
##             print curcol
##             print material
##             print matMem
            
#def firstMaterial( r, g, b ):
#    global mat_r, mat_g, mat_b
#    mat_r = r
#    mat_g = g
#    mat_b = b
#    
#def sameMaterial( r, g, b ):
#    global mat_r, mat_g, mat_b
#    if fabs(mat_r-r) > 0.001 or fabs(mat_g-g) > 0.001 or fabs(mat_b-b) > 0.001:
#        #print r,mat_r,g,mat_g,b,mat_b
#        mat_r = r
#        mat_g = g
#        mat_b = b
#        return 0
#    return 1


def OneColor(color, alpha = 1.0):
    """get a color RGBA Add alpha if missing, return values from 0.0 to 1.0"""

    if max(color) > 1.0:
        color = map( lambda x: x/255., color)
    if len(color) == 3:
	return list(color)+[alpha]
    elif len(color)==4:
	return color
    else:
	raise ValueError( 'Color has to be a 3 or 4 tuple (RGBA)' )


def TkColor(col):
    """col should be a rgb triplet of int 0-255 or 0-1"""
    if max(col)<=1.0: col = map( lambda x: x*255, col)
    return '#%02X%02X%02X' % (col[0],col[1],col[2])


def ToRGB(hsv):
    """Convert an HSV triplet into an RGB triplet.
    Values range from 0.0 to 1.0. Alpha values are optional"""

    l = len(hsv)
    assert l in (3,4)
    assert max(hsv) <= 1.0
    assert min(hsv) >= 0.0
    v = hsv[2]
    if v == 0.0:
        if l==3: return (0.0, 0.0, 0.0)
        else: return (0.0, 0.0, 0.0,hsv[3])
    s = hsv[1]
    if s == 0.0:
        if l==3: return (v, v, v)
        else: return (v, v, v,hsv[3])
    h = hsv[0]*6.0
    if h>=6.0: h = 0.0
    i = int(h)
    f = h - i
    p = v*(1.0 - s)
    q = v*(1.0-(s*f))
    t = v*(1.0-s*(1.0-f))

    if i==0:
	if l==3: return (v,t,p)
	else: return (v,t,p,hsv[3])
    elif i==1:
	if l==3: return (q,v,p)
	else: return (q,v,p,hsv[3])
    elif i==2:
        if l==3: return (p,v,t)
	else: return (p,v,t,hsv[3])
    elif i==3:
        if l==3: return (p,q,v)
	else: return (p,q,v,hsv[3])
    elif i==4:
        if l==3: return (t,p,v)
	else: return (t,p,v,hsv[3])
    elif i==5:
        if l==3: return (v,p,q)
	else: return (v,p,q,hsv[3])
    else:
        print "botch in hsv_to_rgb"


def ToHSV(rgb):
    """Convert an RGB triplet into an HSV triplet.
       Values range from 0.0 to 1.0. Alpha values are optional"""
    l = len(rgb)
    assert l in (3,4)
    assert max(rgb) <= 1.0
    assert min(rgb) >= 0.0
    r = rgb[0];
    g = rgb[1];
    b = rgb[2];
    maxi = max(rgb[:3])
    mini = min(rgb[:3])

    if maxi > 0.0001: s = (maxi - mini)/maxi
    else: s = 0.0
    if s < 0.0001: h = 0.0
    else:
	delta = maxi - mini
	if r == maxi: h = (g - b)/delta
	elif g == maxi: h = 2.0 + (b - r)/delta
	elif b == maxi: h = 4.0 + (r - g)/delta
	h = h/6.0
	if h < 0.0: h = h + 1.0

    if l==3: return (h,s,maxi)
    else: return (h,s,maxi,rgb[3])


def Map(values, colorMap, mini=None, maxi=None):
    """Get colors corresponding to values in a colormap"""

    values = numpy.array(values)
    if len(values.shape)==2 and values.shape[1]==1:
	values.shape = ( values.shape[0], )
    elif len(values.shape) > 1:
	print 'ERROR: values array has bad shape'
	return None

    cmap = numpy.array(colorMap)
    if len(cmap.shape) !=2 or cmap.shape[1] not in (3,4):
	print 'ERROR: colorMap array has bad shape'
	return None

    if mini is None: mini = min(values)
    else: values = numpy.maximum(values, mini)
    if maxi is None: maxi = max(values)
    else: values = numpy.minimum(values, maxi)
    valrange = maxi-mini
    if valrange < 0.0001:
	ind = numpy.ones( values.shape )
    else:
	colrange = cmap.shape[0]-1
	ind = ((values-mini) * colrange) / valrange
    col = numpy.take(colorMap, ind.astype(viewerConst.IPRECISION), axis=0)
    return col


def array2DToImage( array2D, cmap, width=None, 
	height=None, numComponents=4,texture=1, maxi=None, mini=None):
    #build the image:
    if width is None:
        width=array2D.shape[0]
    if height is None:
        height=array2D.shape[1]
    #texture sets the image dimensions to the smallest power of two
    if texture:
        dim1 = dim2=1
        while dim1< width: dim1=dim1<<1
        while dim2< height: dim2=dim2<<1

    colors = Map(array2D.ravel(), cmap,mini=mini, maxi=maxi )

    ###7/19: because these are numpy arrays???!!!!
    ###7/19colors.shape = (width, height,3)
    colors.shape = (height, width,3)
    colors = colors*255
    colors = colors.astype('B')
    tex2Dimage = numpy.ones((dim2,dim1,numComponents), 'B')
    ###7/19: because these are numpy arrays???!!!!
    ###7/19tex2Dimage = numpy.ones((dim1,dim2,numComponents), 'B')
    ###7/19 tex2Dimage[:width,:height,:3] = colors
    tex2Dimage[:height,:width,:3] = colors
    return tex2Dimage

def HSVRamp(size=256, upperValue=.6666666666666667):
    """Generate an HSV color ramp, values range from 0.0 to 1.0
"""
    assert size > 0
    h = []
    v = upperValue
    if size > 1:
        step = v/(size-1)
    else:
        step = v
    for i in range(size):
        h.append(v)
        v -= step
        if v < 0:
            v = 0
    h = numpy.array(h,'f')
    h.shape = (size, 1)
    sv = numpy.ones ( (size, 2), viewerConst.FPRECISION )
    hsv = numpy.concatenate( (h,sv), 1 )
    return hsv


def RGBRamp(size=256, upperValue=.6666666666666667):
    """Generate an RGB color ramp, values range from 0.0 to 1.0"""

    assert size > 0
    hsv = HSVRamp(size, upperValue)
    rgb = numpy.zeros( (hsv.shape[0], 3), viewerConst.FPRECISION )
    for i in xrange(hsv.shape[0]):
        rgb[i] = ToRGB(hsv[i])
    return rgb


def RGBARamp(size=256, upperValue=.6666666666666667):
    """Generate an RGBA color ramp, values range from 0.0 to 1.0"""

    assert size > 0
    hsv = HSVRamp(size, upperValue)
    rgb = numpy.zeros( (hsv.shape[0], 4), viewerConst.FPRECISION )
    for i in xrange(hsv.shape[0]):
        rgb[i][:3] = ToRGB(hsv[i])
        rgb[i][3] = 1.0
    return rgb


def RedWhiteBlueRamp(size=256):
    ramp = numpy.ones( (size, 3), 'f')
    mid = size/2
    incr = 1./(mid-1)
    for i in xrange(mid):
        ramp[i][1] = i*incr
        ramp[i][2] = i*incr
    for i in xrange(mid):
        ramp[mid+i][0] = 1.-(i*incr)
        ramp[mid+i][1] = 1.-(i*incr)
    return ramp

def RedWhiteBlueARamp(size=256):
    ramp = numpy.ones( (size, 4), 'f')
    mid = size/2
    incr = 1./(mid-1)
    for i in xrange(mid):
        ramp[i][1] = i*incr
        ramp[i][2] = i*incr
    for i in xrange(mid):
        ramp[mid+i][0] = 1.-(i*incr)
        ramp[mid+i][1] = 1.-(i*incr)
    return ramp

def colorsForValueRanges(rangeList, size=256):
    """
    create an array of size times RGBS colors giver a list of value ranges
    and associated colors.
    rangList should be a list of 4-tuples in the form of:
        [(mini, maxi, color1, color2=None), ...] where colors are RGBA
    """
    mini = min([x[0] for x in rangeList])
    maxi = max([x[1] for x in rangeList])
    linVal = (maxi-mini)/float(size)
    ramp = numpy.zeros( (size, 4), 'f')

    for m, M, c1, c2 in rangeList:
        fstLine = int(round((m-mini)/linVal))
        lstLine = int(round((M-mini)/linVal))
        h1,s1,v1 = ToHSV(c1[:3])
        if c2:
            h2,s2,v2 = ToHSV(c2[:3])
            hstep = (h2-h1)/(lstLine-fstLine)
            sstep = (s2-s1)/(lstLine-fstLine)
            vstep = (v2-v1)/(lstLine-fstLine)
            for i in range(lstLine-fstLine):
                print i, ToRGB( (h1 + i*hstep, s1 + i*sstep, v1 +i*vstep, 1) )
                ramp[i] = ToRGB( (h1 + i*hstep, s1 + i*sstep, v1 +i*vstep, 1) )
        else:
            for i in range(fstLine, lstLine):
                print i, [c1[0], c1[1], c1[2]]
                ramp[i] = [c1[0], c1[1], c1[2], 1.0]
    return ramp

## had to re-implement this ramp because the version below made python crash
## aftter the ramp was created !
## FIXME all following ramp may be problematic

##  def RedWhiteBlueRamp(size=256):
##      ramp = numpy.ones( (size, 3), 'd')
##      mid = size/2
##      ramp[:mid,1] = numpy.arange(mid) / float(mid-1)
##      ramp[:mid,2] = numpy.arange(mid) / float(mid-1)
##      ramp[size:mid:-1,0] = numpy.arange(mid) / float(mid-1)
##      ramp[size:mid:-1,1] = numpy.arange(mid) / float(mid-1)
##      return ramp.astype('f')

## def RedWhiteBlueARamp(size=256):
##     ramp = numpy.ones( (size, 4), 'd')
##     mid = size/2
##     ramp[:mid,1] = numpy.arange(mid) / float(mid-1)
##     ramp[:mid,2] = numpy.arange(mid) / float(mid-1)
##     ramp[size:mid:-1,0] = numpy.arange(mid) / float(mid-1)
##     ramp[size:mid:-1,1] = numpy.arange(mid) / float(mid-1)
##     return ramp.astype('f')

def RedWhiteRamp(size=256):
    ramp = numpy.ones( (size, 3), 'd')
    ramp[:,1] = numpy.arange(size) / float(size-1)
    ramp[:,2] = numpy.arange(size) / float(size-1)
    return ramp.astype('f')

def RedWhiteARamp(size=256):
    ramp = numpy.ones( (size, 4), 'd')
    ramp[:,1] = numpy.arange(size) / float(size-1)
    ramp[:,2] = numpy.arange(size) / float(size-1)
    return ramp.astype('f')

def WhiteBlueRamp(size=256):
    ramp = numpy.ones( (size, 3), 'd')
    ramp[::-1,0] = numpy.arange(size) / float(size-1)
    ramp[::-1,1] = numpy.arange(size) / float(size-1)
    return ramp.astype('f')
        
def WhiteBlueARamp(size=256):
    ramp = numpy.ones( (size, 4), 'd')
    ramp[::-1,0] = numpy.arange(size) / float(size-1)
    ramp[::-1,1] = numpy.arange(size) / float(size-1)
    return ramp.astype('f')

def GreyscaleRamp(size=256):
    ramp = numpy.ones( (size, 3), 'd')
    ramp[:,0] = numpy.arange(size) / float(size-1)
    ramp[:,1] = numpy.arange(size) / float(size-1)
    ramp[:,2] = numpy.arange(size) / float(size-1)
    return ramp.astype('f')

def GreyscaleARamp(size=256):
    ramp = numpy.ones( (size, 4), 'd')
    ramp[:,0] = numpy.arange(size) / float(size-1)
    ramp[:,1] = numpy.arange(size) / float(size-1)
    ramp[:,2] = numpy.arange(size) / float(size-1)
    return ramp.astype('f')
    
    

    
if __name__ == '__main__':

    def test():

	print '10000 hsv->rgb'
	for i in range(100):
		for j in range(100):
			d = rgb( (i/100, j/100, (i+j)/100) )

	print '10000 rgb->hsv'
	for i in range(100):
		for j in range(100):
			d = hsv( (i/100, j/100, (i+j)/100) )

	cmap = RGBRamp()
