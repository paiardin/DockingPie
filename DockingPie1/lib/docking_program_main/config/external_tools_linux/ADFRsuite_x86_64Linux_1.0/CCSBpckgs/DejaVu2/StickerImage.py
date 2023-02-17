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
# Date: October 2006 Authors: Guillaume Vareille, Michel Sanner
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
# $Header: /mnt/raid/services/cvs/DejaVu2/StickerImage.py,v 1.1.1.1.4.2 2017/11/08 23:27:21 annao Exp $
#
# $Id: StickerImage.py,v 1.1.1.1.4.2 2017/11/08 23:27:21 annao Exp $
#

import os
from PIL import Image
from copy import deepcopy

from opengltk.extent import _gllib
from opengltk.OpenGL import GL
from DejaVu2.Insert2d import Insert2d

class StickerQuad(Insert2d):

    keywords = Insert2d.keywords + [
        'color',
        ]

    def __init__(self, name='StickerPylgon', check=1, **kw):
        #print "StickerImage.__init__"

        self.color = [1,1,1,1]

        apply(Insert2d.__init__, (self, name, check), kw)

        self.needsRedoDpyListOnResize = True


    def Set(self, check=1, redo=1, updateOwnGui=True, **kw):
        """set data for this object
check=1 : verify that all the keywords present can be handle by this func 
redo=1 : append self to viewer.objectsNeedingRedo
updateOwnGui=True : allow to update owngui at the end this func
"""
        #print "Insert2d.Set"
        redoFlags = apply( Insert2d.Set, (self, check, 0), kw)

        val = kw.pop( 'transparent', None)
        if not val is None:
            redoFlags |= self._setTransparent(val)

        anchor = kw.get( 'anchor')
        if anchor \
           and (anchor[0] >= 0.) and (anchor[0] <= 1.) \
           and (anchor[1] >= 0.) and (anchor[1] <= 1.):
            self.anchor = anchor
            redoFlags |= self._redoFlags['redoDisplayListFlag']

        position = kw.get( 'position')
        if position:
            self.position = position
            redoFlags |= self._redoFlags['redoDisplayListFlag']

        size = kw.get( 'size')
        if size:
            self.size = size
            redoFlags |= self._redoFlags['redoDisplayListFlag']

        color = kw.get( 'color')
        if color:
            print color
            self.color = color
            print self.color
            redoFlags |= self._redoFlags['redoDisplayListFlag']

        return self.redoNow(redo, updateOwnGui, redoFlags)


    def Draw(self):

        GL.glMatrixMode(GL.GL_PROJECTION)
        GL.glPushMatrix()
        GL.glLoadIdentity()
        Insert2d.Draw(self)
        GL.glMatrixMode(GL.GL_MODELVIEW)
        GL.glPushMatrix()
        GL.glLoadIdentity()

        GL.glDisable(GL.GL_DEPTH_TEST)
        GL.glDisable(GL.GL_LIGHTING); 

        width = self.size[0]
        height = self.size[1]

        fullWidth = self.viewer.currentCamera.width
        fullHeight = self.viewer.currentCamera.height

        # we want the anchor of the image to be at the given position
        posxFromLeft = self.position[0] * fullWidth - self.anchor[0] * width
        posyFrombottom = (1.-self.position[1]) * fullHeight - (1.-self.anchor[1]) * height
        #print "posxFromLeft, posyFrombottom", posxFromLeft, posyFrombottom

        GL.glColor4fv(self.color)
        
        xo, yo = self.position
        dx, dy = self.size
        GL.glBegin(GL.GL_QUADS);
        GL.glVertex2f(xo, yo)
        GL.glVertex2f(xo+dx, yo)
        GL.glVertex2f(xo+dx, yo+dy)
        GL.glVertex2f(xo, yo+dy)
        GL.glEnd()

       # used for picking
       # self.polygonContour = [ (posxFromLeft, posyFrombottom),
       #                         (posxFromLeft+width, posyFrombottom),
       #                         (posxFromLeft+width, posyFrombottom+height),
       #                         (posxFromLeft, posyFrombottom+height)
       #                       ]

        GL.glEnable(GL.GL_DEPTH_TEST)
        GL.glEnable(GL.GL_LIGHTING); 

        GL.glMatrixMode(GL.GL_PROJECTION)
        GL.glPopMatrix()
        GL.glMatrixMode(GL.GL_MODELVIEW)
        GL.glPopMatrix()

        return 1


class StickerImage(Insert2d):

    keywords = Insert2d.keywords + [
        'image',
        ]

    def __init__(self, name='StickerImage', check=1, **kw):
        #print "StickerImage.__init__"

        self.image =  None

        apply(Insert2d.__init__, (self, name, check), kw)

        self.needsRedoDpyListOnResize = True


    def Set(self, check=1, redo=1, updateOwnGui=True, **kw):
        """set data for this object:
check=1 : verify that all the keywords present can be handle by this func 
redo=1 : append self to viewer.objectsNeedingRedo
updateOwnGui=True : allow to update owngui at the end this func
"""
        #print "StickerImage.Set"
        redoFlags = apply( Insert2d.Set, (self, check, 0), kw)

        image = kw.get('image')
        if image is not None:
            kw.pop('image')
            self.image = image
            self.size[0] = self.image.size[0]
            self.size[1] = self.image.size[1]
            redoFlags |= self._redoFlags['redoDisplayListFlag']

        return self.redoNow(redo, updateOwnGui, redoFlags)


    def Draw(self):
        #print "StickerImage.Draw", self

        if self.image is None:
            return

        GL.glMatrixMode(GL.GL_PROJECTION)
        GL.glPushMatrix()
        GL.glLoadIdentity()
        Insert2d.Draw(self)
        GL.glMatrixMode(GL.GL_MODELVIEW)
        GL.glPushMatrix()
        GL.glLoadIdentity()

        GL.glDisable(GL.GL_DEPTH_TEST)
        GL.glDisable(GL.GL_LIGHTING); 

        width = self.size[0]
        height = self.size[1]

        fullWidth = self.viewer.currentCamera.width
        fullHeight = self.viewer.currentCamera.height

        # we want the anchor of the image to be at the given position
        posxFromLeft = self.position[0] * fullWidth - self.anchor[0] * width
        posyFrombottom = (1.-self.position[1]) * fullHeight - (1.-self.anchor[1]) * height
        #print "posxFromLeft, posyFrombottom", posxFromLeft, posyFrombottom

        # used for picking
        self.polygonContour = [ (posxFromLeft, posyFrombottom),
                                (posxFromLeft+width, posyFrombottom),
                                (posxFromLeft+width, posyFrombottom+height),
                                (posxFromLeft, posyFrombottom+height)
                              ]

        # this accept negative values were GL.glRasterPos2f(x,y) doesn't
        GL.glRasterPos2f(0, 0)
        _gllib.glBitmap(0, 0, 0, 0, posxFromLeft, posyFrombottom, 0)

        GL.glPixelStorei(GL.GL_UNPACK_ALIGNMENT, 1)
        if self.image.mode == 'RGBA':
                _gllib.glDrawPixels(self.image.size[0], self.image.size[1], 
                                    GL.GL_RGBA, GL.GL_UNSIGNED_BYTE, 
                                    self.image.tostring() )
        elif self.image.mode == 'RGB':
                _gllib.glDrawPixels(self.image.size[0], self.image.size[1], 
                                    GL.GL_RGB, GL.GL_UNSIGNED_BYTE, 
                                    self.image.tostring() )
        elif self.image.mode == 'L':
                _gllib.glDrawPixels(self.image.size[0], self.image.size[1], 
                                    GL.GL_LUMINANCE, GL.GL_UNSIGNED_BYTE, 
                                    self.image.tostring() )

        GL.glMatrixMode(GL.GL_PROJECTION)
        GL.glPopMatrix()
        GL.glMatrixMode(GL.GL_MODELVIEW)
        GL.glPopMatrix()

        return 1


    def setSize(self, event, redo=1):
        """the trackball transmit the translation info
"""
        pass
