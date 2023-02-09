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
# Date: August 2006 Authors: Guillaume Vareille, Michel Sanner
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
# $Header: /mnt/raid/services/cvs/DejaVu2/glfSticker.py,v 1.1.1.1.4.1 2017/07/13 22:28:33 annao Exp $
#
# $Id: glfSticker.py,v 1.1.1.1.4.1 2017/07/13 22:28:33 annao Exp $
#

import os
import Tkinter

from math import floor, ceil

from opengltk.OpenGL import GL

from pyglf import glf
from DejaVu2.Insert2d import Insert2d
from DejaVu2.viewerFns import checkKeywords
from DejaVu2.ColorChooser import ColorChooser


class Sticker(Insert2d):  
   
    keywords = Insert2d.keywords + [
        'label',
        'fontName',
        'fontSpacing',
        'fontColor',
        'fontScales',
        'wireFont',
        'frameColor',
        'framePolygonMode',
        'frameSpace',
        ]
    
    glfVectorFontList = [
        'arial1.glf',
        'courier1.glf',
        'crystal1.glf',
        'techno0.glf',
        'techno1.glf',
        'times_new1.glf',
        'aksent1.glf',
        'alpine1.glf',
        'broadway1.glf',
        'chicago1.glf',
        'compact1.glf',
        'cricket1.glf',
        'garamond1.glf',
        'gothic1.glf',
        'penta1.glf',
        'present_script1.glf'
    ]
    
    framePolygonModeDict = { 'none' : GL.GL_NONE,
                             'line' : GL.GL_LINE,
                             'fill' : GL.GL_FILL }

    framePolygonModeDictRev = { GL.GL_NONE :'none',
                                GL.GL_LINE :'line',
                                GL.GL_FILL :'fill'}

    
    def __init__(self, name='sticker', check=1, **kw):
        
        # GLF FONTS Initialisations
        glf.glfInit()
        glf.glfEnable(glf.GLF_CONSOLE_MESSAGES)
        lGlfModulePath = os.path.split(glf.__file__)[-2]
        lPathToFonts = lGlfModulePath+os.sep+'fonts'+os.sep
        self.glfVectorFontLoadIdDict = {}
        for font in self.glfVectorFontList:
            self.glfVectorFontLoadIdDict[font] = glf.glfLoadFont(lPathToFonts+font)   

        self.label = 'No Label Yet'
        self.fontName = 'arial1.glf'
        self.glfFontID = self.glfVectorFontLoadIdDict[self.fontName]
        self.fontSpacing = .2
        self.fontColor = (1, 1, 1, .8)
        self.fontScales = (10, 10)
        self.wireFont = False
        self.frameColor = (0, 1, 1, .5)        
        self.framePolygonMode = GL.GL_LINE #GL.GL_FILL # GL.GL_NONE for no frame
        self.frameSpace = (.2, .2)

        apply( Insert2d.__init__, (self, name), kw)

        self.needsRedoDpyListOnResize = True


    def Set(self, check=1, redo=1, updateOwnGui=True, **kw):
        """set data for this object:
check=1 : verify that all the keywords present can be handle by this func 
redo=1 : append self to viewer.objectsNeedingRedo
updateOwnGui=True : allow to update owngui at the end this func
"""
        #print "Sticker.Set"
        redoFlags = apply( Insert2d.Set, (self, check, 0), kw)
        
        label = kw.get('label')
        if label:
            kw.pop('label')
            self.label = label
            redoFlags |= self._redoFlags['redoDisplayListFlag']

        fontName = kw.get('fontName')
        if fontName:
            kw.pop('fontName')
            self.glfFontID = self.glfVectorFontLoadIdDict[fontName]
            redoFlags |= self._redoFlags['redoDisplayListFlag']

        fontSpacing = kw.get('fontSpacing')
        if fontSpacing:
            kw.pop('fontSpacing')
            self.fontSpacing = fontSpacing
            redoFlags |= self._redoFlags['redoDisplayListFlag']

        position = kw.get('position')
        if position:
            kw.pop('position')
            self.position = position
            redoFlags |= self._redoFlags['redoDisplayListFlag']

        fontScales = kw.get('fontScales')
        if fontScales:
            kw.pop('fontScales')
            self.fontScales = fontScales
            redoFlags |= self._redoFlags['redoDisplayListFlag']

        wireFont = kw.get('wireFont')
        if wireFont is not None:
            kw.pop('wireFont')
            self.wireFont = wireFont
            redoFlags |= self._redoFlags['redoDisplayListFlag']

        fontColor = kw.get('fontColor')
        if fontColor:
            kw.pop('fontColor')
            self.fontColor = fontColor
            redoFlags |= self._redoFlags['redoDisplayListFlag']

        frameColor = kw.get('frameColor')
        if frameColor:
            kw.pop('frameColor')
            self.frameColor = frameColor
            redoFlags |= self._redoFlags['redoDisplayListFlag']

        framePolygonMode = kw.get('framePolygonMode')
        if framePolygonMode is not None:
            kw.pop('framePolygonMode')
            self.framePolygonMode = framePolygonMode
            redoFlags |= self._redoFlags['redoDisplayListFlag']

        frameSpace = kw.get('frameSpace')
        if frameSpace:
            kw.pop('frameSpace')
            self.frameSpace = frameSpace
            redoFlags |= self._redoFlags['redoDisplayListFlag']

        return self.redoNow(redo, updateOwnGui, redoFlags)


    def calculateSize(self):
        #print "Sticker.calculateSize", self
        glf.glfSetCurrentFont(self.glfFontID)
        glf.glfSetSymbolSpace(self.fontSpacing)

        lLabelbounds = glf.glfGetStringBounds(self.label)

        # this is the correct way to use glfGetStringBounds to obtain the exact frame around the label
        # (it still needs more correction when there are spaces in the label)
        lFramebounds = ( -1,
                         lLabelbounds[1],
                         -1+(lLabelbounds[2]-lLabelbounds[0])-glf.glfGetSymbolSpace(),
                         lLabelbounds[3] )

        # we want the frame to be a bit larger than the label
        lFramebounds = ( lFramebounds[0]-self.frameSpace[0],
                         lFramebounds[1]-self.frameSpace[1],
                         lFramebounds[2]+self.frameSpace[0],
                         lFramebounds[3]+self.frameSpace[1])
        
        self.size[0] = (lFramebounds[2] - lFramebounds[0]) * self.fontScales[0]
        self.size[1] = (lFramebounds[3] - lFramebounds[1]) * self.fontScales[1]
        self.size[0] = int(ceil( self.size[0] ))
        self.size[1] = int(ceil( self.size[1] ))
        #print "self.size", self.size      

        return lLabelbounds


    def drawSticker(self, lLabelbounds):   

        if self.framePolygonMode != GL.GL_NONE:
            lenFrameColor = len(self.frameColor)
            if lenFrameColor == 4:
                GL.glColor4fv(self.frameColor)
            elif lenFrameColor == 3:
                GL.glColor3fv(self.frameColor)
            GL.glPolygonMode(GL.GL_FRONT, self.framePolygonMode)
            GL.glBegin(GL.GL_QUADS)
            GL.glVertex2f(0, 0)
            GL.glVertex2f(float(self.size[0]), 0)
            GL.glVertex2f(float(self.size[0]), float(self.size[1]))
            GL.glVertex2f(0, float(self.size[1]))
            GL.glEnd()

        GL.glScalef(float(self.fontScales[0]),
                    float(self.fontScales[1]),
                    0)
        GL.glTranslatef(float(self.frameSpace[0]),
                        float(self.frameSpace[1]),
                        0)

        # this corrects the glf draw function to start the label at the proper position
        GL.glTranslatef(1, float(-lLabelbounds[1]), 0)

        lenFontColor = len(self.fontColor)
        if lenFontColor == 4:
            GL.glColor4fv(self.fontColor)
        elif lenFontColor == 3:
            GL.glColor3fv(self.fontColor)
        GL.glPolygonMode(GL.GL_FRONT, GL.GL_FILL)
        if self.wireFont in [0, False]:
            glf.glfDrawSolidString(self.label)
        else:
            glf.glfDrawWiredString(self.label)


    def Draw(self):
        #print "Sticker.Draw", self
        
        GL.glMatrixMode(GL.GL_PROJECTION)
        GL.glPushMatrix()
        GL.glLoadIdentity()
        Insert2d.Draw(self)
        GL.glMatrixMode(GL.GL_MODELVIEW)
        GL.glPushMatrix()
        GL.glLoadIdentity()

        GL.glDisable(GL.GL_DEPTH_TEST)
        GL.glDisable(GL.GL_LIGHTING); 

        lLabelbounds = self.calculateSize()

        width = self.size[0]
        height = self.size[1]

        fullWidth = self.viewer.currentCamera.width
        fullHeight = self.viewer.currentCamera.height

        # we want the anchor of the image to be at the given position
        posxFromLeft = self.position[0] * fullWidth - self.anchor[0] * width
        posyFrombottom = (1.-self.position[1]) * fullHeight - (1.-self.anchor[1]) * height
        posxFromLeft = int(floor( posxFromLeft ) )
        posyFrombottom = int(floor( posyFrombottom ) )
        
        if (self.position[1] == 0.):
                posyFrombottom -= 1
        #print "posxFromLeft, posyFrombottom", posxFromLeft, posyFrombottom

        # used for picking
        self.polygonContour = [ (posxFromLeft, posyFrombottom),
                                (posxFromLeft+width, posyFrombottom),
                                (posxFromLeft+width, posyFrombottom+height),
                                (posxFromLeft, posyFrombottom+height)
                              ]

        GL.glTranslatef(float(posxFromLeft), float(posyFrombottom), 0)

        self.drawSticker(lLabelbounds)

        GL.glEnable(GL.GL_DEPTH_TEST)
        GL.glEnable(GL.GL_LIGHTING); 

        GL.glMatrixMode(GL.GL_PROJECTION)
        GL.glPopMatrix()
        GL.glMatrixMode(GL.GL_MODELVIEW)
        GL.glPopMatrix()

        return 1


    def setFontColor(self, color):
        """
"""
        self.Set(fontColor=color)


    def setFrameColor(self, color):
        """
"""
        self.Set(frameColor=color)


    def respondToDoubleClick(self, event):
        """
"""
        self.showOwnGui()


    def createOwnGui(self):
        #print "GlfLabels.createOwnGui", self
        pass
