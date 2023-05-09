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
# Date: Febuary 2006 Authors: Guillaume Vareille, Michel Sanner
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
# Revision:
#
#########################################################################
#
# $Header: /mnt/raid/services/cvs/DejaVu2/Insert2d.py,v 1.1.1.1.4.1 2017/07/13 22:28:32 annao Exp $
#
# $Id: Insert2d.py,v 1.1.1.1.4.1 2017/07/13 22:28:32 annao Exp $
#

import numpy
from copy import deepcopy
from opengltk.OpenGL import GL
from DejaVu2.Common2d3dObject import Common2d3dObject


class Insert2d(Common2d3dObject):
    """Base class for 2d inserts that can be displayed in a camera
"""
    keywords = Common2d3dObject.keywords + [
        'transparent',                                            
        'anchor',
        'position',
        'size',
    ]

    def __init__(self, name='insert2d', check=1, **kw):
        #print "Insert2d.__init__"

        self.anchor = [1, 0] # or 0 - left/top; .5 - center; 1 - right/bottom
        self.position = [1, 0] # or 0 - left/top; .5 - center; 1 - right/bottom

        self.size = [40, 30] # in pixels

        self.polygonContour = None # the list of 2d points (pixels coordinates from bottom left)
                          # that defines the 2d polygon drawned by Draw()

        self.lastPickEventTime = 0
        self.isMoving = None # tells what is the current behavior 
                             # the legend can be move (True) or resize (False) or nothing (None)

        apply( Common2d3dObject.__init__, (self, name, check), kw)

        self.initialAnchor = deepcopy(self.anchor)
        self.initialPosition = deepcopy(self.position)
        self.initialSize = deepcopy(self.size)


    def Set(self, check=1, redo=1, updateOwnGui=True, **kw):
        """set data for this object
check=1 : verify that all the keywords present can be handle by this func 
redo=1 : append self to viewer.objectsNeedingRedo
updateOwnGui=True : allow to update owngui at the end this func
"""
        #print "Insert2d.Set"
        redoFlags = apply( Common2d3dObject.Set, (self, check, 0), kw)

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
        if position \
           and (position[0] >= 0.) and (position[0] <= 1.) \
           and (position[1] >= 0.) and (position[1] <= 1.):
            self.position = position
            redoFlags |= self._redoFlags['redoDisplayListFlag']

        size = kw.get( 'size')
        if size \
           and (size[0] >= 0.) \
           and (size[1] >= 0.):
            self.size = size
            redoFlags |= self._redoFlags['redoDisplayListFlag']

        return self.redoNow(redo, updateOwnGui, redoFlags)


    def _setTransparent(self, val):
        #print "_setTransparent", val
        assert val in (False, True)
        self.transparent = val
        return self._redoFlags['redoDisplayListFlag'] | self._redoFlags['redoViewerDisplayListFlag']


    def getState(self, full=1):
        """return a dictionary describing this object's state
This dictionary can be passed to the Set method to restore the object's state
"""
        state = Common2d3dObject.getState(self).copy()
        state.update( {
                       'transparent':self.transparent,
                       'anchor': self.anchor,
                       'position': self.position,
                       'size': self.size,
                      } ) 
        return state


    def setViewer(self, viewer, buildDisplayList):
        # could be renamed as 'onViewerAddition'
        #print "Insert2d.setViewer"
        Common2d3dObject.setViewer(self, viewer, buildDisplayList)
        self.viewer.AddPickingCallback(self.processHit_cb)


    def getTile(self):
        #print "Insert2d.getTile"
        # compute projection parameters
        tr = self.viewer.tileRenderCtx
        Left = 0
        Right = tr.TileWidth #tr.CurrentTileWidth
        Bottom = 0
        Top = tr.TileHeight #tr.CurrentTileHeight
        left = Left + (Right - Left) \
              * (tr.CurrentColumn * tr.TileWidthNB - tr.TileBorder) / float(tr.ImageWidth)
        right = left + (Right - Left) * tr.TileWidth / float(tr.ImageWidth)
        bottom = Bottom + (Top - Bottom) \
                * (tr.CurrentRow * tr.TileHeightNB - tr.TileBorder) / float(tr.ImageHeight)
        top = bottom + (Top - Bottom) * tr.TileHeight / float(tr.ImageHeight)
        tile = (left, right, bottom, top)
        #print "tile", tile
        return tile


    def Draw(self):
        """set the projection matrices, the subclass must call this in their Draw func,
otherwise tile rendering and stereo won't work
"""
        #print "Insert2d.Draw"
        
        GL.glMatrixMode(GL.GL_PROJECTION)
        if self.viewer.tileRender:
            tile = self.getTile()
            GL.glOrtho(float(tile[0]), float(tile[1]), float(tile[2]), float(tile[3]), -1, 1)
        else:
            GL.glOrtho(0, 
                       float(self.viewer.currentCamera.width), 
                       0, 
                       float(self.viewer.currentCamera.height), 
                       -1, 1)
        GL.glMatrixMode(GL.GL_MODELVIEW)


    def pickDraw(self):
        """called by the picking process to operate the selection
"""
        #print "Insert2d.pickDraw", self
        if self.polygonContour is not None:
            # we draw just flat quad of the insert2d
            GL.glMatrixMode(GL.GL_PROJECTION)
            GL.glPushMatrix()       
            #GL.glLoadIdentity()
            GL.glLoadMatrixf(self.viewer.currentCamera.pickMatrix) 
            GL.glOrtho(0, float(self.viewer.currentCamera.width),
                       0, float(self.viewer.currentCamera.height), -1, 1)
            GL.glMatrixMode(GL.GL_MODELVIEW)
            GL.glPushMatrix()
            GL.glLoadIdentity()
            GL.glPolygonMode(GL.GL_FRONT, GL.GL_FILL)
            GL.glPushName(0)
            GL.glBegin(GL.GL_QUADS)
            GL.glVertex2fv(self.polygonContour[0])
            GL.glVertex2fv(self.polygonContour[1])
            GL.glVertex2fv(self.polygonContour[2])
            GL.glVertex2fv(self.polygonContour[3])
            GL.glEnd()
            GL.glPopName()
        
            GL.glMatrixMode(GL.GL_PROJECTION)
            GL.glPopMatrix()
            GL.glMatrixMode(GL.GL_MODELVIEW)
            GL.glPopMatrix()


    def setSize(self, event, redo=1):
        """the trackball transmit the translation info
"""
        self.size[0] = event.x - self.coord2d[0]        
        self.size[1] = self.coord2d[1] - event.y

        if self.size[0] < 1:
            self.size[0] = 1
        if self.size[1] < 1:
            self.size[1] = 1

        if self.size[0] > self.viewer.currentCamera.width:
            self.size[0] = self.viewer.currentCamera.width
        if self.size[1] > self.viewer.currentCamera.height:
            self.size[1] = self.viewer.currentCamera.height

        if self.viewer:
            self.viewer.objectsNeedingRedo[self] = None


    def setPosition(self, event, redo=1):
        """the trackball transmit the translation info
"""
        #print "Insert2d.setPosition", event.x, event.y
        self.position = [ event.x / float(self.viewer.currentCamera.width),
                          event.y / float(self.viewer.currentCamera.height)
                        ]

        if self.position[0] < 0:
            self.position[0] = 0
        elif self.position[0] > 1:
            self.position[0] = 1

        if self.position[1] < 0:
            self.position[1] = 0
        elif self.position[1] > 1:
            self.position[1] = 1

        self.viewer.objectsNeedingRedo[self] = None


    def calculateAnchorAndPosition(self, event):
        insert2dTopLeftFromScreenTopleft0 =   self.position[0] * self.viewer.currentCamera.width \
                                            - self.anchor[0] * self.size[0]
        insert2dTopLeftFromScreenTopleft1 =   self.position[1] * self.viewer.currentCamera.height \
                                            - self.anchor[1] * self.size[1]
        #print "insert2dTopLeftFromScreenTopleft", insert2dTopLeftFromScreenTopleft0, insert2dTopLeftFromScreenTopleft1
        self.anchor = [ (event.x - insert2dTopLeftFromScreenTopleft0) / float(self.size[0]),
                        (event.y - insert2dTopLeftFromScreenTopleft1) / float(self.size[1])
                      ]
        #print "self.anchor", self.anchor
        self.position = [ event.x / float(self.viewer.currentCamera.width),
                          event.y / float(self.viewer.currentCamera.height) 
                        ]


    def ResetPosition(self):
        self.anchor = deepcopy(self.initialAnchor)
        self.position = deepcopy(self.initialPosition)
        if self.viewer:
            self.viewer.objectsNeedingRedo[self] = None


    def ResetSize(self):
        #print "ResetSize", self.initialSize
        self.size = deepcopy(self.initialSize)
        if self.viewer:
            self.viewer.objectsNeedingRedo[self] = None


    def ResetTransformation(self):
        """Reset the tranformations
"""
        self.ResetPosition()       
        self.ResetSize()       

            
    def respondToMouseMove(self, event, redo=1):
        """the trackball transmit the event info
"""
        if self.isMoving is not None:
            if self.isMoving is True:
                self.setPosition(event)
            else:
                self.setSize(event)
            self.viewer.Redraw()


    def respondToDoubleClick(self, event):
        """to be overidden
"""
        pass


    def processHit_cb(self, pick):
        #print "Insert2d.processHit_cb", self, (self.viewer==None)
        #print "pick",pick
        #print "pick.event",dir(pick)
        #print "pick.type",pick.type
        #print "pick.event",dir(pick.event)
        #print "pick.event",pick.event
        #print "pick.event.type",pick.event.type
        #print "pick.event.state",pick.event.state
        #print "pick.event.time",pick.event.time
        #print "pick.hits",pick.hits

        if ( len(pick.hits) == 1) and  pick.hits.has_key(self):
            if self.viewer.currentObject != self:
                    # if the only hit is the legend, 
                    # it becomes the current object
                    self.viewer.SetCurrentObject(self)
                    self.isMoving = True
                    self.calculateAnchorAndPosition(pick.event)
            elif pick.event.time - self.lastPickEventTime < 200: #double click
                self.viewer.SetCurrentObject(self.viewer.rootObject)
                self.respondToDoubleClick(pick.event)
            elif pick.hits[self][0][0] == 1:
                # the click in inside the resize button
                #print "resize"
                self.isMoving = False
            elif pick.hits[self][0][0] == 0:
                # the click in inside the legend but outside
                # the resize button
                self.isMoving = True
                #print "pick.event", pick.event.x, pick.event.y
                self.calculateAnchorAndPosition(pick.event)
                #print "self.position", self.position
            if self.viewer:
                self.viewer.objectsNeedingRedo[self] = None

        elif self.viewer.currentObject == self:
            #print "the insert2d is selected, but picking is outside"
            self.isMoving = None
            self.viewer.SetCurrentObject(self.viewer.rootObject)

        self.lastPickEventTime = pick.event.time
