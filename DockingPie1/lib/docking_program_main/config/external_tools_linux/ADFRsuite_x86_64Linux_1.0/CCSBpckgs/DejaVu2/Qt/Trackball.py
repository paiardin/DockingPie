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

'''
Created on Jun 20, 2012

@author: brunsc
'''

from rotation import Rotation, Vec3
import rotation
from PySide import QtCore
from PySide.QtCore import QObject
from PySide.QtCore import Qt
from math import sqrt, pi

class Trackball(QObject):
    "Trackball converts mouse events into rotation/translation/scale signals"

    rotation_incremented = QtCore.Signal(Rotation)
    zoom_incremented = QtCore.Signal(float)
    pixel_translated = QtCore.Signal(int, int, int)

    def __init__(self, parent=None):
        QObject.__init__(self, parent)
        self.old_pos = None

    def mouseMoveEvent(self, event, windowSize):
        if self.old_pos is None:
            # Without a previous position, we cannot compute the change
            self.old_pos = event.pos()
            return
        dx = event.pos() - self.old_pos
        x = dx.x()
        y = dx.y()
        self.old_pos = event.pos()
        # Press shift to translate rather than rotate
        doTranslate = False
        if event.modifiers() & Qt.ShiftModifier: # shift-drag to translate
            doTranslate = True
        if event.buttons() & Qt.MidButton: # middle-drag to translate
            doTranslate = True
        if doTranslate:
            # Translate view
            self.pixel_translated.emit(-x, y, 0)
        else:
            # Rotate like a trackball
            w = (windowSize.width() + windowSize.height()) / 2.0
            # Dragging the whole window size would be 360 degrees
            angle = 2.0 * pi * sqrt(x*x + y*y) / w
            if 0.0 == angle:
                return # no rotation
            axis = Vec3([y, x, 0]).unit() # orthogonal to drag direction
            # print axis, angle
            r = Rotation().set_from_angle_about_unit_vector(angle, axis)
            self.rotation_incremented.emit(r)
        
    def mouseDoubleClickEvent(self, event, windowSize):
        "center view when double click occurs"
        dx = event.pos().x() - windowSize.width() / 2.0
        dy = event.pos().y() - windowSize.height() / 2.0
        self.pixel_translated.emit(dx, -dy, 0)

    def mousePressEvent(self, event):
        self.old_pos = None
        
    def mouseReleaseEvent(self, event):
        self.old_pos = None
        
    def wheelEvent(self, event):
        degrees = event.delta() / 8.0
        zoom_ratio = 2.0 ** (-degrees / 200.0)
        self.zoom_incremented.emit(zoom_ratio)
