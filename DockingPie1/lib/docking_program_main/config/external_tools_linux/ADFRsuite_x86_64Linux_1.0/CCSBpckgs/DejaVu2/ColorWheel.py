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
# $Header: /mnt/raid/services/cvs/DejaVu2/ColorWheel.py,v 1.1.1.1.4.1 2017/07/13 22:28:32 annao Exp $
#
# $Id: ColorWheel.py,v 1.1.1.1.4.1 2017/07/13 22:28:32 annao Exp $
#

import types
import Tkinter
import numpy
import math
import colorTool
import Slider
from EventHandler import CallbackFunctions
from DejaVu2.colorTool import ToHSV


class ColorWheel(CallbackFunctions):


    def MakeWheelGeometry(self):
	"""Compute vertex coordinates to draw the wheel's color patches"""

	# makes the master circle
	r = float(self.radius / self.circles)
	angle = 2.0 * math.pi / self.stripes
	c = numpy.zeros( (self.stripes+1,2),'f' )
	c[0][0] = c[self.stripes][0] = r;
	c[0][1] = c[self.stripes][1] = 0.0;
	for i in range(1, self.stripes):
	    c[i][0] = r * math.cos(i*angle);
	    c[i][1] = r * math.sin(i*angle);

	# now build the color wheel geometry

	# center point
	self.centerX = self.width / 2
	self.centerY = self.height / 2
	self.c[0][0] = self.centerX
	self.c[0][1] = self.centerY

	nco = 1
	# strips+1 points on a circle (first point == last point)
	for j in range(1,self.circles+1):
	    for i in range(self.stripes+1):
		self.c[nco][0] = self.centerX + (j * c[i][0])
		self.c[nco][1] = self.centerY + (j * c[i][1])
		nco = nco+1


    def MakeWheelColors(self, intensity):
	"""build the array of colors to draw a wheel"""

	self.hsvColor[2] = intensity
	Hinc = 1.0 / float(self.stripes)
	Sinc = 1.0 / float(self.circles)

	# center color
	self.rgbCol[0][0] = self.rgbCol[0][1] = self.rgbCol[0][2] = intensity
	ncol = 1;

	# now make the wheel colors
	for j in range(1,self.circles+1):
	    s = j * Sinc;
	    for i in range(self.stripes):
		self.rgbCol[ncol] = colorTool.ToRGB( (i*Hinc, s, intensity) )
		ncol = ncol+1


    def DrawWheel(self):
	"""draw the color wheel"""

	c = self.c
	col = self.tkColors

	# draw center circle as triangle fan
	for i in range(1, self.stripes+1):
	    self.canvas.create_polygon(
		c[0][0], c[0][1], c[i][0], c[i][1], c[i+1][0], c[i+1][1],
		fill = col[i], outline = col[i], tag='wheel' )

	s = self.stripes
	ncol = self.stripes+1
	for i in range(1,self.circles):
	    c1 = 1 + (i-1)*s + i-1
	    c2 = c1 + s + 1
	    k = ncol
	    for j in range(self.stripes):
		self.canvas.create_polygon(
		    c[c1+1+j][0], c[c1+1+j][1], c[c1+j][0], c[c1+j][1],
		    c[c2+j][0], c[c2+j][1], c[c2+j+1][0], c[c2+j+1][1], 
		    fill = col[ncol], outline = col[ncol], tag='wheel' )
		ncol = ncol+1


    def DrawCursor(self):
	"""Compute cursor position from hsv value and draw it"""

	rad = self.hsvColor[1] * self.radius
	angle = 2.0 * math.pi * self.hsvColor[0]
	x = self.centerX + int(rad * math.cos(angle))
	y = self.centerY + int(rad * math.sin(angle))
	dx = x-self.cursorX
	dy = y-self.cursorY
	self.canvas.move( 'cursor', dx, dy)
	self.cursorX = self.cursorX + dx
	self.cursorY = self.cursorY + dy


    def _MoveCursor(self, x, y):
	# find the saturation based on distance
	s = math.sqrt(x*x + y*y) / self.radius
	if s > 1.0:
	    s = 1.0
	    x = x / s
	    y = y / s

	# now find the hue based on the angle 
	if x or y:
	    angle = math.atan2(y, x)
	    if angle < 0.0: angle = angle + (2.0 * math.pi)
	    h = angle / (2.0 * math.pi)
	else:
	    h = 0

	# check if redraw and callback are needed
	if self.hsvColor[0] != h or self.hsvColor[1] != s:
            if type(self.hsvColor) == types.TupleType:
                self.hsvColor = list(self.hsvColor)
	    self.hsvColor[0] = h;
	    self.hsvColor[1] = s;
	    self.DrawCursor()
	    if self.immediate: self.Callbacks()


    def MoveCursor(self, event=None):
	"""Compute the new color for the cursor position, redraw cursor"""

	x = event.x - self.width*0.5
	y = event.y - self.height*0.5
	self._MoveCursor(x,y)


    def Callbacks(self):
	"""Implement callback functions
"""
	rgb = colorTool.ToRGB(self.hsvColor)
	for f in self.callbacks:
	    f( colorTool.OneColor(rgb, self.alpha) )


    def Get(self, mode='HSV'):
	"""Get the current color"""
	if mode == 'RGB':
	    rgb = colorTool.ToRGB(self.hsvColor)
	    return colorTool.OneColor(rgb)
	elif mode == 'HSV':
	    return colorTool.OneColor(self.hsvColor)
	elif mode == 'TkRGB':
	    col = numpy.array(colorTool.ToRGB(self.hsvColor[:3]), 'f') * 255
	    return colorTool.TkColor(col)


    def Set(self, color, mode='HSV', run=True):
	"""Set the current color
"""
	assert len(color) in (3,4)
	color2 = colorTool.OneColor(color)
	if mode=='RGB':
            color2 = colorTool.ToHSV(color[:3])
        if len(color)==4:
            self.alpha = color[3]
        else:
            self.alpha = 1.0
	if color2[2] != self.hsvColor[2] and self.wysiwyg:
	    self.BuildWheel(color2[2])
	self.hsvColor = color2
	self.DrawCursor()
	if self.immediate and run: self.Callbacks()


    def BuildWheel(self, value):
	"""Compute the colors, convert them to Tk, draw wheel, draw cursor"""

	self.canvas.delete('wheel')
	self.canvas.delete('cursor')
	self.MakeWheelColors(value)
	self.tkColors = map( colorTool.TkColor, self.rgbCol*255 )
	self.DrawWheel()
	cx = self.cursorX
	cy = self.cursorY
	cs = self.cursorSize
	self.canvas.create_line( cx-cs, cy-cs, cx+cs, cy-cs, cx+cs, cy+cs,
				 cx-cs, cy+cs, cx-cs, cy-cs,
				 width=2, fill='black', tag = 'cursor')


    def setWysiwyg(self, onOff):
	"""Toggle Wysiwyg mode. When on, wheel colors are recomputed"""

	assert onOff in (0,1)
	self.wysiwyg = onOff
	if onOff: self.BuildWheel(self.hsvColor[2])
	else: self.BuildWheel(1.0)


    def setImmediate(self, onOff):
	"""Toggle Wysiwyg mode. When on, wheel colors are recomputed"""

	assert onOff in (0,1)
	self.immediate = onOff


    def __init__(self, parent, width=100, height=100, radius=None,
		 circles=5, stripes=10, immediate=1, wysiwyg=0, wheelPad=10,
                 startingColorRGB=(1,1,1,1)):

	CallbackFunctions.__init__(self)
	self.width = width
	self.height = height
        self.wheelPad = wheelPad
	if radius is None:
	    self.radius = max(width,height) / 2
	else:
	    self.radius = min( max(width,height) / 2, self.radius )
	self.circles = circles
	self.stripes = stripes
	self.wysiwyg = wysiwyg
	self.immediate = immediate
	self.ncol = 1 + self.circles * (self.stripes + 1)
	self.rgbCol = numpy.zeros( (self.ncol,3), 'f' )
	self.c = numpy.zeros( (self.ncol,2), 'f' )
	self.MakeWheelGeometry()
	self.hsvColor = numpy.array(ToHSV(startingColorRGB))
	self.cursorX = self.centerX
	self.cursorY = self.centerY
	self.cursorSize = 2
        self.alpha = 1.0

	self.canvas = Tkinter.Canvas(parent, relief=Tkinter.SUNKEN,
				     width = self.width, height = self.height)
	self.BuildWheel(1.0)
	self.canvas.pack(padx=wheelPad, pady=wheelPad, side=Tkinter.RIGHT,
                         expand='no')

	self.canvas.bind( "<1>", self.MoveCursor)
	self.canvas.bind( "<B1-Motion>", self.MoveCursor)

	if not immediate:
	    self.canvas.bind( "<ButtonRelease-1>", self.MouseUp)
	self.callbacks = []

        self.DrawCursor()



if __name__ == '__main__':

    def MyCallback(color):
	print color

    def MyCallback2(color):
	print 'hello'

    root = Tkinter.Tk()
    cw = ColorWheel(root)
    cw.AddCallback(MyCallback)

    cw1 = ColorWheel(root, immediate=0)
    cw1.AddCallback(MyCallback2)

    cw2 = ColorWheel(root, immediate=0)
    cw1.AddCallback(MyCallback)
    cw2.AddCallback(MyCallback2)

