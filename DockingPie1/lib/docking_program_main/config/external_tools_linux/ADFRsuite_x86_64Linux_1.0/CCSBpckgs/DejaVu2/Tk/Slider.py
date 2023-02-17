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

""" Slider Module

This module implements a class for Sliders.

constructor:
    def __init__(self, master=None, label=None, minval=0.0, maxval=100.0,
		 incr=None, init=None, width=150, height=20, withValue=1,
		 immediate=1, left=10, right=10 ):

Usage:
    After instanciating a Slider 's' one needs to pack() it s.frame.pack().
    This is not done automatically in order to provide the user with some
    flexibility with respect to the slider's placement.
    One or more callback function(s) also need to be registered. A callback
    function takes one argument that will hold the slider's value.

Methods:
    slider.Set(value): set the slider programatically to a value
    slider.Get(): get the current value
    slider.AddCallback(func): register a  callback function
    
Examples:
    def MyCallback1(color):
	print 'in MyCallback1: ',color

    def MyCallback2(color):
	print 'in MyCallback1: ',color

    sl1 = Slider(None, label='Continuous slider', height=40,
		 minval=0.1, maxval = 10.0, init=5.0)
    sl1.frame.pack()
    sl1.AddCallback(MyCallback1)

    sl2 = Slider(None, label='multiple of 1.5 slider',
		 minval=0.0, maxval = 10.0, incr=1.5)
    sl2.frame.pack()
    sl2.AddCallback(MyCallback2)

    sl3 = Slider(None, label='No value slider', minval=0.1,
		 maxval = 1.0, init=1.0, withValue=0)
    sl3.frame.pack()

    sl4 = Slider(None, label='Non immediate slider', minval=5.,
		 maxval = 25.0, init=10.0, immediate=0, left=50, right=5)
    sl4.frame.pack()
    sl4.AddCallback(MyCallback1)
    sl4.AddCallback(MyCallback2)
    sl1.Set(1.0)
    sl2.Set(5.2)
    sl3.Set(100.0)
"""
#
# $Header: /mnt/raid/services/cvs/DejaVu2/Tk/Slider.py,v 1.1.1.1.4.1 2017/07/13 22:22:07 annao Exp $
#
# $Id: Slider.py,v 1.1.1.1.4.1 2017/07/13 22:22:07 annao Exp $
#

#from Tkinter import *
import Tkinter
from DejaVu2.EventHandler import CallbackFunctions
import types

class Slider(CallbackFunctions):
    """Class for a simple slider"""


    def Callbacks(self):
	"""Implement call to all callbacks"""

	for f in self.callbacks:
            if self.lookup:
                f(self.lookup[int(round(self.val))])
            else:
                f(self.val)


    def MoveCursor(self, event):
	"""Callback function for left mouse button"""

	x = event.x - self.left
	self._MoveCursor(x)


    def DrawCursor(self):
	"""Update graphics representatin of the cursor"""

	# compute position
	x = (self.val-self.min)*self.cst
	deltax = x - self.lastx

	if self.withValue:
            #self.draw.itemconfig(self.valueLabel, text = str(self.val) )
            if self.lookup:
                val = self.lookup[int(round(self.val))]
            else:
                val = self.val
	    self.draw.itemconfig(self.valueLabel,
                                 text = (self.labelformat)%val )
	self.draw.move('cursor', deltax, 0 )
	self.lastx = self.lastx + deltax


    def _MoveCursor(self, x):
	"""Compute value and move the cursor to the new position"""

	# compute the new value
        #minType = type(self.min)
	val = self.min+self.cst1 * x
	if self.incrCanvas:
            #valType = type(val)
	    val = round(val/self.incr)*self.incr
            val = eval(self.cursortype+'('+str(val)+')')
            #if minType is types.IntType:
                #val = int(val)
	if val<self.min: val = self.min
	elif val>self.max: val = self.max

	# check if redraw and callback are needed
	if self.val != val:
	    self.val = val
	    self.DrawCursor()
	    if self.immediate: self.Callbacks()


    def Set(self, val, update=1):
	"""Set the cursor"""
	if val<self.min: val = self.min
	elif val>self.max: val = self.max
	if self.incrCanvas:
            #valType = type(val)
	    val = round(val/self.incr)*self.incr
            val = eval(self.cursortype+'('+str(val)+')')
            #if valType is types.IntType:
                #val = int(val)
	self.val = val
	self.DrawCursor()
	if update: self.Callbacks()
	return self.val

    def SetMin(self, val):
        """Set the minimum value of the slider"""
        if val>self.max: return
        self.min = val
        self._ComputeCst()
        self.DrawCursor()

    def SetMax(self, val):
        """Set the maximum value of the slider"""
        if val<self.min: return
        self.max = val
        self._ComputeCst()
        self.DrawCursor()
        
    def Get(self):
	"""Get the slider's value"""
        if self.lookup:
            return self.lookup[int(round(self.val))]
        else:
            return self.val


    def _ComputeCst(self, event=None):
	"""Conversion factor between values and pixels"""
	self.cst = (self.right-self.left) / (self.max-self.min)
	self.cst1 = 1.0 / self.cst


    def __init__(self, master=None, label=None, minval=0.0, maxval=100.0,
		 incr=None, init=None, width=150, height=20, withValue=1,
		 immediate=1, left=10, right=10 , labelformat = '%4.2f',
                 cursortype='float', lookup=None):

	CallbackFunctions.__init__(self)
	self.frame = Tkinter.Frame(master)
	self.lastx=0
	self.immediate = immediate
        self.labelformat = labelformat
        self.cursortype = cursortype
        
	if not label: label = '' #'slider'
	fnt='-*-helvetica-medium-r-narrow-*-*-120-*-*-*-*-*-*'
	self.label = Tkinter.Label(self.frame, text=label, font=fnt)
	self.label.pack(side=Tkinter.LEFT)
	if withValue: height = max(30, height)
	self.withValue = withValue

	self.draw = Tkinter.Canvas(self.frame, width=width, height=height,
			   relief=Tkinter.SUNKEN)
	self.width = width
	self.height = height
	self.left = left
	self.right = width-right

        # MS May 1st add support for discrete set of values
        self.lookup = lookup # can be set to a list of values indexed by
                             # int(value)
        if lookup is not None:
            # force min and max to provide proper indices
            minval = 0.0
            maxval = len(lookup)-1
            incr = 1
            
	self.max=maxval
	self.min=minval
	self._ComputeCst()
	self.incr = incr
	self.incrCanvas = None

	if incr:
	    self.incrCanvas = round(incr*self.cst)

	if withValue: m = self.height / 2
	else: m = int(self.height * 0.75)
	self.middle = m

	self.draw.create_line( self.left, m, self.right, m,
			       width=2, fill='black')

	y = m-10 # 10 is the cursor's height
	l = self.left
	self.cursor = self.draw.create_polygon( l, m, l+5, y, l-5, y, l, m,
						fill='blue', tag='cursor')

	if withValue:
	    y = self.middle+10
	    self.valueLabel = self.draw.create_text( l, y, text= str(minval),
						     font=fnt, tag='cursor')

	self.draw.pack(side = Tkinter.RIGHT)

	if init is None:
            init = self.min
        if self.lookup:
            init = self.lookup.index(init)
	self.Set(init)
	Tkinter.Widget.bind(self.draw, "<1>", self.MoveCursor)
	Tkinter.Widget.bind(self.draw, "<B1-Motion>", self.MoveCursor)
	if not immediate:
	    Tkinter.Widget.bind(self.draw, "<ButtonRelease-1>", self.MouseUp)


if __name__ == '__main__':
    def MyCallback(color):
	print color

    def MyCallback2(color):
	print 'hello'

    sl1 = Slider(None, label='Continuous slider', height=40,
		 minval=0.1, maxval = 10.0, init=5.0)
    sl1.frame.pack()
    sl1.AddCallback(MyCallback)

    sl2 = Slider(None, label='multiple of 1.5 slider',
		 minval=0.0, maxval = 10.0, incr=1.5)
    sl2.frame.pack()
    sl2.AddCallback(MyCallback2)

    sl3 = Slider(None, label='No value slider', minval=0.1,
		 maxval = 1.0, init=1.0, withValue=0)
    sl3.frame.pack()

    sl4 = Slider(None, label='Non immediate slider', minval=5.,
		 maxval = 25.0, init=10.0, immediate=0, left=50, right=5)
    sl4.frame.pack()
    sl4.AddCallback(MyCallback)
    sl4.AddCallback(MyCallback2)
    #sl1.Set(1.0)
    sl2.Set(5.2)
    #sl3.Set(100.0)

    #slider for discrete values
    sl5 = Slider(None, label='discrete slider',
                 labelformat='%4d', immediate=0,
                 lookup=[10, 15, 25, 46, 78, 99] )
    sl5.frame.pack()

    # slider for non numeric values
    sl6 = Slider(None, label='string slider',
                 labelformat='%s', immediate=0,
                 lookup=['Pierre', 'Paul', 'Jean', 'Jessica'])
    sl6.frame.pack()
    
