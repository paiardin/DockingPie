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

# $Header: /mnt/raid/services/cvs/DejaVu2/PropertyEditor.py,v 1.1.1.1.4.1 2017/07/13 22:28:32 annao Exp $
#
# $Id: PropertyEditor.py,v 1.1.1.1.4.1 2017/07/13 22:28:32 annao Exp $
#

import Tkinter, numpy
import Slider, ColorWheel, ColorChooser, colorTool
from EventHandler import CallbackFunctions


class MaterialEditor(CallbackFunctions):
    """Class for a material editor"""

    property = [ 'ambient', 'diffuse', 'specular', 'emission' ]

    def Callbacks(self, event=None):
	"""Implement callback functions"""

	if type(event) == type(0.0):
	    for f in self.callbacks:
		f('shininess', event)
	else:
	    tkrgb = self.colorChooser.hsWheel.Get(mode='TkRGB')
	    c = self.currentComponent.get()	
	    self.mat[c] = self.colorChooser.hsWheel.Get(mode='RGB')[:3]
	    self.colw[c].config( background = tkrgb )
	    for f in self.callbacks:
		f(self.property[c], self.mat[c])


    def RestoreColor(self):
	"""Set the wheel cursor to the color of the current component"""
	c = self.currentComponent.get()	
	self.colorChooser.Set(self.mat[c], 'RGB')


    def Set(self, ambi=None, diff=None, spec=None, emis=None, shini=None,
	    mode = 'RGB'):
	"""Set the material Editor to a given material"""

	assert mode in ('HSV', 'RGB')
	c = self.currentComponent.get()	
	if ambi:
	    ambi = colorTool.OneColor(ambi)
	    tkrgb = colorTool.TkColor(ambi[:3])
	    self.colw[0].config( background = tkrgb )
	    self.mat[0][:3] = ambi[:3]
	    if c==0: self.colorChooser.Set( ambi, 'RGB' )
	if diff:
	    diff = colorTool.OneColor(diff)
	    tkrgb = colorTool.TkColor(diff[:3])
	    self.colw[1].config( background = tkrgb )
	    self.mat[1][:3] = diff[:3]
	    if c==1: self.colorChooser.Set( diff, 'RGB' )
	if spec:
	    spec = colorTool.OneColor(spec)
	    tkrgb = colorTool.TkColor(spec[:3])
	    self.colw[2].config( background = tkrgb )
	    self.mat[2][:3] = spec[:3]
	    if c==2: self.colorChooser.Set( spec, 'RGB' )
	if emis:
	    emis = colorTool.OneColor(emis)
	    tkrgb = colorTool.TkColor(emis[:3])
	    self.colw[3].config( background = tkrgb )
	    self.mat[3][:3] = emis[:3]
	    if c==3: self.colorChooser.Set( emis, 'RGB' )
	if shini:
	    self.shini.Set(shini)


    def __init__(self, root=None, colorChooser=None):

	CallbackFunctions.__init__(self)

	self.frame = Tkinter.Frame(root, relief=Tkinter.RIDGE, borderwidth=3)
	self.mat = numpy.ones( (4,3), 'f')
	self.currentComponent = Tkinter.IntVar()
	self.currentComponent.set(0)
	self.colw = [0,]*4
	width = 9
	self.colw[0] = Tkinter.Radiobutton(self.frame, text='Ambient',
			    variable = self.currentComponent, value=0,
			    command=self.RestoreColor, relief = Tkinter.SUNKEN,
			    width=width, borderwidth=3,
			    background = '#FFFFFF', foreground = '#000000',  )
	self.colw[0].grid(row=0,column=0)
	self.colw[1] = Tkinter.Radiobutton(self.frame, text='Diffuse',
			    value=1, variable = self.currentComponent,
			    command=self.RestoreColor, borderwidth=3,
			    width=width, relief = Tkinter.SUNKEN,
			    background = '#FFFFFF', foreground = '#000000' )
	self.colw[1].grid(row=1,column=0)
	self.colw[2] = Tkinter.Radiobutton(self.frame, text='Specular',
			    value=2, variable = self.currentComponent,
			    command=self.RestoreColor, relief = Tkinter.SUNKEN,
			    width=width, borderwidth=3,
			    background = '#FFFFFF', foreground = '#000000' )
	self.colw[2].grid(row=0,column=1)
	self.colw[3] = Tkinter.Radiobutton(self.frame, text='Emissive',
		            value=3, variable = self.currentComponent,
			    command=self.RestoreColor, relief = Tkinter.SUNKEN,
			    width=width, borderwidth=3,
			    background = '#FFFFFF', foreground = '#000000')
	self.colw[3].grid(row=1,column=1)

	self.shini = Slider.Slider(self.frame, label='Shininess', immediate=0,
				   minval=0.0, maxval = 128.0, init=30.0)
	self.shini.frame.grid(row=2,column=0,columnspan=2)
	self.shini.AddCallback(self.Callbacks)

	if not colorChooser:
	    self.colorChooser = ColorChooser.ColorChooser(self.frame)
	    self.colorChooser.frame.grid(row=3, column=0, columnspan=2)
	else:
	    self.colorChooser = colorChooser
	self.colorChooser.AddCallback(self.Callbacks)


class LightColorEditor(CallbackFunctions):
    """Class for a light source color editor"""

    property = [ 'ambient', 'diffuse', 'specular' ]

    def Callbacks(self, event=None):
	"""Implement callback functions"""

	tkrgb = self.colorChooser.hsWheel.Get(mode='TkRGB')
	c = self.currentComponent.get()	
	self.mat[c] = self.colorChooser.hsWheel.Get(mode='RGB')[:3]
	self.colw[c].config( background = tkrgb )
	for f in self.callbacks:
	    f(self.property[c], self.mat[c])


    def RestoreColor(self):
	"""Set the wheel cursor to the color of the current component"""
	c = self.currentComponent.get()	
	self.colorChooser.Set(self.mat[c], 'RGB')


    def Set(self, ambi=None, diff=None, spec=None, mode = 'RGB'):
	"""Set the material Editor to a given material"""

	assert mode in ('HSV', 'RGB')
	c = self.currentComponent.get()	
	if ambi:
	    ambi = colorTool.OneColor(ambi)
	    tkrgb = colorTool.TkColor(ambi[:3])
	    self.colw[0].config( background = tkrgb )
	    self.mat[0][:3] = ambi[:3]
	    if c==0: self.colorChooser.Set( ambi, 'RGB' )

	if diff:
	    diff = colorTool.OneColor(diff)
	    tkrgb = colorTool.TkColor(diff[:3])
	    self.colw[1].config( background = tkrgb )
	    self.mat[1][:3] = diff[:3]
	    if c==1: self.colorChooser.Set( diff, 'RGB' )

	if spec:
	    spec = colorTool.OneColor(spec)
	    tkrgb = colorTool.TkColor(spec[:3])
	    self.colw[2].config( background = tkrgb )
	    self.mat[2][:3] = spec[:3]
	    if c==2: self.colorChooser.Set( spec, 'RGB' )


    def __init__(self, root=None, colorChooser=None):

	CallbackFunctions.__init__(self)

	self.frame = Tkinter.Frame(root, relief=Tkinter.RIDGE, borderwidth=3)
	self.mat = numpy.ones( (3,3), 'f')
	self.currentComponent = Tkinter.IntVar()
	self.currentComponent.set(0)
	self.colw = [0,]*4
	width = 9
	self.colw[0] = Tkinter.Radiobutton(self.frame, text='Ambient',
			    value=0, variable = self.currentComponent,
			    command=self.RestoreColor, relief = Tkinter.SUNKEN,
			    width=width, borderwidth=3,
			    background = '#FFFFFF', foreground = '#000000',  )
	self.colw[0].grid(row=0,column=0)
	self.colw[1] = Tkinter.Radiobutton(self.frame, text='Diffuse',
			    value=1, variable = self.currentComponent,
			    command=self.RestoreColor, borderwidth=3,
			    width=width, relief = Tkinter.SUNKEN,
			    background = '#FFFFFF', foreground = '#000000' )
	self.colw[1].grid(row=1,column=0)
	self.colw[2] = Tkinter.Radiobutton(self.frame, text='Specular',
			    value=2, variable = self.currentComponent,
			    command=self.RestoreColor, relief = Tkinter.SUNKEN,
			    width=width, borderwidth=3,
			    background = '#FFFFFF', foreground = '#000000' )
	self.colw[2].grid(row=0,column=1)

	if not colorChooser:
	    self.colorChooser = ColorChooser.ColorChooser(self.frame)
	    self.colorChooser.frame.grid(row=3, column=0, columnspan=2)
	else:
	    self.colorChooser = colorChooser
	self.colorChooser.AddCallback(self.Callbacks)


if __name__ == '__main__':
    root = Tkinter.Tk()
    root.title( 'Material Editor' )
    me = MaterialEditor(root)
    me.frame.pack()

    lce = LightColorEditor(root)
    lce.frame.pack()
