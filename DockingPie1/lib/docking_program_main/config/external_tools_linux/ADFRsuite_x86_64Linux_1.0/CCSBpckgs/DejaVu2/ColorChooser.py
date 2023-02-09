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
# $Header: /mnt/raid/services/cvs/DejaVu2/ColorChooser.py,v 1.1.1.1.4.1 2017/07/13 22:28:32 annao Exp $
#
# $Id: ColorChooser.py,v 1.1.1.1.4.1 2017/07/13 22:28:32 annao Exp $
#

from copy import deepcopy
import Pmw
import Tkinter, numpy
import Slider, ColorWheel

class ColorChooser:
    """Class of a color picker
"""

# example of valid targetDict
#        lTargetDict = {
#                       'ambient light':   
#                           (self.viewer.lightModel.ambient[:3],
#                            'RGB',
#                            self.viewer.LMColor
#                           ),
#                       'background':
#                           (self.viewer.currentCamera.backgroundColor[:3],
#                            'RGB',
#                            self.viewer.CurrentCameraBackgroundColor
#                           ),
#                      }


    def __init__(self, master=None, targetDict={}, 
                 targetKey=None, gridCfg=None):

        self.master = master
        if hasattr(self.master, 'withdraw'):
            self.master.title('Color Chooser')
            self.master.protocol('WM_DELETE_WINDOW', self.master.withdraw )
        self.frame = Tkinter.Frame(self.master, relief=Tkinter.RIDGE, borderwidth=3)
        self.topFrame = Tkinter.Frame(self.frame)

        # target
        self.targetDict = targetDict
        self.targetKey = None
        lTargetDictKeys = targetDict.keys()[:]
        lTargetDictKeys.sort()
        if targetKey is not None:
            self.comboBoxTarget = Pmw.ComboBox(
                self.master, 
                label_text='target:',
                labelpos='w',
                #entryfield_value=,
                scrolledlist_items=lTargetDictKeys,
                selectioncommand=self.setTarget
                )
            self.comboBoxTarget.pack(side='top', fill='x')

        self.saveCol = Tkinter.Frame(self.topFrame)
        self.currentColor = Tkinter.Label(self.saveCol, text='Current',
                 width=10, relief = Tkinter.SUNKEN, background = '#FFFFFF',
                 foreground = '#FF988E', borderwidth=3 )
        self.currentColor.pack()
        f = Tkinter.Frame(self.saveCol)
        self.save = Tkinter.Button(f, text='Save', relief = Tkinter.RAISED,
                                   command=self.SaveColor,width=7)
        self.save.pack(side=Tkinter.TOP)
        self.swap = Tkinter.Button(f, text='Swap', relief = Tkinter.RAISED,
                                      command=self.SwapColor,width=7)
        self.swap.pack(side=Tkinter.TOP)
        self.restore = Tkinter.Button(f, text='Restore', relief = Tkinter.RAISED,
                                      command=self.RestoreColor,width=7)
        self.restore.pack(side=Tkinter.TOP)
        f.pack()
        self.savedColor = Tkinter.Label(self.saveCol, text='Saved', 
                width=10, relief = Tkinter.SUNKEN, background = '#FFFFFF',
                foreground = '#FF988E', borderwidth=3 )
        self.savedColor.pack()
        self.saveCol.pack(side = Tkinter.LEFT)

        # target
        if targetKey is not None:
            startingColorRGB = self.targetDict[targetKey][0]
        else:
	    startingColorRGB = (1,1,1,1)

        self.hsWheel = ColorWheel.ColorWheel(self.topFrame, circles=6, 
                                             stripes=30,
                                             startingColorRGB=startingColorRGB )
        self.hsWheel.AddCallback(self.UpdateCurrentColor)
        self.hsWheel.canvas.pack(side = Tkinter.LEFT)
        self.topFrame.pack()

        self.savedHsvColor = deepcopy(self.hsWheel.Get())
        if len(startingColorRGB) == 4:
            self.savedHsvColor[3] = startingColorRGB[3]

        coltk = self.hsWheel.Get(mode='TkRGB')
        self.currentColor.config( background = coltk )
        self.savedColor.config( background = coltk )

        self.value = Slider.Slider(self.frame, label='Value', immediate=1,
                                   minval=0.0, maxval = 1.0, init=1.0)

        self.value.frame.pack()
        self.value.AddCallback(self.NewValue)

        self.alpha = Slider.Slider(self.frame, label='Alpha', immediate=1,
                                   minval=0.0, maxval = 1.0, init=1.0)
        self.alpha.frame.pack()
        self.alpha.AddCallback(self.NewValue)
        
        if gridCfg:
            self.frame.grid( **gridCfg)
        else:
            self.frame.pack(ipadx=2, ipady=2, padx=2, pady=2)

        # target
        if targetKey is not None:
            self.setTarget(targetKey)


    def setTarget(self, val):
        #print "setTarget", val
        if self.targetKey is not None:
            self.RemoveCallback(self.targetDict[self.targetKey][2])
        self.targetKey = val
        self.comboBoxTarget.setentry(val)
        lTargetTuple = self.targetDict[val]
        self.Set(lTargetTuple[0], mode=lTargetTuple[1])
        self.AddCallback( lTargetTuple[2] )
        

    def AddCallback(self, func):
	"""Add a callback fuction"""

	self.hsWheel.AddCallback(func)


    def RemoveCallback(self, func):
	"""Delete a callback fuction"""

	self.hsWheel.RemoveCallback(func)


    def UpdateCurrentColor(self, color):
	"""Change the color of the current color label"""
	col = self.hsWheel.Get(mode='TkRGB')
	self.currentColor.config( background = col )


    def NewValue(self, event):
	"""Redo the color wheel with a new value, update color boxes"""

	v = self.value.Get()
	hsv = self.hsWheel.Get('HSV')

        if len(hsv) == 4:
            hsv[3] = self.alpha.Get() 

	if v != hsv[2]:
	    hsv[2]=v
	    self.hsWheel.Set(hsv, 'HSV')
	    self.hsWheel.Callbacks()
        elif len(hsv) == 4:
            self.hsWheel.Set(hsv, 'HSV')
            self.hsWheel.Callbacks()


    def SaveColor(self):
	"""Save current color into second color box"""

	self.savedHsvColor = deepcopy(self.hsWheel.Get())
        self.savedHsvColor[3] = self.alpha.Get()
	col = self.hsWheel.Get(mode='TkRGB')
	self.savedColor.config( background = col )


    def RestoreColor(self):
	"""Restore color from second color box
"""
	self.Set(deepcopy(self.savedHsvColor), 'HSV')


    def SwapColor(self):
	"""Exchange colors in first and second box"""

	coltk = self.hsWheel.Get(mode='TkRGB')
	self.savedColor.config( background = coltk )
	saved_col = deepcopy(self.savedHsvColor)
	self.savedHsvColor = deepcopy(self.hsWheel.Get())
        self.savedHsvColor[3] = self.alpha.Get()
	self.Set(saved_col, 'HSV')


    def Set(self, col, mode='RGB', run=True):
	"""Set the color chooser to a given RGB or HSV triplet (0. to 1.)
"""
	assert mode in ('HSV', 'RGB')
	assert len(col) in (3,4)
	self.hsWheel.Set( col, mode, run )
        self.value.Set( self.hsWheel.Get('HSV')[2], update=0 )
        if len(col) == 4:
            self.alpha.Set( col[3], update=0 )
	self.UpdateCurrentColor(None)


    def Wysiwyg(self, onOff):
	"""Toggle Wysiwyg mode for color wheel. When this flag iset,
	the colors of the wheel are recomputed for each new value"""
	assert onOff in (0,1)
	self.hsWheel.Wysiwyg(onOff)


    def showColorChooser(self, target=None):
        """
"""
        if self.master.winfo_ismapped() == 0:
            self.master.deiconify()
        self.master.lift()
        if target is not None:
            self.setTarget(target)


    def hideColorChooser(self, event=None):
        #print "hideColorChooser", self
        if self.master.winfo_ismapped() == 1:
            self.master.withdraw()


    def toggleColorChooser(self, event=None):
        #print "toggleColorChooser", self
        if self.master.winfo_ismapped() == 1:
            self.master.withdraw()
        else:
            self.master.deiconify()
            self.master.lift()


if __name__ == '__main__':
    root = Tkinter.Tk()
    cc = ColorChooser(root)
    def MyCallback(col):
	print 'col', col
    cc.AddCallback(MyCallback)

