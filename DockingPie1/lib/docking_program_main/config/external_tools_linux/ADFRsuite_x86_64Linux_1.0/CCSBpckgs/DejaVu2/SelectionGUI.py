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
# Date: 2000 Author: Michel F. SANNER, Guillaume Vareille
#
#    sanner@scripps.edu
#    vareille@scripps.edu
#
#       The Scripps Research Institute (TSRI)
#       Molecular Graphics Lab
#       La Jolla, CA 92037, USA
#
# Copyright: Michel Sanner and TSRI
#
# revision:
#
#########################################################################
#
# $Header: /mnt/raid/services/cvs/DejaVu2/SelectionGUI.py,v 1.1.1.1.4.1 2017/07/13 22:28:32 annao Exp $
#
# $Id: SelectionGUI.py,v 1.1.1.1.4.1 2017/07/13 22:28:32 annao Exp $
#

import Tkinter
import DejaVu2
import Slider
from mglutil.gui.BasicWidgets.Tk.thumbwheel import ThumbWheel
from mglutil.util.callback import CallBackFunction
from mglutil.util.colorUtil import TkColor


class SelectionGUI:

    def __init__(self, master=None, viewerGui=None):

        if master is None:
            self.master = Tkinter.Toplevel()
            self.master.withdraw()
            self.master.title('Selection Settings')
            self.master.protocol('WM_DELETE_WINDOW', self.master.withdraw )
        else:
            self.master = master

        self.mainFrame = Tkinter.Frame(self.master, relief='ridge', borderwidth=3)

        # CONTOUR WIDTH
        if DejaVu2.enableSelectionContour is True:
            lLabel =Tkinter.Label(self.mainFrame, text='Contour width:')
            lLabel.grid(row=0, column=0, sticky='w')
            self.contourSizeSlider = Slider.Slider(
                                        self.mainFrame,
                                        minval=0,
                                        maxval = 10,
                                        immediate=1,
                                        incr=1,
                                        labelformat = '%d',
                                        cursortype='int'
                                        )
            self.contourSizeSlider.frame.grid(row=0, column=1, columnspan=3)
            def setContourWidthSliderValue_cb(val):
                DejaVu2.selectionContourSize = val
                viewerGui.viewer.Redraw()
            self.contourSizeSlider.AddCallback(setContourWidthSliderValue_cb)
        else:
            lLabel =Tkinter.Label(self.mainFrame, text='Contour selection is disabled')
            lLabel.grid(row=0, column=0, columnspan=3, sticky='w')

        # BACKGROUND COLOR
        lLabel =Tkinter.Label(self.mainFrame, text='Contour color:')
        lLabel.grid(row=1, column=0, sticky='w')
        lFunc = CallBackFunction(viewerGui.colorChooser.showColorChooser,'selection contour')
        self.contourColorButton = Tkinter.Button(self.mainFrame,
                                 #text='Contour color',
                                 width=7,
                                 background=TkColor(DejaVu2.selectionContourColor),
                                 command=lFunc)
        self.contourColorButton.grid(row=1, column=1, columnspan=3)

        # PATTERN SIZE
        lLabel =Tkinter.Label(self.mainFrame, text='Pattern size:')
        lLabel.grid(row=2, column=0, sticky='w')
        
        def setPatternSizeSliderValue_cb(val):
            DejaVu2.selectionPatternSize = val
            viewerGui.viewer.Redraw()

        self.patternSizeThumbwheel = ThumbWheel(
                             self.mainFrame, width=70, height=16, type=int,
                             value=DejaVu2.selectionPatternSize,
                             callback=setPatternSizeSliderValue_cb,
                             continuous=True, oneTurn=10,
                             wheelPad=2, min=0, max=50)
        self.patternSizeThumbwheel.grid(row=2, column=1, columnspan=3)

        self.mainFrame.pack(side='top')

        self.updateWidgets()


    def show(self, event=None):
        #print "showSpinGui", self
        if self.master.winfo_ismapped() == 0:
            self.master.deiconify()
        self.master.lift()


    def hide(self, event=None):
        #print "hideSpinGui", self
        if self.master.winfo_ismapped() == 1:
            self.master.withdraw()


    def toggle(self, event=None):
        #print "toggleSpinGui", self
        if self.master.winfo_ismapped() == 1:
            self.master.withdraw()
        else:
            self.master.deiconify()
            self.master.lift()


    def setWithWidgets(self):
        if DejaVu2.enableSelectionContour is True:
            DejaVu2.selectionContourSize = self.contourSizeSlider.Get()
        DejaVu2.selectionPatternSize = self.patternSizeThumbwheel.Get()


    def updateWidgets(self):
        if DejaVu2.enableSelectionContour is True:
            self.contourSizeSlider.Set(DejaVu2.selectionContourSize)
        self.contourColorButton.configure(background=TkColor(DejaVu2.selectionContourColor))
        self.patternSizeThumbwheel.set(DejaVu2.selectionPatternSize)
