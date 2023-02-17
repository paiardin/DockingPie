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
# Copyright: M. Sanner and TSRI 2018
#
#########################################################################
#

import sys, os, random, numpy, weakref

from PySide import QtGui, QtCore

from mglutil.util.callback import CallbackFunction
from DejaVu2.colorTool import RGBARamp, RedWhiteBlueARamp, RedWhiteARamp, WhiteBlueARamp, GreyscaleARamp

class HLSAColorMapView(QtGui.QWidget):
    """
    class to display a color ramp stored internally in Hue Saturation Luminosity Opacity format
    NOTES:
        The widget can not be shrunk in Y below the bumber of color entries in the map
        This widget is used by the HLSAColorMapEditor class
        
    usage:
        cmv = HLSAColorMapView()
        cmv.setRamp(HSLARamp)
    """
    def __init__(self, editor, parent=None, name='NoName', width=0, height=0, offset=0):
        """constructor """
        #print 'AAA', self, name, parent is None, parent
        super(HLSAColorMapView, self).__init__(parent)
        self._editor = weakref.ref(editor)
        self._ramp = None # reference to the editor's HSLA ramp
        self._name = name
        self._height = height # will be width of drawable (including offset)
        self._width = width  # height of drawable
        self._offset = offset
        self._rectangles = None   # self._rectangles[ramp index] = QRect
        
    def setRamp(self, ramp):
        """Set the reference to the editor's HLSA ramp and enforce size"""
        self._ramp = ramp
        #self._height = len(ramp)
        #print 'resizing', self._width, self._height
        #print 'SIZE1', len(ramp), self.size(), self._name
        if len(ramp) > self._height:
            #print "RESIZING from %d to %d "%(self._height,len(ramp))
            self.resize(QtCore.QSize(self._width, len(ramp)))
        else:
            self.update()
        #print 'SIZEAA', len(ramp), self.size(), self._name, self.minimumSize()
        #if sys.platform!='darwin':
        self.setMinimumSize( QtCore.QSize(30, self._height))
        #print 'SIZE', len(ramp), self.size(), self._name, self.minimumSize()
            
    def resizeEvent(self, e):
        """handle resize events"""
        #if self._ramp is None : return
        w, h = e.size().toTuple()
        if self._ramp is not None:
            if h>len(self._ramp):
                self._height = h
            else:
                self.resize(QtCore.QSize(w, len(self._ramp)))
        self._width = w
        #print 'RESIZE EVENT', self, self._name, w, h, self._height, self.minimumSize()
        self.update()

    def paintEvent(self, e):
        """draw the color ramp"""
        if self._ramp is None: return
        qp = QtGui.QPainter()
        qp.begin(self)
        size = self.size()

        # draw background cross diagonals
        r1 = QtCore.QRect(QtCore.QPoint(0,0),
                          QtCore.QPoint(size.width(),size.height()))
        brush = QtGui.QBrush(QtCore.Qt.gray, QtCore.Qt.DiagCrossPattern)
        qp.fillRect( r1, brush )
        
        # draw the color ramp
        l = len(self._ramp)
        linesPerUnit = size.height()/float(l)
        n = 0
        fst=0
        for h,s,l,a in self._ramp:
            col = QtGui.QColor.fromHsl(h,s,l,a)
            qp.setPen( QtCore.Qt.transparent)
            qp.setBrush( col)
            lst = int((n+1)*linesPerUnit)
            r1 = QtCore.QRect(QtCore.QPoint(0,fst),
                              QtCore.QPoint(size.width(),lst))
            fst = lst+1
            qp.drawRect( r1 )
            n += 1
        qp.end()

class HLSAColorMapEditorBase(HLSAColorMapView):
    """Base class for widgets used to display chanels of of HSLA ramp"""
    mousePress = QtCore.Signal()
    mouseRelease = QtCore.Signal()
    
    def __init__(self, editor, parent=None, name='NoName', offset=10):
        super(HLSAColorMapEditorBase, self).__init__(editor, parent=parent, name=name, offset=offset)
        self._entriesForRect = {} # {QRect: [ramp indices for this rect]}
        self._lastx = None
        self._lasty = None
        self._mouseButtonPressed = False

    def mousePressEvent(self, me):
        self.mousePress.emit()
        self._mouseButtonPressed = True
        self._lastx, self._lasty = me.pos().toTuple()
        self._editor()._mapChanged = False
        
    def mouseReleaseEvent(self, me ):
        self.mouseRelease.emit()
        self._mouseButtonPressed = False

    def getValidCoords(self, xe, ye, startFrom=None):
        """clip event coordinates to x,y values in the drawable and return the list of
        rectangles covered between the current y values and th elast y value to avoid
        skipping lines when moving fast.
        starFrom is used by HLSAColorMapChanelEditor.mouseReleaseEvent"""
        
        if self._mouseButtonPressed:
            if startFrom is None:
                startFrom = self._lasty
            if xe < self._offset: xe = self._offset
            elif xe >= self._width: xe = self._width
            if ye < 0: ye = 0
            elif ye >= self._height: ye = self._height-1
            if ye > self._lasty:
                r1Ind = self._YToRectInd[startFrom]
                r2Ind = self._YToRectInd[ye]
            else:
                r1Ind = self._YToRectInd[ye]
                r2Ind = self._YToRectInd[startFrom]
            #if r1Ind>=self._height or r2Ind >= self.height:
            #    print "XXXXXXX", me.pos(), self._lasty, r1Ind, r2Ind 
            return xe, ye, r1Ind, r2Ind
        else:
            return None, None, None, None
        
    def mouseMoveEvent(self, me):	
        raise
    
class HLSAColorMapChanelEditor(HLSAColorMapEditorBase):
    """Class for widgets used to display chanels of of HSLA ramp
    mouse left click to alter the map and Shift-LEFT-MOVE to draw lines
    """
    chanelChange = QtCore.Signal()
    _channelIndices = {
        'hue': (0, 299),
        'saturation': (1, 255),
        'luminosity': (2,255),
        'opacity': (3,255)
        }
    
    def __init__(self, chanel, editor, parent=None, offset=10):
        super(HLSAColorMapChanelEditor, self).__init__(editor, parent=parent, name=chanel, offset=offset)
        assert chanel in self._channelIndices.keys()
        self._chanel = chanel
        self._chanelInd, self._maxi = self._channelIndices[chanel]
        self._drawLinePos = None # with be (x,y) is SHIFT is pressed

    def mouseReleaseEvent(self, me ):
        if self._drawLinePos is not None:
            # perform interpolation and update
            x , y = self._drawLinePos
            dy = y - self._lasty
            signDy = dy/int(abs(dy))
            dx = (x - self._lastx)/float(abs(dy))
            #print 'RELEASE', self._lastx, self._lasty, 'to', self._drawLinePos, dx, dy,signDy 
            for n in range(int(abs(dy))+1):
                _xe = int(self._lastx + n*dx)
                _ye = self._lasty + n*signDy
                xe, ye, r1Ind, r2Ind = self.getValidCoords(_xe, _ye,
                                                           startFrom=min(_ye, self._height-1))
                #print  _xe, _ye, xe, ye, r1Ind, r2Ind
                if xe is not None:
                    self._updateRamp(xe, ye, r1Ind, r2Ind)
            self.update()
            self._lasty = None
            self._drawLinePos = None
        # this call will reset self._mouseButtonPressed so need to be called after
        super(HLSAColorMapChanelEditor, self).mouseReleaseEvent(me)
                
    def paintEvent(self, e):
        if self._ramp is None: return
        qp = QtGui.QPainter()
        qp.begin(self)
        qp.setRenderHints(qp.Antialiasing|qp.SmoothPixmapTransform)
        size = self.size()

        # draw background cross diagonals
        if self._chanelInd==3:
            r1 = QtCore.QRect(QtCore.QPoint(0,0),
                            QtCore.QPoint(size.width(),size.height()))
            brush = QtGui.QBrush(QtCore.Qt.gray, QtCore.Qt.DiagCrossPattern)
            qp.fillRect( r1, brush )

        # draw the color ramp
        linesPerUnit = size.height()/float(len(self._ramp))
        fst = 0 # first line index in drawable
        self._rectangles = []
        self._entriesForRect = {}
        self._YToRectInd = [-1]*self._height
        #print "in paint:", len(self._ramp)
        for entry, hsla in enumerate(self._ramp):
            lst = int((entry+1)*linesPerUnit) # last line to draw
            w = hsla[self._chanelInd]*self._width/float(self._maxi)
            qp.setPen( QtCore.Qt.transparent)
            h,s,l,a = hsla
            if self._chanelInd==3:
                qp.setBrush(QtGui.QColor.fromHsl(255,255,255,a))
            else:
                qp.setBrush(QtGui.QColor.fromHsl(h,s,l))
            r1 = QtCore.QRect(QtCore.QPoint(0,fst),
                              QtCore.QPoint(self._offset+w,lst))
            qp.drawRect( r1 )
            self._rectangles.append(r1)
            self._entriesForRect[r1] = entry
            for i in range(fst, lst+1):
                if i>=self._height: break
                self._YToRectInd[i]=entry # for a value of y in widget return the color map entry
            fst = lst+1
            #print 'ASSA', self._entriesForRect[r1]
        #print 'FAFA', size.height(), len(self._ramp), linesPerUnit, len(self._entriesForRect), set(self._YToRectInd), len(self._rectangles)
        #print self._YToRectInd
        if self._drawLinePos:
            x, y = self._drawLinePos
            qp.setPen(QtGui.QPen(QtCore.Qt.black, 2, QtCore.Qt.SolidLine))
            qp.drawLine(self._lastx, self._lasty, x, y)

        qp.end()

        if not self._drawLinePos:
            self.chanelChange.emit()
        #print 'FUGU 1 width %d height %d nblines %d maxi %d'%(self._width, self._height, len(self._YToRectInd), self._maxi)
    #
    # Handles mouse move events for the connect widget.
    #
    def _updateRamp(self, xe, ye, r1Ind, r2Ind):
        #print 'height %d width %d x %d y %d'%(self._height, self._width, xe, ye)
        for rInd in [r1Ind]+range(r1Ind,r2Ind):
            newval = int(self._maxi*(xe-self._offset)/(self._width-self._offset))
            #print 'FOGO', xe, x0, dx, self._lastx, r1Ind, r2Ind, newval, newval1
            rect = self._rectangles[rInd]
            self._ramp[self._entriesForRect[rect]][self._chanelInd] = newval
        
    def mouseMoveEvent(self, me):
        modifiers = QtGui.QApplication.keyboardModifiers()
        if modifiers == QtCore.Qt.ShiftModifier: # draw a line
            self._drawLinePos = me.pos().toTuple()
            self.update()
        else: # redraw map for each event
            self._drawLinePos = None
            xe, ye = me.pos().toTuple()
            xe, ye, r1Ind, r2Ind = self.getValidCoords(xe, ye)
            if xe is not None:
                self._updateRamp(xe, ye, r1Ind, r2Ind)
                self.update()
                self._lasty = ye

        self._editor()._mapChanged = True

class LabeledEntry(QtGui.QWidget):
    def __init__(self, parent=None):
        super(LabeledEntry, self).__init__(parent)
        layout = QtGui.QHBoxLayout()
        self._label =  QtGui.QLabel("map size: ")
        layout.addWidget(self._label)
        self._lineEdit = QtGui.QLineEdit()
        self._onlyInt = QtGui.QIntValidator(2, 1024)
        self._lineEdit.setValidator(self._onlyInt)
        layout.addWidget(self._lineEdit)
        self.setLayout(layout)

class HLSAColorMapEditor(QtGui.QWidget):
    """
    Widget for viewing and editing color maps using Hue, Saturation, Luminosity and opacity chanes
    the values in the chanel editor can be modified usuing the mouse, either freely
    or in straight lines when SHIFT is pressed.

    The widget supports undoing and redoing operations.
    The number of colors in the map can be changed using through the edit menu.
    Predefined maps can be set from the Maps menu and maps can be saved and restored o and from file.
    The chanel selector can be radio buttons or a combobox
    
    Usage:
        from DejaVu2.colorTool import RGBARamp
        ramp = RGBARamp(256)
        ramp[:,3]=0.50
        app = QtGui.QApplication(sys.argv)
        editor = HLSAColorMapEditor() # create the editor
        editor.setRampFromRGBA(ramp)  # set the color ramp
        rgbRamp = editor.getRGBRamp() # get the ramp in RGB format    
    """
    mapChange = QtCore.Signal()
    
    def __init__(self, parent=None, chanelSelector='radio'):
        assert chanelSelector in ['radio', 'combo']
        super(HLSAColorMapEditor, self).__init__(parent=parent, f=QtCore.Qt.WindowStaysOnTopHint)
        self._ramps = [] # Will be a list of numpy array  of HSLA values ranging 
                         # from 0 to 300 for H and 0 to 255 for S L and A 
        self._rampIndex = -1 # self._ramps[self._rampIndex] is displayed
        self._currentRamp = None # will be self._ramps[self._rampIndex] 
        self.initUI(chanelSelector)
        self._mapChanged = False # set to True upon mouse motion, used to append maps
                                 # to self._ramp

    def getRGBRamp(self):
        return numpy.array([QtGui.QColor.fromHsl(h,s,l,a).getRgbF() for h,s,l,a in self._currentRamp])

    def loadMap_cb(self):
        filename, selfilter = QtGui.QFileDialog().getOpenFileName(
            self, "Color map file name", "", "Color map files(*.clr);; All files (*)")
        if filename:
            self.loadMap(filename)
    
    def loadMap(self, filename):
         f = open(filename, 'r')
         lines = f.readlines()
         ramp = []
         for line in lines:
             if len(line) and line != "\n":
                ramp.append(numpy.fromstring(line, dtype=float, sep=" "))
         if len(ramp):
             self.setRampFromHSLA(numpy.array(ramp))
                 
    def saveMap_cb(self):
        filename, filtr = QtGui.QFileDialog.getSaveFileName(
            self, "Color map file name", "", "Color map file (*.clr)")
        if filename:
            self.saveMap(os.path.splitext(filename)[0] + ".clr")
            
    def saveMap(self, filename):
        f = open(filename, 'w')
        #nramps = len(self._ramps)-1
        #for i, ramp in enumerate(self._ramps):
        #    numpy.savetxt(f, ramp, fmt="%.4f")
        #    if i < nramps:
        #        f.write("\n")
        numpy.savetxt(f, self._currentRamp, fmt="%.4f")
        f.close()
        

    def resizeMap(self):
        newSize = int(self._resizeActionWidget._lineEdit.text())
        if newSize == len(self._currentRamp): return
        self.setSize(newSize)
        
    def setSize(self, size):
        curSize = len(self._currentRamp)
        if curSize==size: return
        newRamp = numpy.zeros((size, 4), 'int')
        frac = float(curSize)/size
        #if size > curSize: # grow the ramp
        for i in range(size):
            newRamp[i] = self._currentRamp[int(i*frac)]
            newRamp[-1] = self._currentRamp[-1]
        #else: # Shrink the ramp
        #    newRamp = self._currentRamp[:size]
        self.setRampFromHSLA(newRamp)

    def _mkRamp(self, name, length):
        if name =='RGB':
            self.setRampFromRGBA(RGBARamp(length))
        elif name=='Red White Blue':
            self.setRampFromRGBA(RedWhiteBlueARamp(length))
        elif name=='Red White':
            self.setRampFromRGBA(RedWhiteARamp(length))
        elif name=='WhiteBlue':
            self.setRampFromRGBA(WhiteBlueARamp(length))
        elif name=='Greyscale':
            self.setRampFromRGBA(GreyscaleARamp(length))

    def createActions(self):
        self._undoAction = QtGui.QAction(
            "undo", self,
            shortcut="Ctrl+z",
            statusTip="undo last change to color map",
            triggered=self.undo)
        self._undoAction.setDisabled(True)

        self._redoAction = QtGui.QAction(
            "redo", self,
            shortcut="Ctrl+r",
            statusTip="redo last undone change",
            triggered=self.redo)
        self._redoAction.setDisabled(True)
        
        #self._resizeActionWidget = LabeledSpinBox(self)
        self._resizeActionWidget = LabeledEntry()
        self._resizeAction = QtGui.QWidgetAction(self._resizeActionWidget)
        self._resizeAction.setDefaultWidget(self._resizeActionWidget)
        self._resizeActionWidget._lineEdit.returnPressed.connect(self.resizeMap)
        
    def createMenuBar(self):
            self._menu = QtGui.QMenuBar(self)
            self._menu.setNativeMenuBar(False)

            fileMenu = self._menu.addMenu("File")
            #fileMenu.addAction("New")
            fileMenu.addAction("load map", self.loadMap_cb)
            fileMenu.addAction("save map", self.saveMap_cb)

            editMenu = self._menu.addMenu("Edit")
            editMenu.addAction(self._undoAction)
            editMenu.addAction(self._redoAction)
            editMenu.addAction(self._resizeAction)

            editMenu = self._menu.addMenu("Maps")
            for rampName in ['RGB', 'Red White Blue', 'Red White', 'WhiteBlue', 'Greyscale']:
                sub = editMenu.addMenu(rampName)
                for length in [8, 16, 32, 64, 128, 256]:
                    cb = CallbackFunction(self._mkRamp, rampName, length)
                    sub.addAction('%d colors'%length, cb)

            # for some reason using the widget w does not show up in the menu bar
            #w = QtGui.QWidget()
            #l = QtGui.QHBoxLayout()
            #lab = QtGui.QLabel("size:")
            #l.addWidget(lab)

            #w1 = self._mapSizeSpinBox = QtGui.SQpinBox()
            #w1.setRange(2, 1024)
            #w1.valueChanged.connect(self.resizeMap)
            #self._menu.setCornerWidget(w1)
            self._sizeLab = QtGui.QLabel("#colors Na")
            self._menu.setCornerWidget(self._sizeLab)
            #l.addWidget(w1)
            #w.setLayout(l)
            #self._menu.setCornerWidget(w)

            self.mainLayout.addWidget(self._menu)

    def initUI(self, chanelSelector):
        self.setGeometry(300, 300, 300, 300)
        self.setWindowTitle('ColorMap Editor')
        self.mainLayout = QtGui.QVBoxLayout(self)
        self.createActions()
        self.createMenuBar()

        if chanelSelector=='combo':
            w = self.chanelComboBox = QtGui.QComboBox()
            w.addItem("hue")
            w.addItem("saturation")
            w.addItem("luminosity")
            w.addItem("opacity")
            self.mainLayout.addWidget(w)
        else:
            l = QtGui.QHBoxLayout(self)
            bg = QtGui.QButtonGroup()
            for name in ['hue', 'sat.', 'lum.', 'opa.']:
                b = QtGui.QPushButton(name, self)
                bg.addButton(b)
                l.addWidget(b)
            self.mainLayout.addLayout(l)
            
        self.editorsLayout = QtGui.QHBoxLayout()
        w = self.editorsTab = QtGui.QStackedWidget()

        if chanelSelector=='combo':
            self.chanelComboBox.activated.connect(w.setCurrentIndex)
        else:
            bg.buttonReleased.connect(w.setCurrentIndex)

        self.hueEditor = HLSAColorMapChanelEditor('hue', self)
        w.addWidget(self.hueEditor)

        self.satEditor = HLSAColorMapChanelEditor('saturation', self)
        w.addWidget(self.satEditor)

        self.lumEditor = HLSAColorMapChanelEditor('luminosity', self)
        w.addWidget(self.lumEditor)

        self.opaEditor = HLSAColorMapChanelEditor('opacity', self)
        w.addWidget(self.opaEditor)

        self.cmView = HLSAColorMapView(self)

        self.editorsLayout.addWidget(w)
        self.editorsLayout.addWidget(self.cmView, width=30)
        self.mainLayout.addLayout(self.editorsLayout)
        self.setLayout(self.mainLayout)
        self.show()

        # use signals to update the color map view
        self.hueEditor.chanelChange.connect(self.cmView.update)
        self.satEditor.chanelChange.connect(self.cmView.update)
        self.lumEditor.chanelChange.connect(self.cmView.update)
        self.opaEditor.chanelChange.connect(self.cmView.update)

        # use signal to save maps for undo
        self.hueEditor.mouseRelease.connect(self.mouseRelease)
        self.satEditor.mouseRelease.connect(self.mouseRelease)
        self.lumEditor.mouseRelease.connect(self.mouseRelease)
        self.opaEditor.mouseRelease.connect(self.mouseRelease)

    def undo(self):
        if self._rampIndex > 0:
            self._rampIndex -= 1
            self._activateRamp(self._rampIndex)

        if self._rampIndex == 0:
            self._undoAction.setDisabled(True)

        if self._rampIndex<len(self._ramps)-1:
            self._redoAction.setDisabled(False)
            
    def redo(self):
        if self._rampIndex < len(self._ramps):
            self._rampIndex += 1
            self._activateRamp(self._rampIndex)
        if self._rampIndex==len(self._ramps)-1:
            self._redoAction.setDisabled(True)
            
    def mouseRelease(self):
        if self._mapChanged:
            self._addRampTohistory(self._currentRamp)
            self._activateRamp(self._rampIndex)
            self.mapChange.emit()
            #ramp = self.getRGBRamp()
            #print "mouse release:", len(ramp), ramp

    def _addRampTohistory(self, ramp):
        ramp = numpy.array(ramp)
        if min(ramp.flatten())<0:
            print 'WARNING: HLSAColorMapEditor:_addRampTohistory, ramp contains negative values. Clamping at 0'
            #import pdb; pdb.set_trace()
        if max(ramp[:,0])>299:
            print 'WARNING: HLSAColorMapEditor:_addRampTohistory, ramp contains values hue values larger than 299. Clamping at 299'
        if max(ramp[:,1:].flatten())>255:
            print 'WARNING: HLSAColorMapEditor:_addRampTohistory, ramp contains values saturation luminosity and opacity values larger than 255. Clamping at 255'
        ramp[:,0] = numpy.clip(ramp[:,0], 0, 299)
        ramp[:,1:] = numpy.clip(ramp[:,1:], 0, 255)
        #print 'ADDRAMP_A', self._rampIndex, len(self._ramps)
        if self._rampIndex==len(self._ramps)-1:
            # we are pointing to the last map (i.e. no undo was performed)
            # so we append the current map
            self._ramps.append(ramp)
        else:
            # forget what is after _undoIndex
            self._ramps = self._ramps[:self._rampIndex+1] + [self._currentRamp]

        self._rampIndex += 1
        #print 'ADDRAMP_B', self._rampIndex, len(self._ramps)
        if self._rampIndex > 0: # enable undo
            self._undoAction.setDisabled(False)
        #print '_addRampTohistory', self._rampIndex, len(self._ramps)

    def _activateRamp(self, index):
        #print 'ACTIVATE ramp', index
        ramp = self._currentRamp = self._ramps[index].copy()
        self.hueEditor.setRamp(ramp)
        self.satEditor.setRamp(ramp)
        self.lumEditor.setRamp(ramp)
        self.opaEditor.setRamp(ramp)
        self.cmView.setRamp(ramp)
        self._sizeLab.setText("#colors %d"%len(ramp))
        self._resizeActionWidget._lineEdit.setText(str(len(ramp)))
        self.mapChange.emit()
        #if len(ramp) > self._height:
        #    self.resize(self._width, self._height)
        #else:
        #    self.update()
        #self.setMinimumSize(30, self._height)
        
    def setRampFromHSLA(self, hslaRamp):
        assert set([len(x) for x in hslaRamp])==set([4])
        resize = self._currentRamp is None or len(self._currentRamp)!=len(hslaRamp)
        self._addRampTohistory(hslaRamp)
        self._activateRamp(self._rampIndex)
        
    def setRampFromRGBA(self, rgbaRamp):
        assert set([len(x) for x in rgbaRamp])==set([4])
        assert numpy.min(rgbaRamp) >= 0.0
        if numpy.max(rgbaRamp)==1.0:
            ramp = []
            for r,g,b,a, in rgbaRamp:
                h,s,l,a = QtGui.QColor.fromRgbF(r,g,b,a).getHslF()
                if h<0:
                    h += 1.0
                ramp.append( (int(h*299), int(s*255), int(l*255), int(a*255)))
        else:
            ramp = [QtGui.QColor.fromRgbF(r,g,b,a).getHsl() for r,g,b,a, in rgbaRamp]
        self.setRampFromHSLA(ramp)
        
if __name__ == '__main__':
    from DejaVu2.colorTool import RGBARamp
    ramp = RGBARamp(64)
    ramp[:,3]=0.50
    app = QtGui.QApplication(sys.argv)
    editor = HLSAColorMapEditor(chanelSelector='combo')
    editor.setRampFromRGBA(ramp)
    rgbRamp = editor.getRGBRamp()
    #print rgbRamp
    sys.exit(app.exec_())
