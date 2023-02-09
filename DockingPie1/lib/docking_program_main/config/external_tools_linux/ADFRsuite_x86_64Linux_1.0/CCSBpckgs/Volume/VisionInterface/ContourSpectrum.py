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
# Date: Mars 2006 Authors: Guillaume Vareille, Michel Sanner
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
# $Header: /mnt/raid/services/cvs/Volume/VisionInterface/ContourSpectrum.py,v 1.15.12.1 2017/07/28 01:09:21 annao Exp $
#
# $Id: ContourSpectrum.py,v 1.15.12.1 2017/07/28 01:09:21 annao Exp $
#

# third party packages 
import types
import Tkinter
import numpy
#import Pmw
#import os
#import types
#import string
#import warnings

# TSRI MGL packages 
from NetworkEditor.items import NetworkNode
from NetworkEditor.widgets import PortWidget
from mglutil.util.callback import CallbackManager
from mglutil.gui.BasicWidgets.Tk.optionsPanel import OptionsPanel
from mglutil.util.misc import ensureFontCase

# current package

class ContourSpectrumNE(NetworkNode):
    """
Input Ports
    contourspectrum: (bound to contourspectrum widget)
    mini: minimum value (optional)
    maxi: maximum value (optional)

Output Ports
    value: value is of type int or float, depending on the contourspectrum settings
"""

    def __init__(self, name='contourspectrum', **kw):
        
        #import pdb;pdb.set_trace()
        
        kw['name'] = name
        apply( NetworkNode.__init__, (self,), kw )
        self.inNodeWidgetsVisibleByDefault = True

        self.widgetDescr['contourspectrum'] = {
            'class':'NEContourSpectrum', 'master':'node', 'size':50,
            'oneTurn':1, 'type':'float', 'lockedOnPort':True,
            'initialValue':0.0,
            'labelCfg':{'text':''}}
        
        self.inputPortsDescr.append(datatype='float', name='contourspectrum')
        self.inputPortsDescr.append(datatype='float', name='mini',
                                    required=False)
        self.inputPortsDescr.append(datatype='float', name='maxi',
                                    required=False)
        self.inputPortsDescr.append(datatype='Grid3D', name='grid', required=False)
        
        self.outputPortsDescr.append(datatype='float', name='value')

        code = """def doit(self, contourspectrum, mini, maxi, grid):
    if contourspectrum is not None:
        w = self.inputPortByName['contourspectrum'].widget
        if w:
            if grid is not None and self.inputPortByName['grid'].hasNewValidData():
                data = grid.data
                origin = numpy.array(grid.origin).astype('f')
                stepsize = numpy.array(grid.stepSize).astype('f')
                self.newgrid3D = numpy.reshape( numpy.transpose(grid.data),
                                          (1, 1)+tuple(grid.data.shape) )   
                from UTpackages.UTisocontour import isocontour         
                gridData = isocontour.newDatasetRegFloat3D( 
                               self.newgrid3D, origin, stepsize)
                sig = [gridData.getSignature(0, 0, 0),
                       gridData.getSignature(0, 0, 1),
                       gridData.getSignature(0, 0, 2),
                       gridData.getSignature(0, 0, 3)]
                w.widget.setSignatures(sig)
            if mini is not None and self.inputPortByName['mini'].hasNewValidData():
                w.configure(min=mini)
            if maxi is not None and self.inputPortByName['maxi'].hasNewValidData():
                w.configure(max=maxi)
        self.outputData(value=contourspectrum)
"""

        self.setFunction(code)


    def afterAddingToNetwork(self):
        NetworkNode.afterAddingToNetwork(self)
        # run this node so the value is output
        self.run()
        self.inputPorts[0].widget.configure = self.configure_NEContourSpectrum
        self.inputPorts[0].widget.widget.setType = self.setType_NEContourSpectrum
        

    def configure_NEContourSpectrum(self, rebuild=True, **kw):
        """specialized configure method for contourspectrum widget"""
        # overwrite the contourspectrum widget's configure method to set the outputPort
        # data type when the contourspectrum is configured
        w = self.inputPorts[0].widget
        apply( NEContourSpectrum.configure, (w, rebuild), kw)
        dtype = kw.pop('type', None)

        if dtype:
            self.updateDataType(dtype)
            

    def setType_NEContourSpectrum(self, dtype):
        """specialized setTyp method for mglutil contourspectrum object"""
        # overwrite the Dial's setType method to catch type changes through
        # the optionsPanel
        from mglutil.gui.BasicWidgets.Tk.Dial import ContourSpectrum
        contourspectrum = self.inputPorts[0].widget.widget
        apply( ContourSpectrum.setType, (contourspectrum, dtype), {})
        if type(dtype) == types.TypeType:
            dtype = dtype.__name__

        self.updateDataType(dtype)
        self.inputPorts[1].setDataType(dtype, makeOriginal= True)
        self.inputPorts[2].setDataType(dtype, makeOriginal= True)

    def updateDataType(self, dtype):
        port = self.outputPorts[0]
        port.setDataType(dtype, tagModified=False)
        if port.data is not None:
            if type(dtype) == types.TypeType:
                port.data = dtype(port.data)                
            else:
                port.data = eval("%s(port.data)"%dtype)



class NEContourSpectrum(PortWidget):
    """NetworkEditor wrapper for ContourSpectrum widget.
Handles all PortWidget arguments and all ContourSpectrum arguments except for value.
    Name:          default:
    callback       None
    continuous     1
    lockContinuous 0
    lockBMin       0
    lockBMax       0
    lockMin        0
    lockMax        0
    lockPrecision  0
    lockShowLabel  0
    lockType       0
    lockValue      0
    min            None
    max            None
    precision      2
    showLabel      1
    size           50
    type           'float'
"""

    configOpts = PortWidget.configOpts.copy()
    configOpts['initialValue'] = {
        'defaultValue':0.0, 'type':'float',
        }

    ownConfigOpts = {
    'callback': {
        'defaultValue':None, 'type': 'None',
        'description':"???",
        },
    'continuous': {
        'defaultValue':True, 'type':'boolean',
        'description':"",
        },
    'lockContinuous': {
        'defaultValue':False, 'type':'boolean',
        'description':"",
        },
    'lockBMin': {
        'defaultValue':False, 'type':'boolean',
        'description':"",
        },
    'lockBMax': {
        'defaultValue':False, 'type':'boolean',
        'description':"",
        },
    'lockMin': {
        'defaultValue':False, 'type':'boolean',
        'description':"",
        },
    'lockMax': {
        'defaultValue':False, 'type':'boolean',
        'description':"",
        },
    'lockOneTurn': {
        'defaultValue':False, 'type':'boolean',
        'description':"",
        },
    'lockPrecision': {
        'defaultValue':False, 'type':'boolean',
        'description':"",
        },
    'lockShowLabel': {
        'defaultValue':False, 'type':'boolean',
        'description':"",
        },
    'lockType': {
        'defaultValue':False, 'type':'boolean',
        'description':"",
        },
    'lockValue': {
        'defaultValue':False, 'type':'boolean',
        'description':"",
        },
    'min': {
        'defaultValue':None, 'type':'float',
        'description':"",
        },
    'max': {
        'defaultValue':None, 'type':'float',
        'description':"",
        },
    'oneTurn': {
        'defaultValue':360., 'type':'float',
        'description':"",
        },
    'precision': {
        'defaultValue':2, 'type':'int',
        'description':"number of decimals used in label",
        },
    'showLabel': {
        'defaultValue':True, 'type':'boolean',
        'description':"",
        },
    'size':{
        'defaultValue': 50, 'min':20, 'max':500, 'type':'int'
        },
    'type': {
        'defaultValue':'float', 'type':'string',
        'validValues': ['float', 'int'],
        'description':"",
        },
        }
    configOpts.update( ownConfigOpts )


    def __init__(self, port, **kw):
        
##          # create all attributes that will not be created by configure because
##          # they do not appear on kw
##          for key in self.ownConfigOpts.keys():
##              v = kw.get(key, None)
##              if v is None: # self.configure will not do anyting for this key
##                  setattr(self, key, self.ownConfigOpts[key]['defaultValue'])

        # get all arguments handled by NEThumbweel and not by PortWidget
        widgetcfg = {}
        for k in self.ownConfigOpts.keys():
            if k in kw:
                widgetcfg[k] = kw.pop(k)

        # call base class constructor
        apply( PortWidget.__init__, ( self, port), kw)

        # create the Dial widget
        
        #from NetworkEditor.spectrum import ContourSpectrumGUI
        self.widget = apply( ContourSpectrum, (self.widgetFrame,), widgetcfg)
        
        self.widget.callbacks.AddCallback(self.newValueCallback)
        # rename Options Panel to port name
        self.widget.opPanel.setTitle("%s : %s"%(port.node.name, port.name) )

        # overwrite right mouse button click
        self.widget.canvas.bind("<Button-3>", self.postWidgetMenu)

        self.widget.canvas.configure(cursor='cross')

        # add menu entry to open configuration panel
        self.menu.insert_command(0, label='Option Panel', underline=0, 
                                 command=self.toggleOptionsPanel)

        # register new callback for widget's optionsPanel Apply button
        # NOTE: idf.entryByName is at this time not built
        for k in self.widget.opPanel.idf:
            name = k.get('name', None)
            if name and name == 'ApplyButton':
                k['command'] = self.optionsPanelApply_cb
            elif name and name == 'OKButton':
                k['command'] = self.optionsPanelOK_cb

        # first set default value, in case we have a min or max, else the
        # node would run
        if self.initialValue is not None:
            self.set(self.widget.type(self.initialValue), run=0)
            
        # configure without rebuilding to avoid enless loop
        apply( self.configure, (False,), widgetcfg)

        self._setModified(False) # will be set to True by configure method


    def configure(self, rebuild=True, **kw):
        # call base class configure with rebuild=Flase. If rebuilt is needed
        # rebuildDescr will contain w=='rebuild' and rebuildDescr contains
        # modified descr
        action, rebuildDescr = apply( PortWidget.configure, (self, False), kw)
            
        # handle ownConfigOpts entries
        if self.widget is not None:
            widgetOpts = {}
            for k, v in kw.items():
                if k in self.ownConfigOpts:
                    if k =='size':
                        action = 'rebuild'
                        rebuildDescr[k] = v
                    else:
                        widgetOpts[k] = v

            if len(widgetOpts):
                apply( self.widget.configure, (), widgetOpts)

        if action=='rebuild' and rebuild:
            action, rebuildDescr = self.rebuild(rebuildDescr)

        elif action=='resize' and rebuild:
            if self.widget and rebuild: # if widget exists
                action = None

        return action, rebuildDescr


    def set(self, value, run=1):
        #print "set NEContourSpectrum"
        self._setModified(True)
        self.widget.setValue(value)
        self._newdata = True
        if run:
            self.scheduleNode()

        
    def get(self):
        return self.widget.get()

        
    def optionsPanelOK_cb(self, event=None):
        # register this widget to be modified when opPanel is used
        self.widget.opPanel.OK_cb()
        self._setModified(True)

    def optionsPanelApply_cb(self, event=None):
        # register this widget to be modified when opPanel is used
        self.widget.opPanel.Apply_cb()
        self._setModified(True)


    def toggleOptionsPanel(self, event=None):
        # rename the options panel title if the node name or port name has
        # changed.
        self.widget.opPanel.setTitle(
            "%s : %s"%(self.port.node.name, self.port.name) )
        self.widget.toggleOptPanel()
        

    def getDescr(self):
        cfg = PortWidget.getDescr(self)
        for k in self.ownConfigOpts.keys():
            if k == 'type': # type has to be handled separately
                _type = self.widget.type
                if _type == int:
                    _type = 'int'
                else:
                    _type = 'float'
                if _type != self.ownConfigOpts[k]['defaultValue']:
                    cfg[k] = _type
                continue
            val = getattr(self.widget, k)
            if val != self.ownConfigOpts[k]['defaultValue']:
                cfg[k] = val
        return cfg



class ContourSpectrum(Tkinter.Frame):
    """This class implements a ContourSpectrum widget.
it is fully a copy/paste from Dial
"""
    def __init__(self, master=None, type='float',
                 labCfg={'fg':'black','side':'left', 'text':None},
                 min=None, max=None, 
                 showLabel=1, value=0.0, continuous=1, precision=2,
                 callback=None, lockMin=0, lockBMin=0, lockMax=0, lockBMax=0,
                 lockPrecision=0,lockShowLabel=0, lockValue=0,
                 lockType=0, lockContinuous=0,  signatures=None, **kw):

        Tkinter.Frame.__init__(self, master)
        Tkinter.Pack.config(self)

        self.callbacks = CallbackManager() # object to manage callback
                                        # functions. They get called with the
                                        # current value as an argument

        # initialize various attributes with default values
        self.height = 100      # widget height
        self.width = 256      # widget height
        self.widthMinusOne = self.width - 1

        self.min = 0                 # minimum value
        self.max = 1                 # maximum value
        self.range = self.max - self.min

        self.precision = 2              # decimal places
        self.minOld = 0.                # used to store old values 
        self.maxOld = 0.
        self.size = 50                  # defines widget size
        self.offsetValue = 0.           # used to set increment correctly
        self.lab = None                 # label
        self.callback = None            # user specified callback
        self.opPanel = None             # option panel widget
        self.value = 0.0                # current value of widget
        self.oldValue = 0.0             # old value of widget
        self.showLabel = 1              # turn on to display label on
        self.continuous = 1             # set to 1 to call callbacks at 
                                        # each value change, else gets called 
                                        # on button release event
        
        self.labCfg = labCfg            # Tkinter Label options
        self.labelFont = (
            ensureFontCase('helvetica'), 14, 'bold')    # label font
        self.labelColor = 'yellow'      # label color
        self.canvas = None              # the canvas to create the widget in

        self.lockMin = lockMin          # lock<X> vars are used in self.lock()
        self.lockMax = lockMax          # to lock/unlock entries in optionpanel
        self.lockBMin = lockBMin
        self.lockBMax = lockBMax
        self.lockPrecision = 0
        self.lockShowLabel = lockShowLabel
        self.lockValue = lockValue
        self.lockType = lockType
        self.lockContinuous = lockContinuous

        # configure with user-defined values
        self.setCallback(callback)     
        self.setContinuous(continuous)

        self.setType(type)              
        self.setPrecision(precision)
        self.setMin(min)
        self.setMax(max)
        self.setShowLabel(showLabel)
        self.setValue(value)
        self.setLabel(self.labCfg)

        if master is None:
            master = Tkinter.Toplevel()

        self.master = master   # widget master

        self.createCanvas(master)

        Tkinter.Widget.bind(self.canvas, "<ButtonPress-1>", self.mouseDown)
        Tkinter.Widget.bind(self.canvas, "<B1-Motion>", self.mouseMove)
        Tkinter.Widget.bind(self.canvas, "<ButtonRelease-1>", self.mouseUp)

        # create cursor
        self.cursorTk = self.canvas.create_line( 0, 0, 0, 0, tags=['cursor'])

        self.increment = 0.0    
        self.incrementOld = 0.
        self.lockIncrement = 0
        self.lockBIncrement = 0
        self.oneTurn = 360.
        self.lockOneTurn = 0
        self.opPanel = OptionsPanel(master = self, title="Slider graph Options")
                
        self.signatures = None # Signature objects fro isocontouring lib
        self.sigData = []      # list of (x,y) values arrays for each signature
        self.maxFun = []       # max Y value in each signature
        self.minFun = []       # min Y value in each signature
        self.yratios = []      # normalization factors
        self.colors = ['red', 'green', 'blue', 'orange']
        self.tkLines = []   # list of Tkids for lines
        if signatures:
            self.setSignatures(signatures)


    def setCallback(self, cb):
        """Set widget callback. Must be callable function. Callback is called
every time the widget value is set/modified"""

        assert cb is None or callable(cb) or type(cb) is types.ListType,\
               "Illegal callback: must be either None or callable, or list. Got %s"%cb
        if cb is None: return
        elif type(cb) is types.ListType:
            for func in cb:
                assert callable(func), "Illegal callback must be callable. Got %s"%func
                self.callbacks.AddCallback(func)
        else:
            self.callbacks.AddCallback(cb)
        self.callback = cb
        

    def toggleOptPanel(self, event=None):
        if self.opPanel.flag:
           self.opPanel.Dismiss_cb()
        else:
            if not hasattr(self.opPanel, 'optionsForm'):
                self.opPanel.displayPanel(create=1)
            else:
                self.opPanel.displayPanel(create=0)


    def mouseDown(self, event):
        # remember where the mouse went down
        #self.lastx = event.x
        #self.lasty = event.y
        self.setXcursor(event.x)

    def mouseUp(self, event):
    # call callbacks if not in continuous mode
        if not self.continuous:
            self.callbacks.CallCallbacks(self.opPanel.valInput.get())
        if self.showLabel == 2:
            # no widget labels on mouse release
            self.canvas.itemconfigure(self.labelId2, text='')
            self.canvas.itemconfigure(self.labelId, text='')


    def mouseMove(self, event):
    # move the cursor
        self.setXcursor(event.x)
        #self.lastx = event.x


    def printLabel(self):
        if self.canvas is None:
            return
        self.canvas.itemconfigure(self.labelId2,
                                  text=self.labelFormat%self.value)#newVal)
        self.canvas.itemconfigure(self.labelId,
                                  text=self.labelFormat%self.value)#newVal)
        
    
    def drawCursor(self, x):
        if self.canvas:
            self.canvas.coords(self.cursorTk, x, 0, x, self.height)
        
        
    def get(self):
        return self.type(self.value)


    def setXcursor(self, x, update=1, force=0):       
        """ x is a cursor position in pixel between 1 and self.width
"""
        if x < 1: 
            x = 1
        if x > self.width: 
            x = self.width
            
        self.drawCursor(x)
        
        # the mouse return position from 1 to self.width (x=0 is not drawn)
        # we need cursor position from 0 (so last x is self.width-1)
        x = x - 1 
        
        if self.range is not None:
            self.value = self.min + x * self.range / float(self.widthMinusOne)
            
        newVal = self.get()
        
        if self.continuous or force:
            if update and self.oldValue != newVal or force:
                self.oldValue = newVal
                self.callbacks.CallCallbacks(newVal)
            if self.showLabel==2:
                self.printLabel()
        else:
            if self.showLabel==2:
                self.printLabel()
                
        if self.showLabel==1:
            self.printLabel()
        if self.opPanel:
            self.opPanel.valInput.set(self.labelFormat%newVal)


    def set(self, x, update=1, force=0):       
        """ x is a value between self.min and self.max
"""
        #print "set ContourSpectrum"
        if self.range is not None:
            xcursor = (x - self.min) * float(self.widthMinusOne) / self.range  
            xcursor = xcursor + 1
            self.drawCursor(xcursor)
            self.setXcursor(xcursor,update,force)


    def createCanvas(self, master):
        self.frame = Tkinter.Frame(self, borderwidth=3, relief='sunken')
        self.canvas = Tkinter.Canvas(self.frame, width=self.width, height=self.height)
        self.xm = 25
        self.ym = 25
        self.labelId2 = self.canvas.create_text(self.xm+2, self.ym+2,
                                           fill='black',
                                           justify='center', text='',
                                           font = self.labelFont)
        self.labelId = self.canvas.create_text(self.xm, self.ym,
                                          fill=self.labelColor,
                                          justify='center', text='',
                                          font = self.labelFont)

        # pack em up
        self.canvas.pack(side=Tkinter.TOP)
        self.frame.pack(expand=1, fill='x')
        self.toggleWidgetLabel(self.showLabel)


    def toggleWidgetLabel(self, val):
        if val == 0:
            # no widget labels
            self.showLabel=0
            self.canvas.itemconfigure(self.labelId2,
                                      text='')
            self.canvas.itemconfigure(self.labelId,
                                      text='')

        if val == 1:
            # show always widget labels
            self.showLabel=1
            self.printLabel()

        if val == 2:
            # show widget labels only when mouse moves
            self.showLabel=2
            self.canvas.itemconfigure(self.labelId2,
                                      text='')
            self.canvas.itemconfigure(self.labelId,
                                      text='')


    def setValue(self, val):
        #print "setValue"
        assert type(val) in [types.IntType, types.FloatType],\
               "Illegal type for value: expected %s or %s, got %s"%(
                   type(1), type(1.0), type(val) )
        
        # setValue does NOT call a callback!
        if self.min is not None and val < self.min: 
            val = self.min
        if self.max is not None and val > self.max: 
            val = self.max
        self.value = self.type(val)
        self.offsetValue=self.value
        self.oldValue = self.value

        #print "setValue ContourSpectrum"
        if self.range is not None:
            xcursor = (val - self.min) * float(self.widthMinusOne) / self.range  
            xcursor = xcursor + 1
            self.drawCursor(xcursor)

        if self.showLabel == 1:
            self.printLabel()
        if self.opPanel:
            self.opPanel.valInput.set(self.labelFormat%self.value)
        

    def setLabel(self, labCfg):
        self.labCfg = labCfg

        text = labCfg.get('text', None)
        if text is None or text=='':
            return

        d={}
        for k, w in self.labCfg.items():
            if k == 'side': continue
            else: d[k] = w
        if not 'side' in self.labCfg.keys():
            self.labCfg['side'] = 'left'

        if not self.lab:
            self.lab = Tkinter.Label(self, d)
            self.lab.pack(side=self.labCfg['side'])
            self.lab.bind("<Button-3>", self.toggleOptPanel)
        else:
            self.lab.configure(text)

            
 #####################################################################
 # the 'configure' methods:
 #####################################################################

    def configure(self, **kw):
        for key,value in kw.items():
            # the 'set' parameter callbacks
            if key=='labCfg': self.setLabel(value)
            elif key=='type': self.setType(value)
            elif key=='min': self.setMin(value)
            elif key=='max': self.setMax(value)
            elif key=='precision': self.setPrecision(value)
            elif key=='showLabel': self.setShowLabel(value)
            elif key=='continuous': self.setContinuous(value)

            # the 'lock' entries callbacks
            elif key=='lockType': self.lockTypeCB(value)
            elif key=='lockMin': self.lockMinCB(value)
            elif key=='lockBMin': self.lockBMinCB(value)
            elif key=='lockMax': self.lockMaxCB(value)
            elif key=='lockBMax': self.lockBMaxCB(value)
            elif key=='lockPrecision': self.lockPrecisionCB(value)
            elif key=='lockShowLabel': self.lockShowLabelCB(value)
            elif key=='lockValue': self.lockValueCB(value)
            elif key=='lockContinuous': self.lockContinuousCB(value)


    def setType(self, Type):

        assert type(Type) in [types.StringType, types.TypeType],\
               "Illegal type for datatype. Expected %s or %s, got %s"%(
                   type('a'), type(type), type(Type) )
        
        if type(Type) == type(""): # type str
            assert Type in ('int','float'),\
            "Illegal type descriptor. Expected 'int' or 'float', got '%s'"%Type
            self.type = eval(Type)
        else:
            self.type = Type
        
        if self.type == int:
            self.labelFormat = "%d"
            self.int_value = self.value
        else:
            self.labelFormat = "%."+str(self.precision)+"f"

        if hasattr(self.opPanel, 'optionsForm'):
            w = self.opPanel.idf.entryByName['togIntFloat']['widget']
            if self.type == int:
                w.setvalue('int')
            elif self.type == 'float':
                w.setvalue('float')
            
        if self.opPanel:
            self.opPanel.updateDisplay()

        # and update the printed label
        if self.canvas and self.showLabel == 1:
            self.printLabel()


    def setMin(self, min):
        if min is not None:
            assert type(min) in [types.IntType, types.FloatType,
                                 numpy.int, numpy.int8, numpy.int16,
                                 numpy.int32, numpy.int64,
                                 numpy.uint, numpy.uint8, numpy.uint16,
                                 numpy.uint32, numpy.uint64,
                                 numpy.float, numpy.float32, numpy.float64],\
                 "Illegal type for minimum. Expected type %s or %s, got %s"%(
                     type(0), type(0.0), type(min) )
            
            if self.max and min > self.max:
                min = self.max

            self.min = self.type(min)

            if self.showLabel == 1:
                self.printLabel()

            if self.value < self.min:
                self.set(self.min)

            if hasattr(self.opPanel, 'optionsForm'):
                self.opPanel.minInput.set(self.labelFormat%self.min) 

                self.opPanel.toggleMin.set(1)
                self.opPanel.min_entry.configure(state='normal', fg='gray0')
            self.minOld = self.min    

        else:
            self.min = None
            if hasattr(self.opPanel, 'optionsForm'):
                self.opPanel.toggleMin.set(0)
                self.opPanel.min_entry.configure(state='disabled',
                                                 fg='gray40')

        if self.min is not None and self.max is not None:
            self.range = float(self.max - self.min)
        else:
            self.range = None 


    def setMax(self, max):
        if max is not None:
            assert type(max) in [types.IntType, types.FloatType,
                                 numpy.int, numpy.int8, numpy.int16,
                                 numpy.int32, numpy.int64,
                                 numpy.uint, numpy.uint8, numpy.uint16,
                                 numpy.uint32, numpy.uint64,
                                 numpy.float, numpy.float32, numpy.float64],\
                 "Illegal type for maximum. Expected type %s or %s, got %s"%(
                     type(0), type(0.0), type(max) )
             
            if self.min and max < self.min:
                max = self.min

            self.max = self.type(max)

            if self.showLabel == 1:
                self.printLabel()

            if self.value > self.max:
                self.set(self.max)

            if hasattr(self.opPanel, 'optionsForm'):
                self.opPanel.maxInput.set(self.labelFormat%self.max) 

                self.opPanel.toggleMax.set(1)
                self.opPanel.max_entry.configure(state='normal', fg='gray0')
            self.maxOld = self.max
                
        else:
            self.max = None
            if hasattr(self.opPanel, 'optionsForm'):
                self.opPanel.toggleMax.set(0)
                self.opPanel.max_entry.configure(state='disabled', fg='gray40')

        if self.min is not None and self.max is not None:
            self.range = float(self.max - self.min)
        else:
            self.range = None 


    def setPrecision(self, val):
        assert type(val) in [types.IntType, types.FloatType,
                                 numpy.int, numpy.float32],\
               "Illegal type for precision. Expected type %s or %s, got %s"%(
                     type(0), type(0.0), type(val) )
        val = int(val)
        
        if val > 10:
            val = 10
        if val < 1:
            val = 1
        self.precision = val

        if self.type == float:
            self.labelFormat = "%."+str(self.precision)+"f"
        else:
            self.labelFormat = "%d"

        if hasattr(self.opPanel, 'optionsForm'):
            w = self.opPanel.idf.entryByName['selPrec']['widget']
            w.setvalue(val)
            
        if self.opPanel:
            self.opPanel.updateDisplay()

        # and update the printed label
        if self.canvas and self.showLabel == 1:
            self.printLabel()


    def setContinuous(self, cont):
        """ cont can be None, 0 or 1 """

        assert cont in [None, 0, 1],\
             "Illegal value for continuous: expected None, 0 or 1, got %s"%cont

        if cont != 1:
            cont = None
        self.continuous = cont
        if hasattr(self.opPanel, 'optionsForm'):
            w = self.opPanel.idf.entryByName['togCont']['widget']
            if cont:
                w.setvalue('on')#i=1
            else:
                w.setvalue('off')#i=0

        if self.opPanel:
            self.opPanel.updateDisplay()


    def setShowLabel(self, val):
        """Show label can be 0, 1 or 2
0: no label
1: label is always shown
2: show label only when value changes"""
        
        assert val in [0,1,2],\
               "Illegal value for showLabel. Expected 0, 1 or 2, got %s"%val
        
        if val != 0 and val != 1 and val != 2:
            print "Illegal value. Must be 0, 1 or 2"
            return
        self.showLabel = val
        self.toggleWidgetLabel(val)

        if hasattr(self.opPanel, 'optionsForm'):
            w = self.opPanel.idf.entryByName['togLabel']['widget']
            if self.showLabel == 0:
                label = 'never'
            elif self.showLabel == 1:
                label = 'always'
            elif self.showLabel == 2:
                label = 'move'
            w.setvalue(label)

        if self.opPanel:
            self.opPanel.updateDisplay()


 #####################################################################
 # the 'lock' methods:
 #####################################################################


    def lockTypeCB(self, mode):
        if mode != 0: mode = 1
        self.lockType = mode
        if hasattr(self.opPanel, 'optionsForm'):
            self.opPanel.lockUnlockDisplay()
        

    def lockMinCB(self, mode): #min entry field
        if mode != 0: mode = 1
        self.lockMin = mode
        if hasattr(self.opPanel, 'optionsForm'):
            self.opPanel.lockUnlockDisplay()
        

    def lockBMinCB(self, mode): # min checkbutton
        if mode != 0: mode = 1
        self.lockBMin = mode
        if hasattr(self.opPanel, 'optionsForm'):
            self.opPanel.lockUnlockDisplay()


    def lockMaxCB(self, mode): # max entry field
        if mode != 0: mode = 1
        self.lockMax = mode
        if hasattr(self.opPanel, 'optionsForm'):
            self.opPanel.lockUnlockDisplay()


    def lockBMaxCB(self, mode): # max checkbutton
        if mode != 0: mode = 1
        self.lockBMax = mode
        if hasattr(self.opPanel, 'optionsForm'):
            self.opPanel.lockUnlockDisplay()


    def lockPrecisionCB(self, mode):
        if mode != 0: mode = 1
        self.lockPrecision = mode
        if hasattr(self.opPanel, 'optionsForm'):
            self.opPanel.lockUnlockDisplay()
            

    def lockShowLabelCB(self, mode):
        if mode != 0: mode = 1
        self.lockShowLabel = mode
        if hasattr(self.opPanel, 'optionsForm'):
            self.opPanel.lockUnlockDisplay()


    def lockValueCB(self, mode):
        if mode != 0: mode = 1
        self.lockValue = mode
        if hasattr(self.opPanel, 'optionsForm'):
            self.opPanel.lockUnlockDisplay()


    def lockContinuousCB(self, mode):
        if mode != 0: mode = 1
        self.lockContinuous = mode
        if hasattr(self.opPanel, 'optionsForm'):
            self.opPanel.lockUnlockDisplay()


    def setSignatures(self, signatures):
        self.signatures = signatures
        self.sigData = []

        # get the values
        self.maxFun = []
        self.minFun = []
        for s in self.signatures:
            x = numpy.zeros( (s.nval,), 'f')
            s.getFx(x)

            self.minix = mini = min(x)
            if (isinstance(mini, numpy.ndarray)) and (mini.shape == () ):
                mini = mini[0]

            maxi = max(x)
            if (isinstance(maxi, numpy.ndarray)) and (maxi.shape == () ):
                maxi = maxi[0]

            self.rangex = range = maxi-mini
            if range != 0:
                x = (((x-mini)/range)*self.widthMinusOne).astype('i')

            y = numpy.zeros( (s.nval,), 'f')
            s.getFy(y)
            self.sigData.append( (x,y) )
            self.maxFun.append(max(y))
            self.minFun.append(min(y))
            
        self.setMin(mini)
        self.setMax(maxi)

        # iso value with hightest value in first function
        if len(self.sigData):
            ind = list(self.sigData[0][1]).index(max(self.sigData[0][1]))
            self.setXcursor(ind)
        else:
            self.setXcursor(0.0)
        
        self.drawSignatures()


    def drawSignatures(self):
        # compute normalization factors
        self.yratios = []
        maxi = max(self.maxFun)
        for i in range(4):
            h = self.height-1.
            if maxi != 0 and self.maxFun[i] != 0:
                self.yratios.append( (h/maxi)*(maxi/self.maxFun[i]) )
            else:
                self.yratios.append( 0 )

        for l in self.tkLines:
            self.canvas.delete(l)
            
        for i, f in enumerate(self.sigData):
            coords = []
            for x,y in zip (f[0], f[1]):
                coords.append(x)
                coords.append(self.height-y*self.yratios[i])
            self.tkLines.append( apply( self.canvas.create_line, coords,
                                        {'fill':self.colors[i]}) )



