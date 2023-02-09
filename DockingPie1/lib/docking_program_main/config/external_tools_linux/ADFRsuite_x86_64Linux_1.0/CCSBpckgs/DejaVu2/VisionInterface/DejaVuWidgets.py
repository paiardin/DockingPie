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

## Automatically adapted for numpy.oldnumeric Jul 23, 2007 by 

########################################################################
#
# Date:  2003 Authors: Daniel Stoffler, Michel Sanner
#
#    stoffler@scripps.edu
#    sanner@scripps.edu
#
#       The Scripps Research Institute (TSRI)
#       Molecular Graphics Lab
#       La Jolla, CA 92037, USA
#
# Copyright: Daniel Stoffler, Michel Sanner and TSRI
#
# revision: Guillaume Vareille
#
#########################################################################
#
# $Header: /mnt/raid/services/cvs/DejaVu2/VisionInterface/DejaVuWidgets.py,v 1.1.1.1.4.1 2017/07/13 22:20:08 annao Exp $
#
# $Id: DejaVuWidgets.py,v 1.1.1.1.4.1 2017/07/13 22:20:08 annao Exp $
#

import numpy.oldnumeric as Numeric
import Tkinter, Pmw

from NetworkEditor.widgets import TkPortWidget, PortWidget

from DejaVu2.colorMap import ColorMap
from DejaVu2.ColormapGui import ColorMapGUI
from DejaVu2.colorTool import RGBRamp
from DejaVu2.ColorWheel import ColorWheel      
from mglutil.gui.BasicWidgets.Tk.colorWidgets import ColorEditor
from mglutil.util.callback import CallBackFunction


class NEDejaVu2GeomOptions(PortWidget):
    """widget allowing to configure a DejaVu2Geometry
Handle all things that can be set to a finite list of values in a Geometry
"""

    configOpts = PortWidget.configOpts.copy()

    ownConfigOpts = {}
    ownConfigOpts['initialValue'] = {
        'defaultValue':{}, 'type':'dict',
        }

    configOpts.update(ownConfigOpts)

    def __init__(self, port, **kw):
        
        # call base class constructor
        apply( PortWidget.__init__, (self, port), kw)
        
        self.booleanProps = [
                                'disableTexture',
                                'faceNormals',
                                'inheritBackPolyMode',
                                'inheritCulling',
                                'inheritFrontPolyMode',
                                'inheritLighting'
                                'inheritStippleLines',
                                'inheritLineWidth',
                                'inheritMaterial',
                                'inheritPointWidth',
                                'inheritStipplePolygons',
                                'inheritShading',
                                'inheritSharpColorBoundaries',
                                'inheritXform',
                                'invertNormals',
                                'lighting',
                                'protected',
                                'scissor',
                                'sharpColorBoundaries',
                                'vertexNormals', 
                                'visible',
                            ]

        from DejaVu2 import viewerConst
        self.choiceProps = {
            'frontPolyMode':viewerConst.Front_POLYGON_MODES_keys,
            'backPolyMode':viewerConst.Back_POLYGON_MODES_keys,
            'shading':viewerConst.SHADINGS.keys(),
            'culling':viewerConst.CULLINGS.keys(),
        }

        self.frame = Tkinter.Frame(self.widgetFrame, borderwidth=3,
                                   relief = 'ridge')

        self.propWidgets = {} # will hold handle to widgets created
        self.optionsDict = {} # widget's value
        
        import Pmw
        items = self.booleanProps + self.choiceProps.keys()
        w = Pmw.Group(self.frame, tag_text='Add a Property Widget')
        self.chooser = Pmw.ComboBox(
            w.interior(), label_text='', labelpos='w',
            entryfield_entry_width=20, scrolledlist_items=items,
            selectioncommand=self.addProp)
        self.chooser.pack(padx=2, pady=2, expand='yes', fill='both')
 	w.pack(fill = 'x', expand = 1, side='top')

        w = Pmw.Group(self.frame, tag_text='Property Widgets')
        self.propWidgetMaster = w.interior()
        w.pack(fill='both', expand = 1, side='bottom')

        # configure without rebuilding to avoid enless loop
        apply( self.configure, (False,), kw)
        self.frame.pack(expand='yes', fill='both')
        
        if self.initialValue is not None:
            self.set(self.initialValue, run=0)

        self._setModified(False) # will be set to True by configure method


    def addProp(self, prop):
        if self.propWidgets.has_key(prop):
            return

        l = Tkinter.Label(self.propWidgetMaster, text=prop)
        l.grid(padx=2, pady=2, row=len(self.propWidgets)+1, column=0,
               sticky='e')

        if prop in self.booleanProps:
            var = Tkinter.IntVar()
            var.set(0)
            cb = CallBackFunction( self.setBoolean, (prop, var))
            widget = Tkinter.Checkbutton(self.propWidgetMaster,
                                         variable=var, command=cb)
            self.propWidgets[prop] = (widget, var)
            self.setBoolean( (prop, var) )
        else:
            items = self.choiceProps[prop]
            var = None
            cb = CallBackFunction( self.setChoice, (prop,))
            widget = Pmw.ComboBox(
                self.propWidgetMaster, entryfield_entry_width=15,
                scrolledlist_items=items, selectioncommand=cb)
            self.propWidgets[prop] = (widget, var)
            self.setChoice( (prop,), items[0] )

        widget.grid(row=len(self.propWidgets), column=1, sticky='w')


    def setBoolean(self, args):
        prop, var = args
        self.optionsDict[prop] = var.get()
        if self.port.node.paramPanel.immediateTk.get():
            self.scheduleNode()
            

    def setChoice(self, prop, value):
        self.optionsDict[prop[0]] = value
        self.propWidgets[prop[0]][0].selectitem(value)
        if self.port.node.paramPanel.immediateTk.get():
            self.scheduleNode()

        
    def set(self, valueDict, run=1):
        self._setModified(True)
        for k,v in valueDict.items():
            self.addProp(k)
            if k in self.booleanProps:
                self.propWidgets[k][1].set(v)
                self.setBoolean( (k, self.propWidgets[k][1]) )
            else:
                self.setChoice((k,), v)
                
        self._newdata = True
        if run:
            self.scheduleNode()


    def get(self):
        return self.optionsDict

        
    def configure(self, rebuild=True, **kw):
        action, rebuildDescr = apply( PortWidget.configure, (self, 0), kw)
            
        #  this methods just creates a resize action if width changes
        if self.widget is not None:
            
            if 'width' in kw:
                action = 'resize'

        if action=='rebuild' and rebuild:
            action, rebuildDescr = self.rebuild(rebuildDescr)


        if action=='resize' and rebuild:
            self.port.node.autoResize()

        return action, rebuildDescr


class NEColorMap(PortWidget):

    # description of parameters that can only be used with the widget's
    # constructor
    configOpts = PortWidget.configOpts.copy()
    ownConfigOpts = {
        'mini':{'type':'float', 'defaultValue':None},
        'maxi':{'type':'float', 'defaultValue':None},
        #'ramp':{'defaultValue':None},
        'filename':{'type':'string', 'defaultValue':''},
        'viewer':{'defaultValue':None},
        }
    configOpts.update( ownConfigOpts )


    def configure(self, rebuild=True, **kw):
        action, rebuildDescr = apply( PortWidget.configure, (self, 0), kw)

        # handle ownConfigOpts entries
        if self.widget is not None:
            
            widgetOpts = {}

            # handle viewer first so that legend exists if min or max is set
            viewer = kw.get('viewer', None)
            if viewer:
                if self.widget.viewer != viewer:
                    self.widget.SetViewer(viewer)
                
            for k, v in kw.items():
                if k=='viewer':
                    continue
                elif k=='mini' or k=='maxi':
                    apply( self.widget.update, (), {k:v} )
                elif k=='filename':
                    self.widget.read(v)


    def __init__(self, port, **kw):
        
        # create all attributes that will not be created by configure because
        # they do not appear on kw
        for key in self.ownConfigOpts.keys():
            v = kw.get(key, None)
            if v is None: # self.configure will not do anyting for this key
                setattr(self, key, self.ownConfigOpts[key]['defaultValue'])

        # get all arguments handled this widget and not by PortWidget
        widgetcfg = {}
        for k in self.ownConfigOpts.keys():
            if k in kw:
                widgetcfg[k] = kw.pop(k)

        # we build 2 widgets here, let's do the first:
        cmkw = {} # dict for ColorMap keywords
        ramp = widgetcfg.pop('ramp', None)
        #if ramp is None:
            #ramp = RGBRamp(size=16)
            #cmkw['ramp'] = ramp
        cmkw['mini'] = widgetcfg.pop('mini', None)
        cmkw['maxi'] = widgetcfg.pop('maxi', None)

        # create an instance of the ColorMap which is required to instanciate
        # the ColorMapGUI (i.e. the widget)
        cM = apply( ColorMap, ('cmap',), cmkw)
        
        # call base class constructor
        apply( PortWidget.__init__, ( self, port), kw)

        # create the widget
        widgetcfg['master'] = self.widgetFrame
        widgetcfg['modifyMinMax'] = True
        
        #try: #tries to assign viewer when called from Pmv
        #    widgetcfg['viewer'] = self.vEditor.vf.GUI.VIEWER
        #except:
        #    pass

        self.widget = apply( ColorMapGUI, (cM,), widgetcfg)
        
        # we call cmCallback because cm passed the color map as an argument
        # and node.schedule takes no argument
        self.widget.addCallback(self.cmCallback)

        # configure without rebuilding to avoid enless loop
        #apply( self.configure, (False,), widgetcfg)

        if self.initialValue:
            self.set(self.initialValue, run=0)
            
        self.modified = False # will be set to True by configure method

        # and destroy the Apply and Dismiss buttons of the ColorMapGUI
        self.widget.dismiss.forget()
        self.widget.apply.forget()
        self.widget.frame2.forget()

        # overwrite the paramPanel Appply button method to call
        # the ColorMapGUI apply_cb() method
        if hasattr(self.port.node, 'paramPanel'):
            self.port.node.paramPanel.applyFun = self.apply_cb


##      def __init__(self, port=None, master=None, visibleInNodeByDefault=0,
##                   visibleInNode=0, value=None, callback=None, **kw):

##          PortWidget.__init__(self, port, master, visibleInNodeByDefault,
##                              visibleInNode, callback)

##          cmkw = {} # dict for ColorMap keywords
##          cmkw['ramp'] = RGBRamp() #default ramp, can be overwritten below
        
##          for k in kw.keys():
##              if k in ['ramp', 'geoms', 'filename', 'mini', 'maxi']:
##                  cmkw[k] = kw[k]
##                  del kw[k]
##          # create an instance of the ColorMap which is required to instanciate
##          # the ColorMapGUI
##          cM = apply( ColorMap, ('cmap',), cmkw)

##          # instanciate the ColorMapGUI
##          kw['master'] = self.top
##          self.cm = apply( ColorMapGUI, (cM,), kw)
##          # and destroy the Apply and Dismiss buttons of the ColorMapGUI
##          self.cm.dismiss.forget()
##          self.cm.apply.forget()
##          self.cm.frame2.forget()
##          # add callback. Note, we do not call scheduleNode because this would
##          # lead to recursion problem
##          #self.cm.addCallback(self.port.node.schedule)

##          # we call cmCallback because cm passed the color map as an argument
##          # and node.schedule takes no argument
##          self.cm.addCallback(self.cmCallback)
        
##          if value is not None:
##              self.set(value, run=0)

##          # overwrite the paramPanel Appply button method to call
##          # the ColorMapGUI apply_cb() method
##          if hasattr(self.port.node, 'paramPanel'):
##              self.port.node.paramPanel.applyFun = self.apply_cb


    def onUnbind(self):
        legend = self.widget.legend
        if legend is not None:
            viewer = legend.viewer
            if viewer is not None:
                legend.protected = False
                viewer.RemoveObject(legend)
            legend = None


    def onDelete(self):
        self.onUnbind()


    def cmCallback(self, cmap):
        self.port.node.schedule()
        

    def apply_cb(self, event=None):
        if self.widget is not None:
            self._setModified(True)
            # overwrite the paramPanel Appply button method to call
            # the ColorMapGUI apply_cb() method
            self.widget.apply_cb(event)
            self._newdata = 1
            self.port.node.schedule()
        
            
    def set(self, value, run=1):
        #import pdb;pdb.set_trace()
        self._setModified(True)
        
        if isinstance(value, dict):
            apply(ColorMapGUI.configure,(self.widget,),value)
        else:
            self.widget.set(value)
            
        self._newdata = 1
        if run:
            self.port.node.schedule()


    def get(self):
        return self.widget


    def getDataForSaving(self):
        # this method is called when a network is saved and the widget
        # value needs to be saved
        cfg = PortWidget.getDescr(self)
        cfg['name'] = self.widget.name 
        cfg['ramp'] = self.widget.ramp 
        cfg['mini'] = self.widget.mini
        cfg['maxi'] = self.widget.maxi                                             
        return cfg


    def getDescr(self):
        cfg = PortWidget.getDescr(self)
        
        # the whole colormap is the widget value and
        # is returned by self.get() !
        
        #cfg['name'] = self.widget.name 
        #cfg['ramp'] = self.widget.ramp 
        #cfg['mini'] = self.widget.mini
        #cfg['maxi'] = self.widget.maxi                                             
        return cfg



class NEColorWheel(PortWidget):

    # description of parameters that can only be used with the widget's
    # constructor
    configOpts = PortWidget.configOpts.copy()
    ownConfigOpts = {
        'width':{'min':20, 'max':500, 'defaultValue':75, 'type':'int'},
        'height':{'min':20, 'max':500, 'defaultValue':75, 'type':'int'},
        'circles':{'min':2, 'max':20, 'defaultValue':5, 'type':'int'},
        'stripes':{'min':2, 'max':40, 'defaultValue':20, 'type':'int'},
        }
    configOpts.update( ownConfigOpts )


    def __init__(self, port, **kw):
        
        # create all attributes that will not be created by configure because
        # they do not appear on kw
        for key in self.ownConfigOpts.keys():
            v = kw.get(key, None)
            if v is None: # self.configure will not do anyting for this key
                setattr(self, key, self.ownConfigOpts[key]['defaultValue'])

        # get all arguments handled by this widget and not by PortWidget
        widgetcfg = {}
        for k in self.ownConfigOpts.keys():
            if k in kw:
                widgetcfg[k] = kw.pop(k)

        # call base class constructor
        apply( PortWidget.__init__, ( self, port), kw)

        # create the widget
        self.widget = apply( ColorWheel, (self.widgetFrame,), widgetcfg)
        self.widget.AddCallback(self.cwCallback)

        # configure without rebuilding to avoid enless loop
        #apply( self.configure, (False,), widgetcfg)

        if self.initialValue:
            self.set(self.initialValue, run=0)
            
        self.modified = False # will be set to True by configure method

##      def __init__(self, port=None, master=None, visibleInNodeByDefault=0,
##                   visibleInNode=0, value=None, callback=None, **kw):

##          PortWidget.__init__(self, port, master, visibleInNodeByDefault,
##                              visibleInNode, callback)
##          self.cw = apply( ColorWheel, (self.top,), kw)
##          if value is not None:
##              self.set(value, run=0)
##          self.cw.AddCallback(self.scheduleNode)

    def cwCallback(self, color):
        self._newdata = 1
        self.scheduleNode()

    def set(self, value, run=1):
        self.widget.Set(value, mode='RGB')
        self._newdata = 1
        if run:
            self.scheduleNode()
            #self.port.node.schedule()


    def get(self, mode='RGB'):
        value = self.widget.Get(mode)
        value = list( Numeric.array(value) )
        return value

        
    def getDescr(self):
        cfg = PortWidget.getDescr(self)
        cfg['width'] = self.widget.width
        cfg['height'] = self.widget.height
        cfg['circles'] = self.widget.circles
        cfg['stripes'] = self.widget.stripes
        cfg['wheelPad'] = self.widget.wheelPad
        cfg['immediate'] = self.widget.immediate
        cfg['wysiwyg'] = self.widget.wysiwyg
        return cfg

       
##      def configure(self, **kw):
##          if len(kw)==0:
##              cfg = PortWidget.configure(self)
##              cfg['immediate'] = self.cw.immediate
##              cfg['wysiwyg'] = self.cw.wysiwyg
##              return cfg
##          else:
##              if kw.has_key('immediate'):
##                  self.cw.setImmediate(kw['immediate'])
##              elif kw.has_key('wysiwyg'):
##                  self.cw.setWysiwyg(kw['wysiwyg'])


class NEColorEditor(PortWidget):

    # description of parameters that can only be used with the widget's
    # constructor
    configOpts = PortWidget.configOpts.copy()
    ownConfigOpts = {
        'mode':{'defaultValue':'RGB', 'type':'string',
                'validValues':['RGB', 'HSV', 'HEX']},
        'immediate':{'defaultValue':True, 'type':'boolean'},
        }
    
    configOpts.update( ownConfigOpts )


    def __init__(self, port, **kw):
        
        # create all attributes that will not be created by configure because
        # they do not appear on kw
        for key in self.ownConfigOpts.keys():
            v = kw.get(key, None)
            if v is None: # self.configure will not do anyting for this key
                setattr(self, key, self.ownConfigOpts[key]['defaultValue'])

        # get all arguments handled by this widget and not by PortWidget
        widgetcfg = {}
        for k in self.ownConfigOpts.keys():
            if k in kw:
                widgetcfg[k] = kw.pop(k)

        # call base class constructor
        apply( PortWidget.__init__, ( self, port), kw)

        # create the widget
        self.widget = apply( ColorEditor, (self.widgetFrame,), widgetcfg)
        self.widget.cbManager.AddCallback(self.cwCallback)
        self.widget.pack()
        # configure without rebuilding to avoid enless loop
        #apply( self.configure, (False,), widgetcfg)

        if self.initialValue:
            self.set(self.initialValue, run=0)
            
        self.modified = False # will be set to True by configure method

##      def __init__(self, port=None, master=None, visibleInNodeByDefault=0,
##                   visibleInNode=0, value=None, callback=None, **kw):

##          PortWidget.__init__(self, port, master, visibleInNodeByDefault,
##                              visibleInNode, callback)
##          self.cw = apply( ColorWheel, (self.top,), kw)
##          if value is not None:
##              self.set(value, run=0)
##          self.cw.AddCallback(self.scheduleNode)

    def cwCallback(self, color):
        self._newdata = 1
        self.scheduleNode()

    def set(self, value, run=1):
        if value is not None:
            self._setModified(True)
            self.widget.set(value, mode='RGB')
            self._newdata = 1
            if run:
                self.port.node.schedule()


    def get(self, mode='RGB'):
        value = self.widget.get(mode='RGB')
        value = list( Numeric.array(value) )
        return value

        
    def getDescr(self):
        cfg = PortWidget.getDescr(self)
        cfg['mode'] = self.widget.mode
        cfg['immediate'] = self.widget.immediate
        return cfg

