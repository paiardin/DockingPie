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

import os
import Tkinter, Pmw

#from DejaVu2.GeomChooser import GeomChooser
from mglutil.util.callback import CallbackFunction

from DejaVu2.scenarioInterface.animations import FlyInObjectMAA, FlyOutObjectMAA,\
     FadeInObjectMAA, FadeOutObjectMAA, VisibleObjectMAA, ColorObjectMAA,\
     RotationMAAOptionsGUI, RotationMAA, RockMAA, SnapshotMAAGroup
# from DejaVu2.scenarioInterface.animations import FocusMAA

from Scenario2.gui.Tk.clipboard import ClipboardGUI
from Scenario2.sequenceAnimator import SequenceAnimator
from Scenario2.gui.Tk.sequenceAnimator import SequenceAnimatorGUI
from Scenario2 import _clipboard, _MAATargets

from DejaVu2.states import setRendering, getRendering, getOrientation
from DejaVu2.scenarioInterface.animations import OrientationMAA, RenderingTransitionMAA

from mglutil.gui.BasicWidgets.Tk.thumbwheel import ThumbWheel

class MAAEditor:
    """
    Editor for an MAA

    the base class provides mainly the buttons at the bottom of the form

    the .edit(maa) method will save the maa in self.maa, show the form and
    configure it with the paramters gotten from the maa using maa.getValues().

    when OK or Preview are clicked the execute method has to configure
    self.maa and run the maa for Preview
    
    subclass will implement:

    populateForm() to add gui elements to self.dialog.interior()
    getValues() to return a dictionary of parameter names and values
    setValues() to take a dictionary of parameter names and values show them
                in the editor
    execute(name) to configure self.maa and run it for Preview
    """
    def __init__(self, master=None, title='Editor',
                 buttons=['OK', 'Preview', 'Cancel'],
                 defaultButton='OK'):
        """
        base class bur creating MAA editors

        editorObject <- MAAEditor( maa, master=None, title='Editor',
                 buttons=['OK','Preview', 'Cancel'],
                 defaultButton='OK')

        - maa is the MAA instance for which we want to see and modify params

        - master can be used to specify where to place the GUI.
        If master is None:
            if showOptionForm was called as callback from a Tk event,
            the dialog just below the button that was clicked
            else the dialogue will appear at position 100, 100

        - title is the title of the Editor window
        
        - buttons is a list of string specifying which buttons to create
        at the bottom of the dialog

        - defaultButton is a string that appears in buttons and will be
        set as the default button for the form
        """

        self.master = master
        self.maa = None # will be set in self.edit(maa)

        self.exitStatus = None # will be saet to the name of the button

        # save list of desired buttons
        assert len(buttons)>0
        self.buttons = buttons
        assert defaultButton in buttons
        self.defaultButtons = defaultButton

        # save the title
        assert isinstance(title, str)
        self.title = title
        
        # create the dialog
        self.dialog = self.createGUI()
        self.populateForm()


    def createGUI(self):
        """
        Create the form.
        This base class form is a Pmw dialog containing a set of Radiobuttons
        and a counter to select speed of animation(number of frames)
        """
        self.balloon = Pmw.Balloon(self.master)

        # create the Dialog
        dialog = Pmw.Dialog(
            self.master, buttons=self.buttons, defaultbutton='OK',
            title=self.title, command=self.execute)
        dialog.withdraw()
        # create a frame to hold group to force setting orient and rendering
        bframe = Tkinter.Frame(dialog.interior())

        ##
        ## create group to define if action sets orients
        self.actionGrp = grp = Pmw.Group(bframe, tag_text='Action sets')

        frame = grp.interior()
        
        self.forcew = w = Pmw.RadioSelect(
            frame, selectmode='multiple', orient='horizontal',
            buttontype='checkbutton')

        for text in ['Orientation', 'Rendering']:
            w.add(text)
    
        w.pack(side='left', anchor='w', fill='x', expand=1, padx=8, pady=8)

        grp.pack(side = 'left', fill='x', expand=1)

        self.balloon.bind(grp, "When an action is created the current orientation and rendering are saved.\nCheck these buttons to have the action set the orientation\nand/or Rendering on its first frame during playback")

        ##
        ## create a group of buttons to overwrite 
        self.recordGrp = grp = Pmw.Group(bframe, tag_text='Record:')

        frame = grp.interior()

        b1 = Tkinter.Button(frame, text='Orientation',
                            command=self.recordOrient)
        b1.pack(side='left', anchor='w', fill='x', expand=1, padx=8, pady=8)
        b2 = Tkinter.Button(frame, text='Rendering',
                            command=self.recordRendering)
        b2.pack(side='left', anchor='w', fill='x', expand=1, padx=8, pady=8)

        grp.pack(side = 'bottom', fill='x', expand=1)

        self.balloon.bind(grp, "Record current orientation/rendering for this action")
        bframe.pack( fill='x', expand=1)

        return dialog


    def recordOrient(self, event=None):
        self.maa.recordOrient()


    def recordRendering(self, event=None):
        self.maa.recordRendering()

        
    def edit(self, maa, event=None):
        self.maa = maa

        # configure the editor with current values
        self.setValues( **maa.getValues() )
        
        # activate and place the dialog just below the button that was clicked:
        if event:
            x = event.x_root
            y = event.y_root
            if y + self.master.winfo_height() > self.master.winfo_screenheight():
                y = self.master.winfo_screenheight() - self.master.winfo_height()
              #dialog.activate(globalMode = 'grab', geometry = '+%d+%d' % (x, y )) # does not work
            self.dialog.activate(geometry = '+%d+%d' % (x, y ))
        else:
            self.dialog.activate(geometry = '+%d+%d' % (100, 100 ))

        if self.exitStatus=='OK':
            return self.getValues()

        elif self.exitStatus in ('None', 'Cancel'):
            return None
       

    def _hide(self):
        self.dialog.deactivate()


    def execute(self, name):
        # called when buttons are clicked ot windows is closed
        self.exitStatus = name
        if name in ('None', 'Cancel', 'OK'):
            self._hide()


    ##
    ## the following methods should be subclassed
    ##
    def populateForm(self):
        """
        subclassses will place GUI items for various parameters in this method
        """
        return

    def getValues(self):
        """
        return a dictionary of parameterName:parameterValues
        """
        return {}

    def setValues(self, **kw):
        """
        take a dictionary of parameterName:parameterValues set the editor
        to these values
        """
        return


class MAAEditorWithSpeed(MAAEditor):
    """
    Sub class adding basic animation force orient/rendering and speed GUI
    """
    def __init__(self, master=None, title='Orientation Editor',
                 buttons=['OK', 'Preview', 'Cancel'],
                 defaultButton='OK', speedDict =None):
        """
        Provides a dialog form for setting Orientation actions parameters

           speedDict: dictionary of speed anmes and nbframes (constr only)
           
        The following parameters are handled by getValues and setValues:
           nbFrames: int
        """

        if speedDict is None:
            speedDict = {'slow': 50, 'medium': 30, 'fast': 10}
        self.speedDict = speedDict
        
        self.custom_speed_flag = False
        self.speed = 'medium'

        MAAEditor.__init__(self, master=master, title=title,
                           buttons=buttons, defaultButton=defaultButton)
        


    def populateForm(self):
        """
        added radiobuttons for speed dict and custom speed
        """

        # create group for speed
        parent = self.dialog.interior()

        grp = Pmw.Group(parent, tag_text='speed (in frames)')
        frame = grp.interior()
        self.speedw = sp = Pmw.RadioSelect(
            frame, selectmode='single', orient='horizontal',
            buttontype='radiobutton', command=self.setSpeed_cb)

        for text in ['slow', 'medium', 'fast', 'custom']:
            sp.add(text)
            
        sp.pack(side='left', anchor='w', fill='x', expand=1, padx=8, pady=8)

        for key, val in self.speedDict.items():
            self.balloon.bind(sp.button(key), "%d frames"%val)

        # add counter for overall custom speed
        self.speedcr = cr = Pmw.Counter(
            frame, entry_width=6, entryfield_value=100,
            entryfield_validate = {'validator':'integer', 'min':1 })
        
        self.balloon.bind(cr, "Enter overall number of frames")
        cr.pack(side = 'left', anchor = 'w', fill = 'x', expand = 1)
        grp.pack(side='top', fill='x', expand=1 )


    def setSpeed_cb(self, val=None):
        if val == "custom":
           self.custom_speed_flag = True
        else:
            self.custom_speed_flag = False
        self.setCounterState()


    def setCounterState(self):
        cr = self.speedcr
        entry = cr._counterEntry._entryFieldEntry
        up = cr._upArrowBtn
        down = cr._downArrowBtn
        if self.custom_speed_flag:
            # activate counter
            entry.configure(state='normal')
            down.bind('<Button-1>', cr._countDown)
            up.bind('<Button-1>', cr._countUp)
        else:
            # deactivate:
            entry.configure(state='disabled')
            down.unbind('<Button-1>')
            up.unbind('<Button-1>')


    def getValues(self):

        kw = {'forceOrient':False, 'forceRendering':False}
        sp = self.speedw.getvalue()
        if sp == "None":
            return kw
        if sp == "custom":
            nbFrames = int(self.speedcr.get())
        else:
            nbFrames = self.speedDict[sp]
        
        kw['nbFrames'] = nbFrames

        for val in self.forcew.getvalue():
            if val=='Orientation':
                kw['forceOrient'] = True
            elif val=='Rendering':
                kw['forceRendering'] = True

        return kw
    

    def setValues(self, **kw):
        """
        configure the option form with the values provided in **kw
        keys can be: 'nbFrames', 'easeInOut', 'direction'
        """
        forceList = []

        # if kw has kfpos use last frame in list as nbframes
        kfpos = kw.pop('kfpos', None)
        if kfpos:
            kw['nbFrames'] = kfpos[-1]-kfpos[0]

        for k,v in kw.items():
            if k == 'nbFrames':
                assert int(v)
                assert v>0
                # set speed radio button
                found = False
                for name, speed in self.speedDict.items():
                    if speed==v:
                        self.speedw.invoke(name)
                        found = True
                        break
                if not found:
                    self.speedw.invoke('custom')
                    self.speedcr.setvalue(v)
                    
            elif k == 'forceOrient' and v:
                forceList.append('Orientation')

            elif k == 'forceRendering' and v:
                forceList.append('Rendering')

            #else:
            #    print 'WARNING: unknown key:value',k,v

        self.forcew.setvalue(forceList)


    def execute(self, name):
        # called when buttons are clicked ot windows is closed
        self.exitStatus = name
        if name in ('None', 'Cancel', 'OK'):
            self._hide()
        values = self.getValues()
        self.maa.forceOrient = values['forceOrient']
        self.maa.forceRendering = values['forceRendering']
        



class SSp_MAAEditor(MAAEditorWithSpeed):
    """
    Editor providing speed, and sortPoly parameters
    """

    def populateForm(self):
        """
        add radio buttons for direction and easeInOut
        """
        MAAEditorWithSpeed.populateForm(self)

        parent = self.dialog.interior()

        grp = Pmw.Group(parent, tag_text='Zsort Polygons')
        frame = grp.interior()

        ordergrp = Pmw.Group(frame, tag_pyclass = None)
        ordergrp.pack(side = 'left', fill='x', expand=1)
        self.sortOrderw = Pmw.RadioSelect(
            ordergrp.interior(), selectmode='single', orient='horizontal',
            buttontype='radiobutton')
        
        for text in ['+Zsort', '-Zsort']:
            self.sortOrderw.add(text)
        self.sortOrderw.pack(side='left', anchor='w', fill='x', expand=1, padx=8, pady=8)
        self.balloon.bind(self.sortOrderw, "-Z sorts by furthest Z first,\n+Z sorts by closest Z first")

        whengrp = Pmw.Group(frame, tag_pyclass = None)
        whengrp.pack(side = 'left', fill='x', expand=1)
        self.sortPolyw = w = Pmw.RadioSelect(
            whengrp.interior(), selectmode='single', orient='horizontal',
            buttontype='radiobutton')

        for text in ['Never', 'Once', 'Always']:
            w.add(text)
        w.pack(side='left', anchor='w', fill='x', expand=1, padx=8, pady=8)
        self.balloon.bind(w, "Select when to Z-sort polygons for proper trancparency")

        
        grp.pack(side = 'left', fill='x', expand=1)



    def getValues(self):
        """
        return a dictionary of parameterName:parameterValues
        """
        values = MAAEditorWithSpeed.getValues(self)
        values.update({'sortPoly': self.sortPolyw.getvalue(),
                       'sortOrder': self.sortOrderw.getvalue()})
        
        return values


    def setValues(self, **kw):
        """
        take a dictionary of parameterName:parameterValues set the editor
        to these values
        """
        MAAEditorWithSpeed.setValues(self, **kw)
        self.sortPolyw.setvalue(kw['sortPoly'])
        self.sortOrderw.setvalue(kw['sortOrder'])



class OrientationMAAEditor(MAAEditorWithSpeed):

    def __init__(self, master=None, title='Orientation Editor',
                 buttons=['OK', 'Preview', 'Cancel'],
                 defaultButton='OK'):
        """
        Provides a dialog form for setting Orientation actions parameters

        parameters:
           keyframes: [0, a,b,c]
        """
        MAAEditorWithSpeed.__init__( self, master=master, title=title,
                          buttons=buttons, defaultButton=defaultButton)
        

    def populateForm(self):
        """
        add entryfield for 3 integers to specify 3 intervals
        """
        MAAEditorWithSpeed.populateForm(self)

        # add 3 custom intervals entry form
        
        parent = self.dialog.interior()
        grp = Pmw.Group(parent, tag_text='Custom intervals')
        frame = grp.interior()

        self._3intw = Pmw.EntryField(
            frame, labelpos = 'w', value = '10 20 30',
            label_text = 'nb frames for zoom out, rotate, zoom in:',
            validate = self.custom_validate)

        self._3intw.pack(side='top')
        self.balloon.bind(self._3intw, "Enter 3 space-separated integers\ndefining the zoom out, rotate, and zoom in intervals\ne.g. 10 20 30 means zoom out for 10 frames, rotate for 10 and zoom in for 10")

        grp.pack(side = 'top', fill='x', expand=1 )


    def custom_validate(self, text):
        words = text.split()
        if len(words)!=3:
            return -1
        a,b,c = map(int, words)
        if a>0 and b > a and c>b:
            return 1
        else:
            return -1


    def setSpeed_cb(self, val=None):
        # set self._3intw based on the selected speed
        
        val = MAAEditorWithSpeed.getValues(self).get('nbFrames',None)
        if val:
            self._3intw.setentry("%d %d %d"%(val/3, 2*val/3, val))


    def getValues(self):
        """
        return a dictionary of parameterName:parameterValues
        """
        a,b,c = map(int, self._3intw.getvalue().split())
        kw = {'keyframes':[0,a,b,c],
              'forceOrient':False,
              'forceRendering':False}
        for val in self.forcew.getvalue():
            if val=='Orientation':
                kw['forceOrient'] = True
            elif val=='Rendering':
                kw['forceRendering'] = True

        return kw


    def setValues(self, **kw):
        """
        take a dictionary of parameterName:parameterValues set the editor
        to these values
        """
        zero, a, b, c = kw['keyframes']
        found = False
        for name, value in self.speedDict.items():
            if c==value:
                self.speedw.invoke(name)
                found = True
                break

        if not found:
            self.speedw.invoke('custom')
            self.speedcr.setvalue(value)
            
        return


    def execute(self, name):
        # configure the MAA with the values from the editor
        MAAEditorWithSpeed.execute(self, name)
        if name in ('OK', 'Preview'):
            self.maa.setKeyframePositions( self.getValues()['keyframes'] )
            if name == 'Preview':
                self.maa.run()


class SnapshotMAAGroupEditor(SSp_MAAEditor):
    
    def __init__(self, master=None, title='Snapshot Editor',
                 buttons=['OK', 'Preview', 'Cancel'],
                 defaultButton='OK'):
        """
        Provides a dialog form for setting Snapshot actions parameters

        parameters:
           keyframes: [kf1, kf2]
        """
        SSp_MAAEditor.__init__( self, master=master, title=title,
                          buttons=buttons, defaultButton=defaultButton)
        self.actionGrp._tag.configure(text='Interpolate')
        self.balloon.bind(self.actionGrp, "Check these buttons to interpolate orientation/rendering during playback.")
        

    def getValues(self):
        """
        return a dictionary of parameterName:parameterValues
        """

        #print "snapshot getValues"
        kw =  SSp_MAAEditor.getValues(self)
        nbframes = kw.get('nbFrames',None)
        if nbframes is not None:
            kw.pop('nbFrames')
            kw['keyframes'] = [0, nbframes]
        orient = kw.get('forceOrient', None)
        if orient is not None:
            kw['forceOrient'] = False
            kw['interpolateOrient'] = orient
        rend = kw.get('forceRendering', None)
        if rend is not None:
            kw['forceRendering'] = False
            kw['interpolateRendering'] = rend
        return kw


    def setValues(self, **kw):
        """
        take a dictionary of parameterName:parameterValues set the editor
        to these values
        """
        k1, k2 = kw['keyframes']
        found = False
        for name, value in self.speedDict.items():
            if k2==value:
                self.speedw.invoke(name)
                found = True
                break

        if not found:
            self.speedw.invoke('custom')
            self.speedcr.setvalue(k2)

        forceList = []
        if kw.get('interpolateOrient'):
            forceList.append('Orientation')
        if kw.get('interpolateRendering'): 
            forceList.append('Rendering')
        self.forcew.setvalue(forceList)
        self.sortPolyw.setvalue(kw['sortPoly'])
        if kw.has_key('sortOrder'):
            self.sortOrderw.setvalue(kw['sortOrder'])


    def execute(self, name):
        #print "snapshot execute"
        # configure the MAA with the values from the editor
        MAAEditorWithSpeed.execute(self, name)
        self.maa.sortPoly = sortPoly = self.sortPolyw.getvalue()
        self.maa.renderMaa.sortPoly = sortPoly
        self.maa.orientMaa.sortPoly = sortPoly

        self.maa.sortOrder = sortOrder = self.sortOrderw.getvalue()
        self.maa.renderMaa.sortOrder = sortOrder
        self.maa.orientMaa.sortOrder = sortOrder
        values = self.getValues()
        if values.has_key('interpolateRendering'):
            self.maa.interpolateRendering = values['interpolateRendering']
        if values.has_key('interpolateOrient'):
            self.maa.interpolateOrient = values['interpolateOrient']
        if name in ('OK', 'Preview'):
            self.maa.setKeyframePositions( self.getValues()['keyframes'] )
            if name == 'Preview':
                self.maa.run()


class SE_MAAEditor(MAAEditorWithSpeed):
    """
    Editor providing speed and easeInOut parameters
    """

    def __init__(self, master=None, title='Speed EaseInOut Editor',
                 buttons=['OK', 'Preview', 'Cancel'],
                 defaultButton='OK', speedDict =None, ):

        self.easeInEnd = 0.3
        self.easeOutStart = 0.7
        MAAEditorWithSpeed.__init__(self, master=master, title=title,
                                    buttons=buttons, defaultButton=defaultButton,
                                    speedDict=speedDict)
        
        
    def populateForm(self):
        """
        add radio buttons for direction and easeInOut
        """
        MAAEditorWithSpeed.populateForm(self)

        parent = self.dialog.interior()
        grp = Pmw.Group(parent, tag_text='Ease In/Out')
        frame = grp.interior()
        
        self.easew = w = Pmw.RadioSelect(
            frame, selectmode='single', orient='horizontal',
            buttontype='radiobutton', command=self.setEase_cb)

        for text in ['none', 'ease in', 'ease out', 'ease in and out']:
            w.add(text)

        w.pack(side='left', anchor='w', fill='x', expand=1, padx=8, pady=8)
        self.balloon.bind(w, "Ease in will slow down the start of the action,\n ease out will slow down the end of the action")

        self.easeInEndTw =  ThumbWheel(
                frame, width=70, height=14, type=float,oneTurn=1.0,
                value=self.easeInEnd, continuous=False, wheelPad=1,
                min=0.0, max = 1.0, callback=self.setEaseInEnd_cb,
                labCfg = {'text': 'easeInEnd', 'side':'top'}, increment=0.1,
                showLabel=0, lockShowLabel=1)
        self.easeInEndTw.pack(side = 'left', anchor = 'w', fill = 'x', expand = 1)
        self.balloon.bind(self.easeInEndTw, """ 'ease in' occurs from 0 to EaseInEnd (in range [0, 1])""")
        self.easeOutStartTw =  ThumbWheel(
                frame, width=70, height=14, type=float,oneTurn=1.0,
                value=self.easeOutStart, continuous=False, wheelPad=1,
                min=0.0, max = 1.0, callback=self.setEaseOutStart_cb,
                labCfg = {'text': 'easeOutStart', 'side':'top'}, increment=0.1,
                showLabel=0, lockShowLabel=1)
        self.easeOutStartTw.pack(side = 'left', anchor = 'w', fill = 'x', expand = 1)
        self.balloon.bind(self.easeOutStartTw, """ 'ease out' occurs from easeOutStart to 1(in range [0, 1])""")
        grp.pack(side = 'top', fill='x', expand=1)


    def getValues(self):
        """
        return a dictionary of parameterName:parameterValues
        """
        values = MAAEditorWithSpeed.getValues(self)
        nbFrames = values.pop('nbFrames')

        self.easeInEnd = self.easeInEndTw.get()
        self.easeOutStart = self.easeOutStartTw.get()
        easeInOut = self.easew.getvalue()
        if easeInOut == "ease in and out":
            if self.easeOutStart < self.easeInEnd:
                self.easeOutStart = self.easeInEnd
        values.update( {'kfpos': [0,nbFrames],
                        'easeInOut': self.easew.getvalue(),
                        'easeInEnd': self.easeInEnd,
                        'easeOutStart': self.easeOutStart} )
        #print "SE_MAAEditor, getValues:", values
        return values
    

    def setValues(self, **kw):
        """
        take a dictionary of p <arameterName:parameterValues set the editor
        to these values
        """
        MAAEditorWithSpeed.setValues(self, **kw)
        self.easew.setvalue(kw['easeInOut'])
        easeInEnd = kw.get('easeInEnd')
        if easeInEnd is not None:
            self.easeInEnd = easeInEnd
        easeOutStart = kw.get('easeOutStart')
        if easeOutStart is not None:
            self.easeOutStart = easeOutStart
        self.setEase_cb(kw['easeInOut'])
        #print "SE_MAAEditor, setValues:", kw


    def setEase_cb(self, val=None):
        #print "in setEase_cb, val:", val, "easeInOut:", self.easew.getvalue()
        if val == "ease in":
            self.enableThumbWheel(self.easeInEndTw, self.easeInEnd)
            self.disableThumbWheel(self.easeOutStartTw)
        elif val == "ease out":
            self.enableThumbWheel(self.easeOutStartTw, self.easeOutStart)
            self.disableThumbWheel(self.easeInEndTw)
        elif val == "ease in and out":
            if self.easeInEnd > self.easeOutStart:
                self.easeOutStart = self.easeInEnd
            self.enableThumbWheel(self.easeInEndTw, self.easeInEnd)
            self.enableThumbWheel(self.easeOutStartTw, self.easeOutStart)
        else:
            self.disableThumbWheel(self.easeInEndTw)
            self.disableThumbWheel(self.easeOutStartTw)

    def setEaseInEnd_cb(self, val):
        self.easeInEnd = val


    def setEaseOutStart_cb(self, val):
        self.easeOutStart=val
            

    def disableThumbWheel(self, tw):
        """disables a thumbwheel widgets used to specify easeInEnd and easeOutStart"""
        def foo(val):
            pass
        tw.configure(showLabel=0)
        tw.canvas.bind("<ButtonPress-1>", foo)
	tw.canvas.bind("<ButtonRelease-1>", foo)
	tw.canvas.bind("<B1-Motion>", foo)
        tw.canvas.bind("<Button-3>", foo)
        
        
    def enableThumbWheel(self, tw, val=None):
        """enables a thumbwheel widgets used to specify easeInEnd and easeOutStart"""
        tw.canvas.bind("<ButtonPress-1>", tw.mouseDown)
	tw.canvas.bind("<ButtonRelease-1>", tw.mouseUp)
	tw.canvas.bind("<B1-Motion>", tw.mouseMove)
        tw.canvas.bind("<Button-3>", tw.toggleOptPanel)
        tw.configure(showLabel=1)
        if val:
            tw.set(val, update=0)


    def execute(self, name):
        # configure the MAA with the values from the editor
        MAAEditorWithSpeed.execute(self, name)
        if name in ('OK', 'Preview'):
            self.maa.configure( **self.getValues() )
            if name == 'Preview':
                self.maa.run()


class SED_MAAEditor(SE_MAAEditor):
    """
    Editor providing speed, easeInOut and direction parameters
    """
    def __init__(self, master=None, title='Editor',
                 directions=['left', 'right'],
                 buttons=['OK', 'Preview', 'Cancel'],
                 defaultButton='OK', speedDict=None, addEaseInOut=True):
        """
        Provides a dialog form for setting Fly actions parameters

        parameters:
           kfpos: [a,b] 2 integer key frame positions
           direction: can be 'left', 'right', 'top', or 'bottom'
           easeInOut: can be 'none', 'ease in', 'ease out', 'ease in and out'
           speedDict: dictionary of speed anmes and nbframes (constr. only)
        """
        self.directions = directions
        self.addEaseInOut = addEaseInOut
        if addEaseInOut:
            SE_MAAEditor.__init__(
                self, master=master, title=title, buttons=buttons,
                defaultButton=defaultButton, speedDict=speedDict)
        else:
            MAAEditorWithSpeed.__init__(self, master=master, title=title,
                                    buttons=buttons, defaultButton=defaultButton,
                                    speedDict=speedDict)


    def populateForm(self):
        """
        add radio buttons for direction and easeInOut
        """
        if self.addEaseInOut:
            SE_MAAEditor.populateForm(self)
        else:
            MAAEditorWithSpeed.populateForm(self)

        parent = self.dialog.interior()
        grp = Pmw.Group(parent, tag_text='Direction')
        frame = grp.interior()

        self.directionw = w = Pmw.RadioSelect(
            frame, selectmode='single', orient='horizontal',
            buttontype='radiobutton')

        for text in self.directions:
            w.add(text)

        w.pack(side='left', anchor='w', fill='x', expand=1, padx=8, pady=8)
        grp.pack(side = 'top', fill='x', expand=1)
        self.balloon.bind(w, "Select the action direction")


    def getValues(self):
        """
        return a dictionary of parameterName:parameterValues
        """
        if self.addEaseInOut:
            values = SE_MAAEditor.getValues(self)
        else:
            values = MAAEditorWithSpeed.getValues(self)
            nbFrames = values.pop('nbFrames')
            values.update( {'kfpos': [0,nbFrames]})
        values.update({'direction': self.directionw.getvalue()})
        #print 'SED_MAAEditor, getValues:', values
        return values


    def setValues(self, **kw):
        """
        take a dictionary of parameterName:parameterValues set the editor
        to these values
        """
        if self.addEaseInOut:
            SE_MAAEditor.setValues(self, **kw)
        else:
            MAAEditorWithSpeed.setValues(self, **kw)
        self.directionw.setvalue(kw['direction'])
        #print 'SED_MAAEditor, setValues:', kw



class SESp_MAAEditor(SE_MAAEditor):
    """
    Editor providing speed, easeInOut and sortPoly parameters
    """

    def populateForm(self):
        """
        add radio buttons for direction and easeInOut
        """
        SE_MAAEditor.populateForm(self)

        parent = self.dialog.interior()

        grp = Pmw.Group(parent, tag_text='Zsort Polygons')
        frame = grp.interior()

        ordergrp = Pmw.Group(frame, tag_pyclass = None)
        ordergrp.pack(side = 'left', fill='x', expand=1)
        self.sortOrderw = Pmw.RadioSelect(
            ordergrp.interior(), selectmode='single', orient='horizontal',
            buttontype='radiobutton')
        
        for text in ['+Zsort', '-Zsort']:
            self.sortOrderw.add(text)
        self.sortOrderw.pack(side='top', anchor='w', fill='x', expand=1, padx=8, pady=8)
        self.balloon.bind(self.sortOrderw, "-Zsort sorts by furthest z first,\n+Zsort sorts by closest z first")

        whengrp = Pmw.Group(frame, tag_pyclass = None)
        whengrp.pack(side = 'left', fill='x', expand=1)
        self.sortPolyw = w = Pmw.RadioSelect(
            whengrp.interior(), selectmode='single', orient='horizontal',
            buttontype='radiobutton')

        for text in ['Never', 'Once', 'Always']:
            w.add(text)
        w.pack(side='left', anchor='w', fill='x', expand=1, padx=8, pady=8)
        self.balloon.bind(w, "Select when to Z-sort polygons for proper trancparency")
        
        grp.pack(side = 'top', fill='x', expand=1)
        


    def getValues(self):
        """
        return a dictionary of parameterName:parameterValues
        """

        values = SE_MAAEditor.getValues(self)
        values.update({'sortPoly': self.sortPolyw.getvalue(),
                       'sortOrder': self.sortOrderw.getvalue()})
        #print "SESp_MAAEditor getValues:", values
        return values


    def setValues(self, **kw):
        """
        take a dictionary of parameterName:parameterValues set the editor
        to these values
        """

        SE_MAAEditor.setValues(self, **kw)
        self.sortPolyw.setvalue(kw['sortPoly'])
        self.sortOrderw.setvalue(kw['sortOrder'])
        #print "SESp_MAAEditor setValues:", kw


from mglutil.gui.BasicWidgets.Tk.colorWidgets import ColorChooser

class SECol_MAAEditor(SE_MAAEditor):
    
    """
    Editor providing speed, easeInOut parameters and ColorChooser widget
    """
    def __init__(self, master=None, title='Colors Editor',
                 buttons=['OK', 'Preview', 'Cancel'],
                 defaultButton='OK', speedDict =None, choosecolor = False):
        
        self.choosecolor = choosecolor
        self.color = None
        SE_MAAEditor.__init__(
            self, master=master, title=title, buttons=buttons,
            defaultButton=defaultButton, speedDict=speedDict)


    def populateForm(self):
        """
        add radio buttons for direction and easeInOut
        """
        SE_MAAEditor.populateForm(self)
        if self.choosecolor:
            parent = self.dialog.interior()
            self.colorChooser = ColorChooser(master = parent,
                                             commands = self.setColor_cb)
            self.colorChooser.pack(side = "top")


    def setColor_cb(self, colors):
        self.color = colors


    def getValues(self):
        """
        return a dictionary of parameterName:parameterValues
        """

        values = SE_MAAEditor.getValues(self)
        if self.choosecolor:
            if self.color:
                values['color'] = [self.color,]
        #print "SECol_MAAEditor getValues", values
        return values
    

class Rotation_MAAEditor(SED_MAAEditor):
    """
    Editor providing Rotation animation parameters
    """
    def __init__(self, master=None, title='Rotation Editor',
                 directions=['counter clockwise', 'clockwise'],
                 buttons=['OK', 'Preview', 'Cancel'],
                 defaultButton='OK', speedDict=None, addEaseInOut=True):
        """
        Provides a dialog form for setting Rotation action parameters

        parameters:
           nbFrames: number of keyframes
           angle: rotation angular amplitude in degrees
           vector: rotation axis
           direction: can be 'counter clockwise' or 'clockwise'
           easeInOut: can be 'none', 'ease in', 'ease out', 'ease in and out'
        """

        SED_MAAEditor.__init__(
            self, master=master, title=title, buttons=buttons,
            defaultButton=defaultButton, directions=directions,
            speedDict=speedDict, addEaseInOut=addEaseInOut)


    def populateForm(self):
        """
        add counter for angle and GUI for rotation vector
        """
        SED_MAAEditor.populateForm(self)

        parent = self.dialog.interior()

        # Create and pack the dropdown ComboBox for axis
        grp = Pmw.Group(parent, tag_text='rotation axis')
        axes = self.axes = ('X', 'Y', 'Z', 'XY', 'XZ', 'YZ', 'XYZ')
        self.axisDropdownw = w = Pmw.ComboBox(
            grp.interior(), scrolledlist_items=axes)
        w.pack(side='left', anchor='n', fill='x', expand=1, padx=8, pady=8)
        w.selectitem('Y')
        grp.pack(side = 'top', fill='x', expand=1)

        grp = Pmw.Group(parent, tag_text='rotation angle')
        self.anglew = w = Pmw.Counter(
            grp.interior(), entry_width=6, entryfield_value = 360,
            entryfield_validate = {'validator':'real', 'min' : 1 } )

        w.pack(side='left', anchor='w', fill='x', expand=1, padx=8, pady=8)
        grp.pack(side = 'top', fill='x', expand=1)

        #self.balloon.bind(w, "Select when to Z-sort polygons for proper trancparency")

    def getVector(self):
        """
        return (x, y, z) vector according to axis widget value
        """
        
        axis = self.axisDropdownw.getvalue()[0]
        if axis=='X':
            vector = (1., 0., 0.)
        elif axis=='Y':
            vector = (0., 1., 0.)
        elif axis=='Z':
            vector = (0., 0., 1.)
        elif axis=='XY':
            vector = (1., 1., 0.)
        elif axis=='XZ':
            vector = (1., 0., 1.)
        elif axis=='YZ':
            vector = (0., 1., 1.)
        elif axis=='XYZ':
            vector = (1., 1., 1.)
        return vector


    def getAxis(self, vector):
        """
        return axis anme based in (x, y, z) vector
        """
        x,y,z = vector
        if x==1 and y==0 and z==0:
            return 'X'
        elif x==0 and y==1 and z==0:
            return 'Y'
        elif x==0 and y==0 and z==1:
            return 'Z'
        elif x==1 and y==1 and z==0:
            return 'XY'
        elif x==0 and y==1 and z==1:
            return 'YZ'
        elif x==1 and y==0 and z==1:
            return 'XZ'
        elif x==1 and y==1 and z==1:
            return 'XYZ'
        else:
            return 'custom'
        
    
    def getValues(self):
        """
        return a dictionary of parameterName:parameterValues
        """
        values = SED_MAAEditor.getValues(self)
        nbFrames = values.pop('kfpos')[-1]
        values.update({'nbFrames':nbFrames,
                       'angle': float(self.anglew.get()),
                       'vector': self.getVector(),
                       } )
        #print "Rotation_MAAEditor, getValues:", values
        return values


    def setValues(self, **kw):
        """
        take a dictionary of parameterName:parameterValues set the editor
        to these values
        """
        SED_MAAEditor.setValues(self, **kw)
        self.axisDropdownw.selectitem(self.getAxis(kw['vector']))
        self.anglew.setvalue(kw['angle'])
        #print "Rotation_MAAEditor, setValues:", kw



class Rock_MAAEditor(Rotation_MAAEditor):
    """
    Editor providing Rock animation paramters
    """
    def __init__(self, master=None, title='Rock Editor',
                 directions=['counter clockwise', 'clockwise'],
                 buttons=['OK', 'Preview', 'Cancel'],
                 defaultButton='OK', addEaseInOut=False):
        """
        Provides a dialog form for setting Rotation action parameters

        parameters:
           nbFrames: number of keyframes
           angle: rotation angular amplitude in degrees
           vector: rotation axis
           direction: can be 'counter clockwise' or 'clockwise'
           easeInOut: can be 'none', 'ease in', 'ease out', 'ease in and out'
        """

        Rotation_MAAEditor.__init__(
            self, master=master, title=title, buttons=buttons,
            defaultButton=defaultButton, directions=directions,
            speedDict={'slow': 60, 'medium': 30, 'fast': 10}, addEaseInOut=addEaseInOut)



## def expandGeoms(geoms):
##     newgeoms = []
##     for g in geoms:
##         for child in g.AllObjects():
##             if len(child.getVertices())>0 and child.visible:
##                 newgoms.append( child )
##     return newgeoms


    
class orientationGUI:
    """
    Scrolled frame holding orientations (snapshots)
    """
    def __init__(self, viewer, viewerName, master=None):
        """
        orientationGUI constructor

        orientationGUIObject <- orientationGUI(viewer, viewerName, master=None)
        """
        self.viewer = viewer
        self.viewerName = viewerName
        #self.orientations = {}  
        self.snapshots = {}
        self.nbOrients = 0
        self.master = master
	self.row = 0
	self.col = 0
        self.speedDict = {'slow': 50, 'medium': 30, 'fast': 10}

        #self.editor = OrientationMAAEditor(master=master)
        self.editor = SnapshotMAAGroupEditor(master=master)
        self.modifyOrient = Tkinter.IntVar()
        self.modifyOrient.set(1)
        self.modifyRendering = Tkinter.IntVar()
        self.modifyRendering.set(1)
        self.forceRendering = Tkinter.BooleanVar()
        self.forceRendering.set(False)
        self.createGUI()
        self._animNB = None # will become a reference to the AnimationNotebook instance
        self.nButtonsRow = 5


    def createGUI(self):
        """
        Create a ScrolledFrame to old orientations entries
        """
        if self.master is None:
            self.master = master = Tkinter.Toplevel()
            self.ownsMaster = True
        else:
            self.ownsMaster = False

        self.balloon = Pmw.Balloon(self.master)

        # create a group with a button to record an orientation
        w = self.orientsContainer = Pmw.Group(
            self.master, tag_pyclass = Tkinter.Button,
            tag_text='Record Snapshot')
        w.configure(tag_command = self.recordOrient)

        # create a scrolled frame to display recorded orientation
        
        w1 = self.MAAContainer = Pmw.ScrolledFrame(
            w.interior(), usehullsize=0, hull_width=40, hull_height=200,
            vscrollmode='dynamic', hscrollmode='none')

        w1.pack(padx=5, pady=3, fill='both', expand=1)
       

        w.pack(fill='both', expand=1, padx = 6, pady = 6)

        # bind right button to show option form editor
        button = w.component('tag')
        button.bind('<Button-3>', self.startEditor)


    def startEditor(self, event=None):
        objString = self.viewerName+'.rootObject'
        #orient = getOrientation(object)
        orient = None
        rendering = getRendering(self.viewer, checkAnimatable=True)
        orientMaa = OrientationMAA(self.viewer.rootObject, 'temp', orient, rendering,
                                   objectFromString=objString)
        kfpos = [orientMaa.firstPosition, orientMaa.lastPosition]
        renderMaa = RenderingTransitionMAA(self.viewer, rendering,
                                           kfpos=kfpos, startFlag = "with previous")
        maa = SnapshotMAAGroup(orientMaa, renderMaa,"snapshot%d"% (self.nbOrients+1, ) )

        values = self.editor.edit(maa)
        if values:
            self.nbOrients += 1
            self.saveMAA(maa)
            

    def recordOrient(self, event=None):
        """
        build default orientation transition (left clicks)
        """
        self.nbOrients += 1
        object = self.viewer.rootObject
        #orient = getOrientation(object)
        orient = None
        rendering = getRendering(self.viewer, checkAnimatable=True)
        #maa1 = OrientationMAA(object, 'orient%d'% self.nbOrients, orient, rendering,
        #                     objectFromString=self.viewerName+'.rootObject')
        orientMaa = OrientationMAA(object, 'temp', orient, rendering,
                             objectFromString=self.viewerName+'.rootObject')
        kfpos = [orientMaa.firstPosition, orientMaa.lastPosition]
        renderMaa = RenderingTransitionMAA(self.viewer, rendering,
                                           kfpos=kfpos, startFlag = "with previous")
        maa = SnapshotMAAGroup(orientMaa, renderMaa,"snapshot%d"%self.nbOrients )
        self.saveMAA(maa)


    def saveMAA(self, maagroup):
        """
        adds MAA to the list and adds a button for it in the panel
        """

        assert isinstance(maagroup, SnapshotMAAGroup)
        if not maagroup.name:
            maagroup.name = "snapshot%d"%self.nbOrients
        snName = self.checkName(maagroup.name)
        if maagroup.name != snName: maagroup.name = snName
        
        orientMaa = maagroup.orientMaa
        renderMaa = maagroup.renderMaa
        renderMaa._animNB = self._animNB
        orientMaa.name = snName+"orient"
        renderMaa.name = snName+"rendering"
        self.snapshots[snName] = maagroup
        self.addOrientButton(maagroup)
        
    def checkName(self, name):
        """check if the name exists in the self.snapshots or in the sequence player.
        If exists - create unique name"""
        allnames = self.snapshots.keys()
##         if self._animNB:
##             for maa, pos in  self._animNB().seqAnim.maas:
##                 if maa.name not in allnames:
##                     allnames.append(maa.name)
        if name in allnames:
            i = 1
            while(name in allnames):
                name = name+"_%d"%i
                i = i+1
        return name
        
        
    def addOrientButton(self, maa):
        master = self.MAAContainer.interior()
        if hasattr(maa, 'ims') and maa.ims is not None:
            from PIL import ImageTk
            photo = ImageTk.PhotoImage(maa.ims)
        else:
            self.viewer.master.lift()
            self.viewer.master.master.lift()
            self.viewer.OneRedraw()
            photo, ims = self.viewer.GUI.getButtonIcon()
            maa.ims = ims
        b = Tkinter.Button(master=master ,compound='left', image=photo,
                           command=CallbackFunction(self.runMaa, maa))
        b.photo = photo

        b.name = maa.name
        b.grid(row = self.row, column = self.col, sticky = 'nsew')
        b.bind('<Button-3>', CallbackFunction( self.showOrientMenu_cb, maa))

        self.balloon.bind(b, maa.name)
        if self.col == self.nButtonsRow-1:
	    self.col = 0
	    self.row = self.row + 1
	else:
	    self.col = self.col + 1

    def runMaa(self, maagroup):
        orient = maagroup.orientMaa
        render = maagroup.renderMaa
        #print "run maa:", maagroup.name, 'force rendering:',  orient.forceRendering       
        if orient.forceRendering:
            setRendering(orient.viewer, orient.rendering)
            orient.run()
        else:
            #modify (morph) rendering
            #render.setValuesAt(0)
            #render.run()
            #orient.run()
            maagroup.run()


    def editMaa_cb(self, maagroup):
        values = self.editor.edit(maagroup)
        #check if the maa has been added to the sequence animator:
        animNB = self._animNB()
        for i , _maa in enumerate(animNB.seqAnim.maas):
            if _maa[0] == maagroup:
                position = _maa[1]
                animNB.seqAnimGUI.update(i, position)
                return


    def setForceRendering(self, orient, event = None):
        #print "setForceRendering", self.forceRendering.get()
        orient.forceRendering = self.forceRendering.get()

 
    def showOrientMenu_cb(self, maagroup, event=None):
        """
        Create button menu and post it
        """
        # create the button menu
        orient = maagroup.orientMaa
        render = maagroup.renderMaa
        #orient, render = maagroup.maas
        menu = Tkinter.Menu(self.master, title = orient.name)
        
        #cb = CallbackFunction(self.setForceRendering, orient)

        #self.forceRendering.set(orient.forceRendering)
        #menu.add_checkbutton(label="Force Rendering",
        #                     var = self.forceRendering,
        #                     command=cb)
        
        from Scenario2 import addTargetsToMenu
        #addTargetsToMenu(menu, [orient, render])
        #addTargetsToMenu(menu, maagroup)
        cb = CallbackFunction(self.addAsTransition_cb, maagroup)
        #menu.add_command(label="Add to animation as transition", command = cb)
        menu.add_command(label="Add to animation", command = cb)
        #cb = CallbackFunction(self.addAsKeyframe_cb, maagroup)
        #menu.add_command(label="Add to animation as keyframe", command = cb)
        
        cb = CallbackFunction(self.editMaa_cb, maagroup)
        menu.add_command(label="Edit", command = cb)
        
        cb = CallbackFunction(self.renameOrient_cb, maagroup)
        menu.add_command(label="Rename", command = cb)

        cb = CallbackFunction(self.removeOrient_cb, maagroup)
        menu.add_command(label="Delete", command = cb)

        menu.add_command(label="Dismiss")
        menu.post(event.x_root, event.y_root)


    def addToClipboard(self, orient, render=None):
        """
        adds this orientation animation to the clipboard
        """
        _clipboard.addMaa(orient)
        if render is not None:
            _clipboard.addMaa(render)

    def addAsTransition_cb(self, maagroup):
        kf1, kf2 = maagroup.kfpos
        #if kf2 - kf1 <=1:
        #    values = self.editor.edit(maagroup)
        self._animNB().seqAnim.addMAA(maagroup)

    def addAsKeyframe_cb(self, maagroup):
        maagroup.setKeyframePositions([0, 1])
        self._animNB().seqAnim.addMAA(maagroup)
        


    def renameOrient_cb(self, maagroup):
        name = maagroup.name
        container = self.MAAContainer.interior()
        from tkSimpleDialog import askstring
        newname = askstring("Rename %s"%name, "Enter new name:", initialvalue = name,
                            parent = container)
        if newname != None and newname != name:
            if self.snapshots.has_key(newname):
                from tkMessageBox import showwarning
                showwarning("Warning", "Name %s already exists"%newname,parent = self.master)
                return
            #find cooresponding button, rename it and update the bindings:
            self.snapshots.pop(name)
            self.snapshots[newname] = maagroup
            maagroup.name = newname
            orient = maagroup.orientMaa
            render = maagroup.renderMaa
            #orient, render = maagroup.maas
            orient.name = newname+"orient"
            render.name = newname+"rendering"
            for b in container.grid_slaves():
                if hasattr(b, "name"):
                    if b.name == name:
                       b.name = newname
                       self.balloon.bind(b, newname)
                       break
            seqmaas = [maa[0] for maa in self._animNB().seqAnim.maas]
            if maagroup in seqmaas:
                self._animNB().seqAnimGUI.refreshGUI()


    def removeOrient_cb(self, maagroup):
        orientB = None
        name = maagroup.name
        frame = self.MAAContainer.interior()
        for b in frame.grid_slaves():
            if hasattr(b, "name"):
                if b.name == name:
                    orientB = b
                    break
        if orientB:
            # check if this maa has been added to the sequence animator:
            seqAnim = self._animNB().seqAnim
            seqmaas = [maa[0] for maa in seqAnim.maas]
            ind = None
            if maagroup in seqmaas:
                
                # ask the user i f he really wants to delete this snapshot
                # since it will be removed from the sequence animator too:
                import tkMessageBox
                ok = tkMessageBox.askokcancel("Delete %s?"%name,"%s is in Sequence Anim.\nThis will also remove the snapshot from Sequence Anim."%name)
                if not ok:
                    return
                ind = seqmaas.index(maagroup)
            orientB.destroy()
            # regrid the buttons to fill the space freed by the removed button :
            buttons = frame.grid_slaves() # the widgets in this list
            # seem to be stored in "last created, first in the list" order
            buttons.reverse()
            col = 0
            row = 0
            for i, b in enumerate(buttons):
                b.grid(row=row, column= col, sticky='nsew')
                if col == self.nButtonsRow-1:
                    col = 0
                    row = row + 1
                else:
                    col = col + 1
            self.col = col
            self.row = row

            # remove the orient entry from self.orientations
            self.snapshots.pop(name)
            if ind != None:
                # remove from sequence anim.
                seqAnim.removeMAA(maagroup, seqAnim.maas[ind][1])


    def getSavedMaas(self):
        maas = []
        for b in reversed(self.MAAContainer.interior().grid_slaves()):
            if hasattr(b, "name"):
                name = b.name
                if self.snapshots.has_key(name):
                    maas.append(self.snapshots[name])
        return maas



import Pmw, Tkinter
from Scenario2 import  _clipboard

class RotationMAAOptionsGUI:

    def __init__(self, master, object, objectName):

        self.object = object
        self.objectName = objectName
        if not master:
            self.parent = parent = Tkinter.Toplevel()
        else:
            self.parent = parent = master
        self.axis = 'Y'
        self.direction = 'clockwise'
        self.maa = None
	self.dialog = Pmw.Dialog(parent,
	    buttons = ('OK', 'Preview', 'Cancel', 'Copy To Clipboard'),
	    defaultbutton = 'Preview',
	    title = 'Full Rotation Options:',
	    command = self.execute)
	self.dialog.withdraw()

        parent = self.dialog.interior()
        
        self.buttonBox = Pmw.ButtonBox(
            parent, labelpos = 'nw', label_text = 'Full Rotation Options:',
            frame_borderwidth = 2, frame_relief = 'groove')

        # animation name
        self.namew = Pmw.EntryField(parent, labelpos = 'w', label_text = 'Name:',
                               value = 'Full Rotation Y',)

        # file name for movie
        self.filenamew = Pmw.EntryField(parent, labelpos = 'w',
                                   label_text = 'Filename:',
                                   value='pmvMovie.mpeg')
        
        self.namew.pack(side='top', fill='x', expand=1, padx=10, pady=5)
        self.filenamew.pack(side='top', fill='x', expand=1, padx=10, pady=5)

        # Create and pack the dropdown ComboBox for axis
        axes = ('X', 'Y', 'Z', 'XY', 'XZ', 'YZ', 'XYZ')
        self.axisDropdownw = Pmw.ComboBox(
            parent, label_text = 'Axis:', labelpos = 'w',
            selectioncommand = self.changeAxis, scrolledlist_items = axes)
        self.axisDropdownw.pack(side = 'top', anchor = 'n',
                                fill = 'x', expand = 1, padx = 8, pady = 8)
        self.axisDropdownw.selectitem(self.axis)

        self.directionw = Pmw.ComboBox(
            parent, label_text = 'direction:', labelpos = 'w',
            selectioncommand = self.changeDirection,
            scrolledlist_items = ['clockwise', 'counter clockwise'])
        self.directionw.pack(side = 'top', anchor = 'n',
                             fill = 'x', expand = 1, padx = 8, pady = 8)
        self.directionw.selectitem(self.direction)

        # frame number
        self.nbframesw = Pmw.Counter(parent,
                                labelpos = 'w',
                                label_text = 'Number of Frames:',
                                entry_width = 6,
                                entryfield_value = 180,
                                entryfield_validate = {'validator' : 'integer',
                                                       'min' : 1 }
                                )
        self.nbframesw.pack(side='top', padx=10, pady=5)

        self.anglew = Pmw.Counter(parent,
                                  labelpos = 'w',
                                  label_text = 'Angle (degrees):',
                                  entry_width = 6,
                                  entryfield_value = 360,
                                  entryfield_validate = {'validator':'integer',
                                                       'min' : 1 }
                                  )
        self.anglew.pack(side='top', padx=10, pady=5)

        # record checkbutton
        self.recordVar = Tkinter.IntVar()

        c = Tkinter.Checkbutton(parent, text="Record", variable=self.recordVar)
        c.pack(side='top', padx=10, pady=5)


    def changeAxis(self, axis):
        self.axis = axis
        self.namew.setvalue('Full Rotation '+self.axis)


    def changeDirection(self, direction):
        self.direction = direction
        self.namew.setvalue('Full Rotation '+self.axis)

    def _processReturnKey(self, event):
	self.buttonBox.invoke()


    def execute(self, result):
        if result == 'Preview' or result=='OK':
	    self.preview_cb()

        if result == 'OK' or result=='Cancel':
            self.hide()

        if result == 'Copy To Clipboard':
            if self. maa is None:
                return
            
            from tkSimpleDialog import askstring

            angle = float(self.anglew.get())
            nbFrames = int(self.nbframesw.get())
            name = 'Rotation of %d degrees in %d steps about %s'%(
                angle, nbFrames, str(self.axis))

            d = self.buildMAA()
            _clipboard.addMaa(d)


    def getName(self, result):
	if result is None or result == 'Cancel':
	    self.clipboardName = None
	    self.dialog.deactivate(result)
	else:
	    if result == 'OK':
		print 'Password entered ' + self.dialog.get()
		self.dialog.deactivate()

    def buildMAA(self):
        angle = float(self.anglew.get())
        nbFrames = int(self.nbframesw.get())
        
        name = 'Rotation of %d degrees in %d steps about %s'%(
            angle, nbFrames, str(self.axis))
        
        self.maa = d = RotationMAA(
            self.object, self.objectName, angle=angle, axis=self.axis,
            nbFrames=nbFrames, name=self.namew.getvalue(),
            direction=self.directionw.get())

        return d

    
    def preview_cb(self):
        d = self.buildMAA()
        if self.recordVar.get():
            d.redrawActor.startRecording(self.filenamew.getvalue())
        d.run()
        if self.recordVar.get():
            d.redrawActor.stopRecording()


    def hide(self):
        self.dialog.deactivate()

    def show(self, event=None):
        self.dialog.activate(globalMode = 'nograb')


