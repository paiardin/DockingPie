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

##
## GUI panels for creation Actions or MAA for example used in
## PMV/scenarioInterface/animationGui.py AnimationNotebook
##
import Tkinter, Pmw
from Scenario2.multipleActorsActions import MultipleActorsActions
from mglutil.util.callback import CallbackFunction
from Scenario2 import  _clipboard
from DejaVu2.states import setRendering, setOrientation

class AnimationPanel:
    """
    base class for creating a panel in which .

    This class create a frame (self.frame) containing:
        a group for geometry chooser self.geomChooserG
        a group for widget to create Actions self.makeActionsG
        a group for created actions: self.createdActionsG

    the created MAAs are in self.maa list
    """

    def __init__(self, viewer, viewerName, master=None):
        """
        Creates a GUI form containing a groupGeometry Chooser widget and a Pmw.Group widget
        that can be populated with animation effects buttons.
        
        viewer - an instance of molecular viewer
        """
        self.viewer = viewer
        self.viewerName = viewerName
        self.master = master

        self.createGUI()
        self.maas = [] # created MAAs
        self.maaButtons = [] # buttons representing created MAAs
        self.maaButtonNames = {} # for each short name count how many we have
        self._animNB = None # will become a reference to the AnimationNotebook instance


    def createGUI(self):

        # set or create master 
        if self.master is None:
            self.master = master = Tkinter.Toplevel()
            self.ownsMaster = True
        else:
            self.ownsMaster = False

        # create balloon
        self.balloon = Pmw.Balloon(self.master)
        # create top frame
        self.frame = fm = Tkinter.Frame(self.master)

        # create group for geometry chooser
	self.geomChooserG = w = Pmw.Group(fm, tag_text='Animateable objects')
        w.pack(side='left', fill='both', expand=1, padx=6, pady=6)

        # add a group for animation buttons (clips)
        self.makeActionsG = w = Pmw.Group(fm, tag_text='Create Action')
	w.pack(side='top', fill='both', expand=1, padx=6, pady=6)
        self.effectsContainer= w.interior()

        # add a group for saved animation clips
        self.createdActionsG = w = Pmw.Group(fm, tag_text="Created Actions")
        w.pack(side='top', fill='both', expand=1, padx=6, pady=6)

        self.frame.pack(side='left', fill='both', expand=1)


    def makeMAA(self, maaClass, args, kw):
        """
        create an MAA maaClass( *args, **kw) ans add it to self.maas
        """
        maa = maaClass( *args, **kw )
        if len(maa.actors) > 0:
            self.makeMAAButton(maa)


    def runMaa(self, maa):
        if hasattr(maa, "geomOrients") and maa.setOrient:
            for o, val in maa.geomOrients.items():
                #print "setting orient for:", o.name
                setOrientation(o, val)
        elif maa.forceOrient:
            setOrientation(maa.viewer.rootObject, maa.orient)
        if maa.forceRendering:
            setRendering(maa.viewer, maa.rendering)
        maa.run()
        

    def makeMAAButton(self, maa):
        """
        Places a button into created Actions group
        """

        assert isinstance( maa, MultipleActorsActions)
        # add the maa
        self.maas.append(maa)
        
        master = self.createdActionsG.interior()

        cb = CallbackFunction( self.runMaa, maa)
        b = Tkinter.Button(master=master ,compound='left',
                           command=cb, width=0, height=0)
        self.maaButtons.append(b)
        
        if hasattr(maa, "shortName"):
            shortName = maa.shortName
        else:
            shortName = maa.__class__.__name__
        
        if not self.maaButtonNames.has_key(shortName):
            self.maaButtonNames[shortName] = 1

        count = self.maaButtonNames[shortName]
        self.maaButtonNames[shortName] += 1
        
        text = "%s%d" % (maa.shortName, count)
        b.configure(text = text)
            
        nbuttons = len(self.maas)-1
        row = nbuttons/4
        col = nbuttons%4

        b.grid(row = row, column = col, sticky = 'w')
        #self.balloon.bind(b, maa.name)
        # the balloon will display a string returned by maa.__repr__()  
        self.balloon.bind(b, maa)
        
        b.bind('<Button-3>', CallbackFunction( self.maaMenu_cb, maa, b))


    def maaMenu_cb(self, maa, button, event=None):
        """
        Opens the saved clip menu (right mouse click on the clip button)
        """
        menu = Tkinter.Menu(self.master, title=maa.name)

##         cb = CallbackFunction(_clipboard.addMaa, maa)
##         menu.add_command(label="Add to Clipboard", command=cb)

        from Scenario2 import addTargetsToMenu
        addTargetsToMenu(menu, maa)

        if maa.editorClass:
            cb = CallbackFunction(self.editMAA_cb, maa, button)
            menu.add_command(label="Edit ...", command=cb)

        menu.add_command(label = "Delete",
                         command=CallbackFunction(self.deleteMAA_cb, maa, button) )
        menu.add_command(label="Dismiss menu")

        menu.post(event.x_root, event.y_root)


    def editMAA_cb(self, maa, button=None, args=(), kw={}, event=None):
        kw = maa.editorKw
        args = maa.editorArgs
        kw['master'] = self.master
        editor =  maa.editorClass( *args, **kw)
        if editor == None: return "Cancel" 
        editor.edit(maa)
        # the maa name could have been modified:
        #if button:
        #    self.balloon.bind(button, maa.name)
            
        # if this maa is in the sequence animator it should update !
        animNB = self._animNB()
        for i , _maa in enumerate(animNB.seqAnim.maas):
            if _maa[0] == maa:
                position = _maa[1]
                animNB.seqAnimGUI.update(i, position)
                return
        
        return editor.exitStatus
    
        

    def deleteMAA_cb(self, maa, button, event=None):
        """removes the specified maa and its button from Created Actions frame"""
        
        # if this maa is in the sequence animator it will be removed too
        position = None
        seqAnim = self._animNB().seqAnim
        for i , _maa in enumerate(seqAnim.maas):
            if _maa[0] == maa:
                name = maa.name
                import tkMessageBox
                ok = tkMessageBox.askokcancel("Delete %s?"%name,"%s is in Sequence Anim.\nThis will also remove the action from Sequence Anim."%name)
                if not ok: return
                position = _maa[1]
                break
        button.destroy()
        master = self.createdActionsG.interior()
        buttons = master.grid_slaves()
        buttons.reverse()
        nb = 4
        for i, b in enumerate(buttons):
            #print b.maa.name, "row", i/nb, "column", i%nb
            b.grid(row=i/nb, column= i%nb, sticky='w')
        if maa in self.maas:
            self.maas.remove(maa)
        if position != None:
            seqAnim.removeMAA(maa, position)



from DejaVu2.scenarioInterface.animations import \
     FlyInObjectMAA, FlyOutObjectMAA,\
     FadeInObjectMAA, FadeOutObjectMAA, VisibleObjectMAA, ColorObjectMAA,\
     RotationMAAOptionsGUI, RotationMAA, RockMAA


class ShowHideGeomPanel(AnimationPanel):
    """
    Panel providing buttons to create fly, fade, show/hide actions
    """

    def __init__(self, viewer, viewerName, geomChooserClass,
                 geomChooserClassOpt, master=None):
        """

        viewer is a DejaVu2 Viewer instance
        viewerName is a string which when evaluated yields viewer
        geomChooserClass is class defining the geomChooser object to create
        geomChooserClassOpt is a sequence or *args and **kw used to create the
          geomChooser. (viewer,) is automatically inserted at the begining of
          args and 'root', and 'command' are automatically set in kw.
          other possible kw include 'filterFunction', 'refreshButton',
          'showAllButton'
        """
        AnimationPanel.__init__(self, viewer, viewerName, master)

        args, kw = geomChooserClassOpt
        kw['root'] = self.geomChooserG.interior()
        kw['command'] = self.onSelect_cb
        args = (viewer,) + args
        gc = self.geomChooser = geomChooserClass( *args, **kw)
        gc.pack(side='top', fill='both', expand=1, anchor="w")

        self.selectedGeom = None

        # add action creating buttons
        lastRow = 0
        parent = self.makeActionsG.interior()
        
        # fly in button
        cb = CallbackFunction(self.makeMAA, FlyInObjectMAA, (), {})
        w = self.flyinB = Tkinter.Button(parent, text='Fly In', command=cb)
        w.grid(column=0, row=lastRow, sticky='ew')
        w.bind('<Button-3>', CallbackFunction( self.showMaaEditor_cb, FlyInObjectMAA, (), {}))

        # fly out button
        cb = CallbackFunction(self.makeMAA, FlyOutObjectMAA, (), {})
        w = self.flyoutB = Tkinter.Button(parent, text='Fly Out', command=cb)
        w.grid(column=1, row=lastRow, sticky='ew')
        w.bind('<Button-3>', CallbackFunction( self.showMaaEditor_cb, FlyOutObjectMAA , (), {}))
        lastRow += 1

        # fade in button
        cb = CallbackFunction(self.makeMAA, FadeInObjectMAA, (), {})
        w = self.fadeinB = Tkinter.Button(parent, text='Fade In', command=cb)
        w.bind('<Button-3>', CallbackFunction( self.showMaaEditor_cb, FadeInObjectMAA, (), {}))
        w.grid(column=0, row=lastRow, sticky='ew')

        # fade out button
        cb = CallbackFunction(self.makeMAA, FadeOutObjectMAA, (), {})
        w = self.fadeoutB = Tkinter.Button(parent, text='Fade Out', command=cb)
        w.bind('<Button-3>', CallbackFunction( self.showMaaEditor_cb, FadeOutObjectMAA, (), {}))
        w.grid(column=1, row=lastRow, sticky='ew')
        lastRow += 1

        # show button
        cb = CallbackFunction(self.makeMAA, VisibleObjectMAA, (), {'visible':1})
        w = self.showB = Tkinter.Button(parent, text='Show', command=cb)
        w.grid(column=0, row=lastRow, sticky='ew')
        
        # hide button
        cb = CallbackFunction(self.makeMAA, VisibleObjectMAA, (), {'visible':0})
        w = self.hideB = Tkinter.Button(parent, text='Hide', command=cb)
        w.grid(column=1, row=lastRow, sticky='ew')
        lastRow += 1

        # rotation button
        cb = CallbackFunction(self.makeMAA, RotationMAA, (), {})
        w = self.rotationB = Tkinter.Button(parent, text='Rotate', command=cb)
        w.grid(column=0, row=lastRow, sticky='ew')
        w.bind('<Button-3>', CallbackFunction( self.showMaaEditor_cb, RotationMAA, (), {}))

        # rock button
        cb = CallbackFunction(self.makeMAA, RockMAA, (), {})
        w = self.rockB = Tkinter.Button(parent, text='Rock', command=cb)
        w.grid(column=1, row=lastRow, sticky='ew')
        w.bind('<Button-3>', CallbackFunction( self.showMaaEditor_cb, RockMAA, (), {}))
        lastRow += 1

        self.lastRow = lastRow

    def getSelectedGeoms(self):
        gc = self.geomChooser
        geometries = gc.get()
        kw = {}
        if len(geometries):
            # get the name of the currently selected geometry
            en = gc.chooserW.entries
            ind = gc.chooserW.getInd()[0]
            objname= en[ind][0].strip() # remove leading blanks
            # build a name
            gparent = geometries[0].parent
            if gparent is not None and gparent.name != "root":
                objname = gparent.name + "|" + objname
            kw['objectName'] = objname
            if hasattr(gc, "getNodes"):
                kw['nodes'] = gc.getNodes()
        return geometries, kw
            
        
    def makeMAA(self, maaClass, args, kw, event=None):
        """
        callback for action creating buttons
        """
        geometries, geomkw = self.getSelectedGeoms()
        if len(geometries)==0:
            from tkMessageBox import showwarning
            showwarning("Warning", 'No geometry selected',
                        parent = gc.root)
            return
        kw.update(geomkw)
        args = (geometries, )
        maa = maaClass( *args, **kw )
        self.makeMAAButton(maa)


    def showMaaEditor_cb(self, maaClass, args, kw, event=None):
        """
        open maa editor, create maa based on specified options
        """
        geometries, geomkw = self.getSelectedGeoms()
        if len(geometries)==0:
            from tkMessageBox import showwarning
            showwarning("Warning", 'No geometry selected',
                        parent = gc.root)
            return
        kw.update(geomkw)
        args = (geometries, )
        maa = maaClass( *args, **kw )
        st = self.editMAA_cb(maa)
        if st == "OK":
            self.makeMAAButton(maa)


    def onSelect_cb(self, event = None):
        self.selectedGeom = [self.geomChooser.chooserW.getInd(),
                             self.geomChooser.chooserW.get()]

            
        
