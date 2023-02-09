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

"""
This module provides a Geoemtry picker for DejaVu2 Geometry objects
"""

from DejaVu2 import Viewer
from mglutil.gui.BasicWidgets.Tk.customizedWidgets import ObjectChooser
import Tkinter

class GeomChooser:
    """
    The GeomChooser object display a list of geometry object present in a given
    DejaVu2.Viewer object.
    """

    def __init__(self, viewer, showAll=True, filterFunction=None, root=None,
                 command=None, refreshButton=False, showAllButton=False):
        """
        GeomChooser constructor

        GeomChooserObject <- GeomChooser(viewer)

        - viewer is an instance of a DejaVu2.Viewer object
        - showAll is a boolean. When True all objects are shown
        - filterFunction is an optional function that will be called when
        showAll is False. I takes one argument (the viewer) and returns a
        list of (name, [geoms]) pairs. The names will be displayed in the GUI
        for representing the corresponding list of geometry objects.
        - root is the master in which the chooser will be created.
        - command is the fucntion call upon selection in the list
        """

        assert isinstance(viewer, Viewer)
        self.root = root
        self.viewer = viewer
        self.frame = Tkinter.Frame(root)
        self.frame1 = None
        if refreshButton or showAllButton:
            f1 = self.frame1 = Tkinter.Frame(root)
            if refreshButton:
                self.refresh = Tkinter.Button(f1, text='Refresh Lists',
                                              command=self.updateList)
                #self.refresh.pack(side='left', fill="x", expand=1)
                self.refresh.pack(side='top', fill="x", expand=1,pady=3)
            if showAllButton:
                self.showAllVar = Tkinter.IntVar()
                self.showAllVar.set(0)
                self.show = Tkinter.Checkbutton(
                    f1, text='Show all objects', variable=self.showAllVar,
                    command=self.showAll_cb)
                #self.show.pack(side='left', fill="x", expand=1)
                self.show.pack(side='top', fill="x", expand=1, pady=3)
            self.frame1.pack(side = 'top', fill='x', expand=1)
            
        assert isinstance(showAll, bool)
        self.showAll = showAll
        
        if filterFunction is not None:
            assert callable(filterFunction)
        self.filterFunction = filterFunction
        self.geomList = []
        self.root = root
        
        self.chooserW = ObjectChooser(
            self.frame, mode='single', title='Choose Object',
            command=command, lbpackcfg={'fill':'both', 'expand':1} )
        #self.chooserW.pack(side='left', fill='both', expand=1)
        self.chooserW.pack(side='top', fill='both', expand=1)
        self.buildListOfGeoms()
        self.setList()
        self.chooserW.widget.configure(exportselection = False) # this is done to prevent disappearing
#of the ListBox selection  when a text is selected in another widget or window.

    def onSelect(self, event):
        pass
    

    def buildListOfGeoms(self):
        self.geomList = []
        if self.showAll or self.filterFunction is None:
            for g in self.viewer.rootObject.AllObjects():
                self.geomList.append( (g.fullName, [g]) )

        else:
            self.geomList = self.filterFunction(self.viewer)

        #for n, g in self.geomList:
        #    print n

            
    def setList(self):
        self.chooserW.clear()
        self.chooserW.setObjects(self.geomList)


    def updateList(self):
        lb = self.chooserW.lb
        selection = lb.curselection()
        self.buildListOfGeoms()
        self.setList()
        for i in selection:
            lb.selection_set(i)
            

    #def get(self):
    #    return self.chooserW.get()
    def get(self, event = None):
        res = []
        sel = self.chooserW.lb.curselection()
        if not sel:
            return res

        for ent in map( int, sel ):
            obj= self.geomList[ent][1]
            res.append(obj)
        if self.chooserW.mode=='single' and res:
            res = res[0]

        return res
    


    def showAll_cb(self, even=None):
        val = self.showAllVar.get()
        self.showAll = val
        self.updateList()


    def grid(self,**kw):
        """Applies the specified keywords (kw) to the grid method of the frame
        containing geometry chooser"""
        self.frame.grid(kw)
        if not kw.has_key('row'):
            row = 0
        else:
            row = kw['row']
        if not kw.has_key('column'):
            column=0
        else:
            column=kw['column']
        if kw.has_key('weight'):
            weight = kw['weight']
        else:
            weight = 1
        self.frame.rowconfigure(row, weight=weight)
        self.frame.columnconfigure(column, weight=weight)


    def grid_forget(self):
        """Calls grid_forget() method of the frame containing
        geometry chooser"""
        self.frame.grid_forget()


    def pack(self, **kw):
        """Applies the specified keywords (kw) to the pack method of the frame
        containing geometry chooser"""
        self.frame.pack(kw)

        
    def pack_forget(self):
        """Calls pack_forget() method of the frame containing
        geometry chooser"""
        self.frame.pack_forget()




