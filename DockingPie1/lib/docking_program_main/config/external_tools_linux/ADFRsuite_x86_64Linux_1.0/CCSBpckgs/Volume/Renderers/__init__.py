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

#__init__.py
import Tkinter
from mglutil.gui.BasicWidgets.Tk.customizedWidgets import ListChooser
istest = 0
renderer = None

def setvals(val1, val2):
    global istest 
    global renderer
    istest = val1
    renderer = val2

class ListChooserDialog:
    """A dialog widget consisting of a ListChooser, an 'OK' button and
    an optional 'Cancel' button."""

    def __init__(self, master, title = '', text = '',
                 entries = (('',None),('',None)) , cancel = None,
                 mode = 'single', list_width = None, list_height=None,
                 list_font = None, list_command = None):
        
        assert mode in ['single', 'browse', 'multiple', 'extended' ]
        self.root = Tkinter.Toplevel(master)
        if title:
            self.root.title(title)
            self.root.iconname(title)
        self.entry = None
        self.entries = map(lambda x: x[0], entries)
        self.frame = Tkinter.Frame(self.root)
        self.frame.pack()
        self.root.bind('<Return>', self.return_event)
        list_cfg = {}
        if list_height:
            list_cfg['height']=list_height
        if list_width:
            list_cfg['width']=list_width
        
        self.listchooser = ListChooser(self.frame,
                                       title=text,
                                       entries = entries,
                                       lbwcfg=list_cfg,
                                       command=list_command)
        self.listchooser.pack(fill = Tkinter.BOTH, expand=1, padx=5, pady=5)
        ok_button = Tkinter.Button(self.frame, text="OK",
                       command=(lambda self=self, num=0: self.done(num)))
        ok_button.pack(side=Tkinter.LEFT, fill=Tkinter.BOTH, expand=1)
        if cancel :
            cancel_button = Tkinter.Button(self.frame, text="Cancel",
                       command=(lambda self=self, num=1: self.done(num)))
            cancel_button.pack(side=Tkinter.LEFT, fill=Tkinter.BOTH,
                               expand=1)
        if list_font:
            self.set_newfont(list_font)
        self.root.protocol('WM_DELETE_WINDOW', self.wm_delete_window)

    def go(self):
        self.root.grab_set()
        self.root.mainloop()
        self.root.destroy()
        return self.entry

    def return_event(self, event):
        if self.default is None:
            self.root.bell()
        else:
            self.done(self.default)

    def wm_delete_window(self):
        if self.cancel is None:
            self.root.bell()
        else:
            self.done(self.cancel)

    def choseEntry(self, entry):
        self.entry = entry
        
    def done(self, num):
        if num == 0:
            self.entry = self.listchooser.get()
        elif num == 1:
            self.entry = []
        self.root.quit()

    def setentry(self, index):
        ent = self.entries[index]
        self.listchooser.set(ent)

    def set_newfont(self, newfont):
        self.changefont(self.root, newfont)

    def changefont(self, wid, newfont):
        try:
            wid.config(font=newfont)
        except :
            pass
        if len(wid.children)==0: return
        for item in wid.children.values():
            self.changefont(item, newfont)

            
def ChooseRenderer(renderer = None):
    import Tkinter
    liblist = []
    geomdict ={}
    dialog = None
    ans = None
    
    print "************ Volume Renderer Info ************************"
    if not renderer or renderer=="vli":
        print "Trying to import VolumePro VLI library ..."
        try:
            from VLI import vli
            from VLI.DejaVu.VLIGeom import VLIGeom
            vliGeom = VLIGeom('vli', protected=1)
            liblist.append('vli')
            geomdict['vli']= vliGeom
            print "VLI library imported."
        except:
            print "could not import VLI library."

    if not renderer or renderer == 'utvolren':
        print "Trying to import UTVolumeLibrary..."    
        try:
            from UTVolumeLibrary import UTVolumeLibrary
            from UTVolumeLibrary.DejaVu.UTVolRenGeom import UTVolRenGeom
            utvolGeom = UTVolRenGeom('utvolren')
            liblist.append('utvolren')
            geomdict['utvolren'] = utvolGeom
            print "UTVolumeLibrary imported."
        except:
            print "could not import UTVolumeLibrary."
    print "**********************************************"

    flag = None
    
    if len(liblist)==1:
        flag = liblist[0]

    elif len(liblist) > 1:
        # the user has to choose a library
        root = Tkinter.Tk()
        root.withdraw()
        entries = (('VolumePro VLI',None), ('OpenGL based VolumeLibrary',None))
        dialog = ListChooserDialog(root, title = 'Choose Volume Library',
                            text = 'Choose a library', entries = entries,
                                   list_width = 27, list_font = ('Arial',14) )
        dialog.setentry(0)
        ans = dialog.go()

        #print "ans: ", ans
        if ans:
            #print "Lib:", ans
            flag = None
            if ans[0] == entries[0][0]:
                flag = 'vli'
            elif ans[0] == entries[1][0]:
                flag = 'utvolren'

    else:
        print "No volume rendering library imported."
    return (flag, geomdict)



