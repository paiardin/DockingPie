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


import Tkinter
import Pmw
import tkFileDialog


class fromFileForm:
    """ A class to build a GUI form to get the user's input for creating an  actor from file data"""

    def __init__(self, scenario, propnames, object, command=None):

        master = Tkinter.Toplevel(scenario.application.master)
        master.title("Create Actor from File")
        self.master = master
        #master.protocol('WM_DELETE_WINDOW', self.dismissCreateActor_cb )
        frame = Tkinter.Frame(master, borderwidth=1, relief='ridge')
        l = Tkinter.Label(frame, text = "Object: %s" % object.fullName)
        l.grid(column=0, row=0, sticky='w')
        self.chooser =Pmw.ComboBox(frame, labelpos ='n',
                              label_text='Actors: ',
                              entryfield_value=propnames[0],
                              scrolledlist_items=propnames,
                              fliparrow=1,
                              selectioncommand=self.updateFuncName)
        self.chooser.grid(column=0, row=1)#, rowspan=3)

        self.startEF = Pmw.EntryField(frame, labelpos = 'n',
                            label_text = 'start frame: ',
                            entry_width = 8,
                            value = 0,
                            validate = {'validator': 'numeric'})
        self.startEF.grid(column = 1, row = 0)

        self.endEF = Pmw.EntryField(frame, labelpos = 'n',
                            label_text = 'end frame: ',
                            entry_width = 8, )
                            #validate = {'validator': 'numeric'})
        self.endEF.grid(column = 2, row = 0)
        self.funcnameEF = Pmw.EntryField(frame, labelpos = 'n',
                            label_text = 'function name: ',
                            value = "set_"+propnames[0],)

        self.funcnameEF.grid(column = 1, row=1, columnspan = 2)
        b = self.createB = Tkinter.Button(frame, text='Create From File ...',
                           command=self.getValues_cb)

        b.grid(column=1, row=3, columnspan = 2)
        frame.pack()
        self.frame = frame
        self.command = command
        self.object = object
        self.file = None


    def updateFuncName(self, val):
        
        self.funcnameEF.setvalue("set_" + val)


    def getValues_cb(self):

        actor = self.chooser.get()
        #print "actor:", actor
        try:
            kf1 = int(self.startEF.get())
        except:
            kf1 = 0
        #print "start frame:", kf1
        try:
            kf2 = int(self.endEF.get())
        except:
            kf2 = -1
            
        #print "end frame:", kf2
        functionName = self.funcnameEF.get()
        if not functionName:
           functionName = "func" 
        
        if actor is not None:
            file = tkFileDialog.askopenfilename(parent = self.master,
                                  initialdir = '.', title='Actor file',
                                  filetypes=[('', '*.py'), ('all', '*')] ,
                                                initialfile = self.file)
            if file:
                #print "open from file :",  file
                if self.command is not None:
                    self.command(self.object, actor, "DejaVu2Scenario",
                                 start = kf1, end = kf2,
                                 file = file ,functionName = functionName)
                    self.file = file
                    

