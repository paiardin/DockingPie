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
import os
import weakref
from mglutil.gui.InputForm.Tk.gui import InputFormDescr,InputForm,evalString
from mglutil.gui.BasicWidgets.Tk.thumbwheel import ThumbWheel
from mglutil.util.packageFilePath import findFilePath
import tkMessageBox
import tkFileDialog
from SimpleDialog import SimpleDialog
from os import path

class Recorder(Tkinter.Frame):
    """Implements GUI for recording movie clips( *.mpg)"""

    def __init__(self, master=None,
                 height=80,width=100, title="Video recorder",
                 icondir= None, filetypes=[("MPG", ".mpg")],
                 fileName = None, camera = None, 
                 gui = True):
        self.ifd = None
        self.form = None
        self.paramForm = None
        self.ifd2 = None
        self.fileName = fileName
        
        self.autoPauseDelay = 1 # auto pause after 1 second
        self.pauseLength = 30
        self.cameraw = 0
        self.camerah = 0
        self.camera = weakref.ref(camera)
        if self.camera():
            self.cameraw = camera.width
            self.camerah = camera.height
        
        if not icondir:
            icondir = findFilePath('icons', 'mglutil.gui.BasicWidgets.Tk')
        self.stopIcon = Tkinter.PhotoImage(file=os.path.join(icondir, "stop3.gif"), master=master)
        self.pauseIcon = Tkinter.PhotoImage(file=os.path.join(icondir, "stop2.gif"), master=master)
        self.recordIcon = Tkinter.PhotoImage(file=os.path.join(icondir, "record.gif"), master=master)
        self.record1Icon = Tkinter.PhotoImage(file=os.path.join(icondir, "record1.gif"), master=master)
        self.chmodIcon =  Tkinter.PhotoImage(file=os.path.join(icondir, 'chmod.gif'), master=master)
        self.closeIcon = Tkinter.PhotoImage(file=os.path.join(icondir, 'close.gif'), master=master)
        self.fileTypes = filetypes
        if gui:
            self.RecVar =  Tkinter.IntVar()
            self.RecVar.set(0)
            self.PauseVar =  Tkinter.IntVar()
            self.PauseVar.set(0)
            form = self.buildForm(master=master,
                           title=title,height=height,
                           width=width)
            form.deiconify()
            
            
    def buildForm(self, master=None,  width=100, height=80,
                  title=None):
        if self.form:
            self.form.deiconify()
            self.createKeyBindings()
            return
        self.master = master
        ifd = self.ifd = InputFormDescr(title=title)
        ifd.append({'name':'fileopen',
                    'widgetType':Tkinter.Button,
                    'tooltip': "Opens file browser",
                    'wcfg':{'text':"Save As:",
                            'command': self.browseFile,
                            'width': 0, 'height': 0,},
                    'gridcfg':{'sticky':'w', 'column':0}
                    })
        ifd.append({'name':'filename',
                    'widgetType':Pmw.EntryField,
                    'tooltip': "type filename",
                    'gridcfg':{'sticky':'w',
                               #'columnspan': 3},
                               'columnspan': 2, 'row': -1},
                    'wcfg':{'command':self.getFileName,
                            #'label_text':'Save As',
                            'entry_width':12,
                            'value':self.fileName,
                            #'labelpos':'w'}})
                            }})
        
        ## ifd.append({'name':'fileopen',
##                     'widgetType':Tkinter.Button,
##                     'tooltip': "Open file browser",
##                     'wcfg':{'text':"...",
##                             'command': self.browseFile,
##                             'width': 0, 'height': 0,},
##                     'gridcfg':{'sticky':'w', 'column':2, 'row':-1}
##                     })
        ifd.append({'name': 'recordB',
                    'widgetType': Tkinter.Checkbutton,
                    'tooltip':'start/stop recording',
                    'wcfg': {'variable': self.RecVar,
                             'bd':2,
                             'image':self.record1Icon,
                             'width':self.record1Icon.width(),
                             'height':self.record1Icon.height(),
                             'indicatoron':0,
                             },
                    'gridcfg':{'sticky':'nesw','row':-1},
                    'command':self.record_cb
                    })
        
        ifd.append({'name': 'pauseB',
                    'widgetType': Tkinter.Checkbutton,
                    'tooltip':'pause/start recording',
                    'wcfg':{'variable': self.PauseVar,
                            'bd':2,
                            'image':self.pauseIcon,
                            'width':self.pauseIcon.width(),
                            'height':self.pauseIcon.height(),
                            'indicatoron':0,
                            },
                    'gridcfg':{'sticky':'nesw', 'row':-1},
                    'command':self.pause_cb
                    })
        
        ifd.append({'name': 'stopB',
                    'widgetType': Tkinter.Button,
                    'tooltip':'stop recording',
                    'wcfg':{'bd':2,
                            'image':self.stopIcon,
                            'width':self.stopIcon.width(),
                            'height':self.stopIcon.height(),
                            
                            },
                    'gridcfg':{'sticky':'nesw', 'row':-1},
                    'command':self.stop_cb
                    })
        
##         ifd.append({'name': 'recordB',
##                     'widgetType': Tkinter.Button,
##                     'tooltip':'start/pause recording',
##                     'wcfg':{'bd':4,
##                             'image':self.recordIcon,
##                             'width':self.recordIcon.width(),
##                             'height':self.recordIcon.height()
##                             },
##                     'gridcfg':{'sticky':'nesw','row':-1},
##                     'command':self.record_cb})
        
##         ifd.append({'name': 'pauseB',
##                     'widgetType': Tkinter.Button,
##                     'tooltip':'pause recording',
##                     'wcfg':{'bd':4,
##                             'image':self.pauseIcon,
##                             'width':self.pauseIcon.width(),
##                             'height':self.pauseIcon.height()
##                             },
##                     'gridcfg':{'sticky':'nesw', 'row':-1},
##                     'command':self.pause_cb})
        
        
##         ifd.append({'name': 'modeB',
##                     'widgetType': Tkinter.Button,
##                     'text':'Change Mode',
##                     'tooltip':'opens panel to change video parameters',
##                     'wcfg':{'bd':4,
##                             'image':self.chmodIcon,
##                             'width':self.chmodIcon.width(),
##                             'height':self.chmodIcon.height()
##                             },
##                     'gridcfg':{'sticky':'nesw','row':-1},
##                     'command':self.setparams_cb })
        
##         ifd.append({'name': 'closeB',
##             'widgetType': Tkinter.Button,
##             'text':'Close',
##             'tooltip':'closes video recorder',
##             'wcfg':{'bd':2,
##                     'image':self.closeIcon,
##                     'width':self.closeIcon.width(),
##                     'height':self.closeIcon.height(),
##                     },
##             'gridcfg':{'sticky':'nesw','row':-1},
##             'command':self.close_cb})
        form = self.form = InputForm(self.master, None, descr = ifd,
                                     modal = 0, blocking = 0,
                                     closeWithWindow=1, onDestroy = self.close_cb)
        self.createKeyBindings()
        return form
    

    def getFileName(self):
        """ Get file name from the input form's entry field."""
        if self.ifd:
            name = self.ifd.entryByName['filename']['widget'].get()
            if name:
                if self.fileName != name:
                    self.fileName = name
                
                
    def browseFile(self):
        """Opens file browser."""
        
        fileDir = None
        if self.fileName:
            if path.exists(self.fileName):
                fileDir = path.dirname(path.realpath(self.fileName))
            else:
                fileDir = os.getcwd()
        file = tkFileDialog.asksaveasfilename( filetypes=self.fileTypes,
                                               initialdir=fileDir,
                                               initialfile=None,
                                               title="Save file")
        if file:
            if self.fileName != file:
                self.fileName = file
                self.ifd.entryByName['filename']['widget'].setentry(file)
     
    def record_cb(self):
        
        #get the value of the Record checkbutton
        val = self.RecVar.get()
        if val: #checked
            self.getFileName()
            res = 1
            if self.camera:
                if self.camera().videoRecordingStatus == "stopped":
                    #kw = {'filename':self.fileName, 'width':self.cameraw, 'hight':self.camerah,
                    #      'autoPauseDelay': self.autoPauseDelay,'pauseLength':self.pauseLength}
                    #kw = {'filename':self.fileName }
                    #apply(self.setVideoParams, (), kw)
                    if os.path.exists(self.fileName):
                        # ask the user if he wants to overwrite the existing file
                        from Dialog import Dialog
                        d = Dialog(None, {'title': 'File exisits',
                                          'text':
                                          'File "%s" already exists.'
                                          ' Do you want to overwrite it ?'%self.fileName,
                                          'bitmap': 'warning',
                                          'default': 1,
                                          'strings': ('Yes', 'No')})
                        ans = d.num
                        if ans == 1: #no
                            self.RecVar.set(0)
                            return
                    res = self.setVideoParams(filename = self.fileName)
                if res:
                    #print "record_cb: recording"
                    self.camera().start()
            if self.ifd:
                b = self.ifd.entryByName['recordB']['widget']	
                if res:
                    b.config(image=self.recordIcon)
                    pauseval = self.PauseVar.get()
                    if pauseval:
                        self.PauseVar.set(0) #uncheck the Pause checkbutton
        else: #unchecked
            if self.camera:
                #print "record_cb: stop recording"
                self.stop_cb()
                
    def stop_cb (self):
        #print "stop_cb"
        if self.camera:
            if not self.camera().videoRecordingStatus == 'stopped':
                #print "stop recording"
                self.camera().stop()
	    if self.ifd:
                recval = self.RecVar.get()
                if recval:
                   self.RecVar.set(0) #uncheck the Record checkbutton
                pauseval = self.PauseVar.get()
                if pauseval:
                    self.PauseVar.set(0) #uncheck the Pause checkbutton
                b = self.ifd.entryByName['recordB']['widget']
                b.config(image=self.record1Icon)
                
        
    def pause_cb(self):
        #print "pause_cb"
        val = self.PauseVar.get()
        if self.camera:
            camerastatus = self.camera().videoRecordingStatus
            if val: # the Pause button is checked
                if camerastatus == 'recording':
                    #print "pause recording"
                    self.camera().pause()
                    #if self.ifd:
                        #b = self.ifd.entryByName['recordB']['widget']
                        #b.config(command = self.record_cb, image=self.recordIcon)
            else: # the Pause button is unchecked
                if camerastatus == 'paused':
                    #print "start recording after pause"
                    self.record_cb()

    def spacePress_cb(self, event):
        if self.camera:
           camerastatus = self.camera().videoRecordingStatus
           if camerastatus == "stopped" or camerastatus == "paused":
               self.RecVar.set(1)
               self.record_cb()
           elif camerastatus == "recording":
               self.PauseVar.set(1)
               self.pause_cb()
            
    def setparams_cb(self):
        """Opens a panel to set Video Parameters"""
        if self.paramForm:
            self.paramForm.deiconify()
            return
        self.ifd2 = ifd2 = InputFormDescr(title = "Set video options")
        ifd2.append({'widgetType':Pmw.EntryField,
                     'tooltip': 'Set camera width',
                     'name':'cameraw',
                     'gridcfg':{'sticky':'w', 'columnspan':2},
                     'wcfg':{#'command': self.setCameraWidth_cb,
                             'label_text':'width:',
                             'entry_width':10,
                             'validate': {'validator':'real', 'min': 0},
                             'value':str(self.cameraw),
                             'labelpos':'w'}})
        
        ifd2.append({'widgetType':Pmw.EntryField,
                     'name':'camerah',
                     'tooltip':'Set camera height',
                     'gridcfg':{'sticky':'w' , 'columnspan':2},
                     'wcfg':{#'command': self.setCameraHeight_cb,
                             'label_text':'height',
                             'entry_width':10,
                             'validate': {'validator':'real', 'min': 0},
                             'value':str(self.camerah),
                             'labelpos':'w'}})
        ifd2.append({'name': 'autoPause',
                     'wtype':ThumbWheel,
                     'widgetType':ThumbWheel,
                     'tooltip':'set auto pause delay (seconds)',
                     'wcfg':{'labCfg':{'fg':'black', 'side':'left', 'text':'AutoPause Delay'},
                             'showLabel':1, 'width':100,
                             'min':0,
                             'value':self.autoPauseDelay,
                             'oneTurn':100,
                             'type':'int',
                             'increment':1,
                             #'callback':self.setAutoPauseDelay_cb,
                             'canvascfg':{'bg':'red'},
                             'continuous':0, 'wheelPad':1, 'height':15},
                     'gridcfg':{'sticky':'nesw', 'columnspan':2}})
        ifd2.append({'name': 'pauseLength',
                     'wtype':ThumbWheel,
                     'widgetType':ThumbWheel,
                     'tooltip':'set number of frames to be added when\nrecording resumes after autopause',
                     'wcfg':{'labCfg':{'fg':'black', 'side':'left', 'text':'AutoPause Length'},
                             'showLabel':1, 'width':100,
                             'min':0,
                             'value':self.pauseLength,
                             'oneTurn':100,
                             'type':'int',
                             'increment':1,
                             #'callback':self.setPauseLength_cb,
                             'canvascfg':{'bg':'red'},
                             'continuous':0, 'wheelPad':1, 'height':15},
                     'gridcfg':{'sticky':'nesw','columnspan':2}})
        ifd2.append({'name':'okB',
                     'widgetType': Tkinter.Button,
                     'wcfg':{'text': 'Apply',
                             'command': self.apply_cb,},
                     'gridcfg':{'sticky':'nesw'}})
        ifd2.append({'name':'cancelB',
                     'widgetType': Tkinter.Button,
                     'wcfg':{'text': 'Cancel',
                             'command': self.cancel_cb,},
                     'gridcfg':{'sticky':'nesw', 'row':-1, 'column':1}})
        self.paramForm =  InputForm(self.master, None,
                                    descr = ifd2,
                                    modal = 0, blocking = 0)
        self.paramForm.deiconify()
        
        return self.paramForm 
        
        
    def setCameraWidth_cb(self):
        """ 'Width' entry field callback of the video parameters form ."""  
        w = int(self.ifd2.entryByName['cameraw']['widget'].get())
        if w != self.cameraw:
            self.cameraw = w

    def setCameraHeight_cb(self):
        """'Height' entry field callback of the video parameters form ."""
        h = int(self.ifd2.entryByName['camerah']['widget'].get())
        if h != self.camerah:
            self.camerah = h

    def setAutoPauseDelay_cb(self, val):
        """Callback of the autoPause thumbwheel widget (video parameters input form). """
        #print val
        d = int(val)
        if d != self.autoPauseDelay:
            self.autoPauseDelay = d

    def setPauseLength_cb(self, val):
        """Callback of the pauseLength thumbwheel widget (video parameters input form). """
        #print val
        l = int(val)
        if l != self.pauseLength:
           self.pauseLength = l
           
    def apply_cb(self):
        """Apply button callback of the video parameters input form."""
        kw = {}
        ebn = self.ifd2.entryByName
        w = float(ebn['cameraw']['widget'].get())
        if w != self.cameraw:
            self.cameraw = w
        kw ['width']= w
        h = float(ebn['camerah']['widget'].get())
        if h != self.camerah:
            self.camerah = h
        kw['height']=h
            
        d = ebn['autoPause']['widget'].get()
        if d != self.autoPauseDelay:
            self.autoPauseDelay = d
        kw["autoPauseDelay"]=d

        l = ebn['pauseLength']['widget'].get()
        if l != self.pauseLength:
            self.pauseLength = l
        kw['pauseLength']=l

        #apply(self.setVideoParams, (), kw)
        self.cancel_cb()
        

    def setVideoParams(self, **kw):
        """ Uses Recordable Camera methods to set file name and other video parameters"""
        #print kw
        if self.camera:
            filename = kw.get('filename')
            if filename:
	        try:	      
                    self.camera().setVideoOutputFile(filename)
   		except IOError:
		    button = tkMessageBox.showerror("IOError", message="Could not open file %s\nfor writing a movie."%filename, parent = self.master)
                    return False
            ## params = {}
##             w = kw.get('width')
##             if w is not None:
##                 params['width'] = w
##             h = kw.get('height')
##             if h is not None:
##                 params['height']= h
##             pl = kw.get('pauseLength')
##             if pl is not None:
##                 params['pauseLength'] = pl
##             apd = kw.get('autoPauseDelay')
##             if apd is not None:
##                 params['autoPauseDelay'] = apd
##             if len(params):
##                 apply(self.camera().setVideoParameters, (), params)
                
            #use default parameters for now
            self.camera().setVideoParameters()
	    return True
        
    def createKeyBindings(self):
        self.master.bind('<KeyPress-space>', self.spacePress_cb)
        
    def removeKeyBindings(self):
        self.master.unbind('<KeyPress-space>')
    
    def cancel_cb(self):
        self.paramForm.withdraw()

    def close_cb(self, event=None):
        if hasattr(self,'form'):
            self.form.withdraw()
        self.removeKeyBindings()


