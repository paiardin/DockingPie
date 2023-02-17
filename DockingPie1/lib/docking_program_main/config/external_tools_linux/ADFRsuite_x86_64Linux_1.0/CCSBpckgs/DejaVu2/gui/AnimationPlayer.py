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
# Author: Alexandre T. GILLET
#
# TSRI 2003
#
#############################################################################
#
#
# $Header: /mnt/raid/services/cvs/DejaVu2/gui/AnimationPlayer.py,v 1.1.1.1.4.1 2017/07/13 22:20:53 annao Exp $
#
# $Id: AnimationPlayer.py,v 1.1.1.1.4.1 2017/07/13 22:20:53 annao Exp $
#
#
#
#
import string, math,time
from mglutil.gui.BasicWidgets.Tk.player import Player
from mglutil.gui.InputForm.Tk.gui import InputFormDescr,InputForm,CallBackFunction
from mglutil.gui.BasicWidgets.Tk.customizedWidgets import ListChooser
from DejaVu2.gui.geomsChooser import geomsChooser
import Pmw,Tkinter
import types

class AnimPlayer(Player,geomsChooser):
    """ provide a player to display/undisplay DejaVu2 object
    in a cylce making a animation.

    viewer: DejaVu2 viewer
    height,width are the property of the player gui
    startFrame: frame to begin
    endFrame: last frame to display
    stepSize: step between new frame
    playMode:                  #playMode options:
                               #   0   play once and stop
                               #   1   play continuously in 1 direction
                               #   2   play once in 2 directions
                               #   3   play continuously in 2 directions
    titleStr: title of the gui windows
    counter: gui with a counter on or not
    gui: #1 show the player gui
         #0 player gui not display
    framerate =15.             # number of frame per second to be display
    """


    def __init__(self,viewer, master=None, root=None,
                 height=80,width=200,
                 currentFrameIndex=0,
                 startFrame=0, 
                 endFrame=0,
                 maxFrame=0,
                 stepSize=1, 
                 playMode=1,
                 titleStr='AnimPlayer',
                 counter = 1,
                 gui=1,framerate=15.):

        if master is None:
            master = Tkinter.Toplevel()

        # viewer from which we get the geom to animate
        self.viewer = viewer
        # frame list, is a list of geom .
        # each frame is made of n DejaVu2 object to be display
        # it a list of tuple, each tuple is a {} and []
        # {} is: {'objname':['geom1','geom2'],'objname':['geom1','geom2']}
        # [] is list of DejaVu2 object per frame
        self.framelist=[]

        Player.__init__(self, master=master,root=root,height=height,
                        width=width,currentFrameIndex=currentFrameIndex,
                        startFrame=startFrame,endFrame=endFrame,
                        maxFrame=maxFrame,stepSize=stepSize,
                        playMode=playMode,titleStr=titleStr,counter=counter,
                        gui=gui,framerate=framerate)

    def nextFrame(self, id):
        """ id: index number of the frame to display.
        undisplay object from previous frame, display object
        of new frame.
        """
        if not self.framelist:return
        id = int(id)
        if id >=len(self.framelist):return

        i = self.currentFrameIndex
        ## undisplay previous object
        frame = self.framelist[i]
        for obj in frame[1]:
            obj.Set(visible=0,redo=0)
            
        if self.hasCounter and self.gui:
            self.form.ent2.delete(0,'end')
            self.form.ent2.insert(0, str(id))

        ## display new frame
        currentframe = self.framelist[id]
        for obj in currentframe[1]:
            obj.Set(visible=1,redo=0)
        
        self.viewer.deleteOpenglList()
	self.viewer.Redraw()
        self.currentFrameIndex = id

#####################################################################
#### Animation Setting the frame

    def SetAnim_cb(self):
        self.showForm()
        
    def showForm(self):
        """create formdescr for setAnim
        form to set the animation:
        each frame is a list of geom to display
        """

        entryFrame = []
        for i in range(len(self.framelist)):
            val = 'frame_'+str(i)
            entryFrame.append((val,))
        
        if not hasattr(self,'form_Setanim'):
            ifd = InputFormDescr(title="SetAnimation")
            ifd.append({'widgetType':ListChooser,
                        'name':'availableGeom',
                        'tooltip':'geom available in viewer',
                        'wcfg':{'mode':'extended',
                                'lbwcfg':{'exportselection':1},
                                'command':[(self.toggleExpansion,'<Double-Button-1>')],
                                'commandEvent':None,
                                'title':'availableGeom'},
                        
                        'gridcfg':{'row':0,'column':0,
                                   'sticky':'wens','rowspan':3}})
            ifd.append({'name':'newframe',
                        'widgetType':Tkinter.Button,
                        'tooltip':""" Add an empty frame to the animation""",
                        'wcfg':{'text':'NewFrame','command':self.addframe_cb},
                        'gridcfg':{'row':0,'column':1,'rowspan':1 }})
            
            ifd.append({'name':'add',
                        'widgetType':Tkinter.Button,
                        'tooltip':""" Add the selected geom to selected frame""",
                        'wcfg':{'text':'AddGeom','command':self.add_cb},
                        'gridcfg':{'row':1,'column':1,'rowspan':1 }})
            
            ifd.append({'name':'geomtoload',
                        'widgetType':ListChooser,
                        'tooltip':"""list of frame  the user chose to
                        apply to the pattern""",
                        'wcfg':{'entries':entryFrame,
                                'mode':'extended',
                                'lbwcfg':{'exportselection':0},
                                'title':'Frame(geom) to be display'},
                        'gridcfg':{'sticky':'we', 
                                   'row':0, 'column':2,'rowspan':3}})
            ifd.append({'name':'remove',
                        'widgetType':Tkinter.Button,
                        'tooltip':""" Remove the selected entry from the
                        commands to be applied to the object when loaded in the application""",
                        'wcfg':{'text':'REMOVE','width':10,
                                'command':self.remove_cb},
                        'gridcfg':{'sticky':'we','row':0, 'column':3}})
            
            ifd.append({'name':'oneup',
                        'widgetType':Tkinter.Button,
                        'tooltip':"""Move the selected entry up one entry""",
                        'wcfg':{'text':'Move up','width':10,
                                'command':self.moveup_cb},
                        'gridcfg':{'sticky':'we','row':1,'column':3}})
            
            ifd.append({'name':'onedown',
                        'widgetType':Tkinter.Button,
                        'tooltip':"""Move the selected entry down one entry""",
                        'wcfg':{'text':'Move down','width':10,
                                'command':self.movedown_cb},
                        'gridcfg':{'sticky':'we','row':2,'column':3}})

            self.form_Setanim = InputForm(self.master,None,descr=ifd,
                                  scrolledFrame=0)
            
            
            self.lc = self.form_Setanim.descr.entryByName['availableGeom']['widget']
            self.lc2 = self.form_Setanim.descr.entryByName['geomtoload']['widget']
            self.addObject(self.viewer.rootObject,None)
            val = self.form_Setanim.go()
            if val:
                self.assign()
                self.form_Setanim.withdraw()
        else:
            self.form_Setanim.deiconify()
            val = self.form_Setanim.go()
            if val:
                self.assign()
                self.form_Setanim.withdraw()
                
    def deleteAnim(self):
        """ Delete all frame from the player list."""

        self.Stop_cb()
        self.startFrame=0
        self.endFrame = -1
        self.maxFrame = 0
        self.targetFrame = self.endFrame
        self.target = self.endFrame
        self.currentFrameIndex =0


    def assign(self,list=None):
        """ assign the animation.
        framelist ; list of {}
        return a list of obj for each frame
        """
        self.deleteAnim()
        AllObjects = self.viewer.rootObject.AllObjects()
        if list :
            self.framelist=[]
            for f in list:
                self.framelist.append((f,[]))
        print "assign"
        print "self.framelist",self.framelist
        for frame in self.framelist:
            geomDic = frame[0]
            Namelist=[]
            for parents in geomDic.keys():
                parentlist = filter(lambda x, name=parents:x.name == name,AllObjects)
                obj = parentlist[0]
                childrens = geomDic[obj.name]
                if len(childrens) > 0:
                    for namegeom in childrens:
                        child = filter(lambda x,
                                       name=namegeom:x.name == name,
                                       obj.children)[0]
                        if child.children != []:
                            for c in child.children:
                                frame[1].append(c)
                        else:
                            frame[1].append(child)
                else:
                    frame[1].append(obj)

            self.endFrame = self.endFrame +1
            self.maxFrame = self.maxFrame +1
            self.targetFrame = self.endFrame
            self.target = self.endFrame


    def movedown_cb(self):
        """ move entry one down """
        sel = self.lc2.get()
        index = self.lc2.getInd()
        if not sel: return
        sel = sel[0]
        if string.find(sel,"frame") < 0: return
        # get frame Index
        frameInd = int(string.split(sel,'_')[1])
        currentframe = self.framelist[frameInd]
        if (frameInd + 1) >= len(self.framelist):return
        # select next Frame.
        nextframe = self.framelist[frameInd+1]
        # move current frame one down in list
        self.framelist[frameInd] =nextframe
        self.framelist[frameInd+1] =currentframe
        self.updatelistchooser(self.framelist)

    def moveup_cb(self):
        """ move entry one up """
        sel = self.lc2.get()
        index = self.lc2.getInd()
        if not sel: return
        sel = sel[0]
        if string.find(sel,"frame") < 0: return
        # get frame Index
        frameInd = int(string.split(sel,'_')[1])
        currentframe=self.framelist[frameInd]
        if (frameInd - 1) <  0 :return
        # select previous Frame
        prevframe = self.framelist[frameInd-1]
        # move current frame one up in list
        self.framelist[frameInd] =prevframe
        self.framelist[frameInd-1] =currentframe
        self.updatelistchooser(self.framelist)
        
    def addframe_cb(self):
        """ add a new frame entry"""
        frame = ({},[])
        nb = len(self.framelist)
        #frame.number = nbframe 
        value = 'frame_'+ str(nb)
        self.framelist.append(frame)
        self.lc2.deselect(0,last='end')
        self.lc2.insert('end',value)
        self.lc2.select('end')

    def add_cb(self):
        listgeom=[]
        # get the frame index to add new geom
        if len(self.lc2.entries)>0:
            ind = self.lc2.getInd()
            if len(ind) ==0:return # no entry selected so create a new frame
            # get geom name to load
            o = map(int,self.lc.lb.curselection())
            for Ind in o:
                    fullName = self.getFullName(Ind)
                    listgeom.append(fullName)
            # get frame instance
            for index in ind:
                value = self.lc2.entries[index][0]
                if string.find(value,"frame") < 0: return
                frameInd = int(string.split(value,'_')[1])
                frameGeom = self.framelist[frameInd][0]
                for geom in listgeom:
                    # strip the root
                    l = string.split(geom,'|')
                    if not frameGeom.has_key(l[1]):
                        frameGeom[l[1]]=[]
                    for i in l[2:]:
                        if not( i in frameGeom[l[1]]):
                            frameGeom[l[1]].append(i)

            self.updatelistchooser(self.framelist)
        else:
            return
        
    def updatelistchooser(self,entry):
        """ update what is display in the list chooser """
        # save current selection
        sel  = self.lc2.get()
        # remove current entry
        self.lc2.clear()
        prefix = '~'
        prefix2 = '~~'
        nb= 0
        for e in entry:
            v = 'frame_'+str(nb)
            nb = nb+1
            # listchooser entry has to be a tuple(value,comment)
            self.lc2.add((v,))
            for mol in e[0].keys():
                self.lc2.add((prefix+mol,))
                for geom in e[0][mol]:
                    self.lc2.add((prefix2+geom,))
        # select the entry save as sel
        # keep the entry selected after the update
        for i in sel:
            self.lc2.select((i,))


    def findFramefromGeomInd(self,index):
        """ return the frame number from  which the geom was selected """
        for i in range(index+1):
            value = self.lc2.entries[index-i][0]
            if string.find(value,"frame") >= 0:
                frameInd =  int(string.split(value,'_')[1])
                lcInd = index -i
                return (frameInd,lcInd)
        return 0
    
    def findParentGeom(self,index):
        """ find the parent value name from a geom ~~
        go up to find first value with onley ~ """
        parent = None
        for i in range(index+1):
            parent = self.lc2.entries[index -i][0]
            if parent[0] == '~' and parent[1] !='~':
                parent = parent[1:]
                return parent
        return  None
    
    def removeEntryGeom(self,value,index):
        """ the geom entry in the listchooser"""
        # first find from which frame we need to remove geom
        val = self.findFramefromGeomInd(index)
        if val:
            frameInd =val[0]
            lcInd = val[1]
            if value[1] == '~':
                parent = self.findParentGeom(index)
                if parent:
                    listgeom = self.framelist[frameInd][0][parent]
                    rindex = listgeom.index(value[2:])
                    del listgeom[rindex]
                
            elif self.framelist[frameInd][0].has_key(value[1:]):
                del self.framelist[frameInd][0][value[1:]]
            
            self.updatelistchooser(self.framelist)
            
    def removeEntryFrame(self,value,index):
        """ remove the frame entry in the listchooser"""
        # delete frame from framelist
        frameInd = int(string.split(value,'_')[1])
        del self.framelist[frameInd]
        # Update list chooser
        self.updatelistchooser(self.framelist)
        
    def remove_cb(self):
        """ remove entry in litschooser (geomtoload)"""
        selindex = self.lc2.getInd()
        for index in selindex:
            value = self.lc2.entries[index][0]
            if string.find(value,"frame") >= 0:
                self.removeEntryFrame(value,index)
            elif value[0] == '~':
                self.removeEntryGeom(value,index)
                
################## Animation setup end #################################
########################################################################

    
