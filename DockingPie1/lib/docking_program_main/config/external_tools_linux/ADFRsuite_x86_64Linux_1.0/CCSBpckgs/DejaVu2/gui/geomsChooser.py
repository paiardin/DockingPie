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

#############################################################################
#
# Author: Alexandre T. GILLET
#
# TSRI 2003
#
#############################################################################
#
#
# $Header: /mnt/raid/services/cvs/DejaVu2/gui/geomsChooser.py,v 1.1.1.1.4.1 2017/07/13 22:20:53 annao Exp $
#
# $Id: geomsChooser.py,v 1.1.1.1.4.1 2017/07/13 22:20:53 annao Exp $
#
#
#
#

import string, math
import numpy.oldnumeric as Numeric
from types import IntType
import Tkinter
#from Tkinter import *
from mglutil.gui.InputForm.Tk.gui import InputFormDescr,InputForm,CallBackFunction
from mglutil.gui.BasicWidgets.Tk.customizedWidgets import ListChooser
from opengltk.OpenGL import GL
import Pmw
import types


class geomsChooser:

    def __init__(self,viewer):

        #ViewerGUI(viewer,8,2)
        # the viewe where the geoms are display (for us: pmv or DejaVu2)
        if viewer is None:
            return
        self.viewer = viewer
        geomDic = {}
        dpyList = None
        self.assignGeom =[geomDic,dpyList]

 
    def toggleExpansion(self, event):
        # get a 0-based index into list of names
	o = self.lc.lb.nearest(event.y)
        fullName = self.getFullName(o)
	#obj = self.objectByName(self.olist.get(o))
        obj = self.viewer.FindObjectByName(fullName)
	if obj:
            childGeoms = obj.AllObjects()
            if len(childGeoms)==1:  # this geoemtry has no children
                return
            else: # this geometry has children
                if obj.isExpandedInObjectList: self.collapse(obj)
                else: self.expand(obj)

    def expand(self, object):
        # object is a geometry
        if object.isExpandedInObjectList: return
        object.isExpandedInObjectList = 1
        geoms = object.children
        ind = self.objectIndex(object) + 1
        c = self.countParents(object) + 1
        prefix = '~'*c
        for i in range(len(geoms)):
            g = geoms[i]
            if g==object: continue
            if not g.listed: continue
            if len(g.vertexSet) == 0 and len(g.children)==0:continue
            self.lc.insert(ind, prefix + g.name)
            ind = ind + 1

            
    def collapse(self, object):
        # object is a geometry, we recursively collapse the sub-tree
        object.isExpandedInObjectList = 0

        # delete the names from the bject list widget
        nbChildren = self.countDecendentsInWidget(object)
        ind = self.objectIndex(object) + 1
        for i in range(ind, ind+nbChildren):
            self.lc.lb.delete(ind)
        # toggle isExpandedInObjectList for all descendents
        for child in object.AllObjects():
            if child.listed:
                child.isExpandedInObjectList = 0

    def countDecendentsInWidget(self, object):
        # object is a geometry, we count and return the number of 
        # decendents shown in widget
        ind = self.objectIndex(object)
        allNames = self.lc.lb.get(0,'end')
        nbTild = string.count(allNames[ind],'~')+1
        # count children in widget
        nbChildren = 0
        for i in range(ind+1, len(allNames)):
            nbt = string.count(allNames[i],'~')
            if nbt >= nbTild:
                nbChildren = nbChildren + 1
            else:
                break
        return nbChildren

        
    def getFullName(self, ind):
        # strip the leading ~
        allNames = self.lc.lb.get(0,'end')
        nbTild = string.count(allNames[ind],'~')
        fullName = allNames[ind][nbTild:]
        for i in range(ind-1, -1, -1):
            nbt, name = self.lstripChar(allNames[i], '~')
            if nbt >= nbTild: continue
            nbTild = nbt
            fullName = name + '|' + fullName
        return fullName

    def objectIndex(self, object):
        # object is a geometry and we find this object's index in the list of
        # names displayed in te widget. If the ibecjt is not shown we
        # return -1
        l = self.lc.lb.get(0,'end')
        for i in range(len(l)):
            indent, n = self.lstripChar(l[i], '~')
            if n==object.name: break
        if i==len(l): return -1
        else: return i

    def lstripChar(self, name, char):
	n = string.count(name,'~')
	return n, name[ n : ]

    def countParents(self, object):
	c = 0
	while object.parent:
	    c = c+1
	    object = object.parent
	return c

    def addObject(self, obj, parent):
        if not obj.listed: return
	if not parent:
	    self.lc.insert(0,obj.name)
	    self.lc.select(obj.name)
	else:
            if not parent.isExpandedInObjectList: return
                
            i = self.objectIndex(parent)
	    if i==-1: return
	    c = self.countParents(obj)
            prefix = '~'*c
	    name = prefix + obj.name
            # now we need to skip all children already there
            l = self.lc.lb.get(0,'end')
            while 1:
                i = i + 1
                if i==len(l): break
                if self.lc.get(i)[:c]!=prefix: break

	    self.lc.insert(i, name)

    def getFullNameList(self,name,list,dic):
        """ get the list of fullName from geomDic """

        for p in dic.keys():
            name1 = name +'|'+p
            if dic[p] != {}:
                list = self.getFullNameList(name1,list,dic[p])
            else:
                list.append(name1)
                
        return list
    

    def getDpyList(self,objects):
        """function to obtain a display list for any object in DejaVu2"""
        from opengltk.OpenGL import GL
        # Draw build a display function that contains global coloring,
        # transformation nd all global GL properties of this objects.
        # we cannot use the object's display list that does not contains this
        # info
        #if len(objects) == 1 and objects[0].dpyList:
        #    return objects[0].dpyList
        # transparent object need to be drawn last

        # reorder object so transparent object drawn last
        transparent =[]
        opaque = []
        for obj in objects:
            if obj.transparent: #obj.isTransparent():
                transparent.append(obj)
            else:
                opaque.append(obj)
        objlist = opaque + transparent
        
        lNewList = GL.glGenLists(1)
        #print "lNewList geomsChooser.getDpyList", lNewList, self.name
        dpyList = ( lNewList,
                   self.viewer.currentCamera.tk.call(self.viewer.currentCamera._w, 'contexttag')
                  )

        camera = self.mv.GUI.VIEWER.currentCamera
        camera.Activate()

        GL.glNewList(dpyList, GL.GL_COMPILE)
        #print "geomChooser 208", GL.glGetIntegerv(GL.GL_LIST_INDEX)
        for obj in objlist:
            if obj.immediateRendering: # or obj.hasChildWithImmediateRendering:
                camera.drawMode=5
                camera.Draw(obj)
            else:
                camera.drawMode=2
                camera.drawTransparentObjects = 0
                camera.hasTransparentObjects = 0
                # draw opaque object
                for m in obj.instanceMatricesFortran:
                    GL.glPushMatrix()
                    GL.glMultMatrixf(m)
                    if len(obj.children):
                        map( camera.Draw, obj.children)
                    else:
                        camera.Draw(obj)
                    GL.glPopMatrix()
                
                # draw transparent children of object
                if camera.hasTransparentObjects:
                    camera.drawTransparentObjects = 1
                    for m in obj.instanceMatricesFortran:
                        GL.glPushMatrix()
                        GL.glMultMatrixf(m)
                        map( camera.Draw, obj.children)
                        GL.glPopMatrix()
                # draw transparent object that do not have children
                if obj.isTransparent() and not len(obj.children):
                    camera.drawTransparentObjects =1
                    for m in obj.instanceMatricesFortran:
                        GL.glPushMatrix()
                        GL.glMultMatrixf(m)
                        camera.Draw(obj)
                        GL.glPopMatrix()
                    

        GL.glEndList()    
        return dpyList

    def buildForm(self):
        ifd = InputFormDescr(title=self.ftitle)
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
        
        ifd.append({'name':'add',
                    'widgetType':Tkinter.Button,
                    'tooltip':""" Add the selected geom to selected frame""",
                    'wcfg':{'text':'>>','command':self.add_cb},
                    'gridcfg':{'row':1,'column':1,'rowspan':1 }})
        
        ifd.append({'name':'remove',
                    'widgetType':Tkinter.Button,
                    'tooltip':""" remove the selected geom to selected """,
                    'wcfg':{'text':'<<','command':self.remove_cb},
                    'gridcfg':{'row':2,'column':1,'rowspan':1 }})
            
        ifd.append({'name':'toload',
                    'widgetType':ListChooser,
                    'tooltip':"""list of geom  the user chose to
                    apply to the pattern""",
                    'wcfg':{
            'mode':'extended',
            'lbwcfg':{'exportselection':0},
            'title':self.lc2title},
                    'gridcfg':{'sticky':'we', 
                               'row':0, 'column':2,'rowspan':3}})
        ifd.append({'name':'clear',
                    'widgetType':Tkinter.Button,
                    'tooltip':""" Clear entry """,
                    'wcfg':{'text':'Clear','width':10,
                            'command':self.clear_cb},
                    'gridcfg':{'sticky':'we','row':0, 'column':3}})
        ifd.append({'name':'remove',
                    'widgetType':Tkinter.Button,
                    'tooltip':""" remove the selected geom to selected """,
                    'wcfg':{'text':'<<','command':self.remove_cb},
                    'gridcfg':{'row':2,'column':1,'rowspan':1 }})
        return ifd
        
    def showForm(self,master,geomDic=None,func=None):
        """create formdescr for setGeoms
        geomDic = a dict of Geom 
        """
        # assign a func so we can log it in Pmv
        # func is a Pmv command
        if func:
            self.assign_func = func
            
        if geomDic:
            self.assignGeom[0]=geomDic
            
        if hasattr(self,'form'):
            if self.guiOn:
                return
            else:
                self.update_lc2()
                self.form.deiconify()
                val = self.form.go()
               
                if val:
                    geomDic =self.assignGeom[0] 
                    self.assign_func(geomDic)
                    self.form.withdraw()
        else:
            if not master:
                master = Tkinter.Toplevel()
                master.title('GeomChooser')
            else:
                master = master

            ifd = self.buildForm()
            self.form = InputForm(master,None,descr=ifd,
                                  scrolledFrame=0,
                                  cancelCfg={'command':self.cancel_cb})
                
            self.lc = self.form.descr.entryByName['availableGeom']['widget']
            self.lc2 = self.form.descr.entryByName['toload']['widget']
            self.addObject(self.viewer.rootObject,None)
            self.update_lc2()
            val = self.form.go()
            if val:
                geomDic =self.assignGeom[0] 
                self.assign_func(geomDic)
                self.form.withdraw()
                
    def cancel_cb(self,event = None, func=None):
        """ close setup animation form without setting anything"""
        if hasattr(self,'geomDicCopy'):
            self.patDic = self.patDicCopy.copy()
        
    def add_cb(self):
        lgeom=[]
        geomDic = self.assignGeom[0]
        # get geom name to load
        o = map(int,self.lc.lb.curselection())
        for Ind in o:
            fullName = self.getFullName(Ind)
            obj = self.viewer.FindObjectByName(fullName)
            # strip the root
            l = string.split(fullName,'|')[1:]
            for i in range(len(l)):
                if i==0:
                    if not geomDic.has_key(l[i]):
                        dic = geomDic[l[i]]={}
                    else:
                        dic = geomDic[l[i]]
                else:
                    if not dic.has_key(l[i]):
                        dic[l[i]]={}
                        
                    dic = dic[l[i]]
            
        self.update_lc2()

    def insertValue(self,ind,prefix,dic):
        """ insert Value of a dictionnary in the listchooser2"""
        for key in dic.keys():
            if dic[key] == {}:
                ind = ind+1
                self.lc2.insert(ind, prefix+key)
            else:
                ind =ind +1
                self.lc2.insert(ind, prefix+key)
                p = prefix + '~'
                ind=self.insertValue(ind,p,dic[key])
        return ind
        
    def update_lc2(self):
        """ update listchooser entries display"""

        # save current entry selected
        sel  = self.lc2.get()
        prefix = '~'
        self.lc2.clear()
        ind = 0
        geomDic = self.assignGeom[0]
        ind = self.insertValue(ind,prefix,geomDic)

        # select the entry save as sel
        # keep the entry selected after the update
        for i in sel:
            self.lc2.select((i,))

    def findParentMol(self,index):
        """ find the parent value name from a geom ~~
        go up to find first value with onley ~ """
        parent = None
        for i in range(index+1):
            parent = self.lc2.entries[index -i][0]
            if parent[0] == '~' and parent[1] !='~':
                parent = parent[1:]
                return parent
        return  None
         
    def remove_cb(self):
        """ remove entry in litschooser (geomtoload)"""
        geomDic = self.assignGeom[0]
        selindex = self.lc2.getInd()
        for index in selindex:
            value = self.lc2.entries[index][0]
           
            if value[1] == '~':
                mol = self.findParentMol(index)
                rindex = geomDic[mol].index(value[2:])
                del geomDic[mol][rindex]
            else:
                mol = value[1:]
                del geomDic[mol]
        self.update_lc2()

    def clear_cb(self):
        """ clear all entry """
        self.lc2.clear()
        if  self.assignGeom:
            self.assignGeom[0] = {} # geomDic
            self.assignGeom[1] = None # dpyList
