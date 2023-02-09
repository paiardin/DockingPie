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
# Date: Novembre 2005 Authors: Guillaume Vareille, Michel Sanner
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
# $Header: /mnt/raid/services/cvs/DejaVu2/VisionInterface/GeometryNodes.py,v 1.1.1.1.4.1 2017/07/13 22:20:08 annao Exp $
#
# $Id: GeometryNodes.py,v 1.1.1.1.4.1 2017/07/13 22:20:08 annao Exp $
#

# third party packages 
import os
import sys
from warnings import warn
from weakref import ref 
from numpy.oldnumeric import array, identity
from Tkinter import Menu
from os.path import splitext

# TSRI packages
from mglutil.util.callback import CallBackFunction
from NetworkEditor.items import NetworkNode
from NetworkEditor.macros import MacroInputNode, MacroNetwork, \
                        MacroOutputNode, MacroOutputNode, MacroNode

from DejaVu2.Common2d3dObject import Common2d3dObject
from DejaVu2.Geom import Geom
import DejaVu2


class GeometryNode(NetworkNode):
    """Base Class to hold a single geometry, should be only inherited,
never instanciated. it handles naming and parenting.

The geometry nodes handle a list of geometries. When the name changes, 
a new geom is generated. If an empty name is provided, the current geometry 
is deleted.

The parent and children connections reflect the state of the geometry selected 
in the combo box. If the parent or children geoms are currently selected in 
their own node: the connections are present. If the parent or a child is not 
visible (for the selected geometry) the parent port or/and the output port 
contour lines are white (instead of black). Before parenting you can set the 
behaviour of the 'parent' inputport. Right clicking on the port allows you to 
extend the parenting to the sibling geometries present in the node or to all 
the geometries of the node. This behaviour set by right clicking will only 
affect the next parenting. In the same way, deleting the parent connection, 
reparent to 'root' the selected geometry or the sibling geometries or all the 
geometries in the child node (depending on the right click menu choice).

THINGS TO DO TO INHERIT 
(the best is to follow what is done in class IndexedPolygonsNE bellow)
1- inherit :
        class IndexedPolygonsNE(GeometryNode):
2- at the beginning of init, call original init 
   (the 'indexedPolygons' in quotes will be the name of the output port):
        kw['name'] = name
        apply( GeometryNode.__init__, (self, 'indexedPolygons'), kw ) 
3- rearrange ports order after appending additionnal input and output:
        self.rearrangePorts()
4- at the beginning of doit, call original doit:
        GeometryNode.doit(self, name, instanceMatrices, geomOptions, parent)
5- the doit function should really immitate the class IndexedPolygonsNE example 
   (parameters linked to input ports from the base class must be passed as last)
        def doit(self, coords, indices, vnormals=None, colors=None,
        name=None, instanceMatrices=None, geomOptions=None, parent=None):
6- override appendGeometry()
        def appendGeometry(self, name):
            self.removePreviousGeomWithSameName(name)
            from DejaVu2.IndexedPolygons import IndexedPolygons
            self.geoms.append(IndexedPolygons(name))
            return 1 # num of geoms appended
7- don't forget to use the new datatypes
        ip.append(datatype='coord3(0,3)', name='coords')

GeomsFromFile bellow shows how the filename can be use as a director for
the name. It shows also how to append several geoms from one file.

If you want to read a new type of file the best is to add it inside
the appendGeometry of GeomsFromFile
"""

    def __init__(self, geomOutputPortName, **kw):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
        apply( NetworkNode.__init__, (self,), kw )

        self.geometriesOwner = True # tell if the node has created the geometries
                                    # or is just displaying geometries
        self.selectedGeomIndex = None 
        self.geoms = []

        # creat ports common to all geoemtry producing nodes
        ip = self.inputPortsDescr
        
        # for backward compatibility with old networks
        # name and instanceMatrices must be the first ones in this Base class
        ip.append(datatype='string', required=False, name='name')
        ip.append(datatype='string', required=False, name='geoms')
        ip.append(datatype='instancemat(0)', required=False,
                  name='instanceMatrices')
        ip.append(datatype='dict', required=False, name='geomOptions')
        # even if it is not compulsory, 
        # we try to keep the parent port in the last position
        ip.append(datatype='geomOrInsert2d', required=False, name='parent', 
                  singleConnection=True)

        # used in rearrange Port and in afterConnect and afterDisconnect
        self.numOfBaseClassInputPorts = len(ip)

        op = self.outputPortsDescr
        op.append(datatype='geomOrInsert2d', name=geomOutputPortName)
        op.append(datatype='geomOrInsert2d(0)', name='allGeometries')

        self.widgetDescr['name'] = {
            'class':'NEEntryNotScheduling',
            'master':'node',
            'width':12,
            'labelCfg':{'text':'new geom name:'},
            'initialValue':'',#geomOutputPortName,
            }

        self.widgetDescr['geoms'] = {
            'class':'NEComboBox', 
            'master':'node', 
            'entryfield_entry_width':14,
            'autoList':True,
            'labelCfg':{'text':'current geom:'},
            'selectioncommand':self.selectionCommand_cb,
            }

        codeAfterConnectParent = """def afterConnect(self, conn):
    #print "codeAfterConnectParent", self.name
    GeometryNode.afterConnectParent(self.node, self, conn)           
"""
        # this code runs when we connect something to the parent port
        ip[self.numOfBaseClassInputPorts-1]['afterConnect'] = \
                                                codeAfterConnectParent 

        codeAfterDisconnectParent = """def afterDisconnect(self, p1, p2):
    #print "codeAfterDisconnectParent", self.name
    GeometryNode.afterDisconnectParent(self.node, p1)
"""
        # this code runs when we disconnect something to the parent port
        ip[self.numOfBaseClassInputPorts-1]['afterDisconnect'] = \
                                                codeAfterDisconnectParent 

        codeAfterDisconnectChildren = """def afterDisconnect(self, p1, p2):
    #print "codeAfterDisconnectChildren", self.name
    GeometryNode.afterDisconnectChildren(self.node, p1, p2)
"""
        # this code runs when we disconnect something to the children port
        op[0]['afterDisconnect'] = codeAfterDisconnectChildren


    def selectionCommand_cb(self, value):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
        if value != self.name \
          and len(self.getInputPortByName('geoms').widget.widget.get(0, 'end')) != len(self.geoms):
            self.geoms[self.selectedGeomIndex].Set(name=value)
            self.rename(value)
        else:
            selInd = self.getInputPortByName('geoms').widget.widget.curselection()
            if len(selInd):
                self.selectedGeomIndex = int(selInd[0])
            self.rename(self.inputPortByName['geoms'].getData())
            self.drawParentAndChildrenConnections()


    def afterAddingToNetwork(self):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
        NetworkNode.afterAddingToNetwork(self)
        if self.geometriesOwner is False:
            self.inputPortByName['parent'].cascadeMenuVariable.set('none')


    def doit(self, name=None, geoms=None, instanceMatrices=None, 
             geomOptions=None, texture=None, parent=None):
        #print "GeometryNode.doit", self.name
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()

#        lgeomState = None
#        if self.selectedGeomIndex is not None and name == '':
#            lgeomState = self.geom().getState()

        if   (self.selectedGeomIndex is None) \
          or (self.inputPortByName['geoms'].hasNewValidData() is False) \
          or (self.inputPortByName['geoms'].hasNewValidData() is True
              and self.inputPortByName['geoms'].getData() == self.name):

            if name is not None:
                if name == '':
                    if self.name == 'the list is empty':
                        name = self.getOutputPortByType('geomOrInsert2d').name
                    else:
                        name = self.name
                else:
                    wn = self.inputPortByName['name'].widget
                    if hasattr(wn,'widget'):
                        wn.set('', run=0)

                if self.selectedGeomIndex is None:
                    self.addingEntries(name)
                elif self.geom().name != name:
                    self.addingEntries(name)

                self.reparentGeomType(
                    self.inputPortByName['parent'].cascadeMenuVariable.get(),
                    reparentCurrent=True)
            elif  (self.selectedGeomIndex is not None) \
              and (parent != self.geom().parent):
                #self.parenting(parent)
                self.reparentGeomType(
                    self.inputPortByName['parent'].cascadeMenuVariable.get(),
                    reparentCurrent=True)
            else:
                pass
        else:
#            if lgeomState is not None:
#                lgeomState.pop('name')
            self.naming(geoms)

        # we never want this to be set in the saved file, 
        # the network will do it through the saved name
        self.inputPortByName['geoms'].widget._modified = False 

        self.drawParentAndChildrenConnections()
        self.rebuildComboBoxList()

        if instanceMatrices:
            mat = array(instanceMatrices)
        else:
            mat = [identity(4).astype('f')]

        lOpts = {}
        if mat is not None and self.selectedGeomIndex is not None:
            #print self.geom(), mat
            if geomOptions is not None:
                lOpts = geomOptions.copy()

            if isinstance(self.geom(), Geom):
                lOpts['instanceMatrices'] = mat
            #lOpts['redo'] = 0
            lOpts['tagModified'] = False
            apply( self.geom().Set, (), lOpts )

#        if lgeomState is not None and self.selectedGeomIndex is not None:
#            for key in lOpts.keys():
#                if key in lgeomState.keys():
#                    lgeomState.pop(key)
#            apply(self.geom().Set,(),lgeomState)


    def rearrangePorts(self):
        """reoder port list for backwards compatibility.
(old network were saved with port indices instead of port names)        
the parent and instance matrices ports are created by the base class
and therefore appear in position 1 and 2 but should move after the
ports created by the subclass. the port 'name' used to be index #4.
"""
        #print "rearrangePorts"
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
        ip = self.inputPortsDescr
        if len(ip) > self.numOfBaseClassInputPorts:
           self.inputPortsDescr = ip[self.numOfBaseClassInputPorts:] \
                                  + ip[:self.numOfBaseClassInputPorts]

    
    def geom(self):
        #print "geom"
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
        assert self.selectedGeomIndex is not None
        #if len(self.geoms) == 0:
        #    self.selectedGeomIndex = None
        assert len(self.geoms) > 0
        assert self.selectedGeomIndex >= 0, self.selectedGeomIndex
        assert self.selectedGeomIndex < len(self.geoms), self.selectedGeomIndex
        return self.geoms[self.selectedGeomIndex]


    def removeAndCleanMultipleGeometryies(self, geometryies):
        #print "removeAndCleanMultipleGeometryies", geometryies
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
        if geometryies is None:
            return
        
        if isinstance(geometryies, Common2d3dObject):
            if hasattr(geometryies, 'node'):
                for g in geometryies.node().geoms:
                    self.removeAndCleanMultipleGeometryies(g.parent)
                    self.removeSingleGeom(g)
        else: 
            for g in geometryies:
                self.removeAndCleanMultipleGeometryies(g.parent)
                self.removeSingleGeom(g)


    def removeSingleGeomAndReparentChildren(self, geom):
        #print "removeSingleGeomAndReparentChildren", geom.name
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
        if not isinstance(geom, Common2d3dObject):
            return

        if hasattr(geom, 'viewer'):
            viewer = geom.viewer
        else:
            viewer = None
        
        from copy import copy
        lChildren = copy(geom.children)
        for c in lChildren:
            if hasattr(c, 'node'):
                c.node().parenting(None,c)
            elif viewer is not None:
                viewer.ReparentObject(c, geom.parent)

        if viewer is None:
            return

        if self.geometriesOwner is True:
            geom.viewer.RemoveObject(geom)
            geom.parent = None


    def removeSingleGeom(self, geom):
        #print "removeSingleGeom", geom.name
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
        if not isinstance(geom, Common2d3dObject):
            return
        
        if not hasattr(geom, 'viewer'):
            return
        
        if geom.viewer is None:
            return

        if self.geometriesOwner is True:
            lViewer = geom.viewer
            lViewer.RemoveObject(geom)
            if geom.parent == lViewer.rootObject:
               geom.parent = None 


    def removeViewerConnection(self, doNotRemoveThisConnection = None): 
        #print "removeViewerConnection", doNotRemoveThisConnection.port1, doNotRemoveThisConnection.port2
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
        lOutputPort = self.getOutputPortByType('geomOrInsert2d')
        for c in lOutputPort.connections:
            from DejaVu2.VisionInterface.DejaVu2Nodes import Viewer
            if c != doNotRemoveThisConnection:
                if isinstance(c.port2.node, Viewer):
                    self.network.deleteConnectionsNoCB(c)
                elif isinstance(c.port2.node, MacroOutputNode):
                    m = c.port2.node.macroNode
                    lOutputPortMacro = m.getOutputPortByName(c.port2.name)
                    for cm in lOutputPortMacro.connections:
                        if cm == doNotRemoveThisConnection:
                            break
                        elif isinstance(cm.port2.node, Viewer):
                            cm.network.deleteConnectionsNoCB(cm)
                            break


    def removeViewerConnections(self, doNotRemoveThisConnection = None): 
        #print "removeViewerConnections", doNotRemoveThisConnection
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
        for g in self.geoms:
            p = g.LastParentBeforeRoot()
            for c in p.AllObjects():
                if hasattr(c,'node'):
                    c.node().removeViewerConnection(doNotRemoveThisConnection)


    def ensureNameOfNodeAndDescendants(self, obj):
        #print "ensureNameOfNodeAndDescendants"
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
        if len(self.geoms) > 0:
            for geom in obj.AllObjects():   
                if hasattr(geom,'node'):
                    geom.node().rebuildComboBoxList()


    def getIndexFromName(self, name):
        """returns the index of the geom
returns -1 if entry is not in the list
"""
        #print "getIndexFromName"
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
        
        for i in range(len(self.geoms)):
            if name == self.geoms[i].name:
                return i
        return -1

        
    def removePreviousGeomWithSameName(self, name):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
        lSelfGeomNames = []
        for g in self.geoms:
            lSelfGeomNames.append(g.name)
        if name in lSelfGeomNames:
            self.naming(name)


    def removePreviousGeomsStartingWith(self, name):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
        lSelfGeomNames = []
        for g in self.geoms:
            lSelfGeomNames.append(g.name)
        for lSelfGeomName in lSelfGeomNames:
            if lSelfGeomName.startswith(name):
                self.naming(lSelfGeomName)


    def removePreviousGeomsWithSameName(self, geoms):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
        lSelfGeomNames = []
        for g in self.geoms:
            lSelfGeomNames.append(g.name)
        for g in geoms:
            if g.name in lSelfGeomNames:
                self.naming(g.name)


    def naming(self, name): 
        """
"""     
        #print "Naming old:" , self.name, "new:", name
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
        
        if name is None:
            return
        
        w0 = self.inputPortByName['geoms'].widget
        if hasattr(w0,'widget'):
            w = w0.widget
        else:
            w = None
            self.inputPortByName['geoms']._previousWidgetDescr['choices'] = []

        if name == '' or (self.selectedGeomIndex is not None and name == self.geom().name):  
            #print "deleting entry", self.name
            
            if w.size() == 0: 
                return  
            
            if self.selectedGeomIndex is None:
                #print "the list is empty 1"
#                w1 = self.inputPortByName['name'].widget
#                if hasattr(w1,'widget'):
#                    w1.set(self.name, run=0)
                self.rename("the list is empty")
                if w is not None:
                    w.delete(0,'end')
                return

            assert self.geom().name == self.name
            lRemovedName = self.name

            #remove self.geom in the viewer 
            #and verify children renaming/reparenting
            self.removeSingleGeomAndReparentChildren(self.geom())
            
            if w is not None:
                if w.get(self.selectedGeomIndex) == self.name:
                    w.delete(self.selectedGeomIndex)
                
            self.geoms.remove(self.geom())
            
            if len(self.geoms) > 0:  
                self.selectedGeomIndex = 0
                self.rebuildComboBoxList()
                if w is not None:                
                    w.selectitem(self.selectedGeomIndex,setentry=1) 
                self.rename(self.geom().name)
            else:
                #print "the list is empty 2"
                self.selectedGeomIndex = None
                self.rename("the list is empty")

            #self.inputPortByName['name'].widget.set(lRemovedName, run=0)

            return lRemovedName
                
        index = self.getIndexFromName(name)
        
        if index >= 0:
            # name is found in the geoms list

            #remove corresponding geom from the viewer 
            #and verify children renaming/reparenting
            self.removeSingleGeomAndReparentChildren(self.geoms[index])
            self.geoms.remove(self.geoms[index])

            self.selectedGeomIndex = 0

            self.rebuildComboBoxList()

            if w is not None:                
                w.selectitem(self.selectedGeomIndex,setentry=1) 
            self.rename(self.geom().name)
            
#            #print "selecting entry", name
#            self.selectedGeomIndex = index
            assert self.selectedGeomIndex < len(self.geoms) , len(self.geoms)
#            self.rename(self.geom().name)
            return
        elif self.selectedGeomIndex is not None:    
            #print "renaming entries", name
            self.geom().Set(name=name)
            self.rename(self.geom().name)
            return

 
    def addingEntries(self, name): 
        #print "addingEntries", name
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
        lenMoreGeoms = self.appendGeometry(name)
        self.selectedGeomIndex = len(self.geoms) - 1
        if lenMoreGeoms > 0:
            for i in range(lenMoreGeoms):
                lIndex = self.selectedGeomIndex -lenMoreGeoms + i + 1
                self.geoms[lIndex].node = ref(self)
                self.geoms[lIndex].replace = False
            self.rename(self.geom().name)
        elif len(self.geoms) > 0:
            self.selectedGeomIndex = len(self.geoms) - 1
        else:
            self.selectedGeomIndex = None

        assert self.selectedGeomIndex < len(self.geoms) , self.geoms


    #to silent pychecker
    def appendGeometry(self, name):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
        assert False, name
        pass


    def nodeOnlyGeomRemoval(self, geom):
        #print "cleanNodeForGeomRemoval", geom.name
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()

        assert geom in self.geoms    

        if self.geom() == geom:
            lThatIsTheSelectedOne = True
        else:
            lThatIsTheSelectedOne = False

        self.geoms.remove(geom)

        if lThatIsTheSelectedOne is True:
            if len(self.geoms) > 0:  
                self.selectedGeomIndex = 0
                self.rename(self.geom().name)
                w0 = self.inputPortByName['geoms'].widget
                if hasattr(w0,'widget'):
                    w0.widget.selectitem(self.selectedGeomIndex,setentry=1) 
            else:
                #print "the list is empty 3"
                self.selectedGeomIndex = None
                #self.inputPortByName['name'].widget.set(self.name, run=0)
                self.rename("the list is empty")


    def parenting(self, parent, geom=None):
        """ prepare the parenting and launch it if applicable
it will parent the selectedGeomIndex unless geom is provided
"""
        #print "Parenting", parent , geom
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
        #import traceback;traceback.print_stack()
        #import pdb;pdb.set_trace()
        
        if geom is None:
            if self.selectedGeomIndex is None:
                return 'Go'
            geom = self.geom()   

        if geom.viewer is not None and parent is None:
            parent = geom.viewer.rootObject            
        else:
            if isinstance(parent, Common2d3dObject): 
                if hasattr(parent, 'node'):
                   if parent.node().selectedGeomIndex is not None:
                       parent = parent.node().geom()
                   else:
                       return 'Go'
                #else:
                #    parent = parent
            elif isinstance(parent, list):
                if len(parent) > 0 and hasattr(parent[0], 'node'):
                    #Get the correct parent in the list
                    parent = parent[0].node().geom() 
                    #parent = parent[parent[0].node().selectedGeomIndex] #same thing
                else:
                    parent = None
            else:
                assert parent is None
                geom.parent = None
           
        if parent is not None and parent!=geom.parent:
            #print "le parent" , parent
            #print "le geom.parent" , geom.parent
            if not self.effectiveParenting(parent, geom):
                warn( \
    'the parent is either in another viewer or is a child of the geometry')
                return 'Stop'  # node fails to execute    

        return 'Go'
    

    def effectiveParenting(self, lParent, geom=None):
        """try to reparent geometry, return True upon success
after this call the geometry has a valid name
it will parent the selectedGeomIndex unless geom is provided
"""
        #print "EffectiveParenting"
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()

        if geom is None:
            geom = self.geom()

        if geom.viewer is None: # we cannot do anything yet because we do
            geom.parent = lParent # not know the viewer yet
            geom.fullName = geom.BuildFullName()
            return True     # lParent will be used when geom is added to Viewer

        if lParent is None:
            lParent = geom.viewer.rootObject

        # try to reparent geometry
        if lParent.viewer is None:
            lHighestModifiedParent = \
                geom.viewer.AddObjectAndParents(lParent, parent=lParent.parent,
                                            local=False)
        else:
            lHighestModifiedParent = geom

        # get unique name among new parent's children
        name = geom.viewer.ensureUniqueName(geom.name, lParent, local=False, exclude=geom)
        #name = geom.viewer.ensureUniqueName(name, lParent, local=False, exclude=geom)

        # set the geometry's name NOT using geom.Set() because parent is still
        # the old one and geom.Set(name=..) checks name validity
        geom.name = name

        geom.viewer.ReparentObject(geom, lParent,
            objectRetainsCurrentPosition=self.inputPortByName['parent'].retainPosition.get())

        self.ensureNameOfNodeAndDescendants(lHighestModifiedParent)
        
        #except AssertionError: # if parent is different viewer than geom
        #                       # or ...
        #    return False

        return True

    
    def drawParentAndChildrenConnections(self): 
        """ finds the selected children and parent
and draws the connections if they are selected in their own node
"""     
        #print "DrawParentAndChildrenConnections", self.name
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()

        if self.selectedGeomIndex is None:
            return

        # first we remove the existing parent connection   
        self.removeParentConnection() 
        
        # we remove the existing children connections 
        self.removeChildrenConnections()
        
        # set the parent connection if parent is selected in his own node
        self.setParentConnection()
        
        # set the children connections when the children are selected
        self.setChildrenConnections()


    def setParentConnection(self):
        """ set the parent connection if parent is selected in his own node
"""
        #print "setParentConnection", self
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
        # first we set corectly the parent node children port
        for g in self.geoms:
            p = g.parent
            v = g.viewer
            if p is not None and v is not None and p != v.rootObject:
                if hasattr(p,'node'):
                    p.node().setPortsOfTheChildrenConnections()

        # then we set this node parent port
        lParentPort = self.inputPortByName['parent']
        lParentPort.hasHiddenConnections = False

        lParent = self.geom().parent

        if lParent is None:
            lParentPort.deleteIcon()
            lParentPort.createIcon()
            return

        viewer = self.geom().viewer
        if viewer is not None and lParent == viewer.rootObject:
            lParentPort.deleteIcon()
            lParentPort.createIcon()
            return

        if hasattr(lParent,'node'):
            lParentNode = lParent.node()
            if (lParentNode.selectedGeomIndex is not None) \
               and (self.geom().parent == lParentNode.geom() ):

                lParentNodeOutputPort = lParentNode.getOutputPortByType('geomOrInsert2d')

                if lParentNode.network == self.network:
                    self.network.connectNodes(lParentNode, self,
                                      lParentNodeOutputPort.name,
                                      'parent',
                                      doNotSchedule=True, doNotCb=True)
                elif isinstance(self.network, MacroNetwork):
                    #import pdb;pdb.set_trace()
                    lMacroNode = self.network.macroNode
                    while lParentNode.network != lMacroNode.network:
                        lMacroNode = lMacroNode.network.macroNode
                    lPortNameInMacro = self.name + '_parent'
                    if lMacroNode.inputPortByName.has_key(lPortNameInMacro) is False:
                        for p in lMacroNode.inputPorts:
                             if p.name.endswith('_parent'):
                                 lPortNameInMacro = p.name
                                 break
                        else: #we didn't break
                            assert "weird! we should have found the port _parent"
                    lParentNode.network.connectNodes(
                          lParentNode, 
                          lMacroNode,
                          lParentNodeOutputPort.name,
                          lPortNameInMacro,
                          doNotSchedule=True, doNotCb=True)
                elif isinstance(lParentNode.network, MacroNetwork):
                    #import pdb;pdb.set_trace()
                    lMacroNode = lParentNode.network.macroNode
                    while self.network != lMacroNode.network:
                        lMacroNode = lMacroNode.network.macroNode
                    lPortNameInMacro = lParentNode.name + '_' + lParentNodeOutputPort.name
                    if lMacroNode.outputPortByName.has_key(lPortNameInMacro) is False:
                        for p in lMacroNode.outputPorts:
                             if p.name.endswith('_' + lParentNodeOutputPort.name):
                                 lPortNameInMacro = p.name
                                 break
                        else: #we didn't break
                            assert "weird! we should have found the port _%s"%lParentNodeOutputPort.name
                    lMacroNode.network.connectNodes(
                          lMacroNode, 
                          self,
                          lPortNameInMacro,
                          'parent',
                          doNotSchedule=True, doNotCb=True)
                    # TODO: make sure the white/black outline is set on the ports
                    #lMacroNode.inputPortByName[lPortNameInMacro].deleteIcon()
                    #lMacroNode.inputPortByName[lPortNameInMacro].createIcon()
            else:
                lParentPort.hasHiddenConnections = True
        else:
#            if len(lParentPort.connections) == 1:
#                #if lParent != lParentPort.connections[0].port1.data:
#                if lParent != lParentPort.getData():
#                    lParentPort.hasHiddenConnections = True
#            if len(lParentPort.connections) == 0:
#                if not ( (lParent is None) \
#                   or (    (self.geom().viewer is not None) \
#                       and (lParent == self.geom().viewer.rootObject) ) ):
#                    lParentPort.hasHiddenConnections = True

            if lParent != lParentPort.getData():
                if not ( (lParent is None) \
                   or (    (self.geom().viewer is not None) \
                       and (lParent == self.geom().viewer.rootObject) ) ):
                    lParentPort.hasHiddenConnections = True

        lParentPort.deleteIcon()
        lParentPort.createIcon()


    def removeParentConnection(self): 
        #print "removeParentConnection", self
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()

        lParentPort = self.inputPortByName['parent']
        for c in lParentPort.connections[::-1]: # going backwards
            if isinstance(c.port1.node, GeometryNode) is True:
                self.network.deleteConnectionsNoCB(c)
            elif isinstance(c.port1.node, MacroInputNode) is True:
                if self.geom().parent is None:
                    continue
                if self.geom().viewer is not None \
                   and self.geom().parent == self.geom().viewer.rootObject:
                    continue
                lParentNode = self.geom().parent.node()
                lParentNetwork = lParentNode.network
                lMacroNode = c.port1.node
                lMacroNetwork = lMacroNode.network
                while lParentNetwork != lMacroNetwork:
                    lMacroNode = lMacroNetwork.macroNode
                    lMacroNetwork = lMacroNode.network
                lPortNameInMacro = self.name + '_parent'
                if lMacroNode.inputPortByName.has_key(lPortNameInMacro) is False:
                    for p in lMacroNode.inputPorts:
                         if p.name.endswith('_parent'):
                             lPortNameInMacro = p.name
                             break;
                    else: #we didn't break
                        assert "weird! we should have found the port _parent"
                lConnections = lMacroNode.inputPortByName[lPortNameInMacro].connections
                if len(lConnections) > 0:
                    # there is a unique connection 
                    lMacroNetwork.deleteConnectionsNoCB(lConnections[0])
                else:
                    #needed for macro node
                    self.effectiveParenting(None)
            elif self.selectedGeomIndex is not None:
                if self.geom().parent is None:
                    self.network.deleteConnectionsNoCB(c)
                elif self.geom().viewer is not None \
                   and self.geom().parent == self.geom().viewer.rootObject:
                    self.network.deleteConnectionsNoCB(c)
                elif self.geom().parent != lParentPort.getData():
                    self.network.deleteConnectionsNoCB(c)


    def setChildrenConnections(self):      
        #print "setChildrenConnections"  
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
        for lChild in self.geom().children:
            if hasattr(lChild, 'node'):
                if lChild == lChild.node().geom(): # the child is selected in its node
                    lChild.node().setParentConnection()
        
        self.setPortsOfTheChildrenConnections()
            
            
    def setPortsOfTheChildrenConnections(self):        
        #print "setPortsOfTheChildrenConnections"
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
        lChildPort = self.getOutputPortByType('geomOrInsert2d')
        lChildPort.hasHiddenConnections = False

        for lChild in self.geom().children:
            if hasattr(lChild, 'node'):
                lOtherPort = lChild.node().inputPortByName['parent']
                if lChild == lChild.node().geom(): 
                    # the child is selected in its node
                    lOtherPort.hasHiddenConnections = False
                    lOtherPort.deleteIcon()
                    lOtherPort.createIcon()
                else:
                    lChildPort.hasHiddenConnections = True

        lChildPort.deleteIcon()
        lChildPort.createIcon()
        
            
    def removeChildrenConnections(self): 
        #print "removeChildrenConnections"
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
        for lChild in self.geom().children:
            if hasattr(lChild, 'node'):
                if lChild == lChild.node().geom(): # the child is selected in its node
                    lChild.node().removeParentConnection()
            
            
    def rebuildComboBoxList(self):
        """
"""
        #print "RebuildComboBoxList"
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
        w0 = self.inputPortByName['geoms'].widget
        if not hasattr(w0,'widget'):
            return
        w = w0.widget
        w.delete(0,'end')
        if self.selectedGeomIndex is None:
            w.setentry('')
            return
        for g in self.geoms:
            w.insert('end',g.name)
        w.selectitem(self.selectedGeomIndex,setentry=1)
        lName = self.geom().name
        if self.name != lName:
            self.rename(lName)
            
            
    def reparentGeoms(self, list):
        #print "reparentGeoms", list
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()

        parent = self.inputPortByName['parent'].getData()
        #print "parent", parent
        for g in list:
            self.parenting(parent,g)


    def reparentGeomType(self, type, reparentCurrent = True):
        #print "reparentGeomType", type
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()

        if self.selectedGeomIndex is None:
            return

        lGeomsToReparent = []
        geom = self.geom()
        if type == 'current':
            lGeomsToReparent.append(geom)
        elif type == 'siblings':
            # build list of siblings
            for g in self.geoms:
                if g is not geom or reparentCurrent is True:
                    #print g.parent
                    #print geom.parent
                    if g.parent is geom.parent:
                        lGeomsToReparent.append(g)
        elif type == 'all':
            # build list
            assert type == 'all' , type
            for g in self.geoms:
                if g is not geom or reparentCurrent is True:
                    lGeomsToReparent.append(g)
        else:
            assert type == 'none' , type
        self.reparentGeoms(lGeomsToReparent)


    def getViewerConnection(self, node): 
        #print "getViewerConnection", self.name
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
        lOutputPort = node.getOutputPortByType('geomOrInsert2d')
        if lOutputPort is None:
            return None
        for c in lOutputPort.connections:
            from DejaVu2.VisionInterface.DejaVu2Nodes import Viewer
            if isinstance(c.port2.node, Viewer):
                return c
            elif isinstance(c.port2.node, MacroOutputNode):
                m = c.port2.node.macroNode
                lOutputPortMacro = m.getOutputPortByName(c.port2.name)
                for cm in lOutputPortMacro.connections:
                    if isinstance(cm.port2.node, Viewer):
                        return cm
#            elif isinstance(c.port2.node, MacroNode):
#                lViewerNode = self.getViewerConnection(c.port2.node)
#                if lViewerNode is not None:
#                    return lViewerNode
#            elif hasattr(c.port2.node, 'getViewerConnection'):
#                lViewerNode = c.port2.node.getViewerConnection(c.port2.node)
#                if lViewerNode is not None:
#                    return lViewerNode
            else:
                lViewerNode = self.getViewerConnection(c.port2.node)
                if lViewerNode is not None:
                    return lViewerNode
        else:
            return None


    def afterDisconnectParent(self, p1):
        #print "GeometryNode afterDisconnectParent"
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
        #import pdb;pdb.set_trace()

        lViewerConnection = self.getViewerConnection(self)
        if lViewerConnection is not None:
            lParentNode = p1.node
            if isinstance(lParentNode, GeometryNode):
                lParentOutputPort = lParentNode.getOutputPortByType('geomOrInsert2d')
                lViewerGeomPort = lViewerConnection.port2
                if lParentNode.network == lViewerGeomPort.node.network:
                    self.network.connectNodes(
                                  lParentNode, 
                                  lViewerGeomPort.node,
                                  lParentOutputPort.name,
                                  lViewerGeomPort.name,
                                  doNotSchedule=True, doNotCb=True)

        self.parenting(None)
        self.rebuildComboBoxList()


    def afterConnectParent(self, port, conn):
        #print "GeometryNode afterConnectParent"
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
        if port.getData() == 'Stop':
            port.data = None
        else:
            if self.selectedGeomIndex is not None:
                self.reparentGeomType(conn.port2.cascadeMenuVariable.get(),
                                      reparentCurrent=True)
            self.rebuildComboBoxList()


    def afterDisconnectChildren(self, p1, p2):
        #print "GeometryNode afterDisconnectChildren", self.name
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
        #import pdb;pdb.set_trace()

        if hasattr(p2, 'node'):
            lChidrenNode = p2.node
            if isinstance(lChidrenNode, MacroNode):
                lViewerConnection = self.getViewerConnection(lChidrenNode)
                if lViewerConnection is not None:
                    lParentNode = p1.node
                    lParentOutputPort = lParentNode.getOutputPortByType('geomOrInsert2d')
                    lViewerGeomPort = lViewerConnection.port2
                    if lParentNode.network == lViewerGeomPort.node.network:
                        self.network.connectNodes(
                                          lParentNode, 
                                          lViewerGeomPort.node,
                                          lParentOutputPort.name,
                                          lViewerGeomPort.name,
                                          doNotSchedule=True, doNotCb=True)


    def beforeRemovingFromNetwork(self):
        #print "beforeRemovingFromNetwork"
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
        #self.removeParentConnection()
        if len(self.geoms) > 0:
            #self.removeChildrenConnections()
            for g in self.geoms:
                self.removeSingleGeomAndReparentChildren(g)
                del(g.node)


    def textureManagement(self, image, textureCoordinates=None):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
        if image is not None and self.selectedGeomIndex is not None:
            g = self.geom()       
            from DejaVu2.Texture import Texture
            lTexture = Texture(enable=1, image=image)  
            if textureCoordinates is not None:
                g.Set(texture=lTexture, textureCoords=textureCoordinates)
            else:
                g.Set(texture=lTexture)
            if hasattr(g.vertexSet, 'texCoords'):
                g.vertexSet.texCoords.array = \
                          lTexture.getTextureCoordinatesAfterResizeRatio(
                                              g.vertexSet.texCoords.array)


class IndexedPolygonsNE(GeometryNode):
    """Build a geometry for drawing a set of polygons sharing vertices.

The geometry nodes handle a list of geometries. When the name changes, 
a new geom is generated. If an empty string is provided, the current geometry 
is deleted.

The parent and children connections reflect the state of the geometry selected 
in the combobox. If the parent or children geoms are currently selected in 
their own node: the connections are present. If the parent or a child is not 
visible (for the selected geometry) the parent port or/and the output port 
contour lines are white (instead of black). Before parenting you can set the 
behavior of the 'parent' inputport. Right clicking on the port allows you to 
extend the parenting to the sibling geometries present in the node or to all 
the geometries of the node. This behavior set by right clicking will only 
affect the next parenting. In the same way, deleting the parent connection, 
reparent to 'root' the selected geometry or the sibling geometries or all the 
geometries in the child node (depending on the right click menu choice).

Input Ports common to all GeometryNodes:
    new geom name: name of a new geometry to be created with the new data.
                   when the node runs, this new geometry will be selected
                   in the 'current geom' field and sent to the first output port.
    current geom: if 'new geom name' is empty, the new data coming from the
                  input ports are applied to the current geom.
                  at the end of the run the selected geometry is sent to
                  the first output port.
    instanceMatrices: stream of 4x4 transformation matrices (default identity)
    geomOptions: accept a dictionary of properties with values.
                 can be fed with node "SetGeom Options"
    parent:  a geometry to parent the selected geometry (default root)

Output Ports common to all GeometryNodes:
    first output port: currently selected geometry 
    allGeometries: list of all the geometries in the node

Specific Input Ports:
    coords:  Nx3 floating points values (type coordinates3D)
    indices: nested lists of indices (1 tuple per object)
    vnormals: normal vectors at each vertex
    colors:  list of RGB tuples (optional)
    image: images from PIL can be use as texture
    textureCoordinates: used if their number match coords or indices
"""

    def __init__(self, name='indexedPolygons', **kw):
        
        kw['name'] = name
        apply( GeometryNode.__init__, (self, 'indexedPolygons'), kw )        
        
        ip = self.inputPortsDescr
        ip.append(datatype='coord3(0)', name='coords')
        ip.append(datatype='indice3or4(0)', name='indices')
        ip.append(datatype='normal3(0)', required=False, name='vnormals')
        ip.append(datatype='colorfloat3or4(0)', required=False, name='colors')
        ip.append(datatype='image', required=False, name='image')
        ip.append(datatype='coord2(0)', required=False, name='textureCoordinates')

        # for backward compatibility
        # this make sure that the first ports in the node are the ones 
        # that use to be there before we introduced the GeometryNode class
        # (old network were saved with port indices instead of port names)
        self.rearrangePorts()
 
        code = """def doit(self, coords, indices, vnormals, colors, image,
textureCoordinates, name, geoms, instanceMatrices, geomOptions, parent):
    #print "doit", self.name
    #apply( GeometryNode.doit, (self,)+args )
        
    GeometryNode.doit(self, name=name, geoms=geoms,
                      instanceMatrices=instanceMatrices,
                      geomOptions=geomOptions, parent=parent)

    if self.selectedGeomIndex is not None:
        g = self.geom()
        if colors is not None:
            inheritMaterial = False
        else:
            inheritMaterial = None
        g.Set(vertices=coords, faces=indices, materials=colors, 
              tagModified=False, 
              vnormals=vnormals, inheritMaterial=inheritMaterial )

        GeometryNode.textureManagement(self, 
                                       image=image,
                                       textureCoordinates=textureCoordinates)

        self.outputData(indexedPolygons=self.geom(), allGeometries=self.geoms)
    elif len(self.geoms) > 0:
        self.outputData(allGeometries=self.geoms)
"""
        if code: self.setFunction(code)

        
    def appendGeometry(self, name):
        self.removePreviousGeomWithSameName(name)
        from DejaVu2.IndexedPolygons import IndexedPolygons
        self.geoms.append(IndexedPolygons(name=name))
        return 1 # num of geoms appended



class IndexedPolylinesNE(GeometryNode):
    """Build a geometry for drawing a set of lines sharing vertices.

The geometry nodes handle a list of geometries. When the name changes, 
a new geom is generated. If an empty string is provided, the current geometry 
is deleted.

The parent and children connections reflect the state of the geometry selected 
in the combobox. If the parent or children geoms are currently selected in 
their own node: the connections are present. If the parent or a child is not 
visible (for the selected geometry) the parent port or/and the output port 
contour lines are white (instead of black). Before parenting you can set the 
behavior of the 'parent' inputport. Right clicking on the port allows you to 
extend the parenting to the sibling geometries present in the node or to all 
the geometries of the node. This behavior set by right clicking will only 
affect the next parenting. In the same way, deleting the parent connection, 
reparent to 'root' the selected geometry or the sibling geometries or all the 
geometries in the child node (depending on the right click menu choice).

Input Ports common to all GeometryNodes:
    new geom name: name of a new geometry to be created with the new data.
                   when the node runs, this new geometry will be selected
                   in the 'current geom' field and sent to the first output port.
    current geom: if 'new geom name' is empty, the new data coming from the
                  input ports are applied to the current geom.
                  at the end of the run the selected geometry is sent to
                  the first output port.
    instanceMatrices: stream of 4x4 transformation matrices (default identity)
    geomOptions: accept a dictionary of properties with values.
                 can be fed with node "SetGeom Options"
    parent:  a geometry to parent the selected geometry (default root)

Output Ports common to all GeometryNodes:
    first output port: currently selected geometry 
    allGeometries: list of all the geometries in the node

Specific Input Ports:
    coords:  Nx3 floating points values (type coordinates3D)
    indices: nested lists of indices (1 tuple per polygon)
    vnormals: normal vectors at each vertex
    colors:  list of RGB tuples (optional)
    image: images from PIL can be use as texture
    textureCoordinates: used if their number match coords or indices
"""
    
    def __init__(self, name='indexedPolylines', **kw):
        
        kw['name'] = name
        apply( GeometryNode.__init__, (self,'indexedPolylines'), kw )        
        
        ip = self.inputPortsDescr
        ip.append(datatype='coord3(0)', name='coords')
        ip.append(datatype='indice2+(0)', name='indices')
        ip.append(datatype='normal3(0)', required=False, name='vnormals')
        ip.append(datatype='colorfloat3or4(0)', required=False, name='colors')
        ip.append(datatype='image', required=False, name='image')
        ip.append(datatype='coord2(0)', required=False, name='textureCoordinates')

        # for backward compatibility
        # this make sure that the first ports in the node are the ones 
        # that use to be there before we introduced the GeometryNode class
        # (old network were saved with port indices instead of port names)
        self.rearrangePorts()
 
        code = """def doit(self, coords, indices, vnormals, colors, image,
textureCoordinates, name, geoms, instanceMatrices, geomOptions, parent):
    #print "doit", self.name
    GeometryNode.doit(self, name=name, geoms=geoms,
                      instanceMatrices=instanceMatrices,
                      geomOptions=geomOptions, parent=parent)
    
    if self.selectedGeomIndex is not None:
        g = self.geom()
        if colors is not None:
            inheritMaterial = False
        else:
            inheritMaterial = None
        g.Set(vertices=coords, faces=indices, materials=colors, tagModified=False,
              vnormals=vnormals, inheritMaterial=inheritMaterial)

        GeometryNode.textureManagement(self, 
                                       image=image,
                                       textureCoordinates=textureCoordinates)

        self.outputData(indexedPolylines=g, allGeometries=self.geoms)
    elif len(self.geoms) > 0:
        self.outputData(allGeometries=self.geoms)
"""           
        if code: self.setFunction(code)

        
    def appendGeometry(self, name):
        self.removePreviousGeomWithSameName(name)
        from DejaVu2.IndexedPolylines import IndexedPolylines
        self.geoms.append(IndexedPolylines(name))
        return 1 # num of geoms appended



class Cylinders(GeometryNode):
    """Build a geometry for drawing a set of cylinders.    

The geometry nodes handle a list of geometries. When the name changes, 
a new geom is generated. If an empty string is provided, the current geometry 
is deleted.

The parent and children connections reflect the state of the geometry selected 
in the combobox. If the parent or children geoms are currently selected in 
their own node: the connections are present. If the parent or a child is not 
visible (for the selected geometry) the parent port or/and the output port 
contour lines are white (instead of black). Before parenting you can set the 
behavior of the 'parent' inputport. Right clicking on the port allows you to 
extend the parenting to the sibling geometries present in the node or to all 
the geometries of the node. This behavior set by right clicking will only 
affect the next parenting. In the same way, deleting the parent connection, 
reparent to 'root' the selected geometry or the sibling geometries or all the 
geometries in the child node (depending on the right click menu choice).

Input Ports common to all GeometryNodes:
    new geom name: name of a new geometry to be created with the new data.
                   when the node runs, this new geometry will be selected
                   in the 'current geom' field and sent to the first output port.
    current geom: if 'new geom name' is empty, the new data coming from the
                  input ports are applied to the current geom.
                  at the end of the run the selected geometry is sent to
                  the first output port.
    instanceMatrices: stream of 4x4 transformation matrices (default identity)
    geomOptions: accept a dictionary of properties with values.
                 can be fed with node "SetGeom Options"
    parent:  a geometry to parent the selected geometry (default root)

Output Ports common to all GeometryNodes:
    first output port: currently selected geometry 
    allGeometries: list of all the geometries in the node

Specific Input Ports:
    coords:  Nx3 floating points values (type coordinates3D)
    indices: nested lists of indices (1 tuple per polygon)
    colors:  list of RGB tuples (optional)
    image: images from PIL can be use as texture
    radius:  float or list of floats of same length as coords (default is 1.0)
    quality: cylinder quality, the hight, the more polygons are drawn (default 5)
"""
    
    def __init__(self, name='cylinders', **kw):
        
        kw['name'] = name
        apply( GeometryNode.__init__, (self,'cylinders'), kw )        
        
        ip = self.inputPortsDescr
        ip.append(datatype='coord3(0)', name='coords', required=False, defaultValue=((0,0,0),(10,0,0)) )
        ip.append(datatype='indice2+(0)', name='indices', required=False, defaultValue=((0,1),) )
        ip.append(datatype='colorfloat3or4(0)', name='colors', required=False)
        ip.append(datatype='float(0)', name='radius', required=False, defaultValue=1.)
        ip.append(datatype='int', name='quality', required=False, defaultValue=5)
        ip.append(datatype='image', name='image', required=False)
        
        # for backward compatibility
        # this make sure that the first ports in the node are the ones 
        # that use to be there before we introduced the GeometryNode class
        # (old network were saved with port indices instead of port names)
        self.rearrangePorts()
 
        self.widgetDescr['quality'] = {
            'class':'NEThumbWheel', 'master':'node', 'initialValue':5,
            'labelCfg':{'text':'quality:'}, 'width':80, 'height':20, 
            'type':'int', 'oneTurn':10, 'lockBMin':1, 'min':0, 'wheelPad':2 }

        self.widgetDescr['radius'] = {
            'class':'NEThumbWheel', 'master':'node', 'initialValue':1.0,
            'labelCfg':{'text':'radius:'}, 'width':80, 'height':20,
            'type':'float', 'oneTurn':10., 'lockBMin':1, 'min':0.1,
            'wheelPad':2 }
                     
        code = """def doit(self, coords, indices, colors, radius, quality, image,
name, geoms, instanceMatrices, geomOptions, parent):
    #print "doit", self.name
    GeometryNode.doit(self, name=name, geoms=geoms,
                      instanceMatrices=instanceMatrices,
                      geomOptions=geomOptions, parent=parent)

    if self.selectedGeomIndex is not None:
        g = self.geom()
        if colors is not None:
            inheritMaterial = False
        else:
            inheritMaterial = None
        g.Set(vertices=coords, faces=indices, materials=colors, tagModified=False, 
              radii=radius, quality=quality, inheritMaterial=inheritMaterial)

        GeometryNode.textureManagement(self, image=image)

        self.outputData(cylinders=g, allGeometries=self.geoms)
    elif len(self.geoms) > 0:
        self.outputData(allGeometries=self.geoms)
"""           
        if code: self.setFunction(code)


    def appendGeometry(self, name):
        self.removePreviousGeomWithSameName(name)
        from DejaVu2.Cylinders import Cylinders
        self.geoms.append(Cylinders(name))
        return 1 # num of geoms appended



class AxisNE(GeometryNode):
    """Build a CylinderArrow geometry for drawing an AXIS.   

The geometry nodes handle a list of geometries. When the name changes, 
a new geom is generated. If an empty string is provided, the current geometry 
is deleted.

The parent and children connections reflect the state of the geometry selected 
in the combobox. If the parent or children geoms are currently selected in 
their own node: the connections are present. If the parent or a child is not 
visible (for the selected geometry) the parent port or/and the output port 
contour lines are white (instead of black). Before parenting you can set the 
behavior of the 'parent' inputport. Right clicking on the port allows you to 
extend the parenting to the sibling geometries present in the node or to all 
the geometries of the node. This behavior set by right clicking will only 
affect the next parenting. In the same way, deleting the parent connection, 
reparent to 'root' the selected geometry or the sibling geometries or all the 
geometries in the child node (depending on the right click menu choice).

Input Ports common to all GeometryNodes:
    new geom name: name of a new geometry to be created with the new data.
                   when the node runs, this new geometry will be selected
                   in the 'current geom' field and sent to the first output port.
    current geom: if 'new geom name' is empty, the new data coming from the
                  input ports are applied to the current geom.
                  at the end of the run the selected geometry is sent to
                  the first output port.
    instanceMatrices: stream of 4x4 transformation matrices (default identity)
    geomOptions: accept a dictionary of properties with values.
                 can be fed with node "SetGeom Options"
    parent:  a geometry to parent the selected geometry (default root)

Output Ports common to all GeometryNodes:
    first output port: currently selected geometry 
    allGeometries: list of all the geometries in the node

Specific Input Ports:
    coords:  Nx3 floating points values (type coordinates3D)
    indices: nested lists of indices (1 tuple per polygon)
    colors:  list of RGB tuples (optional)
    image: images from PIL can be use as texture
    radius:  float or list of floats of same length as coords (default is 1.0)
    quality: cylinder quality, the hight, the more polygons are drawn (default 5)
    headLength:
    headRadius:
"""
    
    def __init__(self, name='axis', **kw):
        
        kw['name'] = name
        apply( GeometryNode.__init__, (self,'axis'), kw )        
        
        ip = self.inputPortsDescr
        ip.append(datatype='coord3(0)', name='coords', required=False, defaultValue=((0,0,0),(10,0,0)) )
        ip.append(datatype='indice2+(0)', name='indices', required=False, defaultValue=((0,1),) )
        ip.append(datatype='colorfloat3or4(0)', name='colors', required=False)
        ip.append(datatype='float(0)', name='radius', required=False, defaultValue=.6)
        ip.append(datatype='int', name='quality', required=False, defaultValue=5)
        ip.append(datatype='float', name='headLength', required=False)
        ip.append(datatype='float', name='headRadius', required=False)
        ip.append(datatype='image', name='image', required=False)
        
        # for backward compatibility
        # this make sure that the first ports in the node are the ones 
        # that use to be there before we introduced the GeometryNode class
        # (old network were saved with port indices instead of port names)
        self.rearrangePorts()
 
        self.widgetDescr['quality'] = {
            'class':'NEThumbWheel', 'master':'node', 'initialValue':5,
            'labelCfg':{'text':'quality:'}, 'width':80, 'height':20, 
            'type':'int', 'oneTurn':10, 'lockBMin':1, 'min':1, 'wheelPad':2 }

        self.widgetDescr['radius'] = {
            'class':'NEThumbWheel', 'master':'node', 'initialValue':1.0,
            'labelCfg':{'text':'radius:'}, 'width':80, 'height':20,
            'type':'float', 'oneTurn':10., 'lockBMin':1, 'min':0.1,
            'wheelPad':2 }
                     
        self.widgetDescr['headLength'] = {
            'class':'NEThumbWheel', 'initialValue':1.,
            'labelCfg':{'text':'headLength:'}, 'width':80, 'height':20, 
            'type':'float', 'oneTurn':10, 'lockBMin':1, 'min':0,
            'wheelPad':2 }
            
        self.widgetDescr['headRadius'] = {
            'class':'NEThumbWheel', 'initialValue':2.,
            'labelCfg':{'text':'headRadius:'}, 'width':80, 'height':20, 
            'type':'float', 'oneTurn':10, 'lockBMin':1, 'min':0,
            'wheelPad':2 }

        code = """def doit(self, coords, indices, colors, radius, quality, headLength,
headRadius, image, name, geoms, instanceMatrices, geomOptions, parent):
    #print "doit", self.name
    GeometryNode.doit(self, name=name, geoms=geoms,
                      instanceMatrices=instanceMatrices,
                      geomOptions=geomOptions, parent=parent)
   
    if self.selectedGeomIndex is not None:
        g = self.geom()
        from opengltk.OpenGL import GL
        if colors is not None:
            inheritMaterial = False
        else:
            inheritMaterial = None
        g.Set(vertices=coords, faces=indices, materials=colors, tagModified=False,
              culling=GL.GL_NONE, inheritCulling=False,
              backPolyMode=GL.GL_FILL, inheritBackPolyMode=False,
              quality=quality, radii=radius, inheritMaterial=inheritMaterial)

        GeometryNode.textureManagement(self, image=image)

        self.outputData(axis=g, allGeometries=self.geoms)
    elif len(self.geoms) > 0:
        self.outputData(allGeometries=self.geoms)
"""           
        if code: self.setFunction(code)

        
    def appendGeometry(self, name):
        self.removePreviousGeomWithSameName(name)
        from DejaVu2.Cylinders import CylinderArrows
        self.geoms.append(CylinderArrows(name))
        return 1 # num of geoms appended



class Spheres(GeometryNode):
    """Build a geometry for drawing a set of spheres.    

The geometry nodes handle a list of geometries. When the name changes, 
a new geom is generated. If an empty string is provided, the current geometry 
is deleted.

The parent and children connections reflect the state of the geometry selected 
in the combobox. If the parent or children geoms are currently selected in 
their own node: the connections are present. If the parent or a child is not 
visible (for the selected geometry) the parent port or/and the output port 
contour lines are white (instead of black). Before parenting you can set the 
behavior of the 'parent' inputport. Right clicking on the port allows you to 
extend the parenting to the sibling geometries present in the node or to all 
the geometries of the node. This behavior set by right clicking will only 
affect the next parenting. In the same way, deleting the parent connection, 
reparent to 'root' the selected geometry or the sibling geometries or all the 
geometries in the child node (depending on the right click menu choice).

Output Ports:
    first output port: currently selected geometry 
    allGeometries: list of all the geometries in the node

Input Ports:
    new geom name: name of a new geometry to be created with the new data.
                   when the node runs, this new geometry will be selected
                   in the 'current geom' field and sent to the first output port.
    current geom: if 'new geom name' is empty, the new data coming from the
                  input ports are applied to the current geom.
                  at the end of the run the selected geometry is sent to
                  the first output port.
    instanceMatrices: stream of 4x4 transformation matrices (default identity)
    geomOptions: accept a dictionary of properties with values.
                 can be fed with node "SetGeom Options"
    parent:  a geometry to parent the selected geometry (default root)

    coords:  Nx3 floating points values (type coordinates3D)
    colors:  list of RGB tuples (optional)
    image: images from PIL can be use as texture
    radius:  float or list of floats of same length as coords (default is 1.0)
    quality: cylinder quality, the hight, the more polygons are drawn (default 5)
"""
 
    def __init__(self, name='spheres', **kw):
        
        kw['name'] = name
        apply( GeometryNode.__init__, (self,'spheres'), kw )        
        
        ip = self.inputPortsDescr
        ip.append(datatype='coord3(0)', name='coords', required=False, defaultValue=((0,0,0),))
        ip.append(datatype='colorfloat3or4(0)', name='colors', required=False)
        ip.append(datatype='float(0)', name='radius', required=False, defaultValue=1.)
        ip.append(datatype='int', name='quality', required=False, defaultValue=5)
        ip.append(datatype='image', name='image', required=False)
         
        # for backward compatibility
        # this make sure that the first ports in the node are the ones 
        # that use to be there before we introduced the GeometryNode class
        # (old network were saved with port indices instead of port names)
        self.rearrangePorts()

        self.widgetDescr['quality'] = {
            'class':'NEThumbWheel', 'master':'node', 'initialValue':5,
            'labelCfg':{'text':'quality:'}, 'width':80, 'height':20, 
            'type':'int', 'oneTurn':10, 'lockBMin':1, 'min':0, 'wheelPad':2 }

        self.widgetDescr['radius'] = {
            'class':'NEThumbWheel', 'master':'node', 'initialValue':1.0,
            'labelCfg':{'text':'radius:'}, 'width':80, 'height':20,
            'type':'float', 'oneTurn':10., 'lockBMin':1, 'min':0.0000001,
            'wheelPad':2 }
                     
        code = """def doit(self, coords, colors, radius, quality, image,
name, geoms, instanceMatrices, geomOptions, parent):
    #print "doit", self.name
    GeometryNode.doit(self, name=name, geoms=geoms,
                      instanceMatrices=instanceMatrices,
                      geomOptions=geomOptions, parent=parent)
    
    try:
        len(radius)
        islist = True
    except TypeError:
        islist = False
    if islist:
        assert (len(radius)==1) or (len(radius)==len(coords))

    if self.selectedGeomIndex is not None:
        g = self.geom()
        if colors is not None:
            inheritMaterial = False
        else:
            inheritMaterial = None
        g.Set(vertices=coords, materials=colors, tagModified=False, 
              radii=radius, quality=quality, inheritMaterial=inheritMaterial)
              
        GeometryNode.textureManagement(self, image=image)
              
        self.outputData(spheres=g, allGeometries=self.geoms)
    elif len(self.geoms) > 0:
        self.outputData(allGeometries=self.geoms)
"""           
        if code: self.setFunction(code)

        
    def appendGeometry(self, name):
        self.removePreviousGeomWithSameName(name)
        from DejaVu2.Spheres import Spheres
        self.geoms.append(Spheres(name))
        return 1 # num of geoms appended



class Points(GeometryNode):
    """Build a geometry for drawing a set of 3D Points

The geometry nodes handle a list of geometries. When the name changes, 
a new geom is generated. If an empty string is provided, the current geometry 
is deleted.

The parent and children connections reflect the state of the geometry selected 
in the combobox. If the parent or children geoms are currently selected in 
their own node: the connections are present. If the parent or a child is not 
visible (for the selected geometry) the parent port or/and the output port 
contour lines are white (instead of black). Before parenting you can set the 
behavior of the 'parent' inputport. Right clicking on the port allows you to 
extend the parenting to the sibling geometries present in the node or to all 
the geometries of the node. This behavior set by right clicking will only 
affect the next parenting. In the same way, deleting the parent connection, 
reparent to 'root' the selected geometry or the sibling geometries or all the 
geometries in the child node (depending on the right click menu choice).

Output Ports:
    first output port: currently selected geometry 
    allGeometries: list of all the geometries in the node

Input Ports:
    new geom name: name of a new geometry to be created with the new data.
                   when the node runs, this new geometry will be selected
                   in the 'current geom' field and sent to the first output port.
    current geom: if 'new geom name' is empty, the new data coming from the
                  input ports are applied to the current geom.
                  at the end of the run the selected geometry is sent to
                  the first output port.
    instanceMatrices: stream of 4x4 transformation matrices (default identity)
    geomOptions: accept a dictionary of properties with values.
                 can be fed with node "SetGeom Options"
    parent:  a geometry to parent the selected geometry (default root)

    coords:  Nx3 floating points values (type coordinates3D)
    colors:  list of RGB tuples (optional)
    size:    point size as an integer
"""
 
    def __init__(self, name='points', **kw):
        
        kw['name'] = name
        apply( GeometryNode.__init__, (self,'points'), kw )        
        
        ip = self.inputPortsDescr
        ip.append(datatype='coord3(0)', name='coords', required=False, defaultValue=((0,0,0),))
        ip.append(datatype='colorfloat3or4(0)', name='colors', required=False)
        ip.append(datatype='int', name='size', required=False, defaultValue=3)
         
        # for backward compatibility
        # this make sure that the first ports in the node are the ones 
        # that use to be there before we introduced the GeometryNode class
        # (old network were saved with port indices instead of port names)
        self.rearrangePorts()

        self.widgetDescr['size'] = {
            'class':'NEThumbWheel', 'master':'node', 'initialValue':3,
            'labelCfg':{'text':'size:'}, 'width':80, 'height':20, 
            'type':'int', 'oneTurn':10, 'lockBMin':1, 'min':1, 'wheelPad':2 }

        code = """def doit(self, coords, colors, size,
name, geoms, instanceMatrices, geomOptions, parent):

    GeometryNode.doit(self, name=name, geoms=geoms,
                      instanceMatrices=instanceMatrices,
                      geomOptions=geomOptions, parent=parent)

    if self.selectedGeomIndex is not None:
        g = self.geom()
        if colors is not None:
            inheritMaterial = False
        else:
            inheritMaterial = None

        g.Set(vertices=coords, materials=colors, tagModified=False, 
              pointWidth=size, inheritPointWidth=False,
              inheritMaterial=inheritMaterial)
              
        self.outputData(points=g, allGeometries=self.geoms)
    elif len(self.geoms) > 0:
        self.outputData(allGeometries=self.geoms)
"""           
        if code: self.setFunction(code)

        
    def appendGeometry(self, name):
        self.removePreviousGeomWithSameName(name)
        from DejaVu2.Points import Points
        self.geoms.append(Points(name))
        return 1 # num of geoms appended



class CrossSet(GeometryNode):
    """Build a geometry for drawing a set of 3D crosses.    

The geometry nodes handle a list of geometries. When the name changes, 
a new geom is generated. If an empty string is provided, the current geometry 
is deleted.

The parent and children connections reflect the state of the geometry selected 
in the combobox. If the parent or children geoms are currently selected in 
their own node: the connections are present. If the parent or a child is not 
visible (for the selected geometry) the parent port or/and the output port 
contour lines are white (instead of black). Before parenting you can set the 
behavior of the 'parent' inputport. Right clicking on the port allows you to 
extend the parenting to the sibling geometries present in the node or to all 
the geometries of the node. This behavior set by right clicking will only 
affect the next parenting. In the same way, deleting the parent connection, 
reparent to 'root' the selected geometry or the sibling geometries or all the 
geometries in the child node (depending on the right click menu choice).

Output Ports:
    first output port: currently selected geometry 
    allGeometries: list of all the geometries in the node

Input Ports:
    new geom name: name of a new geometry to be created with the new data.
                   when the node runs, this new geometry will be selected
                   in the 'current geom' field and sent to the first output port.
    current geom: if 'new geom name' is empty, the new data coming from the
                  input ports are applied to the current geom.
                  at the end of the run the selected geometry is sent to
                  the first output port.
    instanceMatrices: stream of 4x4 transformation matrices (default identity)
    geomOptions: accept a dictionary of properties with values.
                 can be fed with node "SetGeom Options"
    parent:  a geometry to parent the selected geometry (default root)

    coords:  Nx3 floating points values (type coordinates3D)
    colors:  list of RGB tuples (optional)
    size:    float, cross size
"""
 
    def __init__(self, name='crossSet', **kw):
        
        kw['name'] = name
        apply( GeometryNode.__init__, (self,'crossSet'), kw )        
        
        ip = self.inputPortsDescr
        ip.append(datatype='coord3(0)', name='coords', required=False, defaultValue=((0,0,0),))
        ip.append(datatype='colorfloat3or4(0)', name='colors', required=False)
        ip.append(datatype='float', name='size', required=False, defaultValue=1.)
         
        # for backward compatibility
        # this make sure that the first ports in the node are the ones 
        # that use to be there before we introduced the GeometryNode class
        # (old network were saved with port indices instead of port names)
        self.rearrangePorts()

        self.widgetDescr['size'] = {
            'class':'NEThumbWheel', 'master':'node', 'initialValue':3.,
            'labelCfg':{'text':'quality:'}, 'width':80, 'height':20, 
            'type':'float', 'oneTurn':10., 'lockBMin':0., 'min':0., 'wheelPad':2 }

        code = """def doit(self, coords, colors, size,
name, geoms, instanceMatrices, geomOptions, parent):

    GeometryNode.doit(self, name=name, geoms=geoms,
                      instanceMatrices=instanceMatrices,
                      geomOptions=geomOptions, parent=parent)
    
    if self.selectedGeomIndex is not None:
        g = self.geom()
        if colors is not None:
            inheritMaterial = False
        else:
            inheritMaterial = None
        g.Set(vertices=coords, materials=colors, tagModified=False, 
              offset=size, inheritMaterial=inheritMaterial)
              
        self.outputData(crossSet=g, allGeometries=self.geoms)
    elif len(self.geoms) > 0:
        self.outputData(allGeometries=self.geoms)
"""           
        if code: self.setFunction(code)

        
    def appendGeometry(self, name):
        self.removePreviousGeomWithSameName(name)
        from DejaVu2.Points import CrossSet
        self.geoms.append(CrossSet(name))
        return 1 # num of geoms appended



class GeomContainer(GeometryNode):
    """Build a geometry that can be used as a container for other geometries.

The geometry nodes handle a list of geometries. When the name changes, 
a new geom is generated. If an empty string is provided, the current geometry 
is deleted.

The parent and children connections reflect the state of the geometry selected 
in the combobox. If the parent or children geoms are currently selected in 
their own node: the connections are present. If the parent or a child is not 
visible (for the selected geometry) the parent port or/and the output port 
contour lines are white (instead of black). Before parenting you can set the 
behavior of the 'parent' inputport. Right clicking on the port allows you to 
extend the parenting to the sibling geometries present in the node or to all 
the geometries of the node. This behavior set by right clicking will only 
affect the next parenting. In the same way, deleting the parent connection, 
reparent to 'root' the selected geometry or the sibling geometries or all the 
geometries in the child node (depending on the right click menu choice).

Input Ports common to all GeometryNodes:
    new geom name: name of a new geometry to be created with the new data.
                   when the node runs, this new geometry will be selected
                   in the 'current geom' field and sent to the first output port.
    current geom: if 'new geom name' is empty, the new data coming from the
                  input ports are applied to the current geom.
                  at the end of the run the selected geometry is sent to
                  the first output port.
    instanceMatrices: stream of 4x4 transformation matrices (default identity)
    geomOptions: accept a dictionary of properties with values.
                 can be fed with node "SetGeom Options"
    parent:  a geometry to parent the selected geometry (default root)

Output Ports common to all GeometryNodes:
    first output port: currently selected geometry 
    allGeometries: list of all the geometries in the node

Specific Input Ports:
"""
    """Set of polygons sharing vertices."""
    
    def __init__(self, name='geomContainer', **kw):
        
        kw['name'] = name
        apply( GeometryNode.__init__, (self,'geomContainer'), kw )        
        
        #ip = self.inputPortsDescr
        
        # for backward compatibility
        # this make sure that the first ports in the node are the ones 
        # that use to be there before we introduced the GeometryNode class
        # (old network were saved with port indices instead of port names)
        # AS THERE IS NO INPUT PORT HERE, WE DON'T NEED TO CALL THIS
        #self.rearrangePorts()                     

        code = """def doit(self, 
name, geoms, instanceMatrices, geomOptions, parent):
    #print "doit", self.name
    GeometryNode.doit(self, name=name, geoms=geoms,
                      instanceMatrices=instanceMatrices,
                      geomOptions=geomOptions, parent=parent)
    
    if self.selectedGeomIndex is not None:
        g = self.geom()
        #g.Set( tagModified=False)
        self.outputData(geomContainer=g, allGeometries=self.geoms)
    elif len(self.geoms) > 0:
        self.outputData(allGeometries=self.geoms)
"""           
        if code: self.setFunction(code)


    def appendGeometry(self, name):
        self.removePreviousGeomWithSameName(name)
        self.geoms.append(Geom(name))
        return 1 # num of geoms appended



class Ellipsoids(GeometryNode):
    """Build a geometry for drawing a set of Ellipsoids.    

The geometry nodes handle a list of geometries. When the name changes, 
a new geom is generated. If an empty string is provided, the current geometry 
is deleted.

The parent and children connections reflect the state of the geometry selected 
in the combobox. If the parent or children geoms are currently selected in 
their own node: the connections are present. If the parent or a child is not 
visible (for the selected geometry) the parent port or/and the output port 
contour lines are white (instead of black). Before parenting you can set the 
behavior of the 'parent' inputport. Right clicking on the port allows you to 
extend the parenting to the sibling geometries present in the node or to all 
the geometries of the node. This behavior set by right clicking will only 
affect the next parenting. In the same way, deleting the parent connection, 
reparent to 'root' the selected geometry or the sibling geometries or all the 
geometries in the child node (depending on the right click menu choice).

Input Ports common to all GeometryNodes:
    new geom name: name of a new geometry to be created with the new data.
                   when the node runs, this new geometry will be selected
                   in the 'current geom' field and sent to the first output port.
    current geom: if 'new geom name' is empty, the new data coming from the
                  input ports are applied to the current geom.
                  at the end of the run the selected geometry is sent to
                  the first output port.
    instanceMatrices: stream of 4x4 transformation matrices (default identity)
    geomOptions: accept a dictionary of properties with values.
                 can be fed with node "SetGeom Options"
    parent:  a geometry to parent the selected geometry (default root)

Output Ports common to all GeometryNodes:
    first output port: currently selected geometry 
    allGeometries: list of all the geometries in the node

Specific Input Ports:
    ellipsoids: list of ellipsoid objects (can be generated using the EllipsoidFit node)
    colors:  list of RGB tuples (optional)
    image: images from PIL can be use as texture
"""
 
    def __init__(self, name='ellipsoids', **kw):
        
        kw['name'] = name
        apply( GeometryNode.__init__, (self,'ellipsoids'), kw )        
        
        ip = self.inputPortsDescr
        ip.append(datatype='list', name='ellipsoids')
        ip.append(datatype='colorfloat3or4(0)', required=False, name='colors')
        ip.append(datatype='image', required=False, name='image')

        # for backward compatibility
        # this make sure that the first ports in the node are the ones 
        # that use to be there before we introduced the GeometryNode class
        # (old network were saved with port indices instead of port names)
        self.rearrangePorts()
                      
        code = """def doit(self, ellipsoids, colors, image,
name, geoms, instanceMatrices, geomOptions, parent):
    #print "doit", self.name
    GeometryNode.doit(self, name=name, geoms=geoms,
                      instanceMatrices=instanceMatrices,
                      geomOptions=geomOptions, parent=parent)

    if ellipsoids is not None:
        centers = []
        size = []
        orient = []
        for ellipse in ellipsoids:
            centers.append( ellipse.getPosition() )
            size.append( ellipse.getAxis().astype('f') )
            orint = identity(4, 'f')
            orint[:3, :3] = ellipse.getOrientation()
            orient.append( orint )
    
    if self.selectedGeomIndex is not None and ellipsoids is not None:
        g = self.geom()
        if colors is not None:
            inheritMaterial = False
        else:
            inheritMaterial = None
        g.Set(centers=centers, scaling=size, orientation=orient, 
              materials=colors, tagModified=False, 
              quality=30, inheritMaterial=inheritMaterial)

        GeometryNode.textureManagement(self, image=image)

        self.outputData(ellipsoids=g, allGeometries=self.geoms)
    elif len(self.geoms) > 0:
        self.outputData(allGeometries=self.geoms)
"""           
        if code: self.setFunction(code)


    def appendGeometry(self, name):
        self.removePreviousGeomWithSameName(name)
        from DejaVu2.Ellipsoids import Ellipsoids
        self.geoms.append(Ellipsoids(name))
        return 1 # num of geoms appended



class Array2DToHeightFieldGeom(GeometryNode):
    """ convert a 2D array of values into a height field

The geometry nodes handle a list of geometries. When the name changes, 
a new geom is generated. If an empty string is provided, the current geometry 
is deleted.

The parent and children connections reflect the state of the geometry selected 
in the combobox. If the parent or children geoms are currently selected in 
their own node: the connections are present. If the parent or a child is not 
visible (for the selected geometry) the parent port or/and the output port 
contour lines are white (instead of black). Before parenting you can set the 
behavior of the 'parent' inputport. Right clicking on the port allows you to 
extend the parenting to the sibling geometries present in the node or to all 
the geometries of the node. This behavior set by right clicking will only 
affect the next parenting. In the same way, deleting the parent connection, 
reparent to 'root' the selected geometry or the sibling geometries or all the 
geometries in the child node (depending on the right click menu choice).

Input Ports common to all GeometryNodes:
    new geom name: name of a new geometry to be created with the new data.
                   when the node runs, this new geometry will be selected
                   in the 'current geom' field and sent to the first output port.
    current geom: if 'new geom name' is empty, the new data coming from the
                  input ports are applied to the current geom.
                  at the end of the run the selected geometry is sent to
                  the first output port.
    instanceMatrices: stream of 4x4 transformation matrices (default identity)
    geomOptions: accept a dictionary of properties with values.
                 can be fed with node "SetGeom Options"
    parent:  a geometry to parent the selected geometry (default root)

Output Ports common to all GeometryNodes:
    first output port: currently selected geometry 
    allGeometries: list of all the geometries in the node

Specific Input Ports:
    data: 2Darray
    scaleZ: Z scaling factor
"""

    def __init__(self, name='heightField', **kw):

        kw['name'] = name
        apply( GeometryNode.__init__, (self,'heightField'), kw)

        ip = self.inputPortsDescr
        ip.append({'name': 'data', 'datatype': '2Darray'})
        ip.append({'name': 'scaleZ', 'datatype': 'float', 'defaultValue':1.})

        self.rearrangePorts()

        self.widgetDescr['scaleZ'] = {
            'class':'NEThumbWheel',
            'master':'node',
            'initialValue':1.0,
            'labelCfg':{'text':'scale Z:'},
            'width':80, 'height':20, 
            'type':'float', 'oneTurn':10., 'lockBMin':1, 'min':0.0,
            'wheelPad':2 }

        code = """def doit(self, data, scaleZ,
name, geoms, instanceMatrices, geomOptions, parent):
    #print "doit", self.name
    GeometryNode.doit(self, name=name, geoms=geoms,
                      instanceMatrices=instanceMatrices,
                      geomOptions=geomOptions, parent=parent)

    if self.selectedGeomIndex is not None and data is not None:
        # normalize the data to range from 0 to 255
        maxi = max(data.ravel())
        mini = min(data.ravel())
        delta = maxi-mini
        data2 = (data-mini)/delta

        verts = []
        quads = []
        sx, sy = data.shape
        for xi in xrange(sx):
            base = sy*xi
            for yi in xrange(sy):
                verts.append( (xi,yi,data[xi,yi]/scaleZ) )
                if xi<sx-1 and yi<sy-1:
                    quads.append( ( base+sy+yi, base+sy+yi+1, base+yi+1,
                                    base+yi ))

        g = self.geom()
        g.Set(vertices=verts, faces=quads, tagModified=False)
        self.outputData(heightField=g, allGeometries=self.geoms)
    elif len(self.geoms) > 0:
        self.outputData(allGeometries=self.geoms)
"""
        self.setFunction(code)


    def appendGeometry(self, name):
        self.removePreviousGeomWithSameName(name)
        from DejaVu2.IndexedPolygons import IndexedPolygons
        self.geoms.append(IndexedPolygons(name))
        return 1 # num of geoms appended



class GeomsFromFile(GeometryNode):
    """inheriting GeometryNode, this node reads several file formats,
STL, OBJ, VERT ... and generates the correcponding geoms.
the geom's names are deducted from the filename.
you can easily add new file types to this node.

The geometry nodes handle a list of geometries. When the name changes, 
a new geom is generated. If an empty string is provided, the current geometry 
is deleted.

The parent and children connections reflect the state of the geometry selected 
in the combobox. If the parent or children geoms are currently selected in 
their own node: the connections are present. If the parent or a child is not 
visible (for the selected geometry) the parent port or/and the output port 
contour lines are white (instead of black). Before parenting you can set the 
behavior of the 'parent' inputport. Right clicking on the port allows you to 
extend the parenting to the sibling geometries present in the node or to all 
the geometries of the node. This behavior set by right clicking will only 
affect the next parenting. In the same way, deleting the parent connection, 
reparent to 'root' the selected geometry or the sibling geometries or all the 
geometries in the child node (depending on the right click menu choice).

Input Ports common to all GeometryNodes:
    new geom name: name of a new geometry to be created with the new data.
                   when the node runs, this new geometry will be selected
                   in the 'current geom' field and sent to the first output port.
    current geom: if 'new geom name' is empty, the new data coming from the
                  input ports are applied to the current geom.
                  at the end of the run the selected geometry is sent to
                  the first output port.
    instanceMatrices: stream of 4x4 transformation matrices (default identity)
    geomOptions: accept a dictionary of properties with values.
                 can be fed with node "SetGeom Options"
    parent:  a geometry to parent the selected geometry (default root)

Output Ports common to all GeometryNodes:
    first output port: currently selected geometry 
    allGeometries: list of all the geometries in the node

Specific Input Ports:
    filename:
    image: images from PIL can be use as texture
"""
    def __init__(self, name='geomsFromFile', **kw):
        kw['name'] = name
        apply( GeometryNode.__init__, (self,'geomsFromFile'), kw )

        ip = self.inputPortsDescr
        ip.append(datatype='string', name='filename')
        ip.append(datatype='image', required=False, name='image')

        self.rearrangePorts()

        self.widgetDescr['filename'] = {
            'class':'NEEntryWithFileBrowser', 'master':'node', 'width':10,
            'initialValue':'',
            'labelCfg':{'text':'Filename: '},
            'filetypes': [('3d object', '*.stl *.obj *.vert *.indpolvert'),
                          ('all', '*')]}
            
        self.filenames = []

        code = """def doit(self, filename, image,
name, geoms, instanceMatrices, geomOptions, parent):
    #print "doit",name, self.name
    lName = os.path.split(os.path.splitext(filename)[0])[-1]  
    if filename not in self.filenames:
        self.filenames.append(filename)
        name = lName
    else:
        for g in self.geoms:
            if (g.name == lName) or g.name.startswith(lName + '+'):
                break
        else: #we didn't break
            if geoms == '':
                self.filenames.remove(filename)
            else:
                name = lName

    GeometryNode.doit(self, name=name, geoms=geoms,
                      instanceMatrices=instanceMatrices,
                      geomOptions=geomOptions, parent=parent)
    
    if self.selectedGeomIndex is not None:
        g = self.geom()
        #g.Set( tagModified=False )

        GeometryNode.textureManagement(self, image=image)

        self.outputData(geomsFromFile=g, allGeometries=self.geoms)
    elif len(self.geoms) > 0:
        self.outputData(allGeometries=self.geoms)
"""
        self.setFunction(code)

    def appendGeometry(self, name):
        # in this, we dont use the parameter "name" at all
        filename = self.inputPortByName['filename'].getData()
        if filename is None or filename == '':
            return

        # read the geoms from the file      
        lExt = os.path.splitext(filename)[-1][1:].lower()
        if lExt == 'stl':
            moreGeoms = self.readFileSTL(filename)
        elif lExt == 'vert':
            moreGeoms = self.readMsmsFile(filename)
        elif lExt == 'indpolvert':
            moreGeoms = self.readFileIndexedPolygons(filename)
        elif lExt == 'obj':
            moreGeoms = self.readFileOBJ(filename)
        else:
            assert False

        # append the geoms
        if moreGeoms is not None:
            self.removePreviousGeomsWithSameName(moreGeoms)
            self.geoms.extend(moreGeoms)
            return len(moreGeoms)
        else:
            return 0


    def nameGeoms(self, listOfGeoms, filename):
        """renames each geom differently
"""
        if len(listOfGeoms) > 0:
            nameSplit = os.path.splitext(filename)
            name = os.path.split(nameSplit[0])[-1]
            listOfGeoms[0].name = name
            for i in range(1, len(listOfGeoms)):
                listOfGeoms[i].name = name + '+' + str(i)


    def readFileSTL(self, filename):
        """from a filename, returns a list of geoms
"""
        from DejaVu2.DataInput import readAnySTL
        moreGeoms = readAnySTL(filename)
        return moreGeoms


    def readMsmsFile(self, filename):
        """from a filename, returns a list of geoms
"""
        try:
            from mslib.msmsParser import MSMSParser
            faceFilename = os.path.splitext(filename)[0] + '.face'
            vertFilename = filename
            print "faceFilename", faceFilename
            from mslib.msmsParser import MSMSParser
            lMsmsParser = MSMSParser()
            from DejaVu2.IndexedPolygons import IndexedPolygons
            vertices, faces  = lMsmsParser.parse(vertFilename, faceFilename)
            singleGeom = IndexedPolygons(vertices=vertices, faces=faces)
            moreGeoms = [ singleGeom ]
            self.nameGeoms(moreGeoms, filename)
            return moreGeoms
	except:
	    return None


    def readFileIndexedPolygons(self, filename):
        """from a filename, returns a list of geoms
"""
        from DejaVu2.IndexedPolygons import IndexedPolygonsFromFile
        singleGeom = IndexedPolygonsFromFile(filename)
        moreGeoms = [ singleGeom ]
        self.nameGeoms(moreGeoms, filename)
        return moreGeoms


    def readFileOBJ(self, filename):
        """from a filename, returns a list of geoms
"""
        from DejaVu2.DataInput import readOBJ
        moreGeoms = readOBJ(filename)
        return moreGeoms




#class GlutLabelsNE(GeometryNode):
#    """Build a geometry for drawing a set of polygons sharing vertices.
#
#The geometry nodes handle a list of geometries. When the name changes, 
#a new geom is generated. If an empty name is provided, the current geometry 
#is deleted.
#
#The parent and children connections reflect the state of the geometry selected 
#in the combo box. If the parent or children geoms are currently selected in 
#their own node: the connections are present. If the parent or a child is not 
#visible (for the selected geometry) the parent port or/and the output port 
#contour lines are white (instead of black). Before parenting you can set the 
#behaviour of the 'parent' inputport. Right clicking on the port allows you to 
#extend the parenting to the sibling geometries present in the node or to all 
#the geometries of the node. This behaviour set by right clicking will only 
#affect the next parenting. In the same way, deleting the parent connection, 
#reparent to 'root' the selected geometry or the sibling geometries or all the 
#geometries in the child node (depending on the right click menu choice).
#
#Input Ports:
#    coords:  Nx3 floating points values (type coordinates3D)
#    labels: list of labels ex:['label1', 'label2', 15.56, 'label4', 4, 5.78]
#            if one label only, it applies to all coords, 
#            otherwise len(coords) should equal len(labels)
#    font:
#    colors:  list of RGB tuples (optional)
#    name:    name of the geometry
#    instanceMatrices: stream of 4x4 transformation matrices (default identity)
#    parent:  a geometry to parent the selected geometry (default root)
#
#Output Ports:
#    indexedPolygons: IndexedPolygons geometry 
#    allGeometries: list of all the geometries in the node
#"""
#
#    def __init__(self, name='glutlabels', **kw):
#        
#        kw['name'] = name
#        apply( GeometryNode.__init__, (self, 'glutlabels'), kw )        
#        
#        ip = self.inputPortsDescr
#        ip.append(datatype='coord3(0)', name='coords')
#        ip.append(datatype='string(0)', name='labels')
#        ip.append(datatype='colorfloat3or4(0)', required=False, name='colors')
#        ip.append(datatype='string', required=False, name='font')
#
#        # for backward compatibility
#        # this make sure that the first ports in the node are the ones 
#        # that use to be there before we introduced the GeometryNode class
#        # (old network were saved with port indices instead of port names)
#        self.rearrangePorts()
#
#        self.widgetDescr['font'] = {
#            'class':'NEComboBox', 
#            'master':'node', 
#            'entryfield_entry_width':14,
#            'choices':[ 'glutTimesRoman10',   
#                        'glut9by15',
#                        'glut8by13',
#                        'glutTimesRoman10',
#                        'glutTimesRoman24',
#                        'glutHelvetica10',
#                        'glutHelvetica12',
#                        'glutHelvetica18'
#                      ],
#            'initialValue':'glutTimesRoman10',
#            'autoList':True,
#            'labelCfg':{'text':'font:'},
#            }
#
#        self.widgetDescr['labels'] = {
#            'class':'NEEntry', 'master':'node', 'width':10,
#            'labelCfg':{'text':'labels:'},
#            }
# 
#        code = """def doit(self, coords, labels, colors, font,
#name, geoms, instanceMatrices, geomOptions, parent):
#    #print "doit", self.name
#    #apply( GeometryNode.doit, (self,)+args )
#    GeometryNode.doit(self, name=name, geoms=geoms,
#                      instanceMatrices=instanceMatrices,
#                      geomOptions=geomOptions, parent=parent)
#    
#    if self.selectedGeomIndex is not None:
#        g = self.geom()
#        if colors is not None:
#            inheritMaterial = False
#        else:
#            inheritMaterial = None
#        import types
#        if type(labels) == types.StringType:
#            if len(labels):
#                labels = eval(labels)
#        if (labels is not None) and (labels != ''):
#            g.Set( inheritMaterial=inheritMaterial,
#                   materials=colors,
#                   tagModified=False,
#                   visible=1,
#                   labels=labels,
#                   vertices=coords,
#                   font=font)
#        else:
#           g.Set(visible=0)
#        self.outputData(glutlabels=g, allGeometries=self.geoms)
#    else:
#        self.outputData(glutlabels=None, allGeometries=self.geoms)
#"""
#        if code: self.setFunction(code)
#
#
#    def appendGeometry(self, name):
#        self.removePreviousGeomWithSameName(name)
#        from DejaVu2.Labels import Labels
#        self.geoms.append( Labels(name=name, shape=(0,3) ) )
#        return 1 # num of geoms appended



class GlfLabelsNE(GeometryNode):
    """Build a geometry for drawing a set of polygons sharing vertices.

The geometry nodes handle a list of geometries. When the name changes, 
a new geom is generated. If an empty string is provided, the current geometry 
is deleted.

The parent and children connections reflect the state of the geometry selected 
in the combobox. If the parent or children geoms are currently selected in 
their own node: the connections are present. If the parent or a child is not 
visible (for the selected geometry) the parent port or/and the output port 
contour lines are white (instead of black). Before parenting you can set the 
behavior of the 'parent' inputport. Right clicking on the port allows you to 
extend the parenting to the sibling geometries present in the node or to all 
the geometries of the node. This behavior set by right clicking will only 
affect the next parenting. In the same way, deleting the parent connection, 
reparent to 'root' the selected geometry or the sibling geometries or all the 
geometries in the child node (depending on the right click menu choice).

Input Ports common to all GeometryNodes:
    new geom name: name of a new geometry to be created with the new data.
                   when the node runs, this new geometry will be selected
                   in the 'current geom' field and sent to the first output port.
    current geom: if 'new geom name' is empty, the new data coming from the
                  input ports are applied to the current geom.
                  at the end of the run the selected geometry is sent to
                  the first output port.
    instanceMatrices: stream of 4x4 transformation matrices (default identity)
    geomOptions: accept a dictionary of properties with values.
                 can be fed with node "SetGeom Options"
    parent:  a geometry to parent the selected geometry (default root)

Output Ports common to all GeometryNodes:
    first output port: currently selected geometry 
    allGeometries: list of all the geometries in the node

Specific Input Ports:
    coords:  Nx3 floating points values (type coordinates3D)
    vnormals: normal vectors at each vertex
    colors:  list of RGB tuples (optional)
    image: images from PIL can be use as texture
    labels: list of labels ex:['label1', 'label2', 15.56, 'label4', 4, 5.78]
            if one label only, it applies to all coords, 
            otherwise len(coords) should equal len(labels)
    font:
"""

    def __init__(self, name='glflabels', **kw):
        
        kw['name'] = name
        apply( GeometryNode.__init__, (self, 'glflabels'), kw )        
        
        ip = self.inputPortsDescr
        ip.append(datatype='coord3(0)', name='coords', required=False, defaultValue=((0,0,0),) )
        ip.append(datatype='string(0)', name='labels', required=False)
        ip.append(datatype='normal3(0)', name='vnormals', required=False)        
        ip.append(datatype='colorfloat3or4(0)', name='colors', required=False)
        ip.append(datatype='boolean', name='billboard', required=False)
        ip.append(datatype='boolean', name='includeCameraRotationInBillboard', required=False)
        ip.append(datatype='string', name='font', required=False)
        ip.append(datatype='string', name='fontstyle', required=False)
        ip.append(datatype='float', name='fontSpacing', required=False)
        ip.append(datatype='float', name='scaleXYZ', required=False)
        ip.append(datatype='float', name='scaleX', required=False)
        ip.append(datatype='float', name='scaleY', required=False)
        ip.append(datatype='float', name='scaleZ', required=False)
        ip.append(datatype='float', name='fontRotateAngleX', required=False)
        ip.append(datatype='float', name='fontRotateAngleY', required=False)
        ip.append(datatype='float', name='fontRotateAngleZ', required=False)
        ip.append(datatype='float', name='fontTranslateX', required=False)
        ip.append(datatype='float', name='fontTranslateY', required=False)
        ip.append(datatype='float', name='fontTranslateZ', required=False)
        ip.append(datatype='image', name='image', required=False)

        # for backward compatibility
        # this make sure that the first ports in the node are the ones 
        # that use to be there before we introduced the GeometryNode class
        # (old network were saved with port indices instead of port names)
        self.rearrangePorts()

        from DejaVu2.glfLabels import GlfLabels
        self.widgetDescr['labels'] = {
            'class':'NEEntry',
            'master':'node',
            'initialValue':['Aa','Bb','Cc'],
            'width':10,
            'labelCfg':{'text':'labels:'},
            }
 
        self.widgetDescr['billboard'] = {
            'class':'NECheckButton', 
            'master':'node', 
            'initialValue':0,
            'labelCfg':{'text':'billboard:'},
            }

        self.widgetDescr['includeCameraRotationInBillboard'] = {
            'class':'NECheckButton', 
            'master':'node', 
            'initialValue':0,
            'labelCfg':{'text':'includeCameraRotationInBillboard:'},
            }

        self.widgetDescr['font'] = {
            'class':'NEComboBox', 
            'master':'node', 
            'entryfield_entry_width':14,
            'choices':GlfLabels.fontList,
            'initialValue':GlfLabels.fontList[0],
            'autoList':True,
            'labelCfg':{'text':'font:'},
            }

        self.widgetDescr['fontstyle'] = {
            'class':'NEComboBox', 
            'master':'node', 
            'entryfield_entry_width':14,
            'choices':GlfLabels.fontStyleList,
            'initialValue':GlfLabels.fontStyleList[3],
            'autoList':True,
            'labelCfg':{'text':'font style:'},
            }

        self.widgetDescr['fontSpacing'] = {
            'class':'NEDial',
            'master':'ParamPanel',
            'size':50,
            'oneTurn':1,
            'min':0.,
            'max':1.,
            'type':'float',
            'showLabel':1,
            'continuous':0,
            'initialValue':.2,
            'labelCfg':{'text':'fontSpacing'},
            }
 
        self.widgetDescr['scaleXYZ'] = {
            'class':'NEDial',
            'master':'ParamPanel',
            'size':50,
            'oneTurn':1,
            'min':0.,
            'max':100.,
            'type':'float',
            'showLabel':1,
            'continuous':0,
            'initialValue':1.,
            'labelCfg':{'text':'scale XYZ'},
            }
 
        self.widgetDescr['scaleX'] = {
            'class':'NEDial',
            'master':'ParamPanel',
            'size':50,
            'oneTurn':1,
            'min':0.,
            'max':100.,
            'type':'float',
            'showLabel':1,
            'continuous':0,
            'initialValue':1.,
            'labelCfg':{'text':'scale X'},
            }
 
        self.widgetDescr['scaleY'] = {
            'class':'NEDial',
            'master':'ParamPanel',
            'size':50,
            'oneTurn':1,
            'min':0.,
            'max':100.,
            'type':'float',
            'showLabel':1,
            'continuous':0,
            'initialValue':1.,
            'labelCfg':{'text':'scale Y'},
            }
 
        self.widgetDescr['scaleZ'] = {
            'class':'NEDial',
            'master':'ParamPanel',
            'size':50,
            'oneTurn':1,
            'min':0.,
            'max':100.,
            'type':'float',
            'showLabel':1,
            'continuous':0,
            'initialValue':1.,
            'labelCfg':{'text':'scale Z'},
            }
 
        self.widgetDescr['fontRotateAngleX'] = {
            'class':'NEDial',
            'master':'ParamPanel',
            'size':50,
            'oneTurn':90,
            'min':-180.,
            'max':180.,
            'type':'float',
            'showLabel':1,
            'continuous':0,
            'initialValue':0.,
            'labelCfg':{'text':'font rotation angle X'},
            }

        self.widgetDescr['fontRotateAngleY'] = {
            'class':'NEDial',
            'master':'ParamPanel',
            'size':50,
            'oneTurn':90,
            'min':-180.,
            'max':180.,
            'type':'float',
            'showLabel':1,
            'continuous':0,
            'initialValue':0.,
            'labelCfg':{'text':'font rotation angle Y'},
            }

        self.widgetDescr['fontRotateAngleZ'] = {
            'class':'NEDial',
            'master':'ParamPanel',
            'size':50,
            'oneTurn':90,
            'min':-180.,
            'max':180.,
            'type':'float',
            'showLabel':1,
            'continuous':0,
            'initialValue':0.,
            'labelCfg':{'text':'font rotation angle Z'},
            }

        self.widgetDescr['fontTranslateX'] = {
            'class':'NEDial',
            'master':'ParamPanel',
            'size':50,
            'oneTurn':10,
            'min':-100.,
            'max':100.,
            'type':'float',
            'showLabel':1,
            'continuous':0,
            'initialValue':0.,
            'labelCfg':{'text':'font translation X'},
            }

        self.widgetDescr['fontTranslateY'] = {
            'class':'NEDial',
            'master':'ParamPanel',
            'size':50,
            'oneTurn':10,
            'min':-100.,
            'max':100.,
            'type':'float',
            'showLabel':1,
            'continuous':0,
            'initialValue':0.,
            'labelCfg':{'text':'font translation Y'},
            }

        self.widgetDescr['fontTranslateZ'] = {
            'class':'NEDial',
            'master':'ParamPanel',
            'size':50,
            'oneTurn':10,
            'min':-100.,
            'max':100.,
            'type':'float',
            'showLabel':1,
            'continuous':0,
            'initialValue':0.,
            'labelCfg':{'text':'font translation Z'},
            }

        code = """def doit(self, coords, labels, vnormals, colors, 
billboard, includeCameraRotationInBillboard, font, fontstyle, 
fontSpacing, scaleXYZ, scaleX, scaleY, scaleZ,
fontRotateAngleX, fontRotateAngleY, fontRotateAngleZ, 
fontTranslateX, fontTranslateY, fontTranslateZ, image,
name, geoms, instanceMatrices, geomOptions, parent):
    #print "doit", self.name
    #apply( GeometryNode.doit, (self,)+args )
    GeometryNode.doit(self, name=name, geoms=geoms,
                      instanceMatrices=instanceMatrices,
                      geomOptions=geomOptions, parent=parent)
    
    if self.selectedGeomIndex is not None:
        g = self.geom()
        if colors is not None:
            inheritMaterial = False
        else:
            inheritMaterial = None
        import types
        if type(labels) == types.StringType:
            if len(labels):
                labels = eval(labels)
        if (labels is not None) and (labels != ''):
            g.Set( inheritMaterial=inheritMaterial,
                   materials=colors,
                   polyFace='frontandback',
                   tagModified=False,
                   visible=1,
                   labels=labels,
                   vertices=coords,
                   vnormals=vnormals,
                   billboard=billboard,
                   includeCameraRotationInBillboard=includeCameraRotationInBillboard,
                   font=font,
                   fontStyle=fontstyle,
                   fontSpacing=fontSpacing,
                   fontScales=(scaleXYZ*scaleX, scaleXYZ*scaleY, scaleXYZ*scaleZ),
                   fontRotateAngles=(fontRotateAngleX, fontRotateAngleY, fontRotateAngleZ),
                   fontTranslation=(fontTranslateX, fontTranslateY, fontTranslateZ),
                  )
        else:
           g.Set(visible=0)

        GeometryNode.textureManagement(self, image=image)

        self.outputData(glflabels=g, allGeometries=self.geoms)
    else:
        self.outputData(glflabels=None, allGeometries=self.geoms)
"""
        if code: self.setFunction(code)


    def appendGeometry(self, name):
        self.removePreviousGeomWithSameName(name)
        from DejaVu2.glfLabels import GlfLabels
        self.geoms.append( GlfLabels(name=name, shape=(0,3) ) )
        return 1 # num of geoms appended




class BoxNE(GeometryNode):
    """Displays a 3-dimensional box whose colors depict the different
axes upon which they lie:
    x-axis: Red
    y-axis: Green
    z-axis: Blue

The geometry nodes handle a list of geometries. When the name changes, 
a new geom is generated. If an empty string is provided, the current geometry 
is deleted.

The parent and children connections reflect the state of the geometry selected 
in the combobox. If the parent or children geoms are currently selected in 
their own node: the connections are present. If the parent or a child is not 
visible (for the selected geometry) the parent port or/and the output port 
contour lines are white (instead of black). Before parenting you can set the 
behavior of the 'parent' inputport. Right clicking on the port allows you to 
extend the parenting to the sibling geometries present in the node or to all 
the geometries of the node. This behavior set by right clicking will only 
affect the next parenting. In the same way, deleting the parent connection, 
reparent to 'root' the selected geometry or the sibling geometries or all the 
geometries in the child node (depending on the right click menu choice).

Input Ports common to all GeometryNodes:
    new geom name: name of a new geometry to be created with the new data.
                   when the node runs, this new geometry will be selected
                   in the 'current geom' field and sent to the first output port.
    current geom: if 'new geom name' is empty, the new data coming from the
                  input ports are applied to the current geom.
                  at the end of the run the selected geometry is sent to
                  the first output port.
    instanceMatrices: stream of 4x4 transformation matrices (default identity)
    geomOptions: accept a dictionary of properties with values.
                 can be fed with node "SetGeom Options"
    parent:  a geometry to parent the selected geometry (default root)

Output Ports common to all GeometryNodes:
    first output port: currently selected geometry 
    allGeometries: list of all the geometries in the node

Specific Input Ports:
    image: images from PIL can be use as texture if textureCoordinates are provided
    textureCoordinates: used if their number match coords or indices
"""

    def __init__(self, name='box', **kw):
        kw['name'] = name
        apply( GeometryNode.__init__, (self, 'box'), kw )        

        ip = self.inputPortsDescr
        ip.append(datatype='float', required=False, name='centerx')
        ip.append(datatype='float', required=False, name='centery')
        ip.append(datatype='float', required=False, name='centerz')
        ip.append(datatype='float', required=False, name='lengthx')
        ip.append(datatype='float', required=False, name='lengthy')
        ip.append(datatype='float', required=False, name='lengthz')
        ip.append(datatype='int', required=False, name='constrProps')
        ip.append(datatype='image', required=False, name='image')
        ip.append(datatype='coord2(0)', required=False, name='textureCoordinates')

        # for backward compatibility
        # this make sure that the first ports in the node are the ones 
        # that use to be there before we introduced the GeometryNode class
        # (old network were saved with port indices instead of port names)
        self.rearrangePorts()
       
        wd = self.widgetDescr
        
        wd['centerx'] = {
            'class':'NEThumbWheel', 'width':60, 'height':10, 'type':'float',
            'showLabel':1, 'oneTurn':50., 'lockMin':True,
            'labelCfg':{'text':'x center'},
            'master':'node',
            'initialValue':0.0
           }

        wd['centery'] = {
            'class':'NEThumbWheel', 'width':60, 'height':10, 'type':'float',
            'showLabel':1, 'oneTurn':50., 'lockMin':True,
            'labelCfg':{'text':'y center'},
            'master':'node',
            'initialValue':0.0
           }

        wd['centerz'] = {
            'class':'NEThumbWheel', 'width':60, 'height':10, 'type':'float',
            'showLabel':1, 'oneTurn':50., 'lockMin':True,
            'labelCfg':{'text':'z center'},
            'master':'node',
            'initialValue':0.0
           }

        wd['constrProps'] = {
            'class':'NECheckButton', 'master':'node',
            'initialValue':1,
            'labelCfg':{'text':'Constrain scale:'}
            }

        wd['lengthx'] = {
            'class':'NEThumbWheel', 'width':60, 'height':10, 'type':'float',
            'showLabel':1, 'oneTurn':50., 'lockMin':True,
            'labelCfg':{'text':'x length'}, 'min':0.01,
            'master':'node',
            'initialValue':12.0
           }
    
        wd['lengthy'] = {
            'class':'NEThumbWheel', 'width':60, 'height':10, 'type':'float',
            'showLabel':1, 'oneTurn':50., 'lockMin':True,
            'labelCfg':{'text':'y length'}, 'min':0.01,
            'master':'node',
            'initialValue':12.0
           }

        wd['lengthz'] = {
            'class':'NEThumbWheel', 'width':60, 'height':10, 'type':'float',
            'showLabel':1, 'oneTurn':50., 'lockMin':True,
            'labelCfg':{'text':'z length'}, 'min':0.01,
            'master':'node',
            'initialValue':12.0
           }

        code = """xyzprop='' #global set so that the set scale proportion is remembered.
def doit(self, centerx, centery, centerz, lengthx, lengthy, lengthz, constrProps,
image, textureCoordinates,
name, geoms, instanceMatrices, geomOptions, parent):
    #print "doit", self.name
    #apply( GeometryNode.doit, (self,)+args )
    GeometryNode.doit(self, name=name, geoms=geoms,
                      instanceMatrices=instanceMatrices,
                      geomOptions=geomOptions, parent=parent)

    global xyzprop
    if self.selectedGeomIndex is not None:
        g = self.geom()

        if constrProps:
            if xyzprop=='':
                xratio=(lengthy/lengthx,lengthz/lengthx)
                yratio=(lengthx/lengthy,lengthz/lengthy)
                zratio=(lengthx/lengthz,lengthy/lengthz)
                xyzprop=(xratio,yratio,zratio)
            if self.inputPorts[3].hasNewValidData(): #x-length changed
                lengthy,lengthz=xyzprop[0][0]*lengthx,xyzprop[0][1]*lengthx
            if self.inputPorts[4].hasNewValidData(): #y-length changed
                lengthx,lengthz=xyzprop[1][0]*lengthy,xyzprop[1][1]*lengthy
            if self.inputPorts[5].hasNewValidData(): #z-length changed
                lengthx,lengthy=xyzprop[2][0]*lengthz,xyzprop[2][1]*lengthz
            if self.inputPorts[3]:
                self.inputPorts[3].widget.set(lengthx, run=0)
            if self.inputPorts[4]:
                self.inputPorts[4].widget.set(lengthy, run=0)
            if self.inputPorts[5]:
                self.inputPorts[5].widget.set(lengthz, run=0)
        else: 
            xyzprop=''
        center = (centerx,centery,centerz)
        side = (lengthx,lengthy,lengthz)
        g.Set(center=center, xside=side[0], yside=side[1], zside=side[2])

        GeometryNode.textureManagement(self, 
                                       image=image,
                                       textureCoordinates=textureCoordinates)

        self.outputData(box=g, allGeometries=self.geoms)
    else:
        self.outputData(box=None, allGeometries=self.geoms)
"""
        if code: self.setFunction(code)


    def appendGeometry(self, name):
        self.removePreviousGeomWithSameName(name)
        from DejaVu2.Box import Box
        self.geoms.append( Box(name=name, center=(0.,0.,0.), side=10.,
                           inheritMaterial=0, frontPolyMode='line') )
        return 1 # num of geoms appended



class StickerTextNE(GeometryNode):
    """Build an insert2d for drawing a Sticker containing a 2d text Label.

The geometry nodes handle a list of geometries. When the name changes, 
a new geom is generated. If an empty string is provided, the current geometry 
is deleted.

The parent and children connections reflect the state of the geometry selected 
in the combobox. If the parent or children geoms are currently selected in 
their own node: the connections are present. If the parent or a child is not 
visible (for the selected geometry) the parent port or/and the output port 
contour lines are white (instead of black). Before parenting you can set the 
behavior of the 'parent' inputport. Right clicking on the port allows you to 
extend the parenting to the sibling geometries present in the node or to all 
the geometries of the node. This behavior set by right clicking will only 
affect the next parenting. In the same way, deleting the parent connection, 
reparent to 'root' the selected geometry or the sibling geometries or all the 
geometries in the child node (depending on the right click menu choice).

Input Ports common to all GeometryNodes:
    new geom name: name of a new geometry to be created with the new data.
                   when the node runs, this new geometry will be selected
                   in the 'current geom' field and sent to the first output port.
    current geom: if 'new geom name' is empty, the new data coming from the
                  input ports are applied to the current geom.
                  at the end of the run the selected geometry is sent to
                  the first output port.
    instanceMatrices: stream of 4x4 transformation matrices (default identity)
    geomOptions: accept a dictionary of properties with values.
                 can be fed with node "SetGeom Options"
    parent:  a geometry to parent the selected geometry (default root)

Output Ports common to all GeometryNodes:
    first output port: currently selected geometry 
    allGeometries: list of all the geometries in the node

Specific Input Ports:
    label:  the string to be written in the sticker
"""

    def __init__(self, name='stickerText', **kw):
        kw['name'] = name
        apply( GeometryNode.__init__, (self, 'stickerText'), kw )

        self.widgetDescr['label'] = {
            'class':'NEEntry',
            'master':'node',
            'width':12,
            'labelCfg':{'text':'label'},
            'initialValue':'stickerText'
            }

        ip = self.inputPortsDescr
        ip.append(datatype='string', required=False, name='label', defaultValue='')

        # for backward compatibility
        # this make sure that the first ports in the node are the ones 
        # that use to be there before we introduced the GeometryNode class
        # (old network were saved with port indices instead of port names)
        self.rearrangePorts()

        code = """def doit(self, label,
name, geoms, instanceMatrices, geomOptions, parent):
    #print "doit", self.name

    lgeomState = None
    if self.selectedGeomIndex is not None and label == '':
        label = self.geom().label       

    GeometryNode.doit(self, name=name, geoms=geoms,
                      instanceMatrices=instanceMatrices,
                      geomOptions=geomOptions, parent=parent)

    if self.selectedGeomIndex is not None:
        g = self.geom()
        g.Set(label=label, transparent=True)
        self.outputData(stickerText=g, allGeometries=self.geoms)
    else:
        self.outputData(stickerText=None, allGeometries=self.geoms)

    wn = self.inputPortByName['label'].widget
    if hasattr(wn,'widget'):
        wn.set('', run=0)
"""
        self.setFunction(code)


    def appendGeometry(self, name):
        self.removePreviousGeomWithSameName(name)
        #print "appendGeometry"
        from DejaVu2.glfSticker import Sticker
        self.geoms.append(Sticker(name))
        return 1 # num of geoms appended



class StickerImageNE(GeometryNode):
    """Build an insert2d for drawing a Sticker containing an Image.

The geometry nodes handle a list of geometries. When the name changes, 
a new geom is generated. If an empty string is provided, the current geometry 
is deleted.

The parent and children connections reflect the state of the geometry selected 
in the combobox. If the parent or children geoms are currently selected in 
their own node: the connections are present. If the parent or a child is not 
visible (for the selected geometry) the parent port or/and the output port 
contour lines are white (instead of black). Before parenting you can set the 
behavior of the 'parent' inputport. Right clicking on the port allows you to 
extend the parenting to the sibling geometries present in the node or to all 
the geometries of the node. This behavior set by right clicking will only 
affect the next parenting. In the same way, deleting the parent connection, 
reparent to 'root' the selected geometry or the sibling geometries or all the 
geometries in the child node (depending on the right click menu choice).

Input Ports common to all GeometryNodes:
    new geom name: name of a new geometry to be created with the new data.
                   when the node runs, this new geometry will be selected
                   in the 'current geom' field and sent to the first output port.
    current geom: if 'new geom name' is empty, the new data coming from the
                  input ports are applied to the current geom.
                  at the end of the run the selected geometry is sent to
                  the first output port.
    instanceMatrices: stream of 4x4 transformation matrices (default identity)
    geomOptions: accept a dictionary of properties with values.
                 can be fed with node "SetGeom Options"
    parent:  a geometry to parent the selected geometry (default root)

Output Ports common to all GeometryNodes:
    first output port: currently selected geometry 
    allGeometries: list of all the geometries in the node

Specific Input Ports:
    image:  the Image to be applied on the sticker
"""

    def __init__(self, name='stickerImage', **kw):
        kw['name'] = name
        apply( GeometryNode.__init__, (self, 'stickerImage'), kw )

        ip = self.inputPortsDescr
        ip.append(datatype='image', required=True, name='image')

        # for backward compatibility
        # this make sure that the first ports in the node are the ones 
        # that use to be there before we introduced the GeometryNode class
        # (old network were saved with port indices instead of port names)
        self.rearrangePorts()

        code = """def doit(self, image,
name, geoms, instanceMatrices, geomOptions, parent):
    #print "doit", self.name
    GeometryNode.doit(self, name=name, geoms=geoms,
                      instanceMatrices=instanceMatrices,
                      geomOptions=geomOptions, parent=parent)

    if self.selectedGeomIndex is not None:
        g = self.geom()
        from Image import FLIP_TOP_BOTTOM
        g.Set(image=image.transpose(FLIP_TOP_BOTTOM), transparent=True)
        self.outputData(stickerImage=g, allGeometries=self.geoms)
    else:
        self.outputData(stickerImage=None, allGeometries=self.geoms)
"""
        self.setFunction(code)


    def appendGeometry(self, name):
        self.removePreviousGeomWithSameName(name)
        #print "appendGeometry"
        from DejaVu2.StickerImage import StickerImage
        self.geoms.append(StickerImage(name))
        return 1 # num of geoms appended



class OneTexturedQuadNE(GeometryNode):
    """Build an quad in 3d for drawing a Sticker containing an Image.

The geometry nodes handle a list of geometries. When the name changes, 
a new geom is generated. If an empty string is provided, the current geometry 
is deleted.

The parent and children connections reflect the state of the geometry selected 
in the combobox. If the parent or children geoms are currently selected in 
their own node: the connections are present. If the parent or a child is not 
visible (for the selected geometry) the parent port or/and the output port 
contour lines are white (instead of black). Before parenting you can set the 
behavior of the 'parent' inputport. Right clicking on the port allows you to 
extend the parenting to the sibling geometries present in the node or to all 
the geometries of the node. This behavior set by right clicking will only 
affect the next parenting. In the same way, deleting the parent connection, 
reparent to 'root' the selected geometry or the sibling geometries or all the 
geometries in the child node (depending on the right click menu choice).

Input Ports common to all GeometryNodes:
    new geom name: name of a new geometry to be created with the new data.
                   when the node runs, this new geometry will be selected
                   in the 'current geom' field and sent to the first output port.
    current geom: if 'new geom name' is empty, the new data coming from the
                  input ports are applied to the current geom.
                  at the end of the run the selected geometry is sent to
                  the first output port.
    instanceMatrices: stream of 4x4 transformation matrices (default identity)
    geomOptions: accept a dictionary of properties with values.
                 can be fed with node "SetGeom Options"
    parent:  a geometry to parent the selected geometry (default root)

Output Ports common to all GeometryNodes:
    first output port: currently selected geometry 
    allGeometries: list of all the geometries in the node

Specific Input Ports:
    image:  the Image to be applied on the sticker
"""

    def __init__(self, name='oneTexturedQuad', **kw):
        kw['name'] = name
        apply( GeometryNode.__init__, (self, 'oneTexturedQuad'), kw )

        ip = self.inputPortsDescr
        ip.append(datatype='image', required=False, name='image')

        # for backward compatibility
        # this make sure that the first ports in the node are the ones 
        # that use to be there before we introduced the GeometryNode class
        # (old network were saved with port indices instead of port names)
        self.rearrangePorts()

        code = """def doit(self, image,
name, geoms, instanceMatrices, geomOptions, parent):
    #print "doit", self.name
    GeometryNode.doit(self, name=name, geoms=geoms,
                      instanceMatrices=instanceMatrices,
                      geomOptions=geomOptions, parent=parent)

    if self.selectedGeomIndex is not None:
        g = self.geom()
        if image is not None:
            from DejaVu2.Texture import Texture
            lTexture = Texture(enable=1, image=image)            
            g.Set(texture=lTexture, textureCoords=( (0, 0), (1, 0), (1, 1), (0, 1)))
            g.vertexSet.texCoords.array = lTexture.getTextureCoordinatesAfterResizeRatio(
                                                         g.vertexSet.texCoords.array)
        self.outputData(oneTexturedQuad=g, allGeometries=self.geoms)
    else:
        self.outputData(oneTexturedQuad=None, allGeometries=self.geoms)
"""
        self.setFunction(code)


    def appendGeometry(self, name):
        self.removePreviousGeomWithSameName(name)
        #print "appendGeometry"
        w = 10
        h = 10
        from DejaVu2.IndexedPolygons import IndexedPolygons
        lIndexedPolygons = IndexedPolygons(
                       name,
                       vertices=[(0, 0, 0), (w, 0, 0), (w, h, 0), (0, h, 0)],
                       faces = ((0,1,2,3), )
                      )
        self.geoms.append(lIndexedPolygons)
        return 1 # num of geoms appended

import types
    
class DecimateGeom(GeometryNode):
    """Uses QSlim library to decimate a surface.
    Creates a QSlim model from the input geometry, decimates it
    and outputs a new geometry with new number of faces, vertices, colors.

The geometry nodes handle a list of geometries. When the name changes, 
a new geom is generated. If an empty string is provided, the current geometry 
is deleted.

The parent and children connections reflect the state of the geometry selected 
in the combobox. If the parent or children geoms are currently selected in 
their own node: the connections are present. If the parent or a child is not 
visible (for the selected geometry) the parent port or/and the output port 
contour lines are white (instead of black). Before parenting you can set the 
behavior of the 'parent' inputport. Right clicking on the port allows you to 
extend the parenting to the sibling geometries present in the node or to all 
the geometries of the node. This behavior set by right clicking will only 
affect the next parenting. In the same way, deleting the parent connection, 
reparent to 'root' the selected geometry or the sibling geometries or all the 
geometries in the child node (depending on the right click menu choice).

Input Ports common to all GeometryNodes:
    new geom name: name of a new geometry to be created with the new data.
                   when the node runs, this new geometry will be selected
                   in the 'current geom' field and sent to the first output port.
    current geom: if 'new geom name' is empty, the new data coming from the
                  input ports are applied to the current geom.
                  at the end of the run the selected geometry is sent to
                  the first output port.
    instanceMatrices: stream of 4x4 transformation matrices (default identity)
    geomOptions: accept a dictionary of properties with values.
                 can be fed with node "SetGeom Options"
    parent:  a geometry to parent the selected geometry (default root)

Output Ports common to all GeometryNodes:
    first output port: currently selected geometry 
    allGeometries: list of all the geometries in the node

Specific Input Ports:
      geometry - any DejaVu2.IndexedPolgons geometry (required);
      percent - percentage of faces to be retained in decimated model (required);
      targetfaces -  corresponding to percent number of faces;
      rebuild - (boolean), when true, rebuilds the QSlim model (takes
               current parameters of selected geometry).
"""

    def __init__(self, name = "decimateGeom", **kw):
        kw['name'] = name
        apply( GeometryNode.__init__, (self,"decimateGeom"), kw )
        self.widgetDescr['name']['initialValue'] = ""
        self.widgetDescr['percent'] = {
            'class':'NEDial', 'master':'ParamPanel', 'size':50,
            'oneTurn':100, 'min':0., 'max':100., 'type':'float',
            'showLabel':1, 'continuous':0,
            'initialValue':50.,
            'labelCfg':{'text':'percent'},
            }
        
        self.widgetDescr['targetfaces'] = {
            'class':'NEEntry', 'master':'ParamPanel', 'width':10,
            'labelCfg':{'text':'num of faces:'},
            }

        self.widgetDescr['newgeom'] = {
                'class':'NECheckButton', 'initialValue':0,
                'master':'ParamPanel', 
                'labelGridCfg':{'sticky':'we'},
                'labelCfg':{'text':'create new geometry'},
                }
        
        self.widgetDescr['rebuild'] = {
            #'class':'NECheckButton',
            'class':'NEButton',
            'master':'ParamPanel',
            'labelGridCfg':{'sticky':'we'},
            'labelCfg':{'text':'rebuild Qslim model'}
            }
        
        ip = self.inputPortsDescr
        ip.append(datatype='geom', name='geometry')
        ip.append(name='percent', datatype='float')
        ip.append(name='targetfaces', datatype = "int", required=False)
        ip.append (name='rebuild', datatype='boolean', required=False)
        ip.append(datatype='None', name='u',  required=False)
        ip.append(datatype='None', name='v',  required=False)

        op = self.outputPortsDescr
        #op.append(datatype='geom', name='geometry')
        op.append(datatype='None', name='u')
        op.append(datatype='None', name='v')
        from DejaVu2.VisionInterface.DejaVu2Nodes import QSlimDecimate
        self.model = QSlimDecimate()
        self.geometry = None      # Current geometry
        
        self.rearrangePorts()
        
        code = """def doit(self, geometry , percent, targetfaces, rebuild, u, v,
name, geoms, instanceMatrices, geomOptions, parent):
    #print 'DecimateGeom: name=', name, 'percent=', percent, 'targetfaces=', targetfaces, 'rebuild=', rebuild, self.inputPortByName['rebuild'].hasNewValidData()
    # DecimateGeom node can only handle 1 geom so far, so print out a message if there
    # are more comming in
    if type(geometry)==types.ListType:
        if len(geometry)>1:
            warnings.warn('Only first geometry is being processed by DecimateGeom')
        geometry = geometry[0]

    # check that the incomming geom can be represented as IndexedPolygons
    if not geometry.asIndexedPolygons(run=0):
       warnings.warn(geometry.name , 'DecimateGeom Node: %s can not be represented as IndexedPolygons', geometry.name)
       return

    if name == '':
        #if self.selectedGeomIndex is not None:
        #    name = self.geom().name
        selInd = self.getInputPortByName('geoms').widget.widget.curselection()
        if len(selInd):
            self.selectedGeomIndex = int(selInd[0])
            name = self.geom().name
            if name != self.name: self.name = name
        else:
            name = geometry.name    
    
    # if there is new data on the first input port (geom) rebuild the model
    #if self.inputPortByName['geometry'].hasNewValidData() or \
       #(self.inputPortByName['rebuild'].hasNewValidData() and rebuild):
    if self.inputPortByName['geometry'].hasNewValidData() or \
        self.inputPortByName['rebuild'].hasNewValidData():
       if len(geometry.vertexSet.vertices):
           print 'DecimateGeom Node BUILDING MODEL'
           self.model.build_model(geometry, u, v)
       else:
           return

    from DejaVu2.IndexedPolygons import IndexedPolygons
    # compute the number of triangles to be kept
    maxFaces = self.model.maxFaces
    if self.inputPortByName['percent'].hasNewValidData(): #percent has changed
        targetfaces = int(maxFaces*percent/100.)
    elif self.inputPortByName['targetfaces'].hasNewValidData() and targetfaces: 
        #targetfaces has changed
        targetfaces = int(targetfaces)
    else: # use the percentage
        targetfaces = int(maxFaces*percent/100.)
        #print 'ELSE', maxFaces*percent, targetfaces

    # make sure targetfaces is OK
    if targetfaces > maxFaces:
        targetfaces = maxFaces
        if self.inputPortByName['percent'].widget:
            self.inputPortByName['percent'].widget.set(100, 0)
        if self.inputPortByName['targetfaces'].widget:
            self.inputPortByName['targetfaces'].widget.set(maxFaces, 0)
            
    if self.model:
        decimVerts, decimFaces, decimNormals, decimColors, decimTex = self.model.decimate(targetfaces)
        #print 'creating geometry:', name 
        self.geometry = IndexedPolygons(name, vertices=decimVerts, faces=decimFaces,
                                        vnormals=decimNormals)
        numFaces = len(decimFaces)
        nvert = len(decimVerts)    
        clen = len(decimFaces)
        if self.model.newcolors is not None:
            if self.model.colorBind == 'vertex':
                clen = nvert
            #print 'len colors:', clen
            if clen > 0:
                self.geometry.Set(materials=decimColors[:clen],
                         inheritMaterial=False, redo=True)
            else:
                self.geometry.Set(inheritMaterial=False, redo=True)

        # set the target face widget to the actual number of faces after
        # decimation
        if self.inputPortByName['targetfaces'].widget:
            self.inputPortByName['targetfaces'].widget.set(numFaces, run=0)

        realpercent = 100.*(float(numFaces)/maxFaces)
        if percent != realpercent and self.inputPortByName['percent'].widget:
            self.inputPortByName['percent'].widget.set(realpercent, run=0)

        u = None
        v = None
        if self.model.newtexcoords is not None:
            if hasattr(geometry.vertexSet, 'texCoords'):
                self.geometry.Set(textureCoords=decimTex[:nvert])
            u = decimTex[:nvert,0]
            v = decimTex[:nvert,1]  

        GeometryNode.doit(self, name=name, geoms=geoms,
                      instanceMatrices=instanceMatrices,
                      geomOptions=geomOptions, parent=parent)
                      
        if self.selectedGeomIndex is not None:
            g  = self.geom()
            self.outputData(decimateGeom=g, allGeometries=self.geoms,
                            u = u, v = v )
        else:
            self.outputData(decimateGeom=None, allGeometries=self.geoms, u=u, v=v)

"""
        self.setFunction(code)

    def appendGeometry(self, name):
        #print "append geometry, name:", name
        #print "self.geoms:", self.geoms
        
        self.removePreviousGeomsStartingWith(name)
        self.geoms.append(self.geometry)
        return 1





class ConnectedComponents(GeometryNode):
    """ Node to split connected surfaces.

The geometry nodes handle a list of geometries. When the name changes, 
a new geom is generated. If an empty string is provided, the current geometry 
is deleted.

The parent and children connections reflect the state of the geometry selected 
in the combobox. If the parent or children geoms are currently selected in 
their own node: the connections are present. If the parent or a child is not 
visible (for the selected geometry) the parent port or/and the output port 
contour lines are white (instead of black). Before parenting you can set the 
behavior of the 'parent' inputport. Right clicking on the port allows you to 
extend the parenting to the sibling geometries present in the node or to all 
the geometries of the node. This behavior set by right clicking will only 
affect the next parenting. In the same way, deleting the parent connection, 
reparent to 'root' the selected geometry or the sibling geometries or all the 
geometries in the child node (depending on the right click menu choice).

Input Ports common to all GeometryNodes:
    new geom name: name of a new geometry to be created with the new data.
                   when the node runs, this new geometry will be selected
                   in the 'current geom' field and sent to the first output port.
    current geom: if 'new geom name' is empty, the new data coming from the
                  input ports are applied to the current geom.
                  at the end of the run the selected geometry is sent to
                  the first output port.
    instanceMatrices: stream of 4x4 transformation matrices (default identity)
    geomOptions: accept a dictionary of properties with values.
                 can be fed with node "SetGeom Options"
    parent:  a geometry to parent the selected geometry (default root)

Output Ports common to all GeometryNodes:
    first output port: currently selected geometry 
    allGeometries: list of all the geometries in the node

Specific Input Ports:
          -geometry or a list of geometries(the first object in the list
            will be used as input)
          - onegeom (check button) - if checked, the largest component is output

          -name of output geometries . Index (0-n) will be attached to the name
"""
    
    def __init__(self, name='ConnectedComponents', **kw):
        
        kw['name'] = name
        apply( GeometryNode.__init__, (self, 'ConnectedComponents'), kw )
        self.geometry = None
        self.outgeoms = []
        
        ip = self.inputPortsDescr
        ip.append(name='geometry', datatype='geom', required=True)
        ip.append(name = 'onegeom', datatype='int', required=False, defaultValue=0)
        
        op = self.outputPortsDescr
        #op.append(name='outgeoms', datatype='geom')
        
        wd = self.widgetDescr
        wd['onegeom'] = {
            'class':'NECheckButton', 'master':'node',
            'initialValue':0,
            'labelCfg':{'text':'Output one component:'}
            }
        self.rearrangePorts()
        
        code = """def doit(self, geometry, onegeom,
name, geoms, instanceMatrices, geomOptions, parent):
        #print 'ConnectedComponents', geometry, 'onegeom=', onegeom , 'name=',  name, 'geoms=', geoms, 'geomOptions=', geomOptions
        if type(geometry)==types.ListType:
            if len(geometry):
                geometry = geometry[0]

        from DejaVu2.IndexedPolygons import IndexedPolygons
        outgeoms = []
        #from time import time
        #t1 = time()
        newfaces , newverts = self.findComponents(geometry)
        #t2 = time()
        #print 'time to find %d connected components : %.2f'%(len(newfaces), t2-t1)
        
        if name == '':
            name = geometry.name
        if onegeom == 1:
            # find the largest set of faces:
            maxn = 0
            maxind = None
            for i, newfs  in enumerate(newfaces):
                lenf = len(newfs)
                if lenf > maxn:
                    maxn = lenf
                    maxind = i
            #print 'len faces of %s is %d'%(name, len(newfaces[maxind]))
            outgeoms.append(IndexedPolygons(name, vertices = newverts[maxind], faces = newfaces[maxind]))
        else:
          for i, newfs  in enumerate(newfaces):
            newname = '%s%d' % (name, i)
            #print 'len faces of %s is %d'%(newname, len(newfs))
            outgeoms.append(IndexedPolygons(newname, vertices = newverts[i], faces = newfs))
                   
        self.outgeoms = outgeoms
        GeometryNode.doit(self, name=name, geoms=geoms, instanceMatrices=instanceMatrices,
                          geomOptions=geomOptions, parent=parent)

        self.geometry = geometry
        if self.selectedGeomIndex is not None:
            g = self.geom()
            #GeometryNode.textureManagement(self, image=image)
            self.outputData(ConnectedComponents=g, allGeometries=self.geoms)
        else:
            self.outputData(ConnectedComponents=None, allGeometries=self.geoms)

"""
        
        self.setFunction(code)

    def appendGeometry(self, name):
        #print "append geometry, name:", name
        
        self.removePreviousGeomsStartingWith(name)
        # append the geoms
        if self.outgeoms is not None:
            self.removePreviousGeomsWithSameName(self.outgeoms)
            self.geoms.extend(self.outgeoms)
            return len(self.outgeoms)
        else:
            return 0



    def findComponents(self, geometry):
        faces = geometry.getFaces()
        verts = geometry.getVertices()
        fdict = {} 
        vdict = {} #dictionary with key - vertex index,
                   #value - list of face indices in which the vertex is found

        flag1 = True; flag2 = True
        newfaces = []; newverts = []
        while flag2:
            for i, fs in enumerate(faces):
                for v in fs:
                    if not vdict.has_key(v):
                        vdict[v] = [i]
                    else:
                        vdict[v].append(i)
                fdict[i] = fs
            Vco = faces[0][:]
            newfaces1 = []; newverts1 = []
            vertinds = {} # keys - vertex indices from the input verts list
                          # values - new vertex indices of current surface
            vcount = 0
            # find a surface
            while flag1:
                _Vco = []
                flag1 = False
                # find all vertices that share the same triangles with the vertices in Vco.
                for vert in Vco:
                    vfs = vdict[vert]
                    for i in vfs:
                        if fdict.has_key(i):
                            flag1 = True
                            fs = fdict.pop(i)

                            fsnew = [] # remapped face (with new vertex idices)
                            for v in fs:
                                if v not in Vco: 
                                    if v not in _Vco:
                                        _Vco.append(v)
                                if not  vertinds.has_key(v):
                                    vertinds[v] = vcount
                                    newverts1.append(verts[v])
                                    fsnew.append(vcount) 
                                    vcount = vcount + 1
                                else:
                                    fsnew.append(vertinds[v])
                            newfaces1.append(fsnew) # add found triangle to the list of triangles of current surface

                Vco  = _Vco
            newfaces.append(newfaces1)
            newverts.append(newverts1)
            if  len(fdict):
                faces = fdict.values()
                fdict = {}
                vdict = {}
                flag1 = True
            else: flag2 = False
        return newfaces, newverts

