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
# Date: July 2002 Authors: Daniel Stoffler, Michel Sanner
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
# $Header: /mnt/raid/services/cvs/DejaVu2/VisionInterface/DejaVuNodes.py,v 1.1.1.1.4.2 2017/11/08 23:27:21 annao Exp $
#
# $Id: DejaVuNodes.py,v 1.1.1.1.4.2 2017/11/08 23:27:21 annao Exp $
#

import types
import warnings
import string
import numpy.oldnumeric as Numeric, sys, os
from math import sqrt, atan, pi

from NetworkEditor.items import NetworkNode
from NetworkEditor.macros import MacroNode
from NetworkEditor.macros import MacroOutputNode
from Vision import UserLibBuild

from DejaVu2.VisionInterface.DejaVu2Widgets import NEColorMap, NEColorWheel, NEColorEditor, NEDejaVu2GeomOptions
from DejaVu2.Geom import Geom
from DejaVu2.Insert2d import Insert2d
from DejaVu2.colorMap import ColorMap


def importSymLib(net):
    try:
        from symserv.VisionInterface.SymservNodes import symlib
        net.editor.addLibraryInstance(
            symlib, 'symserv.VisionInterface.SymservNodes', 'symlib')
    except:
        warnings.warn(
            'Warning! Could not import symlib from symserv.VisionInterface.SymservNodes.py', stacklevel=2)


def importImageLib(net):
    try:
        from Vision.PILNodes import imagelib
        net.editor.addLibraryInstance(imagelib, 'Vision.PILNodes', 'imagelib')
    except:
        warnings.warn(
            'Warning! Could not import imagelib from Vision.PILNodes.py', stacklevel=2)
        

def importVolLib(net):
    try:
        from Volume.VisionInterface.VolumeNodes import vollib
        net.editor.addLibraryInstance(
            vollib, 'Volume.VisionInterface.VolumeNodes', 'vollib')
    except:
        warnings.warn(
            'Warning! Could not import vollib from Volume.VisionInterface.VolumeNodes.py')


class GeomOptions(NetworkNode):
    def __init__(self, name='Set Geom Options', viewer=None, **kw):
        kw['name'] = name
        apply( NetworkNode.__init__, (self,), kw )

        ip = self.inputPortsDescr
        ip.append(name='geomOptions', datatype='dict')

        self.widgetDescr['geomOptions'] = {
            'class':'NEDejaVu2GeomOptions', 'lockedOnPort':True
            }
        
	op = self.outputPortsDescr
        op.append(datatype='dict', name='geomOptions')

        code = """def doit(self, geomOptions):
    self.outputData(geomOptions=geomOptions)       
"""
		
	self.setFunction(code)


class RestoreState(NetworkNode):
    """Set a Viewer's state from a file

Input Ports:
  Viewer:   (required) the Viewer for which to restore the state
  filename: (required) the file containing the state
  mode:     (optional) can be 'Viewer', 'objects' or 'both'
             defaults to 'both'.  When this value is 'viewer',
             the state is restored only for objects always present in the
             viewer (i.e. root, lights, clipping planes, fog, etc).
             'objects': means restoring the state of geometries in the viewer.
             'both': does both.
           
Output Ports:
"""
    from DejaVu2.Geom import Geom

    def __init__(self, name='RestoreState', viewer=None, **kw):
        kw['name'] = name
        apply( NetworkNode.__init__, (self,), kw )

        ip = self.inputPortsDescr
        ip.append(name='viewer', datatype='viewer')
        ip.append(name='filename', datatype='string')
        ip.append(name='mode', datatype='string', required=False, defaultValue='both')

        self.widgetDescr['filename'] = {
            'class':'NEEntryWithFileBrowser', 'master':'node', 'width':10,
            'initialValue':'',
            'filetypes':[('state files','*_state.py'),('all','*')],
            'labelCfg':{'text':'Filename: '}
            }

        self.widgetDescr['mode'] = {
            'class':'NEComboBox', 'master':'node',
            'choices':['viewer', 'objects', 'both'],
            'initialValue': 'both',
            'entryfield_entry_width':8,
            'labelCfg':{'text':'mode:'},
            }

        code = """def doit(self, viewer, filename, mode):
    execfile(filename, {'vi':viewer, 'mode':mode})
    viewer.Redraw()
"""
        self.setFunction(code)


class CenterOnPickedVertex(NetworkNode):
    """Set the rotation center of for the scene to the picked vertex

Input Ports:
  Viewer:   (required) the Viewer for which to restore the state

Output Ports:
"""
    from DejaVu2.Geom import Geom

    def __init__(self, name='CenterOnPick', viewer=None, **kw):
        kw['name'] = name
        apply( NetworkNode.__init__, (self,), kw )

        ip = self.inputPortsDescr
        ip.append(name='viewer', datatype='viewer')

        code = """def doit(self, viewer):
    if viewer.lastPick is None:
        return
    hits = viewer.lastPick.hits
    if len(hits)==1:
        g = hits.keys()[0]
        index = hits[g][0][0]
        point = g.vertexSet.vertices.array[index]
        viewer.rootObject.SetPivot(point)"""

        self.setFunction(code)


class CenterOnVertex(NetworkNode):
    """Set the rotation center of for the scene to a vertex

Input Ports:
  Viewer: (required) the Viewer for which to restore the state
  vertex:
  
Output Ports:
"""
    from DejaVu2.Geom import Geom

    def __init__(self, name='CenterOnPick', viewer=None, **kw):
        kw['name'] = name
        apply( NetworkNode.__init__, (self,), kw )

        ip = self.inputPortsDescr
        ip.append(name='viewer', datatype='viewer')
        ip.append(name='point', datatype='coord3')

        code = """def doit(self, viewer, point):
        viewer.rootObject.SetPivot(point)"""

        self.setFunction(code)


class Viewer(NetworkNode):
    """Create an instance of a DejaVu2 Viewer object.
This object provides a full fledged 3D geometry viewer with support for
the OpenGL material and lighting model with multiple light source,
arbitrary clipping planes, a hierarchy of geometries with property
inheritance, a material editor, etc... .

Input Ports:
  Geometries: (required) Accepts any object that is a DejaVu2.Geom instance
              or a list thereof from each parent port and adds the geometry
              objects to the viewer
Output Ports:
  lastPick: outputs a DejaVu2.Camera.PickObject
  DejaVu2: outputs the DejaVu2 Viewer instance
  Redraw: fires at each redraw event
"""
    from DejaVu2.Geom import Geom

    def __init__(self, name='Viewer', viewer=None, **kw):
        kw['name'] = name
        apply( NetworkNode.__init__, (self,), kw )
        #self.readOnly = 1
        self.vi = viewer

        # code for input port callback "beforeDisconnect"
        # this will delete geometries added to the viewer
        codeBeforeDisconnect = """def beforeDisconnect(self, c):
    #print "Viewer Node beforeDisconnect"
    self.node.deleteGeometries(c)
"""        
        # add an entry so the show parma. panel will be enabled
        ip = self.inputPortsDescr
        ip.append(name='geometries', datatype='geom(0)', required=False,
                  singleConnection=False,
                  beforeDisconnect=codeBeforeDisconnect)

        op = self.outputPortsDescr
        op.append(datatype='None', name='lastPick')
        op.append(datatype='viewer', name='dejaVuViewer')
        op.append(datatype='None', name='redraw')

        code = """def doit(self, geometries):
    self.addGeometriesToViewer(geometries)
    
    self.outputData(lastPick=self.vi.lastPick)
    self.outputData(dejaVuViewer=self.vi)
    self.outputData(redraw=0)
"""
        self.setFunction(code)


    def getStateDefinitionCode(self, nodeName, indent=''):
        #print "ViewerNode getStateDefinitionCode"
        return self.vi.getViewerStateDefinitionCode(
                    viewerName=nodeName+'.vi', 
                    indent=indent,
                    withMode=False)


    def addGeometriesToViewer(self, geometries):
        """ this geometries a list, but not a list from a geom node
"""
        if self.vi is None:
            return

        if geometries:
            for g in geometries:
                self.addGeometryiesToViewer(g)


    def addGeometryiesToViewer(self, geometryies):
        """ we don't know if geometryies is a Geom or a list from a geom node
"""                            
        if geometryies and geometryies != [None]:
            if isinstance(geometryies, Geom) or isinstance(geometryies, Insert2d):  
                if hasattr(geometryies, 'node') and \
                   hasattr(geometryies.node(),'geoms'):
                    self.addMultipleGeomsToViewer(geometryies.node().geoms)
                else: 
                    self.addSingleGeomToViewer(geometryies)
                    self.addGeometryiesToViewer(geometryies.parent)
            elif isinstance(geometryies, list):
                self.addMultipleGeomsToViewer(geometryies)
            else:
                assert False
                    
        self.vi.Redraw()


    def addMultipleGeomsToViewer(self, geometries ):
        """ geometries is a list from a geom node
"""
        #print "addMultipleGeomsToViewer", geometries
        for g in geometries:
            if isinstance(g, Geom) or isinstance(g, Insert2d):  
                self.addSingleGeomToViewer(g)
                if g.parent.viewer is None:
                    self.addGeometryiesToViewer(g.parent)
            elif isinstance(g, list):
                self.addMultipleGeomsToViewer(g)


    def addSingleGeomToViewer(self, aGeom ):
        """ geometries is a Geom
"""
        #print "addSingleGeomToViewer", aGeom.name

        #if aGeom is None or aGeom.viewer is not None:
        if aGeom.viewer is not None:
            return

        # add the geom to the viewer
        lHighestModifiedParent = \
          self.vi.AddObjectAndParents(aGeom, parent=aGeom.parent, local=False)

        # rename node geometry in Vision and set input for name
        # FIXME this test has to go after all Vision geoms have .node()
        if hasattr(aGeom, 'node') :
            aGeom.node().ensureNameOfNodeAndDescendants(lHighestModifiedParent)


    def deleteGeometries(self, connection):
        """code for input port callback "beforeDisconnect"
this will delete geometries added to the viewer
"""
        #print "Viewer Node deleteGeometries"
        #import pdb;pdb.set_trace()
        if self.vi is not None:   
            self.vi.SetCurrentObject(self.vi.rootObject)

        if connection.port1.data:
            geoms = connection.port1.data
            # geoms can be packed in a list or not
            import types
            from DejaVu2.Geom import Geom
            if type(geoms) == types.ListType:
                for g in geoms:
                    if self.vi is None:
                        g.viewer = None  
                    else:
                        if hasattr(g,'node'):
                            g.node().removeViewerConnections(connection)
                            g.node().removeAndCleanMultipleGeometryies(g)
            else: # not list
                if self.vi is None:
                    geoms.viewer = None  
                else:
                    if hasattr(geoms,'node'):
                        geoms.node().removeViewerConnections(connection)                                                               
                        geoms.node().removeAndCleanMultipleGeometryies(geoms)       


    def buildIcons(self, canvas, posx, posy):
        apply( NetworkNode.buildIcons, (self, canvas, posx, posy) )
        self.paramPanel.applyButton.forget()
        from DejaVu2 import Viewer
        from DejaVu2.Geom import Geom
        
        if self.vi is None:
            self.vi = Viewer(nogui=1, guiMaster=self.paramPanel.widgetFrame)
            # do not allow user to kill this window
            self.vi.master.protocol('WM_DELETE_WINDOW', self.doNothing)
            self.ownViewer = 1
        else:
            #self.vi = viewer
            self.ownViewer = 0

        camera = self.vi.currentCamera
        # add callback to output picking events
        self.vi.AddPickingCallback(self.handlePick)
        # add callback to output redraw event
        self.vi.afterRedraw = self.afterRedraw
        
        self.outputData(dejaVuViewer=self.vi)


    def doNothing(self):
        # used for overwriting WM_DELETE_WINDOW of Viewer window
        return
    

    def afterRedraw(self):
        if self.vi:
            apply( self.vi.__class__.afterRedraw, (self.vi,))
            # trigger all children of the redraw outputPort
            if len(self.outputPorts):
                self.outputData(redraw=1)
                self.scheduleChildren(portList=[self.outputPorts[2]])
        
        
    def handlePick(self, pick):
        if pick:
            if self.vi:
                self.outputData(lastPick=self.vi.lastPick)
            # schedule children of lastPick
            self.scheduleChildren(portList=[self.outputPorts[0]])


    def beforeRemovingFromNetwork(self):
        NetworkNode.beforeRemovingFromNetwork(self)

        for g in self.vi.rootObject.AllObjects():
            if hasattr(g,'node'):
                g.node().removeSingleGeom(g)

        if self.ownViewer:
            self.vi.Exit()
        self.vi = None
        

    def toggleNodeExpand_cb(self, event=None):
        # overwrite base class method, because this panel is special
        if self.paramPanel.visible:
            self.paramPanel.hide()
        else:
            self.paramPanel.show()



class AccumPickedVertices(NetworkNode):
    """Class to accumulate picked vertices

Input Ports:
    pick:  DejaVu2 pick object instance
    viewer: DejaVu2 viewer object instance
    reset:  resets the list of vertices to empty
    remove: check button to remove vertices form the list

Output Ports:
    pickedVertices: list of 3D vertices
"""
    def __init__(self, **kw):
        apply( NetworkNode.__init__, (self,), kw )
           
        self.vertices = []
        
         # add an entry so the show parma. panel will be enabled
        ip = self.inputPortsDescr
        ip.append(name='pick', datatype='None')
        ip.append(name='viewer', datatype='viewer')
        ip.append(name='reset', datatype='boolean')
        ip.append(name='remove', datatype='boolean')

        self.widgetDescr['reset'] = {
            'class':'NEButton',
            'master':'node',
            'initialValue':False,
            'labelGridCfg':{'sticky':'we'},
            'labelCfg':{'text':'empty list'},
            }

        self.widgetDescr['remove'] = {
            'class':'NECheckButton',
            'master':'node',
            'initialValue':False,
            'labelGridCfg':{'sticky':'we'},
            'labelCfg':{'text':'remove'},
            }

	op = self.outputPortsDescr
        op.append(datatype='list', name='pickedVertices')

        code = """def doit(self, pick, viewer, reset, remove):
    if self.inputPortByName['reset'].hasNewValidData():
        self.vertices = []

    if self.inputPortByName['pick'].hasNewValidData():
        vertices = viewer.transformedCoordinatesWithInstances(pick.hits)
        if remove:
            for v1 in vertices:
                vstr = '%.3f %.3f %.3f'%tuple(v1)
                for v2 in self.vertices:
                    if '%.3f %.3f %.3f'%tuple(v2)==vstr:
                        self.vertices.remove(v2)
        else:
            self.vertices.extend(vertices)

    self.outputData(pickedVertices=self.vertices)
"""

        self.setFunction(code)



class PolyhedronVolumeArea(NetworkNode):
    """compute the numerical volume and surface area of a polyhedron

Input Ports:
  geometry: (required) Accepts any IndexedPolygon geometry from DejaVu2

Output Ports:
  volume:  float
  area:    float
  compactness: area/volume (float)
."""

    def triangleArea(self, p1, p2, p3):
        """Compute the surface area of a triangle.
    """
        x1,y1,z1 = p1
        x2,y2,z2 = p2
        x3,y3,z3 = p3
        dx, dy, dz = x1-x2, y1-y2, z1-z2
        a = sqrt( dx*dx + dy*dy + dz*dz )
        dx, dy, dz = x2-x3, y2-y3, z2-z3
        b = sqrt( dx*dx + dy*dy + dz*dz )
        dx, dy, dz = x1-x3, y1-y3, z1-z3
        c = sqrt( dx*dx + dy*dy + dz*dz )
        s = .5*(a+b+c)
        area = s*(s-a)*(s-b)*(s-c)
        if area <= 0.:
            return 0.
        return sqrt(area)


    def meshVolume(self, verts, tri, norm):
        """Compute the Volume and surface area of a mesh specified by vertices,
    indices of triangular faces and face normals
    """
        assert len(tri)==len(norm)
        volSum = 0.0
        areaSum = 0.0
        oneThird = 1./3.
        for t,n in zip(tri, norm):
            s1 = verts[t[0]]
            s2 = verts[t[1]]
            s3 = verts[t[2]]
            area = self.triangleArea(s1,s2,s3)
            areaSum += area

            g = [ (s1[0]+s2[0]+s3[0])*oneThird,
                  (s1[1]+s2[1]+s3[1])*oneThird,
                  (s1[2]+s2[2]+s3[2])*oneThird ]
            volSum += (g[0]*n[0] + g[1]*n[1] + g[2]*n[2])*area
        return volSum*oneThird, areaSum


    def __init__(self, name='PolyhedronVolumeArea', **kw):
        kw['name'] = name
        apply( NetworkNode.__init__, (self,), kw )
        #self.readOnly = 1

        ip = self.inputPortsDescr
        ip.append(datatype='geom', name='geometry')

        op = self.outputPortsDescr
        op.append(datatype='float', name='volume')
        op.append(datatype='float', name='area')
        op.append(datatype='float', name='compactness')

        code = """def doit(self, geometry):
    vert = geometry.getVertices()
    tri = geometry.getFaces()
    #from opengltk.extent.utillib import glTriangleNormals
    from geomutils.geomalgorithms import  TriangleNormals
    norm = TriangleNormals( vert, tri, 'PER_FACE')
    # FIXME getFNormals() returns old set of normals after sendity has changed
    # on MSMS surface
    # norm = geometry.getFNormals()

    vol, area = self.meshVolume(vert, tri, norm)
    self.outputData(volume=vol, area=area, compactness=area/vol)
"""

        self.setFunction(code)


class writeIndexedPolygon(NetworkNode):
    """write the vertices, normals and faces of an IndexedPolygon to a file.

Input Ports:
  geometry: (required) Accepts any IndexedPolgon geometry from DejaVu2
  filename: (required) string, no extension. filename.vert and filename.face
                       will be created
Output Ports:
."""
    
    def __init__(self, name='writeIndexedPolygon', **kw):
        kw['name'] = name
        apply( NetworkNode.__init__, (self,), kw )
        #self.readOnly = 1

        self.widgetDescr['filename'] = {
            'class':'NEEntryWithFileSaver', 'master':'node', 'width':10,
            'initialValue':'',
            'labelCfg':{'text':'Filename: '}
            }

        ip = self.inputPortsDescr
        ip.append(datatype='geom', name='geometry')
        ip.append(datatype='string', name='filename')

        code = """def doit(self, geometry, filename):

    from DejaVu2.IndexedPolygons import IndexedPolygons
    assert isinstance(geometry, IndexedPolygons)
    geometry.writeToFile(filename)
"""

        self.setFunction(code)


class writeCurvPly(NetworkNode):
    """write the vertices and triangles of an IndexedPolygon to a file in the
format used by the setCurvature program (which ressembles ply2)

Input Ports:
    geometry: (required) Accepts any IndexedPolgon geometry from DejaVu2
    filename: (required) string, no extension. filename.vert and filename.face
                         will be created
Output Ports:
    filename:  outputs the filename after it wrote the file
"""
    def writeFile(self, filename, verts, faces, neighborhoodSize, crestLines):
        f = open(filename, 'w')
        
        # this assert just because the test bellow commented out was weird
        assert crestLines in [0, 1, True, False]
        
        f.write('%d\n%d\n%d\n%d\n'%( len(verts), len(faces), neighborhoodSize,
                                     crestLines ))
        #                             crestLines==True ))
                                     
                                     
                                     
        dum = map(lambda x, f=f: f.write('%f\n'%x), verts)
        dum = map(lambda x, f=f: f.write('%d\n'%x), faces)
        f.close()
        
    def __init__(self, name='writeCurvPly', **kw):
        kw['name'] = name
        apply( NetworkNode.__init__, (self,), kw )
        #self.readOnly = 1

        self.widgetDescr['filename'] = {
            'class':'NEEntryWithFileSaver',
            'master':'node', 'width':10,
            'initialValue':'',
            'labelCfg':{'text':'Filename: '}
            }

        self.widgetDescr['neighborhoodSize'] = {
            'class':'NEThumbWheel',
            'master':'node',
            'initialValue':4,
            'labelCfg':{'text':'neighborhoodSize:'},
            'width':80, 'height':20, 
            'type':'int', 'oneTurn':10, 'lockBMin':1, 'min':1,
            'wheelPad':2 }

        self.widgetDescr['crestLines'] = {
            'class':'NECheckButton',
            'master':'node',
            'initialValue':False,
            'labelGridCfg':{'sticky':'we'},
            'labelCfg':{'text':'crest Lines:'},
            }

        ip = self.inputPortsDescr
        ip.append(datatype='geom', name='geometry')
        ip.append(datatype='string', name='filename')
        ip.append(datatype='int', name='neighborhoodSize')
        ip.append(datatype='boolean', name='crestLines')

        op = self.outputPortsDescr
        op.append(datatype='string', name='filename')

        code = """def doit(self, geometry, filename, neighborhoodSize,
crestLines):

    from DejaVu2.IndexedPolygons import IndexedPolygons
    assert isinstance(geometry, IndexedPolygons)
    verts = geometry.getVertices()
    faces = geometry.getFaces()
    self.writeFile(filename, verts, faces, neighborhoodSize, crestLines)
    geometry.writeToFile(filename)
"""

        self.setFunction(code)


class Curvature(NetworkNode):
    """run the setCurvature program

Input Ports:
    geometry: (required) Accepts any IndexedPolgon geometry from DejaVu2
    neighborhoodSize: size of the neighborhood use to fit quadric
    crestLines: when set to ture crest liens are traced

Output Ports:
    maxCurv:       max curvature for each vertex
    minCurv:       min curvature for each vertex
    vmaxCurv:      max curvature vector 
    vminCurv:      min curvature vector 
    GaussianCurv:  Gaussian curvature vector (Kmax*Kmin)
    meanCurv:      mean curvature vector for (Kmax+Kmin)
    shapeIndex:    shapeIndex -2/pi*arctan(Kmax+Kmin/Kmax-Kmin)
    curvedness:    curvedness sqrt(Kmax**2 + Kmin**2)/2
    normals:       normals to the surface
"""
    def writeFile(self, filename, verts, faces, neighborhoodSize, crestLines):
        f = open(filename, 'w')
        
        # this assert just because the test bellow commented out was weird
        assert crestLines is True or crestLines is False
        
        f.write('%d\n%d\n%d\n%d\n'%( len(verts)/3, len(faces)/3,
                                     neighborhoodSize,
                                     crestLines ))
                                     #crestLines==True ))
                                     
        dum = map(lambda x, f=f: f.write('%f\n'%x), verts)
        dum = map(lambda x, f=f: f.write('%d\n'%x), faces)
        f.close()

        
    def readCurvature(self, filename):
        f = open(filename)
        data = f.readlines()
        f.close()
        nbv = int(data[0])
        nbt = int(data[1])
        import string
        dataf = map(string.split, data[2:])
        curvData = []
        for l in dataf:
            curvData.append(map(float, l))
        curvData = Numeric.array(curvData, 'f')
        return curvData
    
        
    def __init__(self, name='curvature', **kw):
        kw['name'] = name
        apply( NetworkNode.__init__, (self,), kw )

        self.widgetDescr['neighborhoodSize'] = {
            'class':'NEThumbWheel',
            'master':'node',
            'initialValue':4,
            'labelCfg':{'text':'neighborhoodSize:'},
            'width':80, 'height':20, 
            'type':'int', 'oneTurn':10, 'lockBMin':1, 'min':1,
            'wheelPad':2 }

        self.widgetDescr['crestLines'] = {
            'class':'NECheckButton',
            'master':'node',
            'initialValue':False,
            'labelGridCfg':{'sticky':'we'},
            'labelCfg':{'text':'crest Lines:'},
            }

        ip = self.inputPortsDescr
        ip.append(datatype='geom', name='geometry')
        ip.append(datatype='int', name='neighborhoodSize')
        ip.append(datatype='boolean', name='crestLines')

        op = self.outputPortsDescr
        op.append(datatype='list', name='maxCurv')
        op.append(datatype='list', name='minCurv')
        op.append(datatype='list', name='vmaxCurv')
        op.append(datatype='list', name='vminCurv')
        op.append(datatype='list', name='gaussian')
        op.append(datatype='list', name='mean')
        op.append(datatype='list', name='shapeIndex')
        op.append(datatype='list', name='curvedness')
        op.append(datatype='list', name='normals')

        code = """def doit(self, geometry, neighborhoodSize, crestLines):

    from DejaVu2.IndexedPolygons import IndexedPolygons
    assert isinstance(geometry, IndexedPolygons)
    verts = Numeric.reshape(geometry.getVertices(), (-1,))
    faces = Numeric.reshape(geometry.getFaces(), (-1,))
    self.writeFile('curv.txt' , verts, faces, neighborhoodSize, crestLines)
    os.system('setCurvature curv.txt curv.out')
    curvData = self.readCurvature('curv.out')
    #os.system('rm curv.out curv.txt')
    Kmax = curvData[:,0]
    Kmin = curvData[:,1]
    gaussian = Kmax*Kmin
    mean = (Kmax+Kmin)*0.5
    shapeIndex = (-2/pi)*Numeric.arctan((Kmax+Kmin)/(Kmax-Kmin))
    curvedness = Numeric.sqrt(Kmax*Kmax+Kmin*Kmin)*0.5
    self.outputData(maxCurv=Kmax,
                    minCurv=Kmin, 
                    vmaxCurv=curvData[:,2:5].tolist(),
                    vminCurv=curvData[:,5:8].tolist(),
                    gaussian=gaussian, mean=mean,
                    shapeIndex=shapeIndex, curvedness=curvedness,
                    normals=curvData[:,8:].tolist()
                    )
"""
        self.setFunction(code)



class readIndexedPolygon(NetworkNode):
    """ ************* DEPRECATED, USE GeomsFromFile instead **************
    
    read the vertices, normals and faces of an IndexedPolygon from a file.

Input Ports:
  filename: (required) string, no extension. filename.vert and filename.face
                       will be created
  name:     (optional) name for the geoemtry. Defaults to filename

Output Ports:
  geometry: DejaVu2.IndexedPolgons geometry
."""
    
    def __init__(self, name='readIndexedPolygon', **kw):

        #warnings.warn("read Indexed Polygon is deprecated, use 'geoms from file' instead")

        kw['name'] = name
        apply( NetworkNode.__init__, (self,), kw )
        #self.readOnly = 1

        wdescr = self.widgetDescr
        wdescr['filename'] = {
            'class':'NEEntryWithFileBrowser', 'master':'node', 'width':10,
            'initialValue':'',
            'labelCfg':{'text':'Filename: '},
            'filetypes': [('vert','*.vert'), ('all', '*')]}
        wdescr['name'] = {
            'class':'NEEntry', 'master':'node', 'width':10,
            'labelCfg':{'text':'name:'},
            }

        ip = self.inputPortsDescr
        ip.append(datatype='string', name='filename')
        ip.append(datatype='string', required=False, name='name')

        op = self.outputPortsDescr
	op.append(datatype='geom', name='geom')

        code = """def doit(self, filename, name):
    if filename is None or filename == '':
        return
    from DejaVu2.IndexedPolygons import IndexedPolygonFromFile
    from os.path import splitext
    if name == '':
        name=None
    filename = splitext(filename)[0]
    self.outputData(geom=IndexedPolygonFromFile(filename, name))
"""

        self.setFunction(code)


from DejaVu2.utils import RemoveDuplicatedVertices

class RemoveDuplicatedVerticesNE(NetworkNode):
    """Remove duplicated vertices and re-index the polygonal faces such that
they share vertices

input:
    geom -- an IndexedGeom object in which vertices are suplicated
    newGeom -- when true a new geoemtry is created, else the incomming
               geometry is modified

output:
    geom -- a new IndexedPolgons without duplicated vertices
"""
    def __init__(self, name='removeDupVert', **kw):
        kw['name'] = name
        apply(NetworkNode.__init__, (self,), kw)

        ip = self.inputPortsDescr = []
        ip.append({'datatype':'geom', 'name':'geom'})
        ip.append({'datatype':'boolean', 'name':'newGeom'})

        self.widgetDescr['newGeom'] = {
                'class':'NECheckButton',
                'initialValue':False,
                'master':'node', 
                'labelGridCfg':{'sticky':'we'},
                'labelCfg':{'text':'new geometry'},
                }

        op = self.outputPortsDescr = []
        op.append({'datatype':'geom', 'name':'geom'})

        code = """def doit(self, geom, newGeom):

    from DejaVu2.utils import RemoveDuplicatedVertices
    from DejaVu2.IndexedPolygons import IndexedPolygons, Geom
    results = []

    if isinstance(geom, Geom):
        geom = [geom]

    for g in geom:
        vertices = g.vertexSet.vertices.array
        faces = g.faceSet.faces.array
        vertList, faceList  = RemoveDuplicatedVertices(vertices, faces)
        if newGeom:
            results.append( IndexedPolygons('%sMesh'%g.name,
                                   vertices=vertList, faces=faceList))
        else:
            g.Set(vertices=vertList, faces=faceList)
            results.append(g)

    self.outputData(geom=results)
"""
        self.configure(function=code)


class removeDuplicatedVerticesC(NetworkNode):
    """Remove duplicated vertices and re-index the polygonal faces such that
they share vertices (uses the C++ function from opengltk.extent.utillib)

input:
    geom -- an IndexedGeom object in which vertices are suplicated
    newGeom -- when true a new geoemtry is created, else the incomming
               geometry is modified
output:
    geom -- a new IndexedPolygons without duplicated vertices
"""
    def __init__(self, name='removeDupVertC', **kw):
        kw['name'] = name
        apply(NetworkNode.__init__, (self,), kw)

        ip = self.inputPortsDescr = []
        ip.append({'datatype':'geom', 'name':'geom'})
        ip.append({'datatype':'boolean', 'name':'newGeom'})

        self.widgetDescr['newGeom'] = {
                'class':'NECheckButton',
                'initialValue':False,
                'master':'node', 
                'labelGridCfg':{'sticky':'we'},
                'labelCfg':{'text':'new geometry'},
                }

        op = self.outputPortsDescr = []
        op.append({'datatype':'geom', 'name':'geom'})

        code = """def doit(self, geom, newGeom):
    results = []
    from DejaVu2.IndexedPolygons import IndexedPolygons, Geom
    from geomutils.geomalgorithms import removeDuplicatedVertices

    if isinstance(geom, Geom):
        geom = [geom]

    for g in geom:
        vertices = g.vertexSet.vertices.array
        faces = g.faceSet.faces.array
        norms = g.getVNormals()
        if len(vertices) and len(faces):
            vertList, faceList, normList = removeDuplicatedVertices(
                                         vertices, faces, norms  )
            if newGeom:
                results.append( IndexedPolygons('%sMesh'%g.name,
                       vertices=vertList, faces=faceList))#, vnormals=normList))
            else:
                g.Set(vertices=vertList, faces=faceList)#, vnormals=normList)
                results.append(g)
    self.outputData(geom=results)
"""
        self.configure(function=code)

class QslimExt(NetworkNode):
    """Invoked the Qslim algorithm as an external process to decimate a surface.
Input Ports
  geometry: Any DejaVu2.IndexedPolgons geometry. (required)
  persent:  percentage of faces to be retained in decimated model
  
Output Ports:
  geom:     decimated geometry
."""
    
    def __init__(self, name='write SMF', **kw):
        kw['name'] = name
        apply( NetworkNode.__init__, (self,), kw )

        self.widgetDescr['percent'] = {
            'class':'NEDial', 'size':50,
            'oneTurn':100, 'min':0., 'max':100., 'type':'float',
            'showLabel':1, 'continuous':0,
            'initialValue':50.,
            'labelCfg':{'text':'percent'},
            }
        ip = self.inputPortsDescr
        ip.append(name='lines', datatype='list')
        ip.append(name='percent', datatype='float')
        ip.append(name='nbf', datatype='int', required=False)

        op = self.outputPortsDescr
        op.append({'name':'lines', 'datatype':'list'})

        code = """def doit(self, lines, percent, nbf):

    if nbf is None:
        nbf = 0
        for l in lines:
            if l[0]=='f':
                nbf+=1
                
    targetFaces = int(nbf*percent/100.)

    from mglutil import process
    # start propslim process
    p = process.Popen('propslim %d'%targetFaces, 
                       stdin=process.PIPE, stdout=process.PIPE)
    inp, out, err = (p.stdin, p.stdout, p.stderr)

    # send input
    map( lambda x, f=inp: f.write('%s'%x), lines)
    inp.close()

    # FIXME we should select on out and err and handle errors
    # FIXME we should run popslim from the Binaries package
    # FIXME we should have a way to detect if the Binary package contains
    #       popslim and if not not load the node
    # read results
    data = out.readlines()
    print 'propslim reduced from %d  to %d faces'%(nbf,targetFaces)

    self.outputData(lines=data)
"""

        self.setFunction(code)



class GeomtoSMF(NetworkNode):
    """create the source of an SMF file from an DejaVu2.IndexedPolygons
geometry

Input Ports:
  geometry: Any DejaVu2.IndexedPolgons geometry. (required)
  filename: name of the SMF file to be created (required, string)

Output Ports:
  lines:   list of text lines in SMF format
."""
    
    def __init__(self, name='GeomToSMF', **kw):
        kw['name'] = name
        apply( NetworkNode.__init__, (self,), kw )

        ip = self.inputPortsDescr
        ip.append(name='geometry', datatype='geom')

        self.outputPortsDescr.append(name='lines', datatype='list')
        
        code = """def doit(self, geometry):

    from DejaVu2.IndexedPolygons import IndexedPolygons
    assert isinstance(geometry, IndexedPolygons)
    from DejaVu2.DataOutput import IndexedPolgonsAsSMFString
    lines = IndexedPolgonsAsSMFString(geometry)
    self.outputData(lines = lines)
"""

        self.setFunction(code)


class SMFtoGeom(NetworkNode):
    """converts a list of strings describing an SMF file into a
DejaVu2.IndexedPolygons

Input Ports:
  lines: list of lines
  geomname: name for the geoemtry
  parent: parent geometry

Output Ports:
  geometry: DejaVu2.IndexedPolgons geometry
."""
    
    def __init__(self, name='SMFtoGeom', **kw):
        kw['name'] = name
        apply( NetworkNode.__init__, (self,), kw )

        from DejaVu2.IndexedPolygons import IndexedPolygons
        self.geom = IndexedPolygons()

        self.widgetDescr['geomname'] = {
            'class':'NEEntry', 'master':'node', 'width':10,
            'labelCfg':{'text':'name:'},
            'initialValue':'smfModel',
            }

        ip = self.inputPortsDescr
        ip.append(name='lines', datatype='list')
        ip.append(name='geomname', datatype='string', required=False)
        ip.append(name='parent', datatype='geom', required=False)

        op = self.outputPortsDescr
	op.append(name='geom', datatype='geom')

        code = """def doit(self, lines, geomname, parent):
    if len(lines)==0:
        return

    if not geomname: # i.e. None or ''
        geomname = 'smfModel'
    
    from DejaVu2.DataOutput import ParseSMFString
    v,f,n,c,r = ParseSMFString(lines)

    self.geom.Set(name=geomname, vertices=v, faces=f, redo=0)

    if len(n):
        if len(n)==len(v):
            self.geom.Set(vnormals=n, redo=0)
        elif len(n)==len(f):
            geom.Set(fnormals=n, redo=0)
        else:
            errstr = 'number of normals %d does not not match number of vertices (%d) or number of faces(%d)'%(len(v), len(f))
            warnings.warn(errstr)
        
    if len(c):
        self.geom.Set(materials=c, inheritMaterial=False, redo=0)

    if parent is not None:
        self.geom.viewer.ReparentGeom(self.geom. parent)

    self.geom.Set(redo=1)
    
    self.outputData(geom=self.geom)
"""

        self.setFunction(code)


from DejaVu2.IndexedPolygons import IndexedPolygons
import numpy.oldnumeric as Numeric
try:
    from QSlimLib import qlimlib
except:
    pass


class QSlimDecimate:

    def __init__(self):
        
        self.model = None
        self.newverts = None      # Vertex array of new geometry 
        self.newfaces = None      # Face array of new geometry 
        self.newcolors = None     # Color array of new geometry
        self.newnorms = None
        self.numColors = 0        # Number of valid colors and
        self.colorBind = None
        self.maxFaces = 0  # number of faces in the geometry used to build the
                           # QSlim model   


    def build_model(self, geometry, t1=None, t2 = None):
        """Build new QSlim model."""

        verts = geometry.getVertices()
        faces = geometry.getFaces()
        colors = None
        self.colorBind = None
        binding = geometry.materials[1028].binding[1]
        bind_to_face = False
        if binding==11:
            #per vertex color
            bind_to_face = False
            self.colorBind = 'vertex'
        elif binding==12:
            #per face color
            bind_to_face = True
            self.colorBind = 'face'
        if self.colorBind is not None:   
            colors = geometry.materials[1028].prop[1][:,:3].astype('f')
            self.newcolors = Numeric.zeros(colors.shape, 'f')
        else:
            self.newcolors = None
        normals = geometry.getVNormals()
        #print "in build_model: faces %d, normals %d " % (len(faces), len(normals))
        #texcoords
        texcoords = None
        self.newtexcoords = None
        #t1 = Numeric.arange(len(verts)).astype('f')
        if t1 is not None and t2 is not None:
            assert len(t1) == len(t2)
            texcoords = Numeric.zeros((len(t1),2), "f")
            texcoords[:,0] = t1[:]
            texcoords[:,1] = t2[:]
        elif t1 is not None:
            texcoords = Numeric.zeros((len(t1),2), "f")
            texcoords[:,0] = t1[:]
        elif t2 is not None:
            texcoords = Numeric.zeros((len(t2),2), "f")
            texcoords[:,0] = t2[:]
        else:
            if hasattr(geometry.vertexSet, 'texCoords'):
                texcoords = geometry.vertexSet.texCoords.array
        if texcoords is not None:
            #print "len(texcoords): ", len(texcoords)
            #print "len(verts): ", len(verts)
            assert len(texcoords) == len(verts)
            self.newtexcoords = Numeric.zeros(texcoords.shape, 'f') 
        self.model = qslimlib.QSlimModel(verts, faces, bindtoface=bind_to_face,
                                         colors=colors, norms=normals, texcoords=texcoords)
        nverts = self.model.get_vert_count()
        self.maxFaces = nfaces = self.model.get_face_count()
        
        self.newverts = Numeric.zeros((nverts, 3)).astype('f')
        self.newfaces = Numeric.zeros((nfaces, 3)).astype('i')
        self.newnorms = Numeric.zeros((nverts, 3)).astype('f')
        

    def decimate(self, targetfaces): 

        print 'DECIMATING TO', targetfaces
        self.model.slim_to_target(targetfaces)
        # the number of faces we get back can be slightly different 
        # from the taget number of faces: 
        numFaces = self.model.num_valid_faces()
        
        self.model.outmodel(self.newverts, self.newfaces,
                         outcolors=self.newcolors, outnorms=self.newnorms,outtexcoords=self.newtexcoords )
        numVertices = max(self.newfaces.ravel())+1
        # build lookup table of vertices used in faces
        d = {}
        for t in self.newfaces[:numFaces]:
            d[t[0]] = 1
            d[t[1]] = 1
            d[t[2]] = 1
        vl = {}
        decimVerts = Numeric.zeros((numVertices, 3)).astype('f')
        decimFaces = Numeric.zeros((numFaces, 3)).astype('i')
        decimNormals = Numeric.zeros((numVertices, 3)).astype('f')
        decimColors = None
        decimTex = None

        if self.newcolors is not None:
            decimColors = {'vertex': Numeric.zeros((numVertices, 3)).astype('f'),
                            'face':self.newcolors[:numFaces]}.get(self.colorBind)
        if self.newtexcoords is not None:
            decimTex = Numeric.zeros((numVertices, 2)).astype('f')


        nvert = 0
        for j, t in enumerate(self.newfaces[:numFaces]):
            for i in (0,1,2):
                n = t[i]
                if not vl.has_key(n):
                    vl[n] = nvert
                    decimVerts[nvert] = self.newverts[n,:]
                    decimNormals[nvert] = self.newnorms[n,:]

                    if self.newcolors is not None:
                        if self.colorBind == 'vertex':
                            decimColors[nvert] = self.newcolors[n, :]
                    if self.newtexcoords is not None:
                        decimTex[nvert] = self.newtexcoords[n,:]
                    nvert += 1
            decimFaces[j] = (vl[t[0]], vl[t[1]], vl[t[2]]) 

        return (decimVerts[:nvert], decimFaces, decimNormals[:nvert], decimColors, decimTex)


class QSlim(NetworkNode):
    """Uses QSlim library to decimate a surface.
    Creates a QSlim model from the input geometry, decimates it
    and outputs a geometry with new number of faces, vertices, colors.
Input Ports
  geometry: any DejaVu2.IndexedPolgons geometry (required).
  percent:  percentage of faces to be retained in decimated model (required).
  targetfaces: corresponding to persent number of faces.
  newgeom: (checkbutton) if checked - the node outputs a new IndexedPolygons
           geometry, if unchecked - modified input geometry is output.
  rebuild: (boolean), when true, rebuilds the QSlim model (takes
           current parameters of selected geometry).
  
Output Ports:
  geom:     decimated geometry
."""

    def __init__(self, name='QSlim', **kw):
        
        kw['name'] = name
        apply( NetworkNode.__init__, (self,), kw )
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
        op.append(datatype='geom', name='geometry')
        op.append(datatype='None', name='u')
        op.append(datatype='None', name='v')
        
        self.geometry = None      # Current geometry
        self.model = QSlimDecimate()

        code = """def doit(self, geometry , percent, targetfaces, rebuild, u, v):
    #print 'QSlim Node Running QSlim', percent, targetfaces, rebuild, self.inputPortByName['rebuild'].hasNewValidData()
    # QSlim node can only handle 1 geom so far, so print out a message if there
    # are more comming in
    if type(geometry)==types.ListType:
        if len(geometry)>1:
            warnings.warn('Only first geometry is being processed by QSlim node')
        geometry = geometry[0]

    # check that the incomming geom can be represented as IndexedPolygons
    if not geometry.asIndexedPolygons(run=0):
       warnings.warn(geometry.name , 'QSlim Node: %s can not be represented as IndexedPolygons', geometry.name)
       return
    
    # if there is new data on the first input port (geom) rebuild the model
    #if self.inputPortByName['geometry'].hasNewValidData() or \
       #(self.inputPortByName['rebuild'].hasNewValidData() and rebuild):
    if self.inputPortByName['geometry'].hasNewValidData() or \
        self.inputPortByName['rebuild'].hasNewValidData():
       if len(geometry.vertexSet.vertices):
           print 'QSlim Node BUILDING MODEL'
           self.model.build_model(geometry, u, v)
       else:
           return
    maxFaces = self.model.maxFaces       
    # compute the number of triangles to be kept
    if self.inputPortByName['percent'].hasNewValidData(): #percent has changed
        targetfaces = int(maxFaces*percent/100.)
    elif self.inputPortByName['targetfaces'].hasNewValidData() and targetfaces: 
        #targetfaces has changed
        targetfaces = int(targetfaces)
    else: # use the percentage
        targetfaces = int(maxFaces*percent/100.)
        print 'ELSE', maxFaces*percent, targetfaces

    # make sure targetfaces is OK
    if targetfaces > maxFaces:
        targetfaces = maxFaces
        if self.inputPortByName['percent'].widget:
            self.inputPortByName['percent'].widget.set(100, 0)
        if self.inputPortByName['targetfaces'].widget:
            self.inputPortByName['targetfaces'].widget.set(maxFaces, 0)
        
    # decimate
    if self.model:
        decimVerts, decimFaces, decimNormals, decimColors, decimTex = self.model.decimate(targetfaces)
        redo = self.model.newcolors is None
        geometry.Set(vertices=decimVerts, faces=decimFaces, vnormals=decimNormals, redo=redo)
        #print 'len decimated verts = %d, faces = %d, normals = %d' % (len(decimVerts), len(decimFaces), len(decimNormals))
        geometry.vertexSet.vertices.array = decimVerts
        geometry.faceSet.faces.array = decimFaces
        #geometry.Set(faces=decimFaces, vertices=decimVerts, fnormals = decimNormals, redo=redo)
        numFaces = len(decimFaces)
        nvert = len(decimVerts)
        clen = len(decimFaces)
        if self.model.newcolors is not None:
            if self.model.colorBind == 'vertex':
                clen = nvert
            #print 'len colors:', clen
            if clen > 0:
                geometry.Set(materials=decimColors[:clen],
                         inheritMaterial=False, redo=True)
            else:
                geometry.Set(inheritMaterial=False, redo=True)
        else:
            geometry.Set(faces=decimFaces, vertices=decimVerts, fnormals = decimNormals)
        # set the target face widget to the actual number of faces after
        # decimation
        if self.inputPortByName['targetfaces'].widget:
            self.inputPortByName['targetfaces'].widget.set(numFaces, run=0)

        realpercent = 100.*(float(numFaces)/maxFaces)
        if percent != realpercent and self.inputPortByName['percent'].widget:
            self.inputPortByName['percent'].widget.set(realpercent, run=0)
        if self.model.newtexcoords is not None:
            if hasattr(geometry.vertexSet, 'texCoords'):
                geometry.Set(textureCoords=decimTex[:nvert])
            self.outputData(geometry=geometry, u=decimTex[:nvert,0], v = decimTex[:nvert,1] )
        else:
            self.outputData(geometry=geometry)
"""
        self.setFunction(code)


try:
    bhtreelibFound = True
    from bhtree import bhtreelib

    class decimate3DPoints(NetworkNode):
        """loops over a list of 3D points and removes any point to close to
an already seen one

Input Ports:
    points: (required) 3D points to be decimated
    cutoff: distance below which a point in cut

Output Ports:
    points: decimated set of points
."""

        def __init__(self, name='getSurfaceVFN', **kw):
            kw['name'] = name
            apply( NetworkNode.__init__, (self,), kw )
            #self.readOnly = 1

            ip = self.inputPortsDescr
            ip.append(datatype='coordinates3D', name='point')
            ip.append(datatype='float', name='cutoff')

            op = self.outputPortsDescr
            op.append(datatype='coordinates3D', name='points')

            code = """def doit(self, points, cutoff):
    import numpy.oldnumeric as Numeric
    ids = Numeric.arange(len(points)).astype('i')
    bht = bhtreelib.TBHTree( points, ids, 10, 10, 9999.0 )
    result = Numeric.zeros( (len(points),) ).astype('i')
    dist = Numeric.zeros( (len(points),) ).astype('f')
    pts = []
    removed = {}
    for p in points:
        if removed.has_key('%9.6f,%9.6f,%9.6f'%tuple(p)):
            continue
        pts.append(p)
        nb = bht.ClosePointsDist2( tuple(p), cutoff, result, dist )
        for close in result:
            removed['%9.6f,%9.6f,%9.6f'%tuple(points[close])] = True

    self.outputData(points=pts)
"""

            self.setFunction(code)


    class DistFromSphereToGeom(NetworkNode):
        """assign to every sphere the distance to the closest surface vertex.

Input ports:
    centers: a list of sphere centers
    radii: list of sphere radii 
    geom: a polygonal geometry

OutputPort:
    dist: list of distances
"""
        def __init__(self, name='DistFromSphereToGeom', **kw):
            kw['name'] = name
            apply( NetworkNode.__init__, (self,), kw )

            ip = self.inputPortsDescr
            ip.append(datatype='coordinates3D', name='centers')
            ip.append(datatype='None', name='radii')
            ip.append(datatype='geom', name='geom')

            op = self.outputPortsDescr
            op.append(datatype='list', name='distances')


            code = """def doit(self, centers, radii, geom):

    #build a BHtree for geometry's vertices
    vf = geom.getVertices()
    import numpy.oldnumeric as Numeric
    points = Numeric.array(vf, 'f')
    ids = Numeric.arange(len(points)).astype('i')
    result = Numeric.zeros( (len(points),) ).astype('i')
    dist = Numeric.zeros( (len(points),) ).astype('f')
    bht = bhtreelib.TBHTree( points, ids, 10, 10, 9999.0 )

    distances = []
    try:
        len(radii)
    except:
        radii = (radii,)*len(points)

    # compute distance from sphere to closest geometry vertex
    for c,r in zip(centers,radii):
        found = False
        cutoff=5.0
        while not found:
            nb = bht.ClosePointsDist2( tuple(c), cutoff, result, dist)
            if nb:
                distances.append( sqrt(min(dist[:nb]))-r )
                found = True
            else:
                cutoff += 5.0
                
    self.outputData(distances=distances)
"""

            self.setFunction(code)

except:
      bhtreelibFound = False
  

class getSurfaceVFN(NetworkNode):
    """extract vertices, faces and normals from a geometry.

Input Ports:
  geometry: (required) Accepts any object that is a DejaVu2.Geom instance
            and will add it to the viewer only once
Output Ports:
  geom:  either geom that is input or result of asIndexedPolygons
  vertices:
  faces:
  normals:
."""
    
    def __init__(self, name='getSurfaceVFN', **kw):
        kw['name'] = name
        apply( NetworkNode.__init__, (self,), kw )
        #self.readOnly = 1

        ip = self.inputPortsDescr
        ip.append(datatype='geom', name='geometry')

        op = self.outputPortsDescr
	op.append(datatype='geom', name='geom')
	op.append(datatype='coordinates3D', name='vertices')
	op.append(datatype='faceIndices', name='faces')
        op.append(datatype='normals3D', name='normals')
        op.append(datatype='normals3D', name='facenormals')

        code = """def doit(self, geometry):

    if geometry:
        geom = geometry.asIndexedPolygons()
        v = geom.getVertices()
        vn = geom.getVNormals()
        fn = geom.getFNormals()
        f = geom.getFaces()

        self.outputData(geom=geom, vertices=v, faces=f,
                        normals=vn, facenormals=fn)
"""

        self.setFunction(code)

try:
    from opengltk.extent import _gllib
    from PIL import Image
    from PIL import ImageChops

    class ImageViewerNode(NetworkNode):
        """create an instance of an OpenGL-based image viewer
    This object can compute the first and second derivative and the historgram
    of the image
    
    INPUT:
        image: image to be displayed
        operator: can be 'FirstDerivative', 'SecondDerivative' or
                   'Histogram'
    OUTPUT:
        imageViewer: the image viewer instance
        result: the processed image oir the histogram
    """

        def beforeAddingToNetwork(self, net):
            # import imagelib
            importImageLib(net)

            
        def _init__(self, name='ImageViewer', **kw):
            kw['name'] = name
            apply( NetworkNode.__init__, (self,), kw )

            from DejaVu2.imageViewer import ImageViewer
            self.ivi = ImageViewer(name='image renderer')

            self.widgetDescr['operator'] = {
                'class':'NEComboBox', 'master':'node',
                'choices':['FirstDerivative', 'SecondDerivative', 'Histogram'],
                'entryfield_entry_width':8,
                'labelCfg':{'text':'operator:'},
                }

            ip = self.inputPortsDescr
            ip.append(datatype='image', required=False, name='image')
            ip.append(datatype='string', required=False, name='operator')

            op = self.outputPortsDescr
            op.append(datatype='None', name='imageViewer')
            op.append(datatype='image', name='result')

            code = """def doit(self, image, operator):
        if image:
            self.ivi.setImage(image)
            
        if operator == 'FirstDerivative':
            outdata = self.ivi.firstDerivative()
            out = Image.fromstring('RGB', (self.ivi.width, self.ivi.height),
                                   outdata)
            self.outputPorts[1].setDataType('image', tagModified=False)

        elif operator == 'SecondDerivative':
            outata = self.ivi.secondDerivative()
            out = Image.fromstring('RGB', (self.ivi.width, self.ivi.height),
                                    outdata)
            self.outputPorts[1].setDataType('image', tagModified=False)

        elif operator == 'Histogram':
            out = self.ivi.getHistogram()
            self.outputPorts[1].setDataType('list', tagModified=False)

        else:
            out = None

        self.outputData(imageViewer=self.ivi, result=out )
"""

            self.setFunction(code)


    class NPR(NetworkNode):
        """make NPR rendering of the currentCamera of a viewer"""

        def beforeAddingToNetwork(self, net):
            # import imagelib
            importImageLib(net)
            from DejaVu2.imageViewer import ImageViewer
            self.ivi = ImageViewer(name='image renderer')


        def afterRemovingFromNetwork(self):
            NetworkNode.afterRemovingFromNetwork(self)
            if self.ivi:
                self.ivi.master.destroy()
                del self.ivi


        def __init__(self, name='NPR', **kw):
            kw['name'] = name
            apply( NetworkNode.__init__, (self,), kw )
            self.ivi = None  # will be instanciated when node added to network

            self.widgetDescr['d1cut'] = {
                'class':'NEThumbWheel', 'initialValue':0,
                'labelCfg':{'text':'d1cut:'}, 'width':80, 'height':20, 
                'type':'int', 'oneTurn':10, 'lockBMin':1, 'min':0,
                'wheelPad':2 }
            self.widgetDescr['d1scale'] = {
                'class':'NEThumbWheel', 'initialValue':3,
                'labelCfg':{'text':'d1scale:'}, 'width':80, 'height':20, 
                'type':'float', 'oneTurn':10., 'lockBMin':1, 'min':0,
                'wheelPad':2 }
            self.widgetDescr['d1off'] = {
                'class':'NEThumbWheel', 'initialValue':60,
                'labelCfg':{'text':'d1off:'}, 'width':80, 'height':20, 
                'type':'int', 'oneTurn':10, 'lockBMin':1, 'min':0,
                'wheelPad':2 }
            self.widgetDescr['d2cut'] = {
                'class':'NEThumbWheel', 'initialValue':1,
                'labelCfg':{'text':'d2cut:'}, 'width':80, 'height':20, 
                'type':'int', 'oneTurn':10, 'lockBMin':1, 'min':0,
                'wheelPad':2 }
            self.widgetDescr['d2scale'] = {
                'class':'NEThumbWheel', 'initialValue':8,
                'labelCfg':{'text':'d2scale:'}, 'width':80, 'height':20, 
                'type':'float', 'oneTurn':10., 'lockBMin':1, 'min':0,
                'wheelPad':2 }
            self.widgetDescr['d2off'] = {
                'class':'NEThumbWheel', 'initialValue':10,
                'labelCfg':{'text':'d2off:'}, 'width':80, 'height':20, 
                'type':'int', 'oneTurn':10, 'lockBMin':1, 'min':0,
                'wheelPad':2 }

            self.widgetDescr['firstDerivative'] = {
                'class':'NECheckButton', 'initialValue':1,
                'labelGridCfg':{'sticky':'we'},
                'labelCfg':{'text':'first Derivative:'},
                }
            
            self.widgetDescr['secondDerivative'] = {
                'class':'NECheckButton', 'initialValue':1,
                'labelGridCfg':{'sticky':'we'},
                'labelCfg':{'text':'secondtDerivative:'},
                }
            
            self.widgetDescr['smooth'] = {
                'class':'NECheckButton', 'initialValue':1,
                'labelGridCfg':{'sticky':'we'},
                'labelCfg':{'text':'smooth:'},
                }
        
            ip = self.inputPortsDescr
            ip.append(datatype='viewer', name='viewer')
            ip.append(datatype='int', required=False, name='d1cut', defaultValue=0)
            ip.append(datatype='float', required=False, name='d1scale', defaultValue=3.)
            ip.append(datatype='int', required=False, name='d1off', defaultValue=60)
            ip.append(datatype='int', required=False, name='d2cut', defaultValue=1)
            ip.append(datatype='float', required=False, name='d2scale', defaultValue=8.)
            ip.append(datatype='int', required=False, name='d2off', defaultValue=10)
            ip.append(datatype='boolean', required=False,
                      name='firstDerivative', defaultValue=True)
            ip.append(datatype='boolean', required=False,
                      name='secondDerivative', defaultValue=True)
            ip.append(datatype='boolean', required=False, name='smooth', defaultValue=True)

            op = self.outputPortsDescr
            op.append(datatype='image', name='image')
            op.append(datatype='image', name='outline')

            code = """def doit(self, viewer, d1cut, d1scale, d1off, d2cut, d2scale,
d2off, firstDerivative, secondDerivative, smooth):
        camera = viewer.currentCamera
        viewer.stopAutoRedraw()
        camera.frame.master.lift()
        viewer.OneRedraw()
        image = camera.GrabFrontBuffer()
        zbuf = camera.GrabZBuffer()
        camera.viewer.startAutoRedraw()
        width = camera.width
        height = camera.height
        #self.ivi.master.lift()
        self.ivi.setImage(zbuf)
        if firstDerivative:
            deriv = self.ivi.firstDerivative()
            if deriv is None:
                return
            d1strong = Numeric.where(Numeric.greater(deriv,d1cut),
                                     d1scale*(deriv+d1off), 0)

        if secondDerivative:
            deriv2 = self.ivi.secondDerivative()
            d2strong = Numeric.where(Numeric.greater(deriv2,d2cut),
                                     d2scale*(deriv2+d2off), 0)

        if firstDerivative and secondDerivative:
            dstrong = Numeric.maximum(d1strong, d2strong)
        elif firstDerivative:
            dstrong = d1strong
        elif secondDerivative:
            dstrong = d2strong

        mini = min(dstrong)
        maxi = max(dstrong)
        if maxi!=mini:
            dstrong = ((dstrong-mini)/(maxi-mini))*255
           
        derivcompb = Numeric.fabs(dstrong-255).astype('B')
        
##         if smooth:
##             from time import time
##             t1 = time()
##             self.ivi.setImage(derivcompb,width=width, height=height, mode='RGB')
##             derivcompb = self.ivi.smooth()
##             print time()-t1
           
        outline = Image.fromstring('RGB', (width,height), derivcompb)
        if os.name != 'nt': #sys.platform!='win32':
            outline = outline.transpose(Image.FLIP_TOP_BOTTOM)

        if smooth:
            from time import time
            from PIL import ImageFilter, ImageChops
            t1 = time()
            outlinesmooth = outline.filter(ImageFilter.SMOOTH_MORE)
            #outline = ImageChops.multiply(outline, outlinesmooth)
            print time()-t1
            
        from PIL import ImageChops
        contouredImage = ImageChops.multiply(image, outline)
##         mask = dstrong
##         mask.shape = ( width, height, -1)
##         maskr = Numeric.array(mask[:,:,0]).astype('B')
##         maskim = Image.fromstring('L', (width,height), maskr)
##         maskim = maskim.transpose(Image.FLIP_TOP_BOTTOM)
##         image.paste(outline, None, maskim)

        self.outputData(image=contouredImage, outline=outline)
"""

            self.setFunction(code)
except:
    pass



class SelectMultipleGeometry(NetworkNode):
    """Choose one or more geometry using regular expression matching on the name
Input:
    viewer: Viewer in which to search for geoemtry objects
    name:   regular expression
    parent: (optional) if specified, only children geometry objects
            of 'parent' will be searched.  Defaults to the root object.
Output:
    geomList:  list of geometry objects
"""
    
    def __init__(self, name='Choose Geoms', **kw):
        kw['name'] = name
        apply( NetworkNode.__init__, (self,), kw )

        self.widgetDescr['name'] = {
            'class':'NEEntry', 'master':'node',
            'labelCfg':{'text':'name pattern:'},
            'initialValue':'.*',
            }

        ip = self.inputPortsDescr
        ip.append(datatype='viewer', name='viewer')
        ip.append(datatype='string', name='name')
        ip.append(datatype='geom', required=False, name='searchedParent')

        op = self.outputPortsDescr
        op.append(datatype='geom', name='geomList')

        code = """def doit(self, viewer, name, searchedParent):
    result = viewer.findGeomsByName(name, searchedParent)
    self.outputData(geomList=result)
"""

        self.setFunction(code)



class SelectGeometry(NetworkNode):
    """Choose a geometry in the DejaVu2 Viewer."""
    
    def __init__(self, name='Choose Geom', **kw):
        kw['name'] = name
        apply( NetworkNode.__init__, (self,), kw )

        self.widgetDescr['geomName'] = {
            'class':'NEComboBox', 'master':'node',
            'choices':[''],
            'autoList':True,
            'entryfield_entry_width':16,
            'labelCfg':{'text':'geom:'},
            }

        ip = self.inputPortsDescr
        ip.append(datatype='viewer', name='viewer')
        ip.append(datatype='string', required=False, name='geomName')

	op = self.outputPortsDescr
        op.append(datatype='geom', name='geometry')

        code = """def doit(self, viewer, geomName):
    if viewer:
        allObjs = viewer.rootObject.AllObjects()
        allNames = map(lambda x: x.fullName, allObjs)
        w = self.inputPortByName['geomName'].widget
        if w:
            w.configure(choices=allNames)
        if geomName not in allNames:
            geomName = None
            w.widget.setentry('')

    if geomName:
        obj = allObjs[allNames.index(geomName)]
        if obj:
            self.outputData(geometry=obj )
"""

        self.setFunction(code)



class Redraw(NetworkNode):
    """Create a DejaVu2 Redraw event. Trigger can be anything
Input:
    Viewer: one or more viewer
    trigger: 
Output:
"""
    
    def __init__(self, name='Redraw', **kw):
        kw['name'] = name
        apply( NetworkNode.__init__, (self,), kw )

        ip = self.inputPortsDescr
        ip.append(datatype='viewer', singleConnection=False, name='viewer')
        ip.append(name='trigger')

        code = """def doit(self, viewer, trigger):
    for v in viewer:
        v.Redraw()\n"""

        self.setFunction(code)


class StopAutoRedraw(NetworkNode):
    """Stop autoRedraw mode in the viewer
Input:
    Viewer: one or more viewer
    trigger: 
Output:
"""
    
    def __init__(self, name='Stop Auto Redraw', **kw):
        kw['name'] = name
        apply( NetworkNode.__init__, (self,), kw )

        ip = self.inputPortsDescr
        ip.append(datatype='viewer', singleConnection=False, name='viewer')
        ip.append(name='trigger')

        code = """def doit(self, viewer, trigger):
    for v in viewer:
        v.stopAutoRedraw()\n"""

        self.setFunction(code)


class StartAutoRedraw(NetworkNode):
    """Start autoRedraw mode in the viewer
Input:
    Viewer: one or more viewer
    trigger: 
Output:
"""
    
    def __init__(self, name='Start Auto Redraw', **kw):
        kw['name'] = name
        apply( NetworkNode.__init__, (self,), kw )

        ip = self.inputPortsDescr
        ip.append(datatype='viewer', singleConnection=False, name='viewer')
        ip.append(name='trigger')

        code = """def doit(self, viewer, trigger):
    for v in viewer:
        v.startAutoRedraw()\n"""

        self.setFunction(code)


class OneRedraw(NetworkNode):
    """Make one redraw in the viewer. This is used when auto redrwa is stopped
Input:
    Viewer: one or more viewer
    trigger: 
Output:
"""
    
    def __init__(self, name='One Redraw', **kw):
        kw['name'] = name
        apply( NetworkNode.__init__, (self,), kw )

        ip = self.inputPortsDescr
        ip.append(datatype='viewer', singleConnection=False, name='viewer')
        ip.append(name='trigger')

        code = """def doit(self, viewer, trigger):
    for v in viewer:
        v.OneRedraw()\n"""

        self.setFunction(code)



class SetInstances(NetworkNode):
    """Assign instance matrices to a geometry.
"""    
    def __init__(self, name='Set Instances', **kw):
        kw['name'] = name
        apply( NetworkNode.__init__, (self,), kw )

        ip = self.inputPortsDescr
        ip.append(datatype='geom', singleConnection=False, name='geometry')
        ip.append(datatype='instancemat(0)', name='matrices')

        code = """def doit(self, geometry, matrices):
    if geometry and matrices:
        for g in geometry:
            g.Set(instanceMatrices=matrices)
        if g.viewer:
            g.viewer.Redraw()\n"""

        self.setFunction(code)


    def beforeAddingToNetwork(self, net):
        # import symlib
        importSymLib(net)



class ReparentGeom(NetworkNode):
    """******************* DEPRECATED : USE geoms properties INSTEAD
FOR THE GEOMETRY NODES OR THE IMPUT PORT "parent" 
assigns a geometry to a new parent.
This node can handle a single geometry, list of geometry objects or a list
containing both geoemtry objects and lists of them.

Input Ports:
    geom:   geometry to be assigned to a anew parent
    parent: parent geometry

Output Ports:
"""
    def __init__(self, name='Reparent Geom', **kw):
        kw['name'] = name
        apply( NetworkNode.__init__, (self,), kw )
        

        ip = self.inputPortsDescr
        ip.append(datatype='geom(0)', singleConnection=False, name='geoms')
        ip.append(datatype='geom', required=False, name='parent')
        ip.append(datatype='boolean', required=False, name='retainPosition', defaultValue=False)

        self.widgetDescr['retainPosition'] = {
            'class':'NECheckButton', 
            'master':'node',
            'initialValue':0,
            'labelGridCfg':{'sticky':'we'},
            'labelCfg':{'text':'retain position:'},
            }

        code = """def doit(self, geoms, parent, retainPosition):

    def getRidOfInnerLists(alist):
        returnedList = []
        for element in alist:
            if isinstance(element, Geom) or isinstance(element, Insert2d):
                returnedList.append(element)
            elif isinstance(element, list):
                returnedList += getRidOfInnerLists(element)
        return returnedList

    try:
        if parent is not None:
            if len(parent)>1:
                raise TypeError('parent should be a geometry or a list of lenght 1')
            else:
                parent = parent[0]
    except AttributeError:
        pass

    geoms = getRidOfInnerLists(geoms)
    if len(geoms) > 0:
        viewer = geoms[0].viewer
        assert viewer is not None
        root = viewer.rootObject
        for geom in geoms:
            if parent is None:
                parent = root
            viewer.ReparentObject(geom, parent, objectRetainsCurrentPosition=retainPosition)
"""

        if code: self.setFunction(code)


class GeomsProperties(NetworkNode):
    """this node (should not be used with GeometryNode
because it already manages all that)
assigns a geometry to a new parent and set its properties.
This node can handle a single geometry, list of geometry objects or a list
containing both geoemtry objects and lists of them.

Input Ports:
    geom:   geometry to be assigned to a anew parent
    parent: parent geometry

Output Ports:
"""
    def __init__(self, name='geoms properties', **kw):
        kw['name'] = name
        apply( NetworkNode.__init__, (self,), kw )
        

        ip = self.inputPortsDescr
        ip.append(datatype='geomOrInsert2d(0)', singleConnection=False, name='geoms')
        ip.append(datatype='int', required=True, name='noParentReparentToRoot')
        ip.append(datatype='colorfloat3or4(0)', required=False, name='colors')
        ip.append(datatype='instancemat(0)', required=False,
                  name='instanceMatrices')
        ip.append(datatype='dict', required=False, name='geomOptions')
        ip.append(datatype='geomOrInsert2d', required=False, name='parent', 
                  singleConnection=True)
        ip.append(datatype='boolean', required=False, name='retainPosition', defaultValue=False)

        self.widgetDescr['retainPosition'] = {
            'class':'NECheckButton', 
            'master':'node',
            'initialValue':0,
            'labelGridCfg':{'sticky':'we'},
            'labelCfg':{'text':'retain position:'},
            }

        self.widgetDescr['noParentReparentToRoot'] = {
            'class':'NECheckButton',
            'master':'node',
            'initialValue':0,
            'labelGridCfg':{'sticky':'we'},
            'labelCfg':{'text':'no parent reparent to root:'},
            }

        code = """def doit(self, geoms, noParentReparentToRoot, colors,
instanceMatrices, geomOptions, parent, retainPosition):

    def getRidOfInnerLists(alist):
        returnedList = []
        for element in alist:
            if isinstance(element, Geom) or isinstance(element, Insert2d):
                returnedList.append(element)
            elif isinstance(element, list):
                returnedList += getRidOfInnerLists(element)
        return returnedList


    try:
        if parent is not None:
            if len(parent)>1:
                raise TypeError('parent should be a geometry or a list of lenght 1')
            else:
                parent = parent[0]
    except AttributeError:
        pass

    if geomOptions is not None:
        insert2dOpts = geomOptions.copy()
        geomOpts = geomOptions.copy()
    else:
        insert2dOpt = {}
        geomOpts = {}

    if colors is not None:
        geomOpts['inheritMaterial'] = False
        geomOpts['materials'] = colors

    if instanceMatrices:
        geomOpts['instanceMatrices'] = Numeric.array(instanceMatrices)
    else:
        geomOpts['instanceMatrices'] = [Numeric.identity(4).astype('f')]

    geoms = getRidOfInnerLists(geoms)
    if len(geoms) > 0:
        viewer = geoms[0].viewer
        assert viewer is not None
        root = viewer.rootObject
        for geom in geoms:
            #opts['redo'] = 0
            #opts['tagModified'] = False
            if isinstance(geom, Geom):
                apply( geom.Set, (), geomOpts )
            else:
                apply( geom.Set, (), insert2dOpts )

            if parent is None and noParentReparentToRoot == 1: 
                parent = root

            if parent is not None:
                viewer.ReparentObject(geom, parent, objectRetainsCurrentPosition=retainPosition)
"""

        if code: self.setFunction(code)



class QConvexHull(NetworkNode):
    """Compute a convex hull using the external module QConvex.
The convex hull of a set of points is the smallest convex set containing
the points. Based on QHull http://www.qhull.org/

Input:
    coords: 3-D coordinates of vertices
    colors: list of colors. Length of list can be 1, nb. of vertex, nb. of faces
    name: name of the geometry
    parent: parent geometry
    instanceMatrices: list of 4x4 transformation matrices used to draw
                      the object
    tmpDir: directory used for writing temporary files
    
Output:
    geom: indexed polygones geometry object

"""

    def __init__(self, name='Convex Hull', **kw):
        kw['name'] = name
        from mglutil.util.qhullUtils import QConvex
        self.QH = QConvex()

        apply( NetworkNode.__init__, (self,), kw )

        self.widgetDescr['tmpdir'] = {
            'class':'NEEntry', 'master':'node',
            'labelCfg':{'text':'tmp dir:'}}

        ip = self.inputPortsDescr
        ip.append(datatype='coordinates3D', name='coords')
        ip.append(datatype='colorsRGB', required=False, name='colors')
        ip.append(datatype='string', required=False, name='name')
        ip.append(datatype='geom', required=False, name='parent')
        ip.append(datatype='instancemat(0)', required=False,
                  name='instanceMatrices')
        ip.append(datatype='string', required=False, name='tmpdir')

        self.widgetDescr['name'] = {
            'class':'NEEntry', 'master':'node', 'width':8,
            'initialValue':'qhull',
            'labelCfg':{'text':'name: '}
            }

        op = self.outputPortsDescr
        op.append(datatype='geom', name='geom')
       
        code = """def doit(self, coords, colors, name, parent, instanceMatrices, tmpdir):

    if tmpdir is None:
        tmpdir = './'
    self.QH.computeConvex(coords, tmpdir)

    geom = self.QH.getGeom()

    if name is not None and name!=geom.name[:len(name)]:
        geom.name = name
        if geom.viewer:
            vi = geom.viewer
            vi.stopAutoRedraw()
            vi.RemoveObject(geom)
            vi.AddObject(geom, parent=geom.parent)
            vi.startAutoRedraw()

    if parent is not None:
        if geom.viewer:
            vi = geom.viewer
            vi.stopAutoRedraw()
            vi.RemoveObject(geom)
            vi.AddObject(geom, parent=parent)
            vi.startAutoRedraw()
        else:
            geom.parent = parent

    if colors:
        geom.Set(materials=colors)

    if instanceMatrices:
        mat = Numeric.array(instanceMatrices)
        geom.Set(instanceMatrices=mat)
    
    self.outputData(geom=geom)
"""
     
        self.setFunction(code)


class EllipsoidFit(NetworkNode):
    """Fit an ellipsoid into a list of 3-D coordinates.
DejaVu2 geometries can be generated using the Ellipsoid node

Input:  list of 3-D coordinates or list of lists of coordinates. In the latter
        case, multiple ellipsoids are generated, one for each list of coords
Output: list of Ellipsoid Objects"""

    def __init__(self, name='EllipsoidFit', **kw):
        kw['name'] = name
        #from geomutils import geomutilslib
        from geomutils import efitlib

        apply( NetworkNode.__init__, (self,), kw )

        ip = self.inputPortsDescr
        ip.append(datatype='list', name='coords')
        ip.append(datatype='float', required=False, name='cov_scale')
        ip.append(datatype='float', required=False, name='ell_scale')

        op =self.outputPortsDescr 
        op.append(datatype='list', name='ellipsoids')
       
        code = """def doit(self, coords, cov_scale, ell_scale):
    if not coords:
        return

    # we handle list of lists of coordinates and list of coordinates
    test = coords[0][0]
    if not type(test) == types.ListType:
        coords = [coords]

    if cov_scale is None:
        cov_scale = 1.75

    if ell_scale is None:
        ell_scale = 1.0

    out = []
    for coordList in coords:
        # create an ellipsoid structure
        ellipse = efitlib.ellipsoid()
        # create an ellipsoidInfo structure
        ellipseInfo = efitlib.efit_info()
        # compute the ellipsoid
        status = efitlib.fitEllipse(coordList, ell_scale, cov_scale,
                                 ellipseInfo, ellipse)
        #ellipse = geomutilslib.ellipsoid()
        #ellipseInfo = geomutilslib.efit_info()
        #status = geomutilslib.fitEllipse(coordList, ell_scale, cov_scale,
        #                                 ellipseInfo, ellipse)
        if status==0: # OK
            out.append(ellipse)

    if not len(out):
        return

    self.outputData(ellipsoids=out)
"""
     
        self.setFunction(code)


class GyrationSphere(NetworkNode):
    """Calculate the center and radius of the gyration Sphere
Input:
    coords: 3-D coordinates of vertices
Output:
    center:
    radius:
"""

    def __init__(self, name='Gyration Sphere', **kw):
        kw['name'] = name
        apply( NetworkNode.__init__, (self,), kw )

        from DejaVu2.Spheres import Spheres
        self.geom = Spheres('gyration', quality=30)

        ip = self.inputPortsDescr
        ip.append(datatype='coordinates3D', name='coords')

        op =self.outputPortsDescr 
        op.append(datatype='list', name='center')
        op.append(datatype='float', name='radius')
       
        code = """def doit(self, coords):

    def gyrationSphere(points):
        from numpy import sum
        
        if not isinstance(points, Numeric.ArrayType):
            points = Numeric.array(points)

        # compute centroid
        center = sum(points, 0)/len(points)
        
        # compute dist distances (vectorized)
        delta = points-center
        rg = sqrt( sum( sum( delta*delta, 1))/float(len(points)) )
        
        return center, rg

    res, rg = gyrationSphere(coords)
    center = [ res[0], res[1], res[2] ]

    self.outputData(center=center, radius=rg)
"""
     
        self.setFunction(code)


class GroupGeomsNE(NetworkNode):
    """Group DejaVu2 geoms."""

    def __init__(self, name='Group Geoms', **kw):
        kw['name'] = name
        apply( NetworkNode.__init__, (self,), kw )

        self.g = None

        code = """def doit(self, geomsIn, instanceMatrices, parentName):
    from DejaVu2.Geom import Geom
    if self.g is None:
        self.g = Geom(name=parentName)

    if geomsIn is None or len(geomsIn) == 0:
        for g in self.g.children:
           g.viewer = None
           g.parent = None
        self.g.children = []
        return 

    for c in geomsIn:
        if not c in self.g.children:
            c.parent = self.g
            self.g.children.append(c)
     
    if instanceMatrices:
        mat = Numeric.array(instanceMatrices)
    else:
        mat = None

    if self.g is not None and mat is not None:
        self.g.Set(instanceMatrices=mat)    

    self.outputData(geomsOut=self.g)\n"""
            
        if code: self.setFunction(code)

        self.widgetDescr['parentName'] = {
            'class':'NEEntry', 'master':'node',
            'width':12,
            'labelCfg':{'text':'name:'},
            }

        ip = self.inputPortsDescr
        ip.append(datatype='geom', required=False, singleConnection=False,
                  name='geomsIn')
        ip.append(datatype='instancemat(0)', required=False,
                  name='instanceMatrices')
        ip.append(datatype='string', required=False, name='parentName', defaultValue='group')

        op = self.outputPortsDescr
        op.append(datatype='geom', name='geomsOut')


class Decimate(NetworkNode):

    try:
        import multires
    except:
        pass
    
    def getTriangles(self):
        ntri = self.model.CountTriangles()
        fs = Numeric.zeros( (ntri, 3), 'i')
        self.model.GetCurrentTriangles(fs, ntri)
        return fs

    def __init__(self, name='Color By Ramp', **kw):
        kw['name'] = name
        apply( NetworkNode.__init__, (self,), kw )

        self.model = None
        self.resolution = 0
        
        self.widgetDescr['resolution'] = {
            'class':'NEDial', 'size':50,
            'oneTurn':10, 'min':-20, 'max':0, 'increment':0, 'type':'int',
            'showLabel':1,
            'initialValue':0,
            'labelCfg':{'text':'resolution'},
            }
        
        ip = self.inputPortsDescr
        ip.append(name='vertices', datatype='array(-1,3)')
        ip.append(name='indices', datatype='array(-1,3)')
        ip.append(datatype='float', name='resolution')

        op = self.outputPortsDescr
        op.append(datatype='array(-1,3', name='indices')

        code = """def doit(self, vertices, indices, resolution):
    if resolution == 0: # nothing to do, we keep highest resolution
       self.resolution = resolution
       self.outputData(indices=indices)
       return
    if self.model is None or len(vertices)!=self.nbVertices or \
            len(indices)!=self.nbIndices:
        self.nbVertices = len(vertices)
        self.nbIndices = len(indices)
        self.model = self.multires.ModelT(vertices, indices)

    if resolution==self.resolution:
      print 'WARNING ', resolution, self.resolution
      
    if resolution > self.resolution:
        for i in range(resolution-self.resolution):
            self.model.GoFiner()
    else:
        for i in range(self.resolution-resolution):
            self.model.GoCoarser()
    self.resolution = resolution
    tri = self.getTriangles()
    print len(tri), 'triangles'
    if len(tri):
        self.outputData(indices=tri)\n"""

        self.setFunction(code)


class Texture1DCoords(NetworkNode):
    """Compute 1D texture coordinates for a list of values

InputPorts:
    data:  list of values
    mini:  start of value range
    maxi:  end of value range

OutputPorts:
    texCoords: list of texture coordinates
"""

    def __init__(self, name='1D Texture Coordinates', **kw):
        kw['name'] = name
        apply( NetworkNode.__init__, (self,), kw)
            
        ip = self.inputPortsDescr
        ip.append(datatype='list', name='data')
        ip.append(datatype='float', required=False, name='mini')
        ip.append(datatype='float', required=False, name='maxi')

        op = self.outputPortsDescr
        op.append(datatype='list', name='texCoords')
	
        code = """def doit(self, data, mini, maxi):

    if data is None or len(data)==0:
        return

    if mini is None:
        mini = min(data)
    if maxi is None:
        maxi = max(data)

    delta = maxi-mini
    data = Numeric.array(data)
    texc = (data-mini)/delta
    self.outputData(texCoords=texc)
"""
        if code: self.setFunction(code)
        


class Texture1D(NetworkNode):
    """ This node maps a 1-dimensional texture to a geometry.

InputPorts:
    geom: DejaVu2 Geometry
    texture: 1D texture texture (comming from color map node)
    property: list of texture indices

OutputPorts:
    geom: DejaVu2 Geometry with texture

# FIXME ! textture has to be a colormap with mini and maxi
The texture can be of type image (RGB or RGBA, one pixel wide) or of type
Numeric.array (shape (x,1) ).
The property has to be of length coords of the geometry."""

    def __init__(self, name='1D Texture Mapping', **kw):
        kw['name'] = name
        apply( NetworkNode.__init__, (self,), kw)
	
	code = """def doit(self, geom, texture, prop):
    # if we get a new texture, we have to set texture and textureCoords
    
    prop = Numeric.array(prop, 'f')
    prop.shape = (-1,1)

    if self.inputPortByName['texture'].hasNewValidData():
        tex = Numeric.array(Numeric.array(texture)*255).astype('B')
        from DejaVu2.Texture import Texture
        t = Texture()
        t.Set(enable=1, image=tex)
        geom.Set(texture=t, textureCoords=prop, transparent='implicit')

    # if we only get new property values, we only have to update the
    # texture coords
    elif self.inputPortByName['prop'].hasNewValidData():
        geom.Set(textureCoords=prop)

    self.outputData(geom=geom)\n"""
            
        if code: self.setFunction(code)
        
	ip = self.inputPortsDescr
        ip.append(datatype='geom', name='geom')
	ip.append(datatype='texture', name='texture')
	ip.append(datatype='None', name='prop')

        op = self.outputPortsDescr
        op.append(datatype='geom', name='geom')



class Texture2D(NetworkNode):
    """ This node maps a 2-dimensional texture to a geometry.

InputPorts:
    geom: DejaVu2 Geometry
    texture: 2d texture (comming from color map node)

OutputPorts:
    geom: DejaVu2 Geometry with texture
"""

    def __init__(self, name='Texture2D', **kw):
        kw['name'] = name
        apply( NetworkNode.__init__, (self,), kw)

        ip = self.inputPortsDescr
        ip.append(datatype='geom', name='geom')
        ip.append(datatype='image', name='image')
        ip.append(datatype='coord2(0)', required=False, name='textureCoordinates')

        op = self.outputPortsDescr
        op.append(datatype='geom', name='geom')
        
        code = """def doit(self, geom, image, textureCoordinates):
    if self.inputPortByName['image'].hasNewValidData():
        from DejaVu2.Texture import Texture
        lTexture = Texture(enable=1, image=image) 
        geom.Set(texture=lTexture)

    if self.inputPortByName['textureCoordinates'].hasNewValidData():
        geom.Set(textureCoords=textureCoordinates)

    if self.inputPortByName['image'].hasNewValidData() or self.inputPortByName['textureCoordinates'].hasNewValidData():
        if hasattr(geom.vertexSet, 'texCoords'):
            geom.vertexSet.texCoords.array = \
                  lTexture.getTextureCoordinatesAfterResizeRatio(
                                      geom.vertexSet.texCoords.array)

    self.outputData(geom=geom)
"""    
        if code: 
            self.setFunction(code)



class ColorChooser(NetworkNode):
    """********************* DEPRECATED : USE Color Editor INSTEAD ***************************
Interactively choose a color from a color wheel. Double-click the node
to expand the color wheel. Right-click to open the parameter panel. Here you
can change the brightness of the color, turn the alpha value on/off and choose
whether a list of lists or just a list of values is output.

Input Ports
    color:\t\tbound to a color wheel
    brightness:\tbound to a thumbwheel
    alpha:\t\tbound to a checkbutton
    singleCol:\tbound to a checkbutton

Output Ports
    colors: the resulting color"""
    
    def __init__(self, name='Choose Color', **kw):
        kw['name'] = name
        apply( NetworkNode.__init__, (self,), kw )
        #self.readOnly = 1
        self.w_brightness = 255 # variable used to gain some speed in doit()
        self.w_singleCol = 0 # used to switch between datatypes of outp.Port
        
        self.widgetDescr['color'] = {
            'class':'NEColorWheel', 'master':'node',
            'width':75, 'height':75, 'circles':5, 'stripes':20,
            'immediate':1, 'wysiwyg':0, 'wheelPad':0,
            'widgetGridCfg':{'sticky':'we', 'labelSide':'top'},
            'labelCfg':{'text':''}, 
            }

        self.widgetDescr['brightness'] = {
            'class':'NEThumbWheel', 'master':'ParamPanel',
            'width':80, 'height':15, 'type':'int', 'wheelPad':1,
            'min':0, 'max':255,
            'initialValue':255,
            'labelGridCfg':{'sticky':'we'},
            'widgetGridCfg':{'sticky':'we', 'columnspan':2,'labelSide':'top'},
            'labelCfg':{'text':'brightness:'},
            }

        self.widgetDescr['alpha'] = {
            'class':'NECheckButton', 'master':'ParamPanel',
            'initialValue':1,
            'labelGridCfg':{'sticky':'we'},
            'labelCfg':{'text':'alpha val:'},
            }
        
        self.widgetDescr['singleCol'] = {
            'class':'NECheckButton', 'master':'ParamPanel',
            'initialValue':0,
            'labelGridCfg':{'sticky':'we'},
            'labelCfg':{'text':'single col:'},
            }

        ip = self.inputPortsDescr
        ip.append(datatype='None', name='color')
        ip.append(datatype='int', name='brightness')
        ip.append(datatype='int', name='alpha')
        ip.append(datatype='int', name='singleCol')
        
        op = self.outputPortsDescr
        op.append(datatype='colorsRGB', name='colors')

        code = """def doit(self, color, brightness, alpha, \
        singleCol):

    if brightness != self.w_brightness:
        self.w_brightness = brightness
        w = self.inputPorts[0].widget
        w.widget.BuildWheel(brightness/255.0)
        color = w.get()

    if alpha == 0:
        color = color[:-1]

    if singleCol != self.w_singleCol:
        #self.outputPorts[0]._modified = 1
        self.w_singleCol = singleCol
        port = self.outputPorts[0]
        if singleCol == 1:
            port.setDataType('colorRGB')
        else:
            port.setDataType('colorsRGB')

    if self.w_singleCol == 0:
        self.outputData(colors=[color])
    else:
        self.outputData(colors=color)
"""

        self.setFunction(code)


    def afterAddingToNetwork(self):
        NetworkNode.afterAddingToNetwork(self)
        # run this node so the value is output
        self.run()
       

class ColorEditor(NetworkNode):
    """Interactively choose a color from a sophisticated widget exposing
a colorwheel, a value scale, HSV entries, RGB entries and HEX(HexTriplet)
entries. The color is output as a list.

The widget is bound to the parameter panel of this node.

Input Ports
    color:\t\tbound to a color wheel

Output Ports
    color: the resulting color"""
    
    def __init__(self, name='Choose Editor', **kw):
        kw['name'] = name
        apply( NetworkNode.__init__, (self,), kw )
        #self.readOnly = 1
        self.w_brightness = 255 # variable used to gain some speed in doit()
        self.w_singleCol = 0 # used to switch between datatypes of outp.Port
        
        self.widgetDescr['color'] = {
            'class':'NEColorEditor', 'master':'ParamPanel',
            'mode':'RGB', 'immediate':True,
            'widgetGridCfg':{'sticky':'we', 'labelSide':'top',},
            'labelCfg':{'text':''},
            }

        ip = self.inputPortsDescr
        ip.append(datatype='None', name='color')
       
        op = self.outputPortsDescr
        op.append(datatype='colorsRGB', name='color')

        code = """def doit(self, color):
    w = self.inputPorts[0].widget
    color = w.get()
    self.outputData(color=color)
"""

        self.setFunction(code)


    def afterAddingToNetwork(self):
        NetworkNode.afterAddingToNetwork(self)
        # run this node so the value is output
        self.schedule_cb()



class ColorMapNE(NetworkNode):
    """Create a user editable color map
Turn a list of values into a list of colors using a colormap.
    
You can unbind the colormap in the parampanel 
to re-use an already existing colormap.
    
Input:
    colorMap: a ramp of colors to be used for the color map object. If
           omitted a ramp of 256 colors ranging from blue to red will be used
    values: an optionnal list of values
    min: the floating point value to be mapped to the first color. If omitted
         it defaults to the minimum value in the list.
    max: the flaoting point value to be mapped to the last color. If omitted
         it defaults to the maximum value in the list.
Output:
    colorMap
    legend: if connected to a viewer, a legend can be output
"""
    
    def __init__(self, name='Color Map', **kw):
        kw['name'] = name
        apply( NetworkNode.__init__, (self,), kw )
        #self.readOnly = 1
        
        ip = self.inputPortsDescr
        ip.append(datatype='ColorMapType', name='colorMap')
        ip.append(datatype='vector', required=False, name='values')
        ip.append(datatype='float', required=False, name='mini')
        ip.append(datatype='float', required=False, name='maxi')
        ip.append(datatype='string', name='filename')
        
        self.widgetDescr['colorMap'] = {
            'class':'NEColorMap', 
            'master':'ParamPanel', 
            'mini':None, 
            'maxi':None,
            'labelCfg':{'text':'colormap'},
            'widgetGridCfg':{'labelSide':'top'},
            }
        
        self.widgetDescr['filename'] = {
            'class':'NEEntryWithFileBrowser', 
            'master':'ParamPanel', 
            'width':10,
            'initialValue':'',
            'filetypes':[("ColorMap",'*_map.py'), ("any file",'*.*')],
            'labelCfg':{'text':'filename: '},
            'widgetGridCfg':{'labelSide':'top'}
            }
        
        op = self.outputPortsDescr
        op.append(datatype='colorfloat3or4(0)', name='mappedColors')
        op.append(datatype='ColorMapType', name='colorMap')
        op.append(datatype='geom', name='legend')

        code = """def doit(self, colorMap, values, mini, maxi,
                           filename):

    if self.inputPortByName['values'].hasNewValidData() or \
       self.inputPortByName['mini'].hasNewValidData() or \
       self.inputPortByName['maxi'].hasNewValidData():

        if values is not None and len(values)>0:
            if mini is None:
                mini = min(values)
            if maxi is None:
                maxi = max(values)

        colorMap.configure(mini=mini, maxi=maxi)

    p = self.inputPortByName['filename']
    if p.hasNewValidData() and filename:
        colorMap.read(filename=filename)

    for c in self.outputPortByName['legend'].connections:
        from DejaVu2.VisionInterface.DejaVu2Nodes import Viewer
        if isinstance(c.port2.node, Viewer):
            if colorMap.viewer != c.port2.node.vi:
                colorMap.configure(viewer=c.port2.node.vi)
                colorMap.showLegend()
            self.outputData(colorMap=colorMap, legend=colorMap.legend)
            break
        elif isinstance(c.port2.node, MacroOutputNode):
            for macroOutputPort in c.port2.node.macroNode.outputPorts:
                for c2 in macroOutputPort.connections:
                    if isinstance(c2.port2.node, Viewer):
                        if colorMap.viewer != c2.port2.node.vi:
                            colorMap.configure(viewer=c2.port2.node.vi)
                            colorMap.showLegend()
                        self.outputData(colorMap=colorMap, legend=colorMap.legend)
                        break
    else: #we didn't break
        self.outputData(colorMap=colorMap)

    if (values is not None) and (len(values) > 0) :
        lCol = colorMap.Map(values)
        if lCol is not None:
            self.outputData(mappedColors=lCol.tolist())
    elif len(colorMap.ramp) == 1:
        self.outputData(mappedColors=colorMap.ramp[0])
    else:
        self.outputData(mappedColors=colorMap.ramp)
"""
        self.setFunction(code)


class ColorByRamp(NetworkNode):
    """ DEPRECATED : USE A COLORMAP INSTEAD  
    (and unbind the colormap in the colormap parampanel)
    
    Turn a list of values into a list of colors using a colormap.
This node allows you to re-use an already existing colormap.
(To use/create a new colormap, use the node colormap instead)

Input Ports
    values: a list of values
    colormap: a color map (array of rgba values), typically provided by
              a 'Color Map' node, Defaults to an RGB ramp.
    
Output Ports
    colors: a list of colors.
"""

    def __init__(self, name='Color By Ramp', **kw):
        
        kw['name'] = name
        apply( NetworkNode.__init__, (self,), kw )
        #self.readOnly = 1

        from DejaVu2.colorTool import RGBARamp, Map
        self.defaultMap = ColorMap(name='default', ramp=RGBARamp())

        ip = self.inputPortsDescr
        ip.append(datatype='vector', name='values')
        ip.append(datatype='ColorMapType', required=False, name='colormap')

        op = self.outputPortsDescr
        op.append(datatype='colorfloat3or4(0)', name='colors')

        code = """def doit(self, values, colormap):
    # self refers to the node    
    if colormap is None:
        colormap = self.defaultMap
        
    if len(values):
        lCol = colormap.Map(values)
        if lCol:
            self.outputData(colors=lCol.tolist())
"""

        self.setFunction(code)



class GeomsToVRML2(NetworkNode):
    """This node inputs a DejaVu2 viewer object and converts
all visible geoms in vrml2 format. The output is a list of
strings which can be saved using a save text node."""

    def __init__(self, name='Get VRML2', **kw):
        kw['name'] = name
        apply( NetworkNode.__init__, (self,), kw )

        from DejaVu2.DataOutput import OutputVRML2
        self.V = OutputVRML2()

        code ="""def doit(self, viewer):
    # FIXME: expose these fields either as input ports or as GUI
    co = 1 # 0 or 1
    sn = 0 # 0 or 1
    up = 0 # 0 or 1
    sq = 2 # 0, 1, 2, 3, etc
    cq = 10 # 3, 4, 5, 6, etc

    vrml2 = self.V.getVRML2(viewer.rootObject, complete=co,
                            normals=sn, usePROTO=up,
                            sphereQuality=sq, cylinderQuality=cq)

    if vrml2 is not None and len(vrml2) != 0:
        self.outputData(vrml2=vrml2)\n"""
        
        if code: self.setFunction(code)

        self.inputPortsDescr.append(datatype='viewer', name='viewer')
        self.outputPortsDescr.append(datatype='list', name='vrml2')


class TransferFuncEditor(NetworkNode):
    """ Transfer Function Editor widget used to edit color and opacity lookup tables of
    volumetric data. This node can be connected to UT-VolRen node from vollib library.
    Refer to 'Help' menu of the node for more info. """

    def color_cb(self, values):
        if self.utvolGeom: 
            self.utvolGeom.setVolRenColors(values)
        else:
            self.need_update = True

    def alpha_cb(self, values):
        if self.utvolGeom: 
            self.utvolGeom.setVolRenAlpha(values)
        else:
            self.need_update = True
        
    def __init__(self, name='TransferFuncEditor', **kw):
        kw['name'] = name
        apply( NetworkNode.__init__, (self,), kw )

        self.utvolGeom = None
        self.need_update = False
        self.tm = None   # instanciated when node is added to a network
        self.isWithdrawn = 0

        self.mouseAction['<Double-Button-1>'] = self.show_hide
        
        ip = self.inputPortsDescr
        ip.append(datatype='geom', name='utvolGeom')

        op = self.outputPortsDescr
        op.append(datatype='None', name='cmap')

        code = """def doit(self, utvolGeom):
    self.utvolGeom = utvolGeom
    if self.need_update:
        self.tm.applylut_cb()
        self.need_update = False
"""
        self.setFunction(code)
        

    def show_hide(self, event=None):
        if self.isWithdrawn:
            self.tm.master.deiconify()
            self.isWithdrawn = 0
        else:
            self.tm.master.withdraw()
            self.isWithdrawn = 1


    def beforeAddingToNetwork(self, net):
        import Tkinter
        root = Tkinter.Toplevel()
        from mglutil.gui.BasicWidgets.Tk.tablemaker import TableManager
        self.tm = TableManager(None, master=root, ymaxval=255, xminval = 0,
                         alphaCallback=self.alpha_cb,
                         colorCallback=self.color_cb,
                         )
        self.tm.master.protocol('WM_DELETE_WINDOW', self.show_hide)


    def afterRemovingFromNetwork(self):
        NetworkNode.afterRemovingFromNetwork(self)
        self.network.canvas.delete(self.iconTag)
        if self.tm:
            self.tm.master.destroy()
            del self.tm



class ScaleLegend(NetworkNode):
    """This node creates a line geometry for displaying a scale in a
DejaVu2 camera.
        ________________________________________
        |               10.00                   |


Input:
    length       -- the length of the scale line
    tickLength   -- length of the vertical tick at the end of the scale
    labelOffsetx -- horizontal offset of the label from the center of the line
    labelOffsety -- vertical offset of the label from the the line
    labelFont    -- font for the label
    
Output:
    scaleLines -- line geometry
    scaleLabel -- text geometry
"""

    def __init__(self, name='ScaleLegend', **kw):
        kw['name'] = name
        apply( NetworkNode.__init__, (self,), kw )

        self.v = [[0.,0.,0.], [10.,0.,0.], [0.,-1.,0.], [10.,-1.,0.]]
        self.f = ( (0,1), (0,2), (1,3) )

        #from DejaVu2.Labels import Labels
        from DejaVu2.glfLabels import GlfLabels
        vt = [[5., 0, 0]]
        self.textGeom = GlfLabels('ScaleLabel', vertices=vt, labels=['5.0'])
        
        from DejaVu2.IndexedPolylines import IndexedPolylines
        self.lineGeom = IndexedPolylines(
            'scaleLine', vertices=self.v, faces=self.f, lineWidth=2)


        ip = self.inputPortsDescr
        ip.append(datatype='float', name='length')
        ip.append(datatype='float', name='tickLength')
        ip.append(datatype='float', name='labelOffsetx')
        ip.append(datatype='float', name='labelOffsety')
        ip.append(datatype='float', name='labelFont')


        self.widgetDescr['length'] = {
            'class':'NEThumbWheel', 'initialValue':10,
            'labelCfg':{'text':'length:'}, 'width':80, 'height':20, 
            'type':'float', 'oneTurn':10, 'lockBMin':1, 'min':0,
            'wheelPad':2 }

        self.widgetDescr['tickLength'] = {
            'class':'NEThumbWheel', 'initialValue':1,
            'labelCfg':{'text':'length:'}, 'width':80, 'height':20, 
            'type':'float', 'oneTurn':10, 'lockBMin':1, 'min':0,
            'wheelPad':2 }

        self.widgetDescr['labelOffsetx'] = {
            'class':'NEThumbWheel', 'initialValue':0.,
            'labelCfg':{'text':'label Offset x:'}, 'width':80, 'height':20, 
            'type':'float', 'oneTurn':10,
            'wheelPad':2 }

        self.widgetDescr['labelOffsety'] = {
            'class':'NEThumbWheel', 'initialValue':-2.,
            'labelCfg':{'text':'label Offset y:'}, 'width':80, 'height':20, 
            'type':'float', 'oneTurn':10,
            'wheelPad':2 }

        self.widgetDescr['labelFont'] = {
            'class':'NEComboBox',
            'choices':GlfLabels.fontList, #self.textGeom.BMfonts,
            'entryfield_entry_width':8,
            'initialValue':GlfLabels.fontList[0],
            'labelCfg':{'text':'label font:'},
            }

        op = self.outputPortsDescr
        op.append(datatype='geom', name='scaleLines')
        op.append(datatype='geom', name='scaleLabel')

        
        code = """def doit(self, length, tickLength, labelOffsetx,
labelOffsety, labelFont):
    
    self.v[1][0] = length
    self.v[2][1] = -tickLength
    self.v[3][0] = length
    self.v[3][1] = -tickLength
    self.lineGeom.Set(vertices = self.v)

    vt = [ [length*0.5+labelOffsetx, labelOffsety, 0.] ]
    lab = '%.2f'%length
    self.textGeom.Set(vertices=vt, labels=[lab], font=labelFont)

    self.outputData(scaleLines=self.lineGeom, scaleLabel=self.textGeom)
"""
        self.setFunction(code)



class ReadSTL(NetworkNode):
    """ ************* DEPRECATED, USE GeomsFromFile instead **************

Read ASCII STL Files, convert them to IndexedPolygon geometries.
Input:  filename 
Output: a list of geometries
    """

    def __init__(self, name='Read STL ASCII', **kw):

        #warnings.warn("Read STL ASCII is deprecated, use 'geoms from file' instead")

        kw['name'] = name
        apply( NetworkNode.__init__, (self,), kw )

        from DejaVu2.DataInput import ReadASCIISTL
        self.Parser = ReadASCIISTL()
        
        self.widgetDescr['filename'] = {
            'class':'NEEntryWithFileBrowser', 'master':'node', 'width':10,
            'initialValue':'',
            'filetypes':[('stl','*.stl'),('all','*')],
            'labelCfg':{'text':'Filename: '}
            }

        ip = self.inputPortsDescr
        ip.append(datatype='string', name='filename')

        op = self.outputPortsDescr
        op.append(datatype='list', name='geoms')

        code = """def doit(self, filename):
    if not filename:
        return

    lines = self.Parser.read(filename)
    geoms = self.Parser.doit(lines)
    if geoms and len(geoms):
        self.outputData(geoms=geoms)\n"""

        self.setFunction(code)


class ReadSTL_B(NetworkNode):
    """ ************* DEPRECATED, USE GeomsFromFile instead **************

Read Binary STL File, convert to IndexedPolygon geometry.
Input:  filename 
Output: a DejaVu2 IndexedPolygon geometry object
"""

    def __init__(self, name='Read STL Binary', **kw):

        #warnings.warn("Read STL Binary is deprecated, use 'geoms from file' instead")

        kw['name'] = name
        apply( NetworkNode.__init__, (self,), kw )

        from DejaVu2.DataInput import ReadBinarySTL
        self.Parser = ReadBinarySTL()
        
        self.widgetDescr['filename'] = {
            'class':'NEEntryWithFileBrowser', 'master':'node', 'width':10,
            'initialValue':'',
            'filetypes':[('stl','*.stl'),('all','*')],
            'labelCfg':{'text':'Filename: '}
            }

        ip = self.inputPortsDescr
        ip.append(datatype='string', name='filename')

        op = self.outputPortsDescr
        op.append(datatype='geom', name='geom')

        code = """def doit(self, filename):
    if not filename:
        return

    geom = None
    geom = self.Parser.doit(filename)
    if geom and geom is not None:
        self.outputData(geom=geom)\n"""

        self.setFunction(code)


class ComputeRMSD(NetworkNode):
    """Compute RMSD

Input:\tcoords: 3-D coordinates
\trefCoords: 3-D coordinates
Output: rmsd (float)"""

    def __init__(self, name='Compute RMSD', **kw):
        kw['name'] = name
        apply( NetworkNode.__init__, (self,), kw )

        from mglutil.math.rmsd import RMSDCalculator
        self.RMSD = RMSDCalculator()

        ip = self.inputPortsDescr
        ip.append(datatype='coordinates3D', name='coords')
        ip.append(datatype='coordinates3D', required=False, name='refCoords')

        op =self.outputPortsDescr 
        op.append(datatype='float', name='rmsd')
       
        code = """def doit(self, coords, refCoords):
    rmsd = None
    if refCoords:
        self.RMSD.setRefCoords(refCoords)
    if coords and self.RMSD.refCoords:
        rmsd = self.RMSD.computeRMSD(coords)
    if rmsd is not None:
        self.outputData(rmsd=rmsd)\n"""
     
        self.setFunction(code)



class SelectAxis(NetworkNode):
    def __init__(self, constrkw={}, name='SelectAxis', **kw):
        kw['name']=name
        apply( NetworkNode.__init__, (self,), kw )

        ip = self.inputPortsDescr
        ip.append(datatype='boolean', name ='x-axis')
        ip.append(datatype='boolean', name ='y-axis')
        ip.append(datatype='boolean', name ='z-axis')
        
        wd = self.widgetDescr
        wd['x-axis'] = {
            'class':'NECheckButton', 'master':'node',
            'initialValue':0,  'labelGridCfg':{'sticky':'we'},
            'labelCfg':{'text':'x-axis'},
            }
        wd['y-axis'] = {
            'class':'NECheckButton', 'master':'node',
            'initialValue':1,  'labelGridCfg':{'sticky':'we'},
            'labelCfg':{'text':'y-axis:'},
            }
        wd['z-axis'] = {
            'class':'NECheckButton', 'master':'node',
            'initialValue':0,  'labelGridCfg':{'sticky':'we'},
            'labelCfg':{'text':'z-axis:'},
            }

        self.outputPortsDescr.append(datatype='list',name='rotationAxis')

        code = """def doit(self,xaxis,yaxis,zaxis):
    if xaxis: x_rot=1
    else: x_rot=0
    if yaxis: y_rot=1
    else: y_rot=0
    if zaxis: z_rot=1
    else: z_rot=0
    if not xaxis and not yaxis and not zaxis:
        print"no axis select for rotation!  Using y axis by default"
        y_rot=1

    rotation=[x_rot,y_rot,z_rot]
    
    self.outputData(rotationAxis=rotation)
"""
        self.setFunction(code)
        
class DistanceMatrix(NetworkNode):
    """
Compute the NxN distance between 2 sets of 3D points

Input:
    points:  3D coordinates

Output:
    dm:  2D Numeric array of distances
"""
    def __init__(self, name='DM', **kw):
        kw['name'] = name
	apply( NetworkNode.__init__, (self,), kw)
	
	code = """def doit(self, points):
	import numpy.oldnumeric as Numeric
	l = len(points)
	dm = Numeric.zeros( (l,l), 'd')
	points = Numeric.array(points).astype('d')
	for i in range(l):
	    dist = points[i] - points
	    dist = dist*dist
	    dm[i] = Numeric.sqrt(Numeric.sum(dist, 1))
	self.outputData(dm=dm)
"""

	if code: self.setFunction(code)
	
	self.inputPortsDescr.append({'name': 'points',
                                     'datatype': 'coordinates3D'})
	self.outputPortsDescr.append({'name': 'dm', 'datatype': 'None'})


class DDM(NetworkNode):
    """
Compute the difference between two distance matrices

Input:
    points1: sequence of 3D points
    points2: sequence of 3D points

Output
    ddm: 2D Numeric array of differences between distances
"""
    
    def __init__(self, name='DDM', **kw):
        kw['name'] = name
	apply( NetworkNode.__init__, (self,), kw)
	
	code = """def doit(self, points1, points2):
	import numpy.oldnumeric as Numeric
	result = points1-points2
        self.outputData(ddm=result)
"""
		
	if code: self.setFunction(code)
	
	self.inputPortsDescr.append({'name': 'points1',
                                     'datatype': 'None'})
	self.inputPortsDescr.append({'name':'points2', 'datatype':
                                     'None'})
	self.outputPortsDescr.append({'name':'ddm', 'datatype': 'None'})



class ConnectedComponents0(NetworkNode):
    """ Node to split connected surfaces.
    Input ports:
          -geometry or a list of geometries(the first object in the list
            will be used as input)
          -name of output geometries . Index (0-n) will be attached to the name
    Output:
          - list of Indexed Polygons    
    """
    def __init__(self, name='ConnectedComponents', **kw):
        kw['name'] = name
	apply( NetworkNode.__init__, (self,), kw)
        self.widgetDescr['outgeomname'] = {
            'class':'NEEntry', 'master':'node', 'width':10,
            'labelCfg':{'text':'name:'},
            'initialValue':'polygon',
            }

        ip = self.inputPortsDescr
        ip.append(name='ingeoms', datatype='geom', required=True)
        ip.append(name='outgeomname', datatype='str', required=False)
	op = self.outputPortsDescr
        op.append(name='outgeoms', datatype='geom')

        code = """def doit(self, ingeoms, outgeomname):
        try:
            len(ingeoms)
            geometry = ingeoms[0]
        except:
            geometry = ingeoms
        faces = geometry.getFaces()
        verts = geometry.getVertices()
        if not outgeomname:
            outgeomname = 'polygons'
        from time import time
        t1 = time()
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

        t2 = time()
        print 'time to find %d connected components : %.2f'%(len(newfaces), t2-t1)
        from DejaVu2.IndexedPolygons import IndexedPolygons
        outgeoms = []
        for i, newfs  in enumerate(newfaces):
            newname = outgeomname+str(i)
            print 'len faces of %s is %d'%(newname, len(newfs))
            obj = IndexedPolygons(newname, vertices = newverts[i], faces = newfs)
            outgeoms.append(obj)
        self.outputData(outgeoms=outgeoms)
"""
		
	self.setFunction(code)

try:
    bhtreelibFound = True
    from bhtree import bhtreelib
except ImportError:
    bhtreelibFound = False

if bhtreelibFound:
    class ClosestPoints(NetworkNode):
        """
    identify the closest point in 2 sets fo 3D points

    Input:
        points1: sequence of 3D points
        points2: sequence of 3D points
        mini: minimum distance (default 0.0)
        maxi: maximum distance (default 10.0)
        incr: distance increment for searching close points (default 1.0)
        
    Output
        points: (i,j) 0-based indices in set 1 and 2
        distance: distance between the points
    """

        def __init__(self, name='ClosestPoints', **kw):
            kw['name'] = name
            apply( NetworkNode.__init__, (self,), kw)

            code = """def doit(self, points1, points2, mini=0.0, maxi=10., incr=1.0):
    if len(points1) >= len(points2):
        set1 = points1
        set2 = points2
        case = 1
    else:
        set1 = points2
        set2 = points1
        case = 2
    import numpy
    
    bht = bhtreelib.BHtree(numpy.array(set1).astype('f'), None, 10)
    results = numpy.zeros(len(set1), 'i')
    dist2 = numpy.zeros(len(set1), 'f')
    minimum = 99999999.
    pair = [-1,-1]
    for j, pt in enumerate(set2):
        cutOff = mini
        nb = 0
        while nb==0 and cutOff < maxi:
            nb = bht.closePointsDist2(tuple(pt), cutOff, results, dist2)
            if nb:
                mind = numpy.min(dist2[:nb])
                if mind < minimum**2:
                    minimum = sqrt(mind)
                    i = results[numpy.argmin(dist2[:nb])]
                    pair = [i, j]
                break
            cutOff += incr
    if case==2:
        pair = [pair[1], pair[0]]
    self.outputData(points=pair, distance=minimum)
    """

            if code: self.setFunction(code)

            self.widgetDescr['mini'] = {
                'class':'NEThumbWheel',
                'master':'node',
                'initialValue':0.,
                'labelCfg':{'text':'mini:'},
                'width':80, 'height':20, 
                'type':'float', 'oneTurn':10, 'lockBMin':1, 'min':1,
                'wheelPad':2 }
            self.widgetDescr['maxi'] = {
                'class':'NEThumbWheel',
                'master':'node',
                'initialValue':10.,
                'labelCfg':{'text':'maxi.:'},
                'width':80, 'height':20, 
                'type':'float', 'oneTurn':10, 'lockBMin':1, 'min':1,
                'wheelPad':2 }
            self.widgetDescr['incr'] = {
                'class':'NEThumbWheel',
                'master':'node',
                'initialValue':1.,
                'labelCfg':{'text':'incr.:'},
                'width':80, 'height':20, 
                'type':'float', 'oneTurn':10, 'lockBMin':1, 'min':1,
                'wheelPad':2 }
            self.inputPortsDescr.append({'name': 'points1', 'datatype': 'None'})
            self.inputPortsDescr.append({'name':'points2', 'datatype': 'None'})
            self.inputPortsDescr.append({'name':'mini', 'datatype': 'float', 'required':False})
            self.inputPortsDescr.append({'name':'maxi', 'datatype': 'float', 'required':False})
            self.inputPortsDescr.append({'name':'incr', 'datatype': 'float', 'required':False})
            self.outputPortsDescr.append({'name':'points', 'datatype': 'list'})
            self.outputPortsDescr.append({'name':'distance', 'datatype': 'float'})
   

    class PairsCloserThan(NetworkNode):
        """
    identify pairs of 3D points between 2 sets that are closers than maxi

    Input:
        points1: sequence of 3D points
        points2: sequence of 3D points
        maxi: maximum distance between points
        
    Output
        points: [(i,j)] list of 0-based indices in set 1 and 2
        distances: distance between the points
    """

        def __init__(self, name='ClosestPoints', **kw):
            kw['name'] = name
            apply( NetworkNode.__init__, (self,), kw)

            code = """def doit(self, points1, points2, maxi=10.):
    if len(points1) >= len(points2):
        set1 = points1
        set2 = points2
        case = 1
    else:
        set1 = points2
        set2 = points1
        case = 2
    import numpy
    
    bht = bhtreelib.BHtree(numpy.array(set1).astype('f'), None, 10)
    results = numpy.zeros(len(set1), 'i')
    dist2 = numpy.zeros(len(set1), 'f')
    distances = []
    pairs = []
    for j, pt in enumerate(set2):
        nb = bht.closePointsDist2(tuple(pt), maxi, results, dist2)
        for n in range(nb):
            distances.append(sqrt(dist2[n]))
            if case==1:
                pairs.append( (results[n], j) )
            else:
                pairs.append( (j, results[n]) )

    self.outputData(distances=distances, points=pairs)
    """

            if code: self.setFunction(code)

            self.widgetDescr['maxi'] = {
                'class':'NEThumbWheel',
                'master':'node',
                'initialValue':10.,
                'labelCfg':{'text':'maxi.:'},
                'width':80, 'height':20, 
                'type':'float', 'oneTurn':10, 'lockBMin':1, 'min':1,
                'wheelPad':2 }

            self.inputPortsDescr.append({'name': 'points1', 'datatype': 'None'})
            self.inputPortsDescr.append({'name':'points2', 'datatype': 'None'})
            self.inputPortsDescr.append({'name':'maxi', 'datatype': 'float'})
            self.outputPortsDescr.append({'name':'points', 'datatype': 'list'})
            self.outputPortsDescr.append({'name':'distances', 'datatype': 'float'})
   
from Vision.VPE import NodeLibrary
vizlib = NodeLibrary('3D Visualization', '#feec70')

vizlib.addNode(ColorEditor, 'Color Editor', 'Input')
vizlib.addNode(ColorMapNE, 'Color Map', 'Input')

vizlib.addNode(TransferFuncEditor, 'TF Editor', 'Input')

vizlib.addNode(ColorChooser, 'Color Chooser', 'deprecated')
vizlib.addNode(ReadSTL, 'Read STL ASCII', 'deprecated') # deprecated
vizlib.addNode(ReadSTL_B, 'Read STL Binary', 'deprecated') # deprecated
vizlib.addNode(readIndexedPolygon, 'readIndexedPolygon', 'deprecated') # deprecated
vizlib.addNode(ColorByRamp, 'Color', 'deprecated')
vizlib.addNode(ReparentGeom, 'Reparent Geom', 'deprecated')

vizlib.addNode(GeomsProperties, 'geoms properties', 'Mapper')

vizlib.addNode(RestoreState, 'RestoreState', 'Input')
vizlib.addNode(GeomOptions, 'Set Geom Options', 'Input')

vizlib.addNode(Viewer, 'Viewer', 'Output')
vizlib.addNode(Redraw, 'Redraw', 'Output')
vizlib.addNode(StopAutoRedraw, 'Stop Auto Redraw', 'Output')
vizlib.addNode(StartAutoRedraw, 'Start Auto Redraw', 'Output')
vizlib.addNode(OneRedraw, 'One Redraw', 'Output')
vizlib.addNode(writeIndexedPolygon, 'writeIndexedPolygon', 'Output')
vizlib.addNode(writeCurvPly, 'writeCurvPly', 'Output')
vizlib.addNode(ImageViewerNode, 'ImageViewer', 'Output')

from DejaVu2.VisionInterface.GeometryNodes import GeomsFromFile
vizlib.addNode(GeomsFromFile, 'geoms from file', 'Geometry')
from DejaVu2.VisionInterface.GeometryNodes import BoxNE
vizlib.addNode(BoxNE, 'box', 'Geometry')
from DejaVu2.VisionInterface.GeometryNodes import StickerTextNE
vizlib.addNode(StickerTextNE, 'sticker text', 'Geometry')
from DejaVu2.VisionInterface.GeometryNodes import StickerImageNE
vizlib.addNode(StickerImageNE, 'sticker image', 'Geometry')
from DejaVu2.VisionInterface.GeometryNodes import OneTexturedQuadNE
vizlib.addNode(OneTexturedQuadNE, 'oneTexturedQuad', 'Geometry')
from DejaVu2.VisionInterface.GeometryNodes import IndexedPolygonsNE
vizlib.addNode(IndexedPolygonsNE, 'IndPolygons', 'Geometry')
from DejaVu2.VisionInterface.GeometryNodes import IndexedPolylinesNE
vizlib.addNode(IndexedPolylinesNE, 'IndPolylines', 'Geometry')
from DejaVu2.VisionInterface.GeometryNodes import Cylinders
vizlib.addNode(Cylinders, 'Cylinders', 'Geometry')
from DejaVu2.VisionInterface.GeometryNodes import AxisNE
vizlib.addNode(AxisNE, 'axis', 'Geometry')
from DejaVu2.VisionInterface.GeometryNodes import Spheres
vizlib.addNode(Spheres, 'Spheres', 'Geometry')
from DejaVu2.VisionInterface.GeometryNodes import Points
vizlib.addNode(Points, 'points', 'Geometry')
from DejaVu2.VisionInterface.GeometryNodes import CrossSet
vizlib.addNode(CrossSet, 'CrossSet', 'Geometry')
from DejaVu2.VisionInterface.GeometryNodes import GeomContainer
vizlib.addNode(GeomContainer, 'GeomContainer', 'Geometry')
from DejaVu2.VisionInterface.GeometryNodes import Ellipsoids
vizlib.addNode(Ellipsoids, 'Ellipsoids', 'Geometry')
from DejaVu2.VisionInterface.GeometryNodes import Array2DToHeightFieldGeom
from DejaVu2.VisionInterface.GeometryNodes import GlfLabelsNE
vizlib.addNode(GlfLabelsNE, 'GlfLabels', 'Geometry')
#from DejaVu2.VisionInterface.GeometryNodes import GlutLabelsNE
#vizlib.addNode(GlutLabelsNE, 'GlutLabels', 'Geometry')

vizlib.addNode(Array2DToHeightFieldGeom, 'HeightField', 'Mapper')
vizlib.addNode(Curvature, 'Curvature', 'Mapper')
vizlib.addNode(SetInstances, 'Set Instances', 'Mapper')
vizlib.addNode(GroupGeomsNE, 'Group Geoms', 'Mapper')
vizlib.addNode(Texture2D, 'Texture2D', 'Mapper')
vizlib.addNode(Texture1D, '1D Texture Mapping', 'Mapper')
vizlib.addNode(Texture1DCoords, '1D Texture Coordinates', 'Mapper')
vizlib.addNode(GeomsToVRML2, 'Get VRML2', 'Mapper')
vizlib.addNode(GyrationSphere, 'Gyration Sphere', 'Mapper')
vizlib.addNode(ComputeRMSD, 'Compute RMSD', 'Mapper')
vizlib.addNode(SMFtoGeom, 'SMF to Geom', 'Mapper')
vizlib.addNode(GeomtoSMF, 'Geom to SMF', 'Mapper')
vizlib.addNode(ScaleLegend, 'Scale', 'Mapper')
vizlib.addNode(DistanceMatrix, 'DistanceMatrix', 'Mapper')
vizlib.addNode(DDM, 'DDM', 'Mapper')
vizlib.addNode(PolyhedronVolumeArea, 'PolyhedronVolumeArea', 'Mapper')
#vizlib.addNode(AccumPickedVertices, 'AccumPickedVert', 'Mapper')

from DejaVu2.VisionInterface.RotateScene import RotateScene
vizlib.addNode(RotateScene, 'RotateScene', 'Macro')
from DejaVu2.VisionInterface.StereoSep import StereoSep
vizlib.addNode(StereoSep, 'StereoSep', 'Macro')
from DejaVu2.VisionInterface.MapPotOnGeom import MapPotOnGeom
vizlib.addNode(MapPotOnGeom, 'Map Pot On Geom', 'Macro')

vizlib.addNode(getSurfaceVFN, 'getSurfaceVFN', 'Filter')

vizlib.addNode(SelectGeometry, 'Choose Geom', 'Filter')
vizlib.addNode(SelectMultipleGeometry, 'Choose Geoms', 'Filter')


vizlib.addNode(NPR, 'NPR', 'Filter')
vizlib.addNode(QslimExt, 'QslimExt', 'Filter')
vizlib.addNode(RemoveDuplicatedVerticesNE, 'RemoveDupVert','Filter')
vizlib.addNode(CenterOnPickedVertex, 'CenterPick','Filter')
vizlib.addNode(CenterOnVertex, 'CenterVertex','Filter')


try:
    from geomutils.geomalgorithms import removeDuplicatedVertices
    vizlib.addNode(removeDuplicatedVerticesC, 'RemoveDupVert(C++)','Filter')
except:
    pass

try:
    from QSlimLib import qslimlib
    vizlib.addNode(QSlim, 'QSlim', 'Filter')
except:
    pass

try:
    from QSlimLib import qslimlib
    from DejaVu2.VisionInterface.GeometryNodes import DecimateGeom
    vizlib.addNode(DecimateGeom, 'DecimateGeom', 'Filter')
except:
    pass

if bhtreelibFound:
    vizlib.addNode(decimate3DPoints, 'decimate3DPoints', 'Filter')
    vizlib.addNode(DistFromSphereToGeom, 'DistFromSphereToGeom', 'Filter')
    
try:
    import multires
    vizlib.addNode(Decimate, 'Polygon Reduction', 'Filter')
except:
    pass

try:
    #from pyefit import efit
    #from geomutils import geomutilslib
    from geomutils import efitlib
    vizlib.addNode(EllipsoidFit, 'EllipsoidFit', 'Mapper')
except:
    pass

try:
    from mglutil.util.qhullUtils import findQhullModule
    pth = findQhullModule('qconvex')
    if pth:
        vizlib.addNode(QConvexHull, 'Convex Hull', 'Mapper')
except:
    pass

if bhtreelibFound:
    vizlib.addNode(ClosestPoints, 'ClosestPoints', 'Mapper')
    vizlib.addNode(PairsCloserThan, 'PairsCloserThan', 'Mapper')
    
from DejaVu2.VisionInterface.GeometryNodes import ConnectedComponents
vizlib.addNode(ConnectedComponents, "ConnectedComponents", "Filter")

vizlib.addWidget(NEColorMap)
vizlib.addWidget(NEColorWheel)
vizlib.addWidget(NEColorEditor)
vizlib.addWidget(NEDejaVu2GeomOptions)


UserLibBuild.addTypes(vizlib, 'DejaVu2.VisionInterface.DejaVu2Types')

try:
    UserLibBuild.addTypes(vizlib, 'Vision.PILTypes')
except:
    pass
