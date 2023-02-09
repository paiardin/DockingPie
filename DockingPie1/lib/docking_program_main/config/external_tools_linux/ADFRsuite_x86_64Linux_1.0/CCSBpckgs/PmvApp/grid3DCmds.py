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

#
# $Header: /mnt/raid/services/cvs/PmvApp/grid3DCmds.py,v 1.1.4.1 2017/07/13 20:52:52 annao Exp $
#
# $Id: grid3DCmds.py,v 1.1.4.1 2017/07/13 20:52:52 annao Exp $
#
import os
import numpy
from PmvApp.Pmv import MVCommand



class GridMaps:

    def __init__(self):
        self.maps = {}
        self.center = None
        self.spacing = None
        self.flexResStr = None
        self.receptorFilename = None
        self.mapsFolder = None
        self.numGridPoints = None
        
    def mapsFromZip(self, zipFile):
        """
        Unzip a map file in a temporary location , read maps and create grid objects.
        """
        from mglutil.util.packageFilePath import getResourceFolderWithVersion
        tmpFolder = os.path.join(getResourceFolderWithVersion(), 'tmp')
        if not os.path.exists(tmpFolder):
            os.mkdir(tmpFolder)
        assert os.path.isdir(tmpFolder)
        from zipfile import ZipFile
        zf = ZipFile(zipFile)
        zf.extractall(tmpFolder)
        self.mapsFolder = folder = os.path.join(tmpFolder, os.path.splitext(os.path.basename(zipFile))[0])
        #tpointsFilename = os.path.join(folder, "translationPoints.npy")
        # find receptor file name i.e PDBQT file which is NOT receptor.pdbqt
        for filename in zf.namelist():
            if filename.endswith('.pdbqt'):
                filename = os.path.split(filename)[1]
                recname = os.path.splitext(filename)[0]
                if recname != 'receptor':
                    break
        self.receptorFilename = os.path.join(folder, '%s.pdbqt'%recname) 
        for mapFileName in zf.namelist():
            self.readMap(os.path.basename(mapFileName), folder)

    def mapsFromListOfNames(self, mapFiles, receptorFile=None):
        if receptorFile is not None:
            assert receptorFile.endswith('.pdbqt')
            self.receptorFilename = receptorFile
        for item in mapFiles:
            mapFileName = os.path.basename(item)
            folder = os.path.abspath(os.path.dirname(item))
            self.readMap(mapFileName, folder)

    def readMap(self, mapFileName, folder):
        from ADFRcc.adfr import GridMap
        w = mapFileName.split('.')
        if w[-1]=='map':
            _map = GridMap()
            _f, name = os.path.split(mapFileName)
            _map.loadFromMapFile(w[-2], folder.encode('ascii', 'replace'),
                                name.encode('ascii', 'replace'))
            if not (w[-2]=='e' or w[-2]=='d') and self.center is None:# read 'e' map to get center, spacing, and flexres
                self.center = _map.getCenterPy()
                self.spacing = _map.getDistBetweenGridPoints()
                self.flexResStr = _map.getFlexRes()
                self.numGridPoints = _map.getNumGridPointsPy()
                #flexRes = flexResStr2flexRes(flexResStr)
            self.maps[name] = _map


class ReadGrid(MVCommand):
    """ Reads a zip file containig map file, receptor file.
    Saves Grid3D object in self.app().grids3D[gridFile]"""

    def __init__(self):
        MVCommand.__init__(self)
        self.iterateOverFirstArgument = False
    

    def expandArg0(self, obj):
        #if isinstance(obj, list) : return obj[0]
        #else:
        if isinstance(obj, unicode):
            return str(obj)
        return obj

    def checkArguments(self, maps):
        if isinstance(maps, (str, unicode)):
            assert os.path.exists(maps)
            ext = os.path.splitext(maps)[1]
            assert ext in (".zip", ".map")
            if ext == ".map":
                maps = [maps]
        else:
            assert isinstance(maps, (list, tuple))
            for name in maps:
                assert os.path.exists(name)
        return (maps,), {}


    def doit(self, maps):
        self.mapsObj= mapsObj = GridMaps()
        if  isinstance(maps, (str, unicode)):
            mapsObj.mapsFromZip(maps)
        else:
            mapsObj.mapsFromListOfNames(maps)
        receptorFilename = mapsObj.receptorFilename
        #if receptorFilename == None:
        #    receptor = self.app().readMolecules([str(receptorFilename)])[0]
        center = mapsObj.center
        nx, ny, nz = mapsObj.numGridPoints
        spacing = mapsObj.spacing
        sx = spacing*(nx-1)
        sy = spacing*(ny-1)
        sz = spacing*(nz-1)
        stepSize = (spacing, spacing, spacing)
        origin = (center[0]-(nx/2)*spacing, center[1]-(ny/2)*spacing,
                   center[2]-(nz/2)*spacing)
        for name, _map in mapsObj.maps.items():
            data = _map.getGridDataPy()
            grid_file = _map.getMapFilePath()

            header = {"GRID_PARAMETER_FILE":"",
                      "GRID_DATA_FILE":"",
                      "MACROMOLECULE":receptorFilename,
                      "title": "AutoGrid from %s"%grid_file
                      }
            from Volume.Grid3D import Grid3DD
            grid3D = Grid3DD(data, origin, stepSize, header)
            #grid3D.receptor = receptor
            self.addGrid(grid3D, name, (sx,sy,sz), center,)

    def addGrid(self, grid3D, name, sides, center):        
        grid3D.origin_copy = grid3D.origin
        grid3D.stepSize_copy = grid3D.stepSize
        mini, maxi, mean, std = grid3D.stats()
        grid3D.mini = mini
        grid3D.maxi = maxi
        grid3D.mean = mean
        grid3D.std = std
        if name == None:
            name = str(grid3D) 
        if self.app().grids3D.has_key(name):
            name += "_" 
        def returnStringRepr():
            return None, "\"" + name + "\""
        grid3D.returnStringRepr = returnStringRepr
        if not hasattr(grid3D,'geomContainer'):
            grid3D.geomContainer = {}
            from DejaVu2.Geom import Geom
            g = Geom(name)
            grid3D.master_geom = g
            self.app().gui().viewer.AddObject(g)
            grid3D.geomContainer['IsoSurf'] = {}
            grid3D.geomContainer['OrthoSlice'] = {}
            IsoSurf = Geom('IsoSurf')                
            OrthoSlice = Geom('OrthoSlice')
            grid3D.IsoSurf = IsoSurf
            grid3D.OrthoSlice = OrthoSlice
            self.app().gui().viewer.AddObject(IsoSurf,parent=g)
            self.app().gui().viewer.AddObject(OrthoSlice,parent=g)
        self.app().grids3D[name] = grid3D

    def addBox(self, grid3D, center, sides):
        ## from DejaVu2.Box import Box
        ## box = Box('BoundingBox')

        ## pt1 = grid3D.origin
        ## dims = grid3D.data.shape
        ## pt2 = [pt1[0]+(grid3D.stepSize[0]*(dims[0]-1)),
        ##        pt1[1]+(grid3D.stepSize[1]*(dims[1]-1)),
        ##        pt1[2]+(grid3D.stepSize[2]*(dims[2]-1))]
        ## if grid3D.crystal:
        ##     ptList=((pt2[0],pt2[1],pt1[2]),
        ##             (pt1[0],pt2[1],pt1[2]),
        ##             (pt1[0],pt1[1],pt1[2]),
        ##             (pt2[0],pt1[1],pt1[2]),
        ##             (pt2[0],pt2[1],pt2[2]),
        ##             (pt1[0],pt2[1],pt2[2]),
        ##             (pt1[0],pt1[1],pt2[2]),
        ##             (pt2[0],pt1[1],pt2[2]))
        ##     coords = grid3D.crystal.toCartesian(ptList)
        ##     box.Set(vertices=coords)
        ## else:
        ##     coords = (pt1, pt2)
        ##     box.Set(cornerPoints=coords) 
                
        from DejaVu2.Box import NiceBox
        self.boxGeom = box = NiceBox('gridOutline')
        box.addToViewer(self.app().gui().viewer)
        box.setCenter( *center)
        sx, sy, sz = sides
        box.setSides( sx, sy, sz)
        self.setGridVisible(True)
        print "Box, sides:", sides, "center" , center
        grid3D.geomContainer['Box'] = box

    def setGridVisible(self, value):
        # value is 0 for unchecked and 2 for checked for checkbox
        # not(value==0) make it work for 0, 1, 2, False, True
        self.boxGeom.master.Set(visible = not(value==0))
        for c in self.boxGeom.master.children:
            if c.name=='faces':
                c.Set(visible = 0)
            else:
                c.Set(visible = not(value==0))


class Isocontour:
    """Isocontour class calculates and displays isocontours for any given gridname
    """
    def __init__(self, app):
        self.grid = None
        self.iso_data = None
        self.app = app

        
    def calculate(self, grid3D, material = None, isovalue = None , 
                  name = None, invertNormals = False):
        """
        \nRequired Arguments:\n    
        grid3D  : key for self.vf.grids3D object\n
        \nOptional Arguments:\n  
        isovalue : if None given, uses the first element in the Grid3D \n
        material : defaults to (0.,0.,1.0,0.5) - yellow half transparent\n
        name     : the name given to IndexedPolygons that represents isocontour.\n   
        invertNormals : defaults to False """
        if isinstance(grid3D, str):
            assert self.app.grids3D.has_key(grid3D)
            grid3D = self.app.grids3D[grid3D]
        if isovalue == None:
            isovalue = float(grid3D.data[0][0][0])
        from UTpackages.UTisocontour import isocontour
        isocontour.setVerboseLevel(0)
        isoc = isocontour.getContour3d(self.iso_data, 0, 0, isovalue,
                                       isocontour.NO_COLOR_VARIABLE)
        
        vert = numpy.zeros((isoc.nvert,3)).astype('f')
        norm = numpy.zeros((isoc.nvert,3)).astype('f')
        col = numpy.zeros((isoc.nvert)).astype('f')
        tri = numpy.zeros((isoc.ntri,3)).astype('i')

        if invertNormals:
            isocontour.getContour3dData(isoc, vert, norm, col, tri, 1)
        else:
            isocontour.getContour3dData(isoc, vert, norm, col, tri, 0)

        if grid3D.crystal:
            vert = grid3D.crystal.toCartesian(vert)
        if not name:
            name = "Grid3D_Iso_%4.4f"%isovalue
        if name in  grid3D.geomContainer['IsoSurf']:
            g = grid3D.geomContainer['IsoSurf'][name]
        else:
            from DejaVu.IndexedPolygons import IndexedPolygons
            g = IndexedPolygons(name)
            if self.app.userpref['Sharp Color Boundaries for MSMS']['value'] == 'blur':
                g.Set(inheritSharpColorBoundaries=False, sharpColorBoundaries=False,)
            g.Set(culling='none')
            g.Set(vertices=vert,vnormals=norm,faces=tri)
            self.app.gui().viewer.AddObject(g, parent = grid3D.IsoSurf)
            grid3D.geomContainer['IsoSurf'][name] = g

        if material:
            g.inheritMaterial = False
            g.Set(materials=[material,],)

        g.Set(vertices=vert,vnormals=norm,faces=tri)
        if vert is not None:
            g.sortPoly()
        self.app.gui().viewer.Redraw()
        return g
    
    def select(self, grid_name):
        grid = self.vf.grids3D[grid_name]
        if not hasattr(grid,'hist'):
            bound = None #this is used to set the limits on the histogram
            if grid_name.endswith('.map'):
                if grid.mini < 0:
                    bound = [grid.mini,-grid.mini]
            hist = numpy.histogram(grid.data.copy().flat,bins=entry['widget'].width+100)
            grid.hist = hist
        if not hasattr(grid,'isoBarNumber'):
            grid.isoLastX = {}
            grid.isoLastColor = {}
            grid.isoBarNumber = 0
            grid.isoBarTags = []
            
        self.grid = grid
        origin = numpy.array(grid.origin).astype('f')
        stepsize = numpy.array(grid.stepSize).astype('f')
        data = grid.data
        if data.dtype != numpy.float32:
            print 'converting %s from %s to float'%(grid_name,data.dtype)
            data = data.astype('f')
        self.newgrid3D = numpy.ascontiguousarray(numpy.reshape( numpy.transpose(data),
                                              (1, 1)+tuple(data.shape) ) , data.dtype.char)
        if self.iso_data:
            isocontour.delDatasetReg(self.iso_data)
        self.iso_data = isocontour.newDatasetRegFloat3D(self.newgrid3D, origin, 
                                                        stepsize)

        
commandClassFromName = {
    'readGrid' : [ReadGrid,  None]}

def initModule(viewer, gui=True):
    for cmdName, values in commandClassFromName.items():
        cmdClass, guiInstance = values
        viewer.addCommand(cmdClass(), cmdName, guiInstance)
