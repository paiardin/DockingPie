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
# Author: Anna Omelchenko, Michel Sanner
#
# Copyright: M. Sanner TSRI 2003
#
#############################################################################

#
# $Header $
#
# $Id $
#

import string, sys,time
from opengltk.OpenGL import GL
import numpy
from struct import unpack
from Volume.Renderers.UTVolumeLibrary import UTVolumeLibrary
#from VolumeLibrary import VolumeLibrary
#print "UTVolumeLibrary imported from" , UTVolumeLibrary.__file__
from DejaVu.viewerFns import checkKeywords
from Volume.Grid3D import Grid3DUC

from DejaVu.VolumeGeom import VolumeGeom, CropBox

class VolumeCrop:
    # trying to mimic the VolumePro API
    kDisable = 0
    kSubVolume = 1
    k3DCross = 2
    k3DCrossInvert = 3
    k3DCrossInvert = 4
    k3DFence = 5
    k3DFenceInvert = 6
	
    def __init__(self, volgeom):
        self.volgeom = volgeom
	self.data = None
	self.flag = self.kSubVolume
	self.xmin = 0
	self.xmax = 0
	self.ymin = 0
	self.ymax = 0
	self.zmin = 0
	self.zmax = 0
	self.size = (0,0,0)

    def updateData(self):
        self.data = numpy.array(self.volgeom.dataArr)
	nx,ny,nz = self.volgeom.volumeSize
	self.xmin = 0
	self.xmax = nx
	self.ymin = 0
	self.ymax = ny
	self.zmin = 0
	self.zmax = nz
	self.size = (nx,ny,nz)

    def SetXSlab(self, xmin, xmax):
        #if self.data == None: return
	self.xmin = xmin
	self.xmax = xmax+1
	if self.flag == self.kDisable: return
	self.data[self.zmin:self.zmax, self.ymin:self.ymax, xmin:xmax+1] = self.volgeom.dataArr[self.zmin:self.zmax, self.ymin:self.ymax, xmin:xmax+1]
	self.data[:, :, :xmin]=0
	self.data[:, :, xmax+1: ]=0
	nx,ny,nz=self.size
	self.volgeom.volume.uploadColorMappedData(self.data.ravel(), nx,ny,nz)

    def SetYSlab(self, ymin, ymax):
	#if self.data == None: return
	self.ymax = ymax+1
	self.ymin = ymin
	if self.flag == self.kDisable: return
	self.data[self.zmin:self.zmax, ymin:ymax+1, self.xmin:self.xmax] = self.volgeom.dataArr[self.zmin:self.zmax, ymin:ymax+1, self.xmin:self.xmax]
	self.data[:, :ymin, :]=0
	self.data[:, ymax+1:, :]=0
	nx,ny,nz=self.size
	self.volgeom.volume.uploadColorMappedData(self.data.ravel(),
							  nx, ny, nz)

    def SetZSlab(self, zmin, zmax):
	#if self.data == None: return
	self.zmax = zmax+1
	self.zmin = zmin
	if self.flag == self.kDisable: return
	self.data[zmin:zmax+1, self.ymin:self.ymax, self.xmin:self.xmax] = \
		  self.volgeom.dataArr[zmin:zmax+1, self.ymin:self.ymax,
				       self.xmin:self.xmax]
	self.data[:zmin, :, :]=0
	self.data[zmax+1:, :, :]=0
	nx,ny,nz=self.size
	self.volgeom.volume.uploadColorMappedData(self.data.ravel(),
						  nx, ny, nz)

    def SetSlabs(self, xmin, xmax, ymin, ymax, zmin, zmax):
	#if self.data == None: return
	self.data[zmin:zmax+1, ymin:ymax+1, xmin:xmax+1] = \
		  self.volgeom.dataArr[zmin:zmax+1, ymin:ymax+1, xmin:xmax+1]
	nx,ny,nz=self.size
	if ymin > 0:
	    self.data[: , : ymin, :] = 0
	if ymax < ny:
	    self.data[:, ymax+1: , :] = 0
	if xmin > 0:
	    self.data[:xmin, :, :] = 0
	if xmax < nx:
	    self.data[xmax+1:, :, :] = 0
	if zmin > 0:
	    self.data[:, :, :zmin] = 0
        if zmax < nz:
	    self.data[:, :, zmax+1:] = 0
	self.volgeom.volume.uploadColorMappedData(self.data.ravel(),
						  nx, ny, nz)
	self.xmin = xmin
	self.xmax = xmax
	self.ymax = ymax
	self.ymin = ymin
	self.zmax = zmax
	self.zmin = zmin

    def SetFlags(self, flag):
        self.flag = flag
        nx,ny,nz=self.size
	if flag == self.kDisable:
	    #nx,ny,nz=self.size
	    self.xmin = 0
	    self.xmax = nx
	    self.ymin = 0
	    self.ymax = ny
	    self.zmin = 0
	    self.zmax = nz
	self.data[self.zmin:self.zmax+1,
		  self.ymin:self.ymax+1,
		  self.xmin:self.xmax+1] = \
		  self.volgeom.dataArr[self.zmin:self.zmax+1,
				       self.ymin:self.ymax+1,
				       self.xmin:self.xmax+1]
	self.volgeom.volume.uploadColorMappedData(self.data.ravel(),
						  nx, ny,nz)
		

class UTVolRenGeom(VolumeGeom):
	
    def Draw(self):
        """ Display the render volume"""
        if self.dataArr is None:
	    return
        t1 = time.time()
	UTVolumeLibrary.InitTexParameteri() 
	self.volume.renderVolume()
	t2 = time.time()
	#print "time tor render:",t2 -t1
		
    def __init__(self, name=None, **kw):
        
        apply( VolumeGeom.__init__, (self,name,0), kw)
	self.inheritXform=1
        self.Set(frontPolyMode='fill', backPolyMode='fill')
	
	#volume renderer
	self.volume = UTVolumeLibrary.VolumeRenderer()
	self.dataArr =None
	self.byte_map=None
	
	
##  	vi.cameras[0].Activate()
##  	self.InitVolumeRenderer()
	self.volumeSize = (0,0,0)
	self.onaddVolume_list = []
	self.volrenInitialized = 0
	self.updateBoundBox = 1
	self.volume.GetSize = self.GetSize
	self.crop = VolumeCrop(self)
	self.cropBox = CropBox(self)
	self.cropBox.crop = self.crop
	self.masterObject = None
	self.cropOpts = ['off','SubVolume']
	self.firstLoaded = 1

        # rendering quality used when moving
        self.coarseRenderingQuality = 0.5


    def coarseRendering(self, event=None):
        self.volume.setQuality(self.coarseRenderingQuality)


    def fineRendering(self, event=None):
        self.volume.setQuality(1.0)

		
    def InitVolumeRenderer(self):
        """ Initialize the volume rendering"""
	if not self.volume.initRenderer():
	    print "Warning, there was an error initializing the volume renderer\n"
        else:
	    self.volrenInitialized = 1

		    
    def LoadVolume(self,file):
        """ Get the data from binary file,rawiv format"""
	
	#load colormap from file
	myfile = open(file, "rb" )
	#the header
	header=myfile.read(68)
	# unpack header following rawiv format,big endian
	h = unpack('>6f5I6f',header)
	width=int(h[8])
	height=int(h[9])
	depth =int(h[10])
	
	# load the data
	l = myfile.read(width*height*depth)
	self.dataArr = numpy.fromstring(l, numpy.uint8,
					  (width*height*depth))
	self.volumeSize = (width,height,depth)
	#self.dataArr = numpy.reshape(self.dataArr,(width,height,depth))
	myfile.close()
			
	if self.dataArr == None:
	    print " you need to load a volume data"
            return
        if self.firstLoaded:
	    if self.viewer:
	        rootObject = self.viewer.rootObject
		self.minBB=(-0.5, -0.5, -0.5)
		self.maxBB=(0.5, 0.5, 0.5)
		self.viewer.NormalizeCurrentObject()
		self.firstLoaded = 0
	status = self.volume.uploadColorMappedData(self.dataArr,width,height,depth)
	if status !=1:
	    raise RuntimeError("uploadColorMappedData() in LoadVoulume failed. Status %d"% status)
            #print "status uploadColorMappedData: ", status
	self.dataArr = numpy.reshape(self.dataArr,
				       (depth,height,width))
	if self.byte_map == None:
	    self.grayRamp()
	#update crop box
	self.crop.updateData()
	nx,ny,nz = self.volumeSize
	self.cropBox.setVolSize((nx,ny,nz))
	self.cropBox.xmin = 0
	self.cropBox.xmax = nx
	self.cropBox.ymin = 0
	self.cropBox.ymax = ny
	self.cropBox.zmin = 0
	self.cropBox.zmax = nz
	self.cropBox.update()
	for c in self.onaddVolume_list:
	    c.OnAddVolumeToViewer()


    def AddGrid3D(self, grid):
        assert isinstance(grid, Grid3DUC)

	if not self.volrenInitialized:
	    if not self.viewer: # we need an OpenGL context before
                print "self.viewer: ", self.viewer
	        return      # we can initialize the renderer
	    self.InitVolumeRenderer()
            for c in self.viewer.cameras:
                c.addButtonDownCB(self.coarseRendering)
                c.addButtonUpCB(self.fineRendering)

	orig = grid.origin[:]
	step = grid.stepSize

	nx, ny, nz = grid.data.shape
	gridSize = [nx*step[0], ny*step[1], nz*step[2]]
	gridSize2 = [gridSize[0]*0.5, gridSize[1]*0.5, gridSize[2]*0.5]
	maxgrid = [ orig[0]+gridSize[0],
		    orig[1]+gridSize[1],
		    orig[2]+gridSize[2]]
	transl = [ orig[0]+gridSize2[0],
                   orig[1]+gridSize2[1],
                   orig[2]+gridSize2[2]]

        # save scale and tranlation into Matrix which is not affected
	# by reset
        if grid.crystal:
            # compute cartesian voxel sizes
            x, y, z = grid.getStepSizeReal()

            #build crystal object for length with padding
            from mglutil.math.crystal import Crystal
            if hasattr(grid, 'dataDims'):
                # compute ratios of padding along the 3 dimensions
                dx = grid.dimensions[0] - grid.dataDims[0] - 1
                dy = grid.dimensions[1] - grid.dataDims[1] - 1
                dz = grid.dimensions[2] - grid.dataDims[2] - 1
                ry = float(dx)/dy
                rz = float(dx)/dz
            else:
                ry = rz = 1.0

            dims = (grid.dimensions[0]-1, (grid.dimensions[1]-1)*ry,
                    (grid.dimensions[2]-1)*rz ) 
            cryst = Crystal( (dims[0]*x, dims[1]*y, dims[2]*z),
                             grid.crystal.angles)

            matrix = numpy.identity(4, 'f')
            matrix[:3, :3] = cryst.ftoc.astype('f')
            matrix.shape = (16,)

            # cleanup=False because rotation can contain shear which should
            # be kept
            self.MatrixRot, MatrixTransl, self.MatrixScale = self.Decompose4x4(
                matrix, cleanup=False)

            # set utvolgeom's Matrix components
            self.MatrixScale = (self.MatrixScale[0], self.MatrixScale[1]/ry,
                                self.MatrixScale[2]/rz)

            origin = grid.crystal.toCartesian(grid.origin)
            self.MatrixTransl = (origin[0]+MatrixTransl[0]+dims[0]*0.5*x,
                                 origin[1]+MatrixTransl[1]+dims[1]*0.5*y/ry,
                                 origin[2]+MatrixTransl[2]+dims[2]*0.5*z/rz)
            #o = grid.origin
            #t1 = cryst.toFractional(MatrixTransl)
            #t2 = (0.5, 0.5, 0.5)
            #trans = ( o[0]+t1[0]+t2[0], o[1]+t1[1]+t2[1],o[2]+t1[2]+t2[2])
            #self.MatrixTransl = cryst.toCartesian(trans)
            #print o
            #print t1
            #print t2
            #print trans
            #print self.MatrixTransl

            #self.MatrixTransl = (0.,0.,0.)
            #self.MatrixScale = (1.,1.,1.)
            
            RotInv = numpy.transpose(numpy.reshape(self.MatrixRot, (4,4)))
            self.MatrixRotInv = numpy.reshape(RotInv, (16,))
            # 
            self.MatrixRot = self.MatrixRot.astype('f')
            self.MatrixRotInv = self.MatrixRot.astype('f')

        else:
            self.setMatrixComponents(trans=transl)
            self.setMatrixComponents(scale=gridSize, trans=transl)

	if self.firstLoaded:
	    if self.viewer:
	        rootObject = self.viewer.rootObject
		self.minBB = [-0.5, -0.5, -0.5]
		self.maxBB = [ 0.5,  0.5,  0.5]
                #print "all objects:", self.viewer.rootObject.AllObjects()
		self.viewer.NormalizeCurrentObject()
		self.firstLoaded = 0

	# scale and translate volume
	##self.SetScale( gridSize)
	##trans = numpy.array( gridSize2, 'f')
	##self.SetTranslation( trans )
        
	#mat = self.GetMatrix(self)
	#self.SetMatrix(mat)
	#self.ResetTransformation()
        
	arr = numpy.ascontiguousarray( numpy.transpose(grid.data),
			     grid.data.dtype.char)
	upload = self.volume.uploadColorMappedData
	status = upload(arr.ravel(), nx, ny,nz)
	if status !=1:
	    raise RuntimeError("uploadColorMappedData() in AddVolume failed. Status %d"% status)

        self.dataArr = numpy.reshape(arr,(nz,ny,nx))
	self.volumeSize = (nx,ny,nz)
	if self.byte_map == None:
	    self.grayRamp()
	# update cropping box
	self.crop.updateData()
	self.cropBox.setVolSize((nx,ny,nz))
	self.cropBox.xmin = 0
	self.cropBox.xmax = nx
	self.cropBox.ymin = 0
	self.cropBox.ymax = ny
	self.cropBox.zmin = 0
	self.cropBox.zmax = nz
	self.cropBox.update()
	for c in self.onaddVolume_list:
	    c.OnAddVolumeToViewer()

    def AddVolume(self, arr, nx, ny, nz, datatype='numarr'):
	if datatype == 'numarr':
	    assert len(arr.shape) == 1

	if not self.volrenInitialized:
	    if not self.viewer: # we need an OpenGL context before
	        return      # we can initialize the renderer

	    self.InitVolumeRenderer()

	if self.firstLoaded:
	    if self.viewer:
	        rootObject = self.viewer.rootObject
		self.minBB=(-0.5, -0.5, -0.5)
		self.maxBB=(0.5, 0.5, 0.5)
		self.viewer.NormalizeCurrentObject()
		self.firstLoaded = 0
	if datatype == 'numarr':
	    status = self.volume.uploadColorMappedData(arr,nx,ny,nz)
	else:
	    status = self.volume.uploadZeroPaddedData(arr,nx,ny,nz)
	#print "status uploadColorMappedData: ", status
	if status !=1:
	    raise RuntimeError("uploadColorMappedData() in AddVolume failed. Status %d"% status)
	self.dataArr = numpy.reshape(arr,(nz,ny,nx))
	self.volumeSize = (nx,ny,nz)
	if self.byte_map == None:
	    self.grayRamp()
	# update cropping box
	self.crop.updateData()
	self.cropBox.setVolSize((nx,ny,nz))
	self.cropBox.xmin = 0
	self.cropBox.xmax = nx
	self.cropBox.ymin = 0
	self.cropBox.ymax = ny
	self.cropBox.zmin = 0
	self.cropBox.zmax = nz
	self.cropBox.update()
	for c in self.onaddVolume_list:
	    c.OnAddVolumeToViewer()
			
    def LoadColorMap(self,file):
        """ Get the colormap data from binaries files"""
		
	myfile = open(file, "rb" )
	l = myfile.read(256*4)
	self.byte_map =numpy.fromstring(l,(256*4),
					  numpy.uint8) 
	self.byte_map = numpy.reshape(self.byte_map,(256,4))
	myfile.close()
	if self.byte_map == None:
	    print " you need to load a color map"
	    return
        else:
	    self.volume.uploadColorMap(self.byte_map)


    def setVolRenColors(self, values):
	""" update the color from the color map data"""
	#the color editor are in 0.0 .... 1.0 range
	#the color map are in 0 .. 255 range
	val = values[1] *255
	self.byte_map[values[0]:values[0]+len(val),:3] = val.astype(numpy.uint8)
	self.volume.uploadColorMap(self.byte_map)
	self.viewer.Redraw()
		 
    def setVolRenAlpha(self, values):
	""" update the opacity of the colormap"""
	#val = (values[1]/4096.)*255
	val = values[1]
	self.byte_map[values[0]:values[0]+len(val),3] = val.astype(numpy.uint8)
	self.volume.uploadColorMap(self.byte_map)
	self.viewer.Redraw()

    def grayRamp(self):
        self.byte_map = numpy.transpose(numpy.array([range(256),range(256),range(256),range(256)] )).astype(numpy.uint8)
	self.volume.uploadColorMap(self.byte_map)
	
    def Set(self, check=1, redo=1, updateOwnGui=True, **kw):
	redoFlags = apply( VolumeGeom.Set, (self,0,0), kw)
	s = kw.get('volscale')
	tr = kw.get('voltranslate')
	if s:
	    self.SetScale(s)
	if tr:
	    tr = numpy.array(tr)
	    self.SetTranslation(tr)

        return self.redoNow(redo, updateOwnGui, redoFlags)

    def GetSize(self):
        return self.volumeSize
		

if __name__ == '__main__':

    # create a DejaVu Viewer
    from DejaVu import Viewer
    vi = Viewer()
    #Create a volumerender geoms
    g = UTVolRenGeom('volren')
    g.viewer = vi
    
    vi.cameras[0].Activate()
    g.InitVolumeRenderer()
    
    g.LoadVolume("vh4.rawiv")
    g.grayRamp()
    #g.LoadColorMap("/mgl/ms3/python/dev/Volume/SimpleExample/colormap.map")
    #g.InitData()
    
    # add a transparent object without a dpylist
    vi.AddObject(g)
    
    #Color map editor
    
##  	from mglutil.gui.BasicWidgwts.Tk.tablemaker import TableManager
##  	t = TableManager(master=vi.master,
##  			 ymaxval = 255,
##  			 alphaCallback = g.setVolRenAlpha,
##  			 colorCallback = g.setVolRenColors)







