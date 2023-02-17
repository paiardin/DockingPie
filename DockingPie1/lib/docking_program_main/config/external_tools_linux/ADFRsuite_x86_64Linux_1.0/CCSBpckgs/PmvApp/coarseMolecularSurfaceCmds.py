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
# Author: Michel F. SANNER, Ruth Huey
#
# Copyright: M. Sanner TSRI 2015
#
#############################################################################

#
# $Header: /mnt/raid/services/cvs/PmvApp/coarseMolecularSurfaceCmds.py,v 1.1.4.1 2017/07/13 20:52:52 annao Exp $
#
# $Id: coarseMolecularSurfaceCmds.py,v 1.1.4.1 2017/07/13 20:52:52 annao Exp $
#

from PmvApp.Pmv import MVCommand
from PmvApp.Pmv import DeleteAtomsEvent, AddAtomsEvent, EditAtomsEvent
from DejaVu2.IndexedPolygons import IndexedPolygons
from MolKit2.selection import SelectionSet

import numpy

class ComputeCoarseMolecularSurfaceCommand(MVCommand): #, MVAtomICOM):
    """Compute 
    \nPackage : PmvApp
    \nModule  : coarseMolecularSurfaceCmds.py
    \nClass   : ComputeCoarseMolecularSurfaceCommand
    None <--- computeCoarseMolecularSurface(self, molSels, 
                       surfName='coarseMolSurf', perMol=True,
                       gridSize=32, padding=0., resolution=-0.3,
                       isolvalue='fast approximation'.... )

    \nRequired Arguments:\n
    molSels --- any set of MolKit2.Selection describing molecular components\n

    \nOptional Arguments:\n
    lineWidth --- int specifying the width of the lines, dots or doted lines
                 representing the selection. (default = 2)\n
     """

    def handleDeleteAtoms(self, event=None):
        print 'computeCoarseMolecularSurfaceCommand does not yet implement a handler for delete atoms events'

    def handleAddAtoms(self, event=None):
        print 'computeCoarseMolecularSurfaceCommand does not yet implement a handler for add atoms events'

    def handleEditAtoms(self, event=None):
        print 'computeCoarseMolecularSurfaceCommand does not yet implement a handler for edit atoms events'

    def onAddCmdToApp(self):
        self.UTisocontourData = None
        self.app().eventHandler.registerListener(
            DeleteAtomsEvent, self.handleDeleteAtoms)
        self.app().eventHandler.registerListener(
            AddAtomsEvent, self.handleAddAtoms)
        self.app().eventHandler.registerListener(
            EditAtomsEvent, self.handleEditAtoms)

    def initializeMolForCmd(self, mol):
        """
        Creates the coarseMolSurf geometry
        """
        self.initializedFor[mol] = True
        mol._coarseMolSurfParams = {}
            
    def checkArguments(self, molSelSet, surfName='coarseMolSurf', perMol=True,
                       gridSize=32, padding=0., resolution=-0.3, bind_surface_to_molecule=True,
                       isovalue='fast approximation', check_surface_components=False ):
        """
        \nRequired Arguments:\n
        molSel --- MolKit2 selection\n

        \nOptional Arguments:\n
        lineWidth --- int specifying the width of the lines, dots or doted lines
                     representing the selection. (default = 2)\n
        displayBO --- boolean flag specifying whether or not to display the
                     bond order (default = False)\n
        negate --- boolean flag specifying whether or not to negate the
                     current selection. (default = False)\n
        """
        kw = {}
        assert isinstance(molSelSet, SelectionSet)
        assert isinstance(surfName, str)
        assert perMol in [True, False, 0, 1]
        assert isinstance(gridSize, int) and gridSize > 0
        assert isinstance(padding, float)
        assert isinstance(resolution, float) and resolution < 0.0
        assert isinstance(isovalue, float) or isovalue in [
            'fast approximation','precise value']
        kw['surfName'] = surfName
        kw['perMol'] = perMol
        kw['gridSize'] = gridSize
        kw['padding'] = padding
        kw['resolution'] =  resolution
        kw['isovalue'] =  isovalue
        return (molSelSet,), kw

    def UTblur(self, molSel, Xdim, Ydim, Zdim, padding, blobbyness):
        from UTpackages.UTblur import blur
        from Volume.Grid3D import Grid3DF
        coords = molSel.getCoords()
        radii = molSel.getRadii()
        volarr, origin, span = blur.generateBlurmap(
            coords.astype('f'), radii.astype('f'), [Xdim, Ydim, Zdim],
            blobbyness, padding = padding)
        volarr.shape = (Zdim, Ydim, Xdim)
        volarr = numpy.ascontiguousarray(numpy.transpose(volarr), 'f')
        self.UTblurData = (volarr, origin, span)
        #print "volarr", volarr
        h = {}
        grid3D= Grid3DF( volarr, origin, span , h)
        #print "data:", grid3D.data
        h['amin'], h['amax'],h['amean'],h['arms']= grid3D.stats()
        return grid3D

    def UTisocontour(self, grid3D, isovalue, verbosity=None):
        from UTpackages.UTisocontour import isocontour
        if verbosity is not None:
            isocontour.setVerboseLevel(verbosity)

        data = grid3D.data
        #print 'isovalue', isovalue
        origin = numpy.array(grid3D.origin).astype('f')
        stepsize = numpy.array(grid3D.stepSize).astype('f')
        # add 1 dimension for time steps amd 1 for multiple variables
        if data.dtype != numpy.float32:
            #print 'converting from ', data.dtype.char
            data = data.astype('f')

        newgrid3D = numpy.ascontiguousarray(numpy.reshape( numpy.transpose(data),
                                               (1, 1)+tuple(data.shape) ),
                                data.dtype.char)
        # destroy the ConDataset structure
        if self.UTisocontourData:
            isocontour.delDatasetReg(self.UTisocontourData)

        self.UTisocontourData = gridData = isocontour.newDatasetRegFloat3D(\
                newgrid3D, origin, stepsize)
        #print "UTisocontourData", self.UTisocontourData
        sig = [gridData.getSignature(0, 0, 0),
               gridData.getSignature(0, 0, 1),
               gridData.getSignature(0, 0, 2),
               gridData.getSignature(0, 0, 3)]
        if self.UTisocontourData:
            isoc = isocontour.getContour3d(self.UTisocontourData, 0, 0, isovalue,
                                           isocontour.NO_COLOR_VARIABLE)

            vert = numpy.zeros((isoc.nvert,3)).astype('f')
            norm = numpy.zeros((isoc.nvert,3)).astype('f')
            col = numpy.zeros((isoc.nvert)).astype('f')
            tri = numpy.zeros((isoc.ntri,3)).astype('i')
            isocontour.getContour3dData(isoc, vert, norm, col, tri, 0)

            if grid3D.crystal:
                vert = grid3D.crystal.toCartesian(vert)
            self.isoCoords = vert
            self.isoIndices = tri
            self.isoNormals = norm
            return [vert, tri, norm]
        else:
            print "Warning: UTisocontourData is None, using saved coords, indices and normals"
            return [self.isoCoords, self.isoIndices, self.isoNormals]

    def piecewiseLinearInterpOnIsovalue(self, x):
        """Piecewise linear interpretation on isovalue that is a function
        blobbyness.
        """
        import sys
        X = [-3.0, -2.5, -2.0, -1.5, -1.3, -1.1, -0.9, -0.7, -0.5, -0.3, -0.1]
        Y = [0.6565, 0.8000, 1.0018, 1.3345, 1.5703, 1.8554, 2.2705, 2.9382, 4.1485, 7.1852, 26.5335]
        if x<X[0] or x>X[-1]:
            print "WARNING: Fast approximation :blobbyness is out of range [-3.0, -0.1]"
            return None
        i = 0
        while x > X[i]:
            i +=1
        x1 = X[i-1]
        x2 = X[i]
        dx = x2-x1
        y1 = Y[i-1]
        y2 = Y[i]
        dy = y2-y1
        return y1 + ((x-x1)/dx)*dy

    def getIsoValue(self, resolution, isovalue):
        if isinstance(isovalue, float):
            return isovalue
        elif isovalue=='fast approximation':
            return self.piecewiseLinearInterpOnIsovalue(resolution)
        elif isovalue=='precise value':
            raise ValueError("precise value not yet implemented")
        else:
            raise ValueError("ERROR: isovalue can only be a float or 'fast approximation' or 'precise value'. Got %s"%str(isovalue))

    def doit(self, molSel, surfName='coarseMolSurf', perMol=True,
             gridSize=32, padding=0., resolution=-0.3, bind_surface_to_molecule=True,
             isovalue='fast approximation', check_surface_components=False):

        mol = molSel.getAtomGroup().getMolecule()
        if not self.initializedFor.get(mol, False):
            self.initializeMolForCmd(mol)
            
        isovalue = self.getIsoValue(resolution, isovalue)
        mol._coarseMolSurfParams[surfName] = {'perMol':perMol, 'gridSize':gridSize,
                                              'padding':padding, 'resolution':resolution,
                                              'isovalue':isovalue}
        gc = mol.geomContainer
        if not gc.geoms.has_key(surfName):
            geom  = IndexedPolygons(surfName, visible=False)
            gc.addGeom(geom, parent=gc.masterGeom, redo=0 )
        else:
            geom = gc.geoms[surfName]
        grid3D = self.UTblur(molSel, gridSize, gridSize, gridSize,
                             padding, resolution)
        
        coords, indices, normals = self.UTisocontour(grid3D, isovalue)
        if coords is None:
            geom.Set(faces=[], visible=0)
            gc.setAtomsForGeom(surfName, mol.emptySelection())
        else:
            geom.Set(vertices=coords, faces=indices, vnormals=normals,
                     tagModified=False, visible=1, inheritMaterial=None)
            
            gc.setAtomsForGeom(surfName, molSel)
            if bind_surface_to_molecule:
                atoms = self.app().bindGeomToMolecularFragment(geom, molSel)
                if atoms is not None and len(atoms):
                    gc.geomPickToAtoms[geom.name] = self.pickedVerticesToAtoms
                    gc.geomPickToBonds[geom.name] = None
            #if checkComp:
            #    #print "checking for connected components"
            #    newindices, newfaces = self.checkConnectedComponents(surf)
            #    surf.Set(vertices=newfaces, faces=newindices)
    
    def pickedVerticesToAtoms(self, geom, vertInd):
        """Function called to convert picked vertices into atoms"""
        
        # this function gets called when a picking or drag select event has
        # happened. It gets called with a geometry and the list of vertex
        # indices of that geometry that have been selected.
        mol = geom.mol()
        geomC = mol.geomContainer
        atom_inds = geomC.boundGeom[geom.name]['cl_atoms']
        atoms  = geomC.boundGeom[geom.name]['atoms'].getIndices()
        #pickedAtms = []
        #for ind in vertInd:
        #    pickedAtms.append(atoms[atom_inds[ind]])
        #return pickedAtms
        return [atoms[atom_inds[i]] for i in vertInd]

    
commandClassFromName = {
    'computeCoarseMolecularSurface' : [ComputeCoarseMolecularSurfaceCommand, None],
}

def initModule(viewer, gui=True):
    for cmdName, values in commandClassFromName.items():
        cmdClass, guiInstance = values
        viewer.addCommand(cmdClass(), cmdName, guiInstance)

