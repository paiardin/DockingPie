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

# -*- coding: utf-8 -*-
"""
Created on Mon Sep 19 13:31:01 2011

@author: -
"""
import DejaVu2 as DejaVu
from DejaVu2.Spheres import GLUSpheres

try :
    from DejaVu2 import hyperballs
    hyperballsFound = True
except :
    hyperballsFound = False
    
if hyperballsFound:
    from DejaVu2.hyperballs.AtomAndBondGLSL import AtomAndBondGLSL 
    from DejaVu2.hyperballs.ballimproved_frag import ballimproved_frag
    from DejaVu2.hyperballs.ballimproved_vert import ballimproved_vert
    from DejaVu2.hyperballs.stickimproved_frag import stickimproved_frag 
    from DejaVu2.hyperballs.stickimproved_vert import stickimproved_vert 
#    from opengltk.OpenGL.GL import *
#    from opengltk.OpenGL.GLU import *
#    from OpenGL.GLUT import glutSwapBuffers
#    from OpenGL import GL 
#    from OpenGL.GL import *
    
    import numpy
    TOO_BIG = 0
    TOO_SMALL = 1
    SIZE_OK = 2
    
    class hyperBalls(GLUSpheres):

        keywords = list(GLUSpheres.keywords)
        for kw in ['quality', 'stacks', 'slices']:
            keywords.remove(kw)


        def makeTemplate(self):
            if __debug__:
             if hasattr(DejaVu, 'functionName'): DejaVu.functionName()
            pass


        def __init__(self, name=None, check=1, **kw):
            if __debug__:
             if hasattr(DejaVu, 'functionName'): DejaVu.functionName()
            self.r = 1 # default scaling factor for system
            self.pos = numpy.zeros(3) # center of bounding box
#            self.orien = quaternion([0,0,-1,0]) # orientation in space
            self.scaleFactor = 1
            self.idx = None
            self.DirV = []
            self.clipplane = numpy.array([0.,0.,0.,0,], numpy.float32)
            self.excl = numpy.array([], numpy.int32)
            self.sticks = True
            self.shaders = AtomAndBondGLSL()
            #print "ok shaders ", self.shaders
            self.shaders.setAtomFragmentShaderProgramSource(ballimproved_frag);
            self.shaders.setAtomVertexShaderProgramSource(ballimproved_vert);
            self.shaders.setBondFragmentShaderProgramSource(stickimproved_frag);
            self.shaders.setBondVertexShaderProgramSource(stickimproved_vert);
            self.immediateRendering = False
            self.width = kw.get( 'w')
            self.bufferAssigned = False
            self.atompos = None
            self.numbonds = 0
            self.numatoms = 0
            self.bonds = [] 
            self.shrink = 0.01
            self.bScale = 0.1
            self.aScale = 1.0
            self.draw_axes = False
            self.pmvProj = False
            self.pmvTransf = False
            self.lightDir = (0.,0.8,0.6)
            apply( GLUSpheres.__init__, (self, name, check), kw)
            self.dpyList=None
            # shrink, aScale, bScale 
            self.options={"CPK":[0.01,0.2,0.3],
                          "VDW":[0.01,3,0.01],
                          "LIC":[0.01,0.26,0.26],
                          "HBL":[0.3,0.4,0.4]}

        def setOption(self,name):
            self.shrink, self.aScale, self.bScale = self.options[name]
            self.initialise()
            #update ?

        def setAtoms(self,mol):
            self.pmvmol = mol
            self.atoms = mol.allAtoms
            for at in mol.allAtoms:
                at.colors["qutemol"] = tuple( [1.0,.0,.0] )
 
#        def initialise(self,):
#            #se should trigger the wheel event zomm to change he overall cssscle
#            #recompute r
#            self.shaders.initGL()    

        def initialise(self):
            #if not already initialized do it, otherwise update
            self.shaders.bshrinks = self.shrink
            self.shaders.aScales = self.aScale#/10.0
            self.shaders.bScales = self.bScale#/10.0
            if not self.bufferAssigned  :
                self.shaders.setupBuffersAndTextures(self.numatoms,self.numbonds,
                                                     self.atompos,self.bonds,
                                                     self.colors,self.radii);
                self.shaders.initGL() 
                self.bufferAssigned = True
            else :
                self.shaders.updateBuffersAndTextures(self.numatoms,self.numbonds,
                                                      self.atompos,self.bonds,
                                                      self.colors,self.radii);
            
        def Set(self, check=1, redo=0, updateOwnGui=False, **kw):
            """set data for this object check=1 : verify that all the keywords present can be handle by this func redo=1 : append self to viewer.objectsNeedingRedo updateOwnGui=True : allow to update owngui at the end this func """            
            if __debug__:
                if hasattr(DejaVu, 'functionName'): DejaVu.functionName()

            v = kw.get( 'centers')
            if v is None :
                v = kw.get('vertices')
            mat = kw.has_key('materials')
            rad = kw.get( 'radii')
            colors = kw.get( 'colors') 
            bonds = kw.get( "bonds" )
            shrink = kw.get('shrink')
            bScale = kw.get('bScale')
            aScale = kw.get('aScale')
            redoFlags = 0
            ve=False
            bo=False
            ra=False
            if type(v) != type(None)and len(v):
                ve=True
#                self.redoDspLst=1
                self.atompos = coords = numpy.array(v, numpy.float32)
                self.numatoms = len(coords)
                self.radii = numpy.ones(self.numatoms)
                self.shaders.nbAtoms = self.numatoms;
                self.shaders.positions = self.atompos
                #update vbo ?
            if  type(bonds) != type(None)and len(bonds):
                bo=True
                self.bonds = numpy.array(kw["bonds"], int)#.tolist()
                self.numbonds = len(self.bonds)
                self.shaders.nbBonds = self.numbonds;
                self.shaders.links=numpy.array(self.bonds,numpy.uint32)
            if ve or bo :
                self.shaders.updateTextureSize()
                self.shaders.initVertTCoordIndice();#in case nbatoms changed ?
                if not self.bufferAssigned  :
                    self.shaders.initGL() 
                    self.bufferAssigned = True
                self.shaders.updateVBO()
                if ve : 
                    self.shaders.updatePositions(self.atompos)                
            if type(rad) != type(None)and len(rad):
                ra=True
                self.radii = self.vdw_radii = numpy.array(rad, numpy.float32)       
            if ve or ra :
                self.shaders.updateSizes(self.radii)
            if shrink:
                self.shrink = shrink
                #update shrink
                self.shaders.updateShrinks(shrink)
            if bScale:
                self.bScale = bScale
                self.shaders.updateBondScales(bScale)
                #update bScale
            if aScale:
                self.aScale = aScale
                self.shaders.updateAtomScales(aScale)           
            if colors:# colors :
                self.colors = numpy.array(colors, numpy.float32)
            elif mat :
                if len(kw["materials"]) == 1 : #uniq color
                    self.colors = numpy.array(kw["materials"]*self.numatoms,numpy.float32)
                else :
                    self.colors = numpy.array(kw["materials"], numpy.float32)
            if colors or mat :
                self.shaders.updateColors(self.colors) 
#            if self.atompos is not None or 'colors' in kw:
#                self.initialise()
            redoFlags = apply( GLUSpheres.Set, (self, 0, 0), kw)
            redoFlags = 0
            #print "okSET"
#            self.SetSpaceFill()
            return self.redoNow(redo, updateOwnGui, redoFlags)
    
        def updateColors(self,):
            self.colors = numpy.array(self.pmvmol.allAtoms.colors['hpballs'], numpy.float32)

        def doClip(self):
            if self.viewer is not None :
                if self.viewer.currentClip.visible :
                    #recompute the plane ?
                    self.clipplane = self.viewer.currentClip.eqn
                else :
                    self.clipplane = numpy.array([0.,0.,0.,0,], numpy.float32)
                clip = self.viewer.currentClip.clipPlaneNames[self.viewer.currentClip.num]
            else :
                clip = GL_CLIP_PLANE0
            return clip

        def makeTemplate(self):
            pass

        def DisplayFunction(self):
#            print "ok HP DisplayFunction"
            self.Draw()

        def vertexArrayCallback(self):
#            print "ok HP vertexArrayCallback"
            self.Draw()
        
        def Draw(self):
            if __debug__:
                if hasattr(DejaVu, 'functionName'): DejaVu.functionName()
            """Draw function of the geom
    return status 0 or 1
    If you want fast rendering, you need to set self.templateDSPL
    using MakeTemplate.
    """     
            #print "ok HP Draw"
            self.shaders.useAllTextures()
            status = 1
            self.shaders.drawAtoms_vbo()            
            if self.sticks : 
                self.shaders.drawBonds_vbo();
            return status


if __name__=="__main__":
    import numpy
    from DejaVu2.hpBalls import hyperBalls
    hsp = hyperBalls("hpBalls")
    print "ok",hsp
    v=numpy.random.random((100,3))*20.0
    r=numpy.random.random(100)*5.0
    c=numpy.random.random((100,4))
    hsp.Set(vertices=v,radii=r.tolist(),materials=c)#no bons
    hsp.sticks=False
    vi = self.GUI.VIEWER
    vi.useMasterDpyList=0
    vi.AddObject(hsp)
    
