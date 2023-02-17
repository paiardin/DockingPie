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

########################################################################
#
# Date: May 2006 Authors: Guillaume Vareille, Michel Sanner
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
# $Header: /mnt/raid/services/cvs/DejaVu2/glfLabels.py,v 1.1.1.1.4.1 2017/07/13 22:28:33 annao Exp $
#
# $Id: glfLabels.py,v 1.1.1.1.4.1 2017/07/13 22:28:33 annao Exp $
#
import os
import numpy
import types
import math
import warnings

from opengltk.OpenGL import GL
from opengltk.extent.utillib import glCleanRotMat
from mglutil.math import rotax
from mglutil.math.VectorModule import Vector

from DejaVu2.IndexedGeom import IndexedGeom
from DejaVu2.viewerFns import checkKeywords
from DejaVu2 import viewerConst
from DejaVu2.colorTool import glMaterialWithCheck, resetMaterialMemory
from DejaVu2.IndexedPolygons import IndexedPolygons

from pyglf import glf

class GlfLabels(IndexedGeom):
    """Class for sets of 3d labels
"""
    keywords = IndexedGeom.keywords + [
        'billboard',
        'font',
        'fontRotateAngles',
        'fontScales',
        'fontSpacing',
        'fontStyle',
        'fontTranslation',
        'labelTranslation',
        'includeCameraRotationInBillboard',
        'labels',
        ]
    
    bitmapFont3dList = [
#        'arbat.bmf',
#        'arial.bmf',
#        'brushtype.bmf',
#        'chicago.bmf',
#        'courier.bmf',
#        'cricket.bmf',
#        'crystal.bmf',
#        'fixedsys.bmf',
#        'gals.bmf',
#        'greek.bmf',
#        'impact.bmf',
#        'proun.bmf',
#        'techno.bmf',
#        'times_new.bmf'
    ]
    
    vectorFont3dList = [
        'arial1.glf',
        'courier1.glf',
        'crystal1.glf',
        'techno0.glf',
        'techno1.glf',
        'times_new1.glf',
        'aksent1.glf',
        'alpine1.glf',
        'broadway1.glf',
        'chicago1.glf',
        'compact1.glf',
        'cricket1.glf',
        'garamond1.glf',
        'gothic1.glf',
        'penta1.glf',
        'present_script1.glf']

    fontList = vectorFont3dList + bitmapFont3dList

    fontStyleDict = {
        'solid':   glf.glfDrawSolidString,
        'solid3d': glf.glfDraw3DSolidString,
        'wire':   glf.glfDrawWiredString,
        'wire3d': glf.glfDraw3DWiredString,
        }
    fontStyleList = fontStyleDict.keys()


    def __init__(self, name=None, check=1, **kw):
        #print "GlfLabels::__init__"

        self.displayedList = None # can be set to a list as long as vertices
        if not kw.get('shape'):
            kw['shape'] = (0,3)

        #if len(kw.get('labels'))==0:
        #    kw['labels'] = ['Aa']
        if kw.get('billboard') is None:
            kw['billboard'] = True
        if kw.get('includeCameraRotationInBillboard') is None:
            kw['includeCameraRotationInBillboard'] = False
        if not kw.get('font'):
            kw['font'] = self.fontList[0]
        if not kw.get('fontStyle'):
            kw['fontStyle'] = 'solid'
        if not kw.get('fontSpacing'):
            kw['fontSpacing'] = .2
        if kw.get('fontScales', None) is None:
            kw['fontScales'] = (1, 1, 1)
        if kw.get('fontRotateAngles', None) is None:
            kw['fontRotateAngles'] = (0, 0, 0)
        if kw.get('fontTranslation', None) is None:
            kw['fontTranslation'] = (0, 0, 0)
        if kw.get('labelTranslation', None) is None:
            kw['labelTranslation'] = None
        
        # Glf initialisations
        glf.glfInit()
        # usefull glf messages (like "font not found")
        glf.glfEnable(glf.GLF_CONSOLE_MESSAGES)
        # we manage this through the opengl Z scale. (we have to make a call to glScale anyway)
        glf.glfSetSymbolDepth(.8) # .2 is C code default value (it is not deep enough)
        #contouring is far from being nice in the glf 1.4 C code! So we don't use it.
        #glf.glfDisable(glf.GLF_CONTOURING) # C code default value
        #glf.glfDisable(glf.GLF_TEXTURING) # C code default value
        # loading BMF fonts has to be done after the context exists
        # otherwise nothing appears
        self.vectorFonts = {}
        self.bitmapFonts = {}
        self.useGlfDrawLabels = True
        apply( IndexedGeom.__init__, (self, name, check), kw)

        # make ambient to be same as diffuse
        mat = self.materials[GL.GL_FRONT]
        mat.getFrom[mat.ambient] = [mat.diffuse, 1.]

        self.labelTranslation = kw['labelTranslation']


    def Set(self, check=1, redo=1, updateOwnGui=True, **kw):
        """set data for this object
check=1 : verify that all the keywords present can be handle by this func 
redo=1 : append self to viewer.objectsNeedingRedo
updateOwnGui=True : allow to update owngui at the end this func
"""
        #print "Set glfLabels"
        redoFlags = apply( IndexedGeom.Set, (self, check, 0), kw )

        displayed = kw.get('displayed')
        if displayed:
            assert len(displayed)==len(self.vertexSet.vertices)
            self.displayedList = numpy.array(displayed, 'b')

        labels = kw.get('labels')
        if labels is not None and len(labels)!=0:
            kw.pop('labels')
            #self.labels = numpy.array(labels, dtype='|S5')
            #import pdb
            #pdb.set_trace()
            self.labels = labels
            redoFlags |= self._redoFlags['redoDisplayListFlag']
            redoFlags |= self._redoFlags['updateOwnGuiFlag']

        billboard = kw.get('billboard')
        if billboard is not None:
            kw.pop('billboard')
            self.billboard = billboard
            self.immediateRendering = billboard
            redoFlags |= self._redoFlags['redoDisplayListFlag']
            redoFlags |= self._redoFlags['updateOwnGuiFlag']

        includeCameraRotationInBillboard = kw.get('includeCameraRotationInBillboard')
        if includeCameraRotationInBillboard is not None:
            kw.pop('includeCameraRotationInBillboard')
            self.includeCameraRotationInBillboard = includeCameraRotationInBillboard
            redoFlags |= self._redoFlags['updateOwnGuiFlag']

        font = kw.get('font')
        if not font is None:
            # loading BMF fonts has to be done after the context exists
            # otherwise nothing appears
            kw.pop('font')
            lGlfModulePath = os.path.split(glf.__file__)[-2]
            lPathToFonts = lGlfModulePath+os.sep+'fonts'+os.sep
            if font in self.vectorFont3dList:
                self.fontTypeIsVector = True
                if font not in self.vectorFonts:
                    self.vectorFonts[font] = glf.glfLoadFont(lPathToFonts+font)
                self.font = font
                redoFlags |= self._redoFlags['redoDisplayListFlag']
                redoFlags |= self._redoFlags['updateOwnGuiFlag']
            elif font in self.bitmapFont3dList:
                self.fontTypeIsVector = False
                if font not in self.bitmapFonts:
                    self.bitmapFonts[font] = glf.glfLoadBMFFont(lPathToFonts+font)
                self.font = font
                redoFlags |= self._redoFlags['redoDisplayListFlag']
                redoFlags |= self._redoFlags['updateOwnGuiFlag']
            else:
                warnings.warn('not a glf font %s'%font)

        fontStyle = kw.get('fontStyle')
        if fontStyle is not None:
            kw.pop('fontStyle')
            assert fontStyle in self.fontStyleList           
            self.fontStyle = fontStyle
            redoFlags |= self._redoFlags['redoDisplayListFlag']
            redoFlags |= self._redoFlags['updateOwnGuiFlag']

        fontSpacing = kw.get('fontSpacing')
        if fontSpacing is not None:
            # this set the space between letters
            #glf.glfSetSymbolSpace(.2) # C code default value
            kw.pop('fontSpacing')
            self.fontSpacing = fontSpacing
            redoFlags |= self._redoFlags['redoDisplayListFlag']
            redoFlags |= self._redoFlags['updateOwnGuiFlag']

        fontScales = kw.get('fontScales')
        if fontScales is not None:
            kw.pop('fontScales')
            self.fontScales = fontScales
            redoFlags |= self._redoFlags['redoDisplayListFlag']

        fontRotateAngles = kw.get('fontRotateAngles')
        if fontRotateAngles is not None:
            kw.pop('fontRotateAngles')
            self.fontRotateAngles = fontRotateAngles
            redoFlags |= self._redoFlags['redoDisplayListFlag']
            redoFlags |= self._redoFlags['updateOwnGuiFlag']

        fontTranslation = kw.get('fontTranslation')
        if fontTranslation is not None:
            kw.pop('fontTranslation')
            self.fontTranslation = fontTranslation
            redoFlags |= self._redoFlags['redoDisplayListFlag']
            redoFlags |= self._redoFlags['updateOwnGuiFlag']

        labelTranslation = kw.get('labelTranslation')
        if labelTranslation is not None:
            kw.pop('labelTranslation')
            assert labelTranslation.shape == self.vertexSet.vertices.array.shape
            self.labelTranslation = labelTranslation
            redoFlags |= self._redoFlags['redoDisplayListFlag']
            redoFlags |= self._redoFlags['updateOwnGuiFlag']
        return self.redoNow(redo, updateOwnGui, redoFlags)


    def getState(self, full=False):
        """return a dictionary describing this object's state
This dictionary can be passed to the Set method to restore the object's state
"""
        state = IndexedGeom.getState(self, full).copy()
        state.update( self.getSubClassState() )
        return state
            

    def getSubClassState(self):
        """return a dictionary describing this object's state
This dictionary can be passed to the Set method to restore the object's state
"""
        state = {
                'billboard': self.billboard,
                'font': self.font,
                'fontRotateAngles': self.fontRotateAngles,
                'fontScales': self.fontScales,
                'fontSpacing': self.fontSpacing,
                'fontStyle': self.fontStyle,
                'fontTranslation': self.fontTranslation,
                'labelTranslation': self.labelTranslation,
                'includeCameraRotationInBillboard': self.includeCameraRotationInBillboard,
                'labels': self.labels,
            }
        return state
            

    def Add(self, check=1, redo=1, **kw):
        """Add glfLabels
"""
        #print "Add glfLabels"

        if __debug__:
            if check:
                apply( checkKeywords, (self.name,self.keywords), kw)

        v = kw.get('vertices')
        if v:
            self.redoDspLst = 1
        
        labels = kw.get( 'labels')
        if labels:
            #labels is apparently a tuple
            self.labels = self.labels + labels
            self.redoDspLst = 1
        
            IndexedGeom.Add( self, check=0, redo=0,
                      vertices = kw.get( 'vertices'),
                      materials = kw.get( 'materials') 
                    )
        
        if self.viewer and redo:
            if self.redoDspLst:
                self.viewer.objectsNeedingRedo[self] = None


    def drawOne3dTextLine(self, textLine, index):

        c = self.vertexSet.vertices.array
        n = self.vertexSet.normals.array
        
        if self.labelTranslation is not None:
            labelTrans = self.labelTranslation[index]
        else:
            labelTrans = (0,0,0)
        lenn = len(n)
        try:
            GL.glPushMatrix()
            GL.glTranslatef(
                            float(c[index][0]), 
                            float(c[index][1]), 
                            float(c[index][2]),
             )

            if self.billboard:
                GL.glMultMatrixf(self.billboardRotation)

            elif lenn == len(c):
                    lMat = rotax.rotVectToVect(n[index] , (0,0,1) )
                    lMat = [
                                lMat[0][0],
                                lMat[0][1],
                                lMat[0][2],
                                lMat[0][3],
                                lMat[1][0],
                                lMat[1][1],
                                lMat[1][2],
                                lMat[1][3],
                                lMat[2][0],
                                lMat[2][1],
                                lMat[2][2],
                                lMat[2][3],
                                lMat[3][0],
                                lMat[3][1],
                                lMat[3][2],
                                lMat[3][3]
                           ]
                    GL.glMultMatrixf(lMat)
            elif lenn > 0:
                    GL.glMultMatrixf(self.orientation)

            GL.glTranslatef(
                            float(labelTrans[0]), 
                            float(labelTrans[1]), 
                            float(labelTrans[2]),
            )

            GL.glTranslatef(
                            float(self.fontTranslation[0]), 
                            float(self.fontTranslation[1]), 
                            float(self.fontTranslation[2]),
            )
    
            GL.glRotatef(float(self.fontRotateAngles[2]),
                         0.,
                         0.,
                         1.,
                         )
    
            GL.glRotatef(float(self.fontRotateAngles[1]),
                         0.,
                         1.,
                         0.,
                         )
    
            GL.glRotatef(float(self.fontRotateAngles[0]),
                         1.,
                         0.,
                         0.,
                         )
    
            if self.fontTypeIsVector:
                GL.glScalef(float(self.fontScales[0]),
                            float(self.fontScales[1]),
                            float(self.fontScales[2]))
                if textLine is not None:
                    self.fontStyleDict[self.fontStyle](textLine)
            else:
                GL.glScalef(float(self.fontScales[0])*40,
                            float(self.fontScales[1])*40,
                            0) # to be the same size as vector fonts
                if textLine is not None:
                    glf.glfDrawBString(textLine)
        finally:
            if textLine is None:
                lMatrix = GL.glGetDoublev(GL.GL_MODELVIEW_MATRIX)
                GL.glPopMatrix()
                return lMatrix
            else:
                GL.glPopMatrix()


    def prepareBillboardAndNormalForAllTextLines(self):
        if self.billboard:
            m = self.GetMatrix()
            m = numpy.reshape(m, (16,))           
            rot = glCleanRotMat(m) #much faster than self.Decompose4x4(m)

            if self.includeCameraRotationInBillboard:
                # this permit billboarding even if the camera is not in the Z axis
                lCameraTransformation = self.viewer.currentCamera.GetMatrix()
                lCameraTransformation = numpy.reshape(lCameraTransformation, (16,))           
                lCameraRotation = glCleanRotMat(lCameraTransformation) #much faster than self.Decompose4x4(m)
                lCameraRotation = numpy.transpose(lCameraRotation)
                rot = numpy.dot(lCameraRotation, rot)

            rot = numpy.reshape(rot, (16,))
            self.billboardRotation = rot.astype('f')
        else:
            c = self.vertexSet.vertices.array
            n = self.vertexSet.normals.array
            lenn = len(n)
            if lenn > 0 and lenn != len(c):
                lMat = rotax.rotVectToVect(n[0] , (0,0,1) )
                self.orientation = [
                            lMat[0][0],
                            lMat[0][1],
                            lMat[0][2],
                            lMat[0][3],
                            lMat[1][0],
                            lMat[1][1],
                            lMat[1][2],
                            lMat[1][3],
                            lMat[2][0],
                            lMat[2][1],
                            lMat[2][2],
                            lMat[2][3],
                            lMat[3][0],
                            lMat[3][1],
                            lMat[3][2],
                            lMat[3][3]
                       ]

    def Draw(self):
        #print "GlfLabels.Draw", self

        if len(self.faceSet.faces)==0:
            centers = self.vertexSet.vertices.array
            labels = self.labels
        else:
            faces = self.faceSet.faces.array
            centers = self.vertexSet.vertices.array[faces.flatten()]
            labels = self.labels[faces.flatten()]

        if len(centers)==0: 
            return
        #import pdb
        #pdb.set_trace()
        if labels is None or len(labels) == 0:
            return

        elif len(labels) == centers.shape[0]:
            txt = None
        else:
            txt = labels[0]
            if type(txt) != types.StringType:
                txt= str(txt)
        labels = numpy.array(self.labels)
        self.prepareBillboardAndNormalForAllTextLines()

        glf.glfSetSymbolSpace(self.fontSpacing)

        if self.fontTypeIsVector:
            font = self.vectorFonts[self.font]
            glf.glfSetCurrentFont(font)
        else:
            font = self.bitmapFonts[self.font]
            glf.glfSetCurrentBMFFont(font)
            glf.glfStartBitmapDrawing()

        resetMaterialMemory()
        if self.inheritMaterial:
            fp = None
            bp = None
            fpProp = None
            frontAndBack = 0
        else:
            self.InitMaterial()
            colorFront = numpy.array(self.materials[GL.GL_FRONT].prop[1], copy=1)
            fp = self.materials[GL.GL_FRONT]
            if fp.binding[1] == viewerConst.OVERALL:
                GL.glColor4fv(colorFront[0])
                bp = None
                face = GL.GL_FRONT_AND_BACK
                frontAndBack = 1
            else:
                bp = self.materials[GL.GL_BACK]
                face = GL.GL_FRONT
                frontAndBack = 0

        ## if fp: # this block is not used
        ##     for m in (0,1,2,3,4):
        ##         if fp.binding[m] == viewerConst.OVERALL:
        ##             glMaterialWithCheck( face,
        ##                                  viewerConst.propConst[m],
        ##                                  fp.prop[m][0])
        ##     if fp.binding[1] == viewerConst.OVERALL:
        ##         GL.glColor4fv(colorFront[0])
        #import pdb; pdb.set_trace()
        #print self.labels
        if self.useGlfDrawLabels: # use C function from pyglf to draw labels
            #print labels
            #print "type:", labels.dtype, "strides:", labels.strides
            if fp:
                #fp = self.materials[GL.GL_FRONT]
                fpProp = []
                for propInd in range(4):
                    b, p = fp.GetProperty(propInd)
                    fpProp.append(p)
                fpProp.append(fp.prop[4]) 
            fontStyles = {'solid' :0, 'solid3d': 1, 'wire':2, 'wire3d':3 }
            if self.displayedList is None:
                indices = numpy.arange(len(centers), dtype='i4')
                #displayedList = numpy.ones(len(centers), 'b')
            else:
                #displayedList = numpy.ones(self.displayedList, 'b')
                indices = numpy.where(self.displayedList==True)[0].astype('i4')
            labTrans = self.labelTranslation
            if labTrans is None:
                labTrans = [[0.,0.,0.]]*len(centers)
            glf.glfDrawLabels(
                centers, labels, indices, self.billboardRotation,
                self.fontTranslation, self.fontRotateAngles,
                self.fontTypeIsVector, self.fontScales,
                fontStyles[self.fontStyle], fpProp, frontAndBack,
                labTrans)
        else:
            for i in xrange(centers.shape[0]):
                if self.displayedList is not None and not self.displayedList[i]:
                    continue
                if fp:
                    for m in (0,1,2,3,4):
                        if fp.binding[m] == viewerConst.PER_VERTEX:
                            b, p = fp.GetProperty(m)
                            glMaterialWithCheck( face,
                                                 viewerConst.propConst[m],
                                                 p[i] )
                    if fp.binding[1] != viewerConst.OVERALL:
                        GL.glColor4fv(colorFront[i])
                if bp:
                    for m in (0,1,2,3,4):
                        if bp.binding[m] == viewerConst.PER_VERTEX:
                            b, p = bp.GetProperty(m)
                            glMaterialWithCheck( face,
                                                 viewerConst.propConst[m],
                                                 p[i] )

                #GL.glPushName(i)
                if txt is None:
                    txt2 = self.labels[i]
                    if type(txt2) != types.StringType:
                        txt2 = str(txt2)
                    self.drawOne3dTextLine(txt2, i)
                else:       
                    self.drawOne3dTextLine(txt, i)
                #GL.glPopName()

        if self.fontTypeIsVector is False:
            glf.glfStopBitmapDrawing()

        return 1


    def asIndexedPolygons(self, run=1, removeDupVerts=0, **kw):
        """Should return an IndexedPolygons object if this object can be
represented using triangles, else return None. run=0 returns 1
if this geom can be represented as an IndexedPolygon and None if not
run=1 returns the IndexedPolygon object.
"""
        #print "GlfLabels::asIndexedPolygons"
        #import pdb;pdb.set_trace()

        if run == 0:
            return 1 # yes, I can be represented as IndexedPolygons

        self.prepareBillboardAndNormalForAllTextLines()

        lLenLabelsVertices = len(self.vertexSet.vertices.array)
        if lLenLabelsVertices != len(self.labels):
            lSameText = True
            assert len(self.labels) > 0
            if type(self.labels[0]) == types.StringType:
                output = glf.glfGet3DSolidStringTriangles(self.labels[0])
            else:
                output = glf.glfGet3DSolidStringTriangles( str(self.labels[0]) )
            lNumOfTriangleVertices = len(output[0])
            assert lNumOfTriangleVertices == len(output[1])
            #assert is a multiple of 3
        else:
            lSameText = False

        if lLenLabelsVertices != len(self.vertexSet.normals.array):
            if len(self.vertexSet.normals.array) > 0:
                lSameOrientation = True
                lLocalNormal = self.normals[0]
            else:
                lLocalNormal = None

        self.prepareBillboardAndNormalForAllTextLines()

        lOverallTriangles = []
        lOverallVertices = []
        lOverallTrianglesColors = []
        for i in range(lLenLabelsVertices):
            #print "i", i
            if lSameText is False:
                if type(self.labels[0]) == types.StringType:
                    output = glf.glfGet3DSolidStringTriangles(self.labels[i])
                else:
                    output = glf.glfGet3DSolidStringTriangles( str(self.labels[i]) )
                lNumOfTriangleVertices = len(output[0])
                assert lNumOfTriangleVertices == len(output[1])
                #assert is a multiple of 3

            colorFront = numpy.array(self.materials[GL.GL_FRONT].prop[1], copy=1)

            lNumOfOverallVertices = len(lOverallVertices)
            for j in range(
                       lNumOfOverallVertices, 
                       lNumOfOverallVertices+lNumOfTriangleVertices,
                       3):
                lOverallTriangles.append( ( j, j+1, j+2 ) )
                lOverallTrianglesColors.append(colorFront[i])

            lMatrix = self.drawOne3dTextLine(None, i)
            lMatrix.shape = (4, 4)
            for j in range(lNumOfTriangleVertices):
                lVertexBeforeTransform = output[0][j]
                lVertex = []
                lVertex.append( lMatrix[0][0] * lVertexBeforeTransform[0] \
                              + lMatrix[1][0] * lVertexBeforeTransform[1] \
                              + lMatrix[2][0] * lVertexBeforeTransform[2] \
                              + lMatrix[3][0] )
                lVertex.append( lMatrix[0][1] * lVertexBeforeTransform[0] \
                              + lMatrix[1][1] * lVertexBeforeTransform[1] \
                              + lMatrix[2][1] * lVertexBeforeTransform[2] \
                              + lMatrix[3][1] )
                lVertex.append( lMatrix[0][2] * lVertexBeforeTransform[0] \
                              + lMatrix[1][2] * lVertexBeforeTransform[1] \
                              + lMatrix[2][2] * lVertexBeforeTransform[2] \
                              + lMatrix[3][2] )
                lOverallVertices.append(lVertex)
        
        lIndexedPolygons = IndexedPolygons(
                               self.name+'_glfTriangles', 
                               vertices=lOverallVertices, 
                               faces=lOverallTriangles,
                               materials=lOverallTrianglesColors,
                               visible=1,
                               invertNormals=self.invertNormals,
                               )

        return lIndexedPolygons


    def createOwnGui(self):
        #print "GlfLabels.createOwnGui", self
        pass

    def setWithOwnGui(self, event=None):
        """
"""
        #print "setWithOwnGui", event
        lGlobalScale = self.ownGui.guiGlobalScale.get()
        self.Set(labels=eval(self.ownGui.labelsEnt.get()),
                 billboard=self.ownGui.billboardVar.get(),
                 includeCameraRotationInBillboard=self.ownGui.includeCameraRotationInBillboardVar.get(),
                 #lighting=self.ownGui.lightingVar.get(),
                 font=self.ownGui.guiFontComboBox.get(), 
                 fontStyle=self.ownGui.guiFontStyleComboBox.get(),
                 fontSpacing=self.ownGui.guiSpacing.get(),
                 fontScales = ( lGlobalScale*self.ownGui.guiScaleX.get(),
                                lGlobalScale*self.ownGui.guiScaleY.get(),
                                lGlobalScale*self.ownGui.guiScaleZ.get() ),
                 fontTranslation = ( self.ownGui.guiTranslateX.get(),
                                     self.ownGui.guiTranslateY.get(),
                                     self.ownGui.guiTranslateZ.get() ),
                 fontRotateAngles = ( self.ownGui.guiRotateX.get(),
                                      self.ownGui.guiRotateY.get(),
                                      self.ownGui.guiRotateZ.get() ),
                 updateOwnGui=False)


    def updateOwnGui(self):
        self.ownGui.title(self.name)
        self.ownGui.labelsEnt.setvalue(str(self.labels))
        self.ownGui.billboardVar.set(self.billboard)
        self.ownGui.includeCameraRotationInBillboardVar.set(self.includeCameraRotationInBillboard)
        #self.ownGui.lightingVar.set(self.lighting)
        self.ownGui.guiFontComboBox.selectitem(self.font)
        self.ownGui.guiFontStyleComboBox.selectitem(self.fontStyle)
        self.ownGui.guiSpacing.set(self.fontSpacing, update=0)
        self.ownGui.guiGlobalScale.set(1., update=0)
        self.ownGui.guiScaleX.set(self.fontScales[0], update=0)
        self.ownGui.guiScaleY.set(self.fontScales[1], update=0)
        self.ownGui.guiScaleZ.set(self.fontScales[2], update=0)
        self.ownGui.guiTranslateX.set(self.fontTranslation[0], update=0)
        self.ownGui.guiTranslateY.set(self.fontTranslation[1], update=0)
        self.ownGui.guiTranslateZ.set(self.fontTranslation[2], update=0)
        self.ownGui.guiRotateX.set(self.fontRotateAngles[0], update=0)
        self.ownGui.guiRotateY.set(self.fontRotateAngles[1], update=0)
        self.ownGui.guiRotateZ.set(self.fontRotateAngles[2], update=0)



