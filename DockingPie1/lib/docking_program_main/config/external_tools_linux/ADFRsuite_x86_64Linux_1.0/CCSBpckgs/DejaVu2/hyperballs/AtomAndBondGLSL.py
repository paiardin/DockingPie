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

##// $Id: AtomAndBondGLSL.py,v 1.6.4.1 2017/07/13 22:27:37 annao Exp $
#// ****************************************************************** #//
#//                                                                    #//
#// Copyright (C) 2010-2011 by                                         #//
#// Laboratoire de Biochimie Theorique (CNRS),                         #//
#// Laboratoire d'Informatique Fondamentale d'Orleans (Universite      #//
#// d'Orleans),                                                        #//
#// (INRIA) and                                                        #//
#// Departement des Sciences de la Simulation et de l'Information      #//
#// (CEA).                                                             #//
#// ALL RIGHTS RESERVED.                                               #//
#//                                                                    #//
#// contributors :                                                     #//
#// Matthieu Chavent,                                                  #//
#// Antoine Vanel,                                                     #//
#// Alex Tek,                                                          #//
#// Marc Piuzzi,                                                       #//
#// Jean-Denis Lesage,                                                 #//
#// Bruno Levy,                                                        #//
#// Sophie Robert,                                                     #//
#// Sebastien Limet,                                                   #//
#// Bruno Raffin and                                                   #//
#// Marc Baaden                                                        #//
#//                                                                    #//
#// October 2011                                                       #//
#//                                                                    #//
#// Contact: Marc Baaden                                               #//
#// E-mail: baaden@smplinux.de                                         #//
#// Webpage: http:#//hyperballs.sourceforge.net                         #//
#//                                                                    #//
#// This software is a computer program whose purpose is to visualize  #//
#// molecular structures. The source code is part of FlowVRNano, a     #//
#// general purpose library and toolbox for interactive simulations.   #//
#//                                                                    #//
#// This software is governed by the CeCILL-C license under French law #//
#// and abiding by the rules of distribution of free software. You can #//
#// use, modify and/or redistribute the software under the terms of    #//
#// the CeCILL-C license as circulated by CEA, CNRS and INRIA at the   #//
#// following URL "http:#//www.cecill.info".                            #//
#//                                                                    #//
#// As a counterpart to the access to the source code and  rights to   #//
#// copy, modify and redistribute granted by the license, users are    #//
#// provided only with a limited warranty and the software's author,   #//
#// the holder of the economic rights, and the successive licensors    #//
#// have only limited liability.                                       #// 
#//                                                                    #//
#// In this respect, the user's attention is drawn to the risks        #//
#// associated with loading, using, modifying and/or developing or     #//
#// reproducing the software by the user in light of its specific      #//
#// status of free software, that may mean  that it is complicated to  #//
#// manipulate, and that also therefore means that it is reserved for  #//
#// developers and experienced professionals having in-depth computer  #//
#// knowledge. Users are therefore encouraged to load and test the     #//
#// software's suitability as regards their requirements in conditions #//
#// enabling the security of their systems and/or data to be ensured   #//
#// and, more generally, to use and operate it in the same conditions  #//
#// as regards security.                                               #//
#//                                                                    #//
#// The fact that you are presently reading this means that you have   #//
#// had knowledge of the CeCILL-C license and that you accept its      #//
#// terms.                                                             #//
#// ****************************************************************** #//
from math import *
from .AtomAndBondShader import *
import numpy as np
from math import *
useopengltk=False
try :
    import opengltk
    useopengltk=True
except:
    useopengltk=False
if useopengltk :
    from opengltk.OpenGL import GL
    from opengltk.OpenGL import GLU
    from opengltk.extent import glextlib
    from opengltk.extent import _gllib
    GL.glBindBuffer=glextlib.glBindBufferARB#GL.glBindBuffer
    GL.glBufferData=glextlib.glBufferDataARB
    GL.GL_STREAM_DRAW= glextlib.GL_STREAM_DRAW_ARB
    GL.GL_ELEMENT_ARRAY_BUFFER = glextlib.GL_ELEMENT_ARRAY_BUFFER# 0x8893
    GL.GL_ARRAY_BUFFER = glextlib.GL_ARRAY_BUFFER_ARB#0x8892
    GL.glGetUniformLocation = glextlib.glGetUniformLocation
    GL.glGetAttribLocation = glextlib.glGetAttribLocation
    GL.glUseProgram = glextlib.glUseProgram
    GL.glUniform1i = glextlib.glUniform1i
    GL.glDisableVertexAttribArray = glextlib.glDisableVertexAttribArray
    GL.glEnableVertexAttribArray = glextlib.glEnableVertexAttribArray
    GL.glVertexAttribPointer = glextlib.glVertexAttribPointer
    for k in glextlib.__dict__:
        if not hasattr(GL,k):
            setattr(GL,k,glextlib.__dict__[k])
else :
    from OpenGL import GL
    from OpenGL import GLU
#from OpenGL import GL as oGL
#from OpenGL.GLUT import *

import sys

sizeof = sys.getsizeof

GL_TEXTURE_RECTANGLE_NV = 0x84F5
GL_RGB32F_ARB = 0x8815

TEX_POSITIONS=0
TEX_COLORS=1
TEX_SIZES=2
TEX_ASCALES=3
TEX_BSCALES=4
TEX_SHRINKS=5

DEBUG=False

#pragma mark -
#pragma mark Constructor & Destructor
#//-----------------------------------------------------------------------------
#//
#// Constructor
#//
class AtomAndBondGLSL( AtomAndBondShaders ) :
    #pragma mark -
    #pragma mark OpenGL Drawing Methods
    def __init__(self):
        AtomAndBondShaders.__init__(self)
        self.vertex_shader=0
        self.fragment_shader=0
        self.program=0
        self.Atp=0
        self.Atc=0
        self.Atsz=0
        self.Atsc=0
        self.Bvertex_shader=0
        self.Bfragment_shader=0
        self.Bprogram=0
        self.Btp=0
        self.Btc=0
        self.Btsz=0
        self.Btsc=0
        self.Btsh=0
        self.Dvertex_shader=0
        self.Dprogram=0
        self.Ivertex_shader=0
        self.Ifragment_shader=0
        self.Iprogram=0
        self.framebufferID = np.ones(2,int)#[2];                // FBO names
        self.renderbufferID=0                  #// renderbuffer object name
        self.maxRenderbufferSize=0;              #// maximum allowed size for FBO renderbuffer
        self.texCoordOffsets = np.ones(18,float);
        self.windowWidth=0
        self.windowHeight=0
        self.textureWidth=0
        self.textureHeight=0
        
        self.maxDrawBuffers=0;                   #// maximum number of drawbuffers supported
        self.renderTextureID = np.ones(2,int)#[2];              // 1 in 1st pass and 1 in 2nd pass
        
        self.pos_loc=0

    def _drawSimple(self,):
#        oGL.glFrontFace(oGL.GL_CCW)
        print "test draw simpl vbo"
        GL.glUseProgram(self.program);
        print "bind VBO"
        try :
            self.vertice_vbo.bind()
            self.indice_vbo.bind()
            try :
                GL.glEnableVertexAttribArray(self.pos_loc)
                GL.glVertexAttribPointer(self.pos_loc, 3, GL.GL_FLOAT, False, 0, self.vertice_vbo)
#                self.texturecoord_vbo.bind()
#                oGL.glTexCoordPointer(2,GL.GL_FLOAT,0,self.texturecoord_vbo);
                GL.glDrawElements(GL.GL_TRIANGLES,len(self.indice_vbo.data.flatten()), GL.GL_UNSIGNED_INT, self.indice_vbo );#24 integer for indices
            finally:
                GL.glDisableVertexAttribArray(self.pos_loc)    
                self.indice_vbo.unbind()
                self.vertice_vbo.unbind()
        finally:
            GL.glUseProgram(0) 
        
    def drawAtoms_vbo(self,):
        GL.glUseProgram(self.program);
        #print "Draw Atoms vbo"
        try :
            self.useAllTextures()
#            self.vertice_vbo.bind()
#            self.indice_vbo.bind()
#            stride = self.vertice_vbo.data[0].nbytes
            self.bindVBO(self.vertice_vbo,self.ball_data)            
            self.bindVBO(self.indice_vbo,self.indice,target=GL.GL_ELEMENT_ARRAY_BUFFER)
            stride = self.ball_data[0].nbytes
            #            self.texturecoord_vbo.bind()                
            try :
                #//Bindings between textures handles and textures ids
#                self.useTexture(GL.GL_TEXTURE0,TEX_POSITIONS)
                GL.glUniform1i(self.Atp,TEX_POSITIONS);
                GL.glUniform1i(self.Atc,TEX_COLORS);
                GL.glUniform1i(self.Atsz,TEX_SIZES);
                GL.glUniform1i(self.Atsc,TEX_ASCALES);
                GL.glUniform1i(self.ATextureSize,int(self.textureSize))
                
                GL.glEnableVertexAttribArray(self.pos_loc)
                GL.glEnableVertexAttribArray(self.texco_loc)
#                GL.glEnableVertexAttribArray(self.off_loc)#the atom position?
#                oGL.glEnableVertexAttribArray(self.colors_loc)#the atom position?
#                oGL.glVertexAttribPointer(self.pos_loc, 3, GL.GL_FLOAT, False, stride, self.vertice_vbo)#self.vertice_vbo
#                oGL.glVertexAttribPointer(self.texco_loc, 2, GL.GL_FLOAT, False, stride, self.vertice_vbo+(3*4))
#                oGL.glVertexAttribPointer(self.off_loc, 3, GL.GL_FLOAT, False, stride, self.vertice_vbo+(5*4))
                GL.glVertexAttribPointer(self.pos_loc, 3, GL.GL_FLOAT, False, stride, 0)#self.vertice_vbo
                GL.glVertexAttribPointer(self.texco_loc, 2, GL.GL_FLOAT, False, stride, 3*4)
#                GL.glVertexAttribPointer(self.off_loc, 3, GL.GL_FLOAT, False, stride, 5*4)
#                oGL.glVertexAttribPointer(self.colors_loc, 3, GL.GL_FLOAT, False, stride, self.vertice_vbo+(8*4))
#                oGL.glTexCoordPointer(self.texco_loc,GL.GL_FLOAT,0,None);
                _gllib.glDrawElements( GL.GL_TRIANGLES,len(self.indice.flatten()),
                                          GL.GL_UNSIGNED_INT, None)
#                oGL.glDrawElements(GL.GL_TRIANGLES,
#                                   len(self.indice_vbo.data.flatten()), 
#                                   GL.GL_UNSIGNED_INT, self.indice_vbo );#self.indice_vbo#24 integer for indices
#                oGL.glDrawElementsInstanced(GL.GL_TRIANGLES, len(self.indice_vbo.data.flatten()),
#                                            GL.GL_UNSIGNED_INT, self.indice_vbo,self.nbAtoms)
            finally:
                GL.glDisableVertexAttribArray( self.pos_loc )  
                GL.glDisableVertexAttribArray( self.texco_loc )                
#                GL.glDisableVertexAttribArray( self.off_loc )
#                oGL.glDisableVertexAttribArray( self.colors_loc )
#                self.indice_vbo.unbind()
#                self.vertice_vbo.unbind()
                GL.glBindBuffer(GL.GL_ARRAY_BUFFER, 0)
                GL.glBindBuffer(GL.GL_ELEMENT_ARRAY_BUFFER, 0)
#                self.texturecoord_vbo.unbind()
        finally:
            GL.glUseProgram(0) 
        if DEBUG:
            for v in self.vertice:
                GL.glTranslatef(float(v[0]),float(v[1]),float(v[2]));    
                glutSolidSphere(0.1,20,20)
                GL.glTranslatef(float(-v[0]),float(-v[1]),float(-v[2]));    
            
    #//-----------------------------------------------------------------------------
    # should we use pyrex like pyqutemol ?
#    def drawAtoms(self):
#    	
#    	GL.glEnableClientState(GL.GL_VERTEX_ARRAY);
#    	
#    	#//Buffer binding  
#    	oGL.glBindBuffer(GL.GL_ARRAY_BUFFER,int(self.buffer_id[BUF_VERTICE]));
#    	oGL.glVertexPointer(3,GL.GL_FLOAT,0,0);
#    	
#    	oGL.glBindBuffer(GL.GL_ELEMENT_ARRAY_BUFFER,int(self.buffer_id[BUF_INDICE]));
#    	
#    	#//Indices in the texture of atoms
#    	oGL.glClientActiveTexture(GL.GL_TEXTURE0);
#    	oGL.glEnableClientState(GL.GL_TEXTURE_COORD_ARRAY);
#    	oGL.glBindBuffer(GL.GL_ARRAY_BUFFER,int(self.buffer_id[BUF_TCOORD]));
#    	oGL.glTexCoordPointer(2,GL.GL_FLOAT,0,0);
#    	
#    	#//Load Shader 
#    	GL.glUseProgram(self.program);
#    	
#    	print "#//Handles for textures"
#    	self.Atp = GL.glGetUniformLocation(self.program,"texturePosition");
#    	self.Atc = GL.glGetUniformLocation(self.program,"textureColors");
#    	self.Atsz =GL.glGetUniformLocation(self.program,"textureSizes");
#    	self.Atsc =GL.glGetUniformLocation(self.program,"textureScale");
#    	
#      	#//Bindings between textures handles and textures ids
#    	GL.glUniform1i(self.Atp,TEX_POSITIONS);
#    	GL.glUniform1i(self.Atc,TEX_COLORS);
#    	GL.glUniform1i(self.Atsz,TEX_SIZES);
#    	GL.glUniform1i(self.Atsc,TEX_ASCALES);
#
#    
#    	#//MB test :: Bindings between textures handles and textures ids
#    #//	glUniform1i(Atp,0);
#    #//	glUniform1i(Atc,1);
#    #//	glUniform1i(Atsz,2);
#    #//	glUniform1i(Atsc,4);
#    	
#    	#//Draws the atoms
#    	oGL.glDrawElements(GL.GL_TRIANGLES,len(self.indice),#self.nbAtoms*12*3,#
#    				   GL.GL_UNSIGNED_INT,0);
#    	oGL.glDisableClientState(GL.GL_VERTEX_ARRAY);
#    	#//glDisableClientState(GL_TEXTURE_COORD_ARRAY);
#    	
#    	#//Detach shader programs
#    #//MB OpenGL error
#    #//	if (vertex_shader) glDetachShader(program, vertex_shader);
#    #//	if (fragment_shader) glDetachShader(program, fragment_shader);
#    #//	glDeleteShader(vertex_shader);  
#    #//	glDeleteShader(fragment_shader);
#    	GL.glUseProgram(0);	


    def bindBuffer(self,gl_id,buffer_id):
    	print "bindBuffer"    
    	GL.glClientActiveTexture(gl_id);
    	GL.glEnableClientState(GL.GL_TEXTURE_COORD_ARRAY);
    	GL.glBindBuffer(GL.GL_ARRAY_BUFFER,buffer_id);
    	GL.glTexCoordPointer(2,GL.GL_FLOAT,0,0);


    def drawBonds_vbo(self,):
        GL.glUseProgram(self.Bprogram);
        #bind the vbo
        #print "Draw Bonds vbo"
        try :
            self.useAllTextures()
#            self.stick_vertice_vbo.bind()
#            self.stick_indice_vbo.bind()
#            stride = self.stick_vertice_vbo.data[0].nbytes
            self.bindVBO(self.stick_vertice_vbo,self.stick_data)            
            self.bindVBO(self.stick_indice_vbo,self.indice_sticks,target=GL.GL_ELEMENT_ARRAY_BUFFER)
            stride = self.stick_data[0].nbytes
#            self.texturecoord_vbo.bind()                
            try :
                #//Bindings between textures handles and textures ids
                #print "location set"
                #/#/Bindings between textures handles and textures ids
                GL.glUniform1i(self.Btsh,TEX_SHRINKS); #shrink
                GL.glUniform1i(self.Btp,TEX_POSITIONS);  #position
                GL.glUniform1i(self.Btc,TEX_COLORS);  #color
                GL.glUniform1i(self.Btsz,TEX_SIZES); #size
                GL.glUniform1i(self.Btsc,TEX_BSCALES); #scale
                GL.glUniform1i(self.BTextureSize,int(self.textureSize))
                GL.glUniform1i(self.BTextureBondsSize,int(self.textureSizeBonds))
                
                GL.glEnableVertexAttribArray(self.stick_pos_loc)
                GL.glEnableVertexAttribArray(self.texco1_loc)
                GL.glEnableVertexAttribArray(self.texco2_loc)
                GL.glEnableVertexAttribArray(self.texco3_loc)
#                GL.glEnableVertexAttribArray(self.pos_loc1)
#                GL.glEnableVertexAttribArray(self.pos_loc2)
#                oGL.glEnableVertexAttribArray(self.off_loc)
#                oGL.glVertexAttribPointer(self.stick_pos_loc, 3, GL.GL_FLOAT, False, stride, self.stick_vertice_vbo)#self.vertice_vbo
#                oGL.glVertexAttribPointer(self.texco1_loc, 2, GL.GL_FLOAT, False, stride, self.stick_vertice_vbo+(3*4))
#                oGL.glVertexAttribPointer(self.texco2_loc, 2, GL.GL_FLOAT, False, stride, self.stick_vertice_vbo+(5*4))
#                oGL.glVertexAttribPointer(self.texco3_loc, 2, GL.GL_FLOAT, False, stride, self.stick_vertice_vbo+(7*4))
#                oGL.glVertexAttribPointer(self.pos_loc1, 3, GL.GL_FLOAT, False, stride, self.stick_vertice_vbo+(9*4))#self.vertice_vbo
#                oGL.glVertexAttribPointer(self.pos_loc2, 3, GL.GL_FLOAT, False, stride, self.stick_vertice_vbo+(12*4))#self.vertice_vbo
                GL.glVertexAttribPointer(self.stick_pos_loc, 3, GL.GL_FLOAT, False, stride, 0)#self.vertice_vbo
                GL.glVertexAttribPointer(self.texco1_loc, 2, GL.GL_FLOAT, False, stride, (3*4))
                GL.glVertexAttribPointer(self.texco2_loc, 2, GL.GL_FLOAT, False, stride, (5*4))
                GL.glVertexAttribPointer(self.texco3_loc, 2, GL.GL_FLOAT, False, stride, (7*4))
#                GL.glVertexAttribPointer(self.pos_loc1, 3, GL.GL_FLOAT, False, stride, (9*4))#self.vertice_vbo
#                GL.glVertexAttribPointer(self.pos_loc2, 3, GL.GL_FLOAT, False, stride, (12*4))#self.vertice_vbo
#                oGL.glVertexAttribPointer(self.off_loc, 3, GL.GL_FLOAT, False, stride, self.vertice_vbo+(9*4))
#                oGL.glTexCoordPointer(self.texco_loc,GL.GL_FLOAT,0,None);
#                oGL.glDrawElements(GL.GL_TRIANGLES,len(self.stick_indice_vbo.data.flatten()), 
#                                   GL.GL_UNSIGNED_INT, self.stick_indice_vbo );#self.indice_vbo#24 integer for indices
                _gllib.glDrawElements( GL.GL_TRIANGLES,len(self.indice_sticks.flatten()),
                                          GL.GL_UNSIGNED_INT, None)
#                oGL.glDrawElementsInstanced(GL.GL_TRIANGLES, len(self.indice_vbo.data.flatten()),
#                                            GL.GL_UNSIGNED_INT, self.indice_vbo,self.nbAtoms)
            finally:
                GL.glDisableVertexAttribArray(self.stick_pos_loc)  
                GL.glDisableVertexAttribArray( self.texco1_loc )                
                GL.glDisableVertexAttribArray( self.texco2_loc )                
                GL.glDisableVertexAttribArray( self.texco3_loc )                
#                GL.glDisableVertexAttribArray(self.pos_loc1)  
#                GL.glDisableVertexAttribArray(self.pos_loc2)  
#                oGL.glDisableVertexAttribArray( self.off_loc )
#                self.stick_indice_vbo.unbind()
#                self.stick_vertice_vbo.unbind()
                GL.glBindBuffer(GL.GL_ARRAY_BUFFER, 0)
                GL.glBindBuffer(GL.GL_ELEMENT_ARRAY_BUFFER, 0)
#                self.texturecoord_vbo.unbind()
        finally:
            GL.glUseProgram(0) 
#        if DEBUG:
#            for v in self.vertice:
#                GL.glTranslatef(float(v[0]),float(v[1]),float(v[2]));    
#                glutSolidSphere(0.1,20,20)
#                GL.glTranslatef(float(-v[0]),float(-v[1]),float(-v[2]));    
            
    #//-----------------------------------------------------------------------------
#    def drawBonds(self,):
#    	#//Activates use of Vertex and Texture_Coord
#    	oGL.glEnableClientState(GL.GL_VERTEX_ARRAY);
#    	#//glEnableClientState(GL_TEXTURE_COORD_ARRAY);
#    	
#      	#//Buffer binding for vertices
#    	oGL.glBindBuffer(GL.GL_ARRAY_BUFFER,self.buffer_id[BUF_VERTICE]);
#    	oGL.glVertexPointer(3,GL.GL_FLOAT,0,0);
#    	oGL.glBindBuffer(GL.GL_ELEMENT_ARRAY_BUFFER,self.buffer_id[BUF_INDICE]);
#    	
#    	#//Indices in the texture of atom 1 for each bond
#    	self.bindBuffer(GL.GL_TEXTURE0,self.buffer_id[BUF_TCOORD0])
#    	
#    	#//Indices in the texture of atom 2 for each bond
#    	self.bindBuffer(GL.GL_TEXTURE1,self.buffer_id[BUF_TCOORD1])
#    	
#    	#//Indices in the texture of bonds
#    	self.bindBuffer(GL.GL_TEXTURE2,self.buffer_id[BUF_TCOORD2])
#    	
#    	#//Load shader program
#    	GL.glUseProgram(self.Bprogram);
#    	print "okprogram"
#
#    	print "location set"
#    	#//Bindings between textures handles and textures ids
#    	GL.glUniform1i(self.Btsh,5);
#    	GL.glUniform1i(self.Btp,0);
#    	GL.glUniform1i(self.Btc,1);
#    	GL.glUniform1i(self.Btsz,2);
#    	GL.glUniform1i(self.Btsc,4);
#    	
#    	print "#//Drawing of bonds"
#    	oGL.glDrawElements(GL.GL_TRIANGLES,self.nbBonds*12*3,
#    				   GL.GL_UNSIGNED_INT,None);
#    	oGL.glDisableClientState(GL.GL_VERTEX_ARRAY);
#    	#//glDisableClientState(GL_TEXTURE_COORD_ARRAY);
#    	
#    	#//Detach shader programs
#    	#//MB OpenGL error
#    #//	if (Bvertex_shader) glDetachShader(Bprogram, Bvertex_shader);
#    #//	if (Bfragment_shader) glDetachShader(Bprogram, Bfragment_shader);
#    #//	glDeleteShader(Bvertex_shader);  
#    #//	glDeleteShader(Bfragment_shader);
#    	GL.glUseProgram(0);

    
    #pragma mark -
    #pragma mark Load Shader Methods
    
    #//-----------------------------------------------------------------------------
    #//
    #// Load shader source code
    #//
    def LoadSource(self,fn):
    	#read the shader
        if fn is not None :
            sfile = open(fn,"r")
            lines = sfile.readlines()
            sfile.close()
            code=""
            for l in lines :
                code+=l
            return code

#    
    #//-----------------------------------------------------------------------------
    #//
    #// Compile the shader source code given as an argument
    #//
    def LoadShaderFromString(self,shadertype, src):
        logsize = 0;
        compile_status = GL.GL_TRUE;
        log = None;
        lStatus = 0x7FFFFFFF
        shader = glextlib.glCreateShader(shadertype);
        glextlib.glShaderSource(shader, 1, src, lStatus)
        glextlib.glCompileShader(shader)
        lStatus1 = glextlib.glGetShaderiv(shader, glextlib.GL_COMPILE_STATUS,
                                           lStatus )
        if lStatus1 == 0:
            charsWritten  = 0
            shaderInfoLog = '\0' * 2048
            charsWritten, infoLog = glextlib.glGetShaderInfoLog(shader, 
                        len(shaderInfoLog), charsWritten, shaderInfoLog)
            print "shaderInfoLog", shaderInfoLog
            print "shader didn't compile"
            return 0
        return shader;	
    
    #//-----------------------------------------------------------------------------
    #//
    #// Read and compile the Shader source code.
    #//
    def LoadShaderFromFile(self,shadertype, filename):
    	src = self.LoadSource(filename);
    	r = self.LoadShaderFromString(shadertype, src);
    	src = None;
    	return r;
    
    #//-----------------------------------------------------------------------------
    #//
    #// Read and compile the Shader source code.
    #// Equivalent to LoadShaderFromFile.
    #//
    def LoadShader(self,shadertype, filename):
    	src = self.LoadSource(filename);
    	r = self.LoadShaderFromString(shadertype, src);
    	src = None;
    	return r;

    def getLoc(self):
        print "#//Handles for textures in ball program"

        self.Atp = GL.glGetUniformLocation(self.program,"texturePosition");
        self.Atc = GL.glGetUniformLocation(self.program,"textureColors");
        self.Atsz =GL.glGetUniformLocation(self.program,"textureSizes");
        self.Atsc =GL.glGetUniformLocation(self.program,"textureScale");
        self.ATextureSize =GL.glGetUniformLocation(self.program,"textureSize");
        self.pos_loc= GL.glGetAttribLocation(self.program, "position") 
        self.texco_loc= GL.glGetAttribLocation(self.program, "texco") 
#        self.off_loc= GL.glGetAttribLocation(self.program, "posoffset") 

        print "#//Handles for textures in stick program"
        self.Btsh = GL.glGetUniformLocation(self.Bprogram,"textureShrink");
        self.Btp = GL.glGetUniformLocation(self.Bprogram,"texturePosition");
        self.Btc = GL.glGetUniformLocation(self.Bprogram,"textureColors");
        self.Btsz = GL.glGetUniformLocation(self.Bprogram,"textureSizes");
        self.Btsc = GL.glGetUniformLocation(self.Bprogram,"textureScale");
        self.BTextureSize =GL.glGetUniformLocation(self.Bprogram,"textureSize");
        self.BTextureBondsSize =GL.glGetUniformLocation(self.Bprogram,"textureBondsSize");

        self.stick_pos_loc= GL.glGetAttribLocation(self.Bprogram, "position") 
        self.texco1_loc= GL.glGetAttribLocation(self.Bprogram, "texco1") 
        self.texco2_loc= GL.glGetAttribLocation(self.Bprogram, "texco2") 
        self.texco3_loc= GL.glGetAttribLocation(self.Bprogram, "texco3") 
#        self.pos_loc1= GL.glGetAttribLocation(self.Bprogram, "position1") 
#        self.pos_loc2= GL.glGetAttribLocation(self.Bprogram, "position2") 


    #pragma mark -
    #pragma mark Initialization Methods
    #//-----------------------------------------------------------------------------
    #//
    #// OpenGL init for Buffers and Textures and shaders
    #//
    def initGL(self,):
        print "init parent GL"
        self._initGL();
        print "init GLSL"
        self.initGLSL();
        self.getLoc()
        #validate the program
    #//-----------------------------------------------------------------------------
    #//
    #// OpengGl init for GLSL
    #//
    def initGLSL(self,):
        #//
        #// Setup Atom GLSL
        #//
        if self.atomVertexShaderProgramName=="":
            self.vertex_shader = self.LoadShaderFromString(glextlib.GL_VERTEX_SHADER, 
            							 self.atomVertexShaderProgramSource);
            print ":: Using internal atom vertex shaders";
        else :
            self.vertex_shader = self.LoadShader(glextlib.GL_VERTEX_SHADER, self.atomVertexShaderProgramName);
            print ":: Using external atom vertex shaders " + self.atomVertexShaderProgramName;
        
        if(self.vertex_shader == 0) :
            print( "Atoms vertex shader error\n");
            return 0
        
        if (self.atomFragmentShaderProgramName==""):
            self.fragment_shader = self.LoadShaderFromString(glextlib.GL_FRAGMENT_SHADER, 
            								   self.atomFragmentShaderProgramSource);
            print ":: Using internal atom fragment shaders";
        else:
            self.fragment_shader = self.LoadShader(glextlib.GL_FRAGMENT_SHADER, 
            							 self.atomFragmentShaderProgramName.c_str());
            print ":: Using external atom fragment shaders " + self.atomFragmentShaderProgramName;
        
        if(self.fragment_shader == 0) :
            print( "Atoms fragment shader error\n");
            return 0
        
        self.program = glextlib.glCreateProgram();
        if (self.vertex_shader) :glextlib.glAttachShader(self.program, self.vertex_shader);  
        if (self.fragment_shader): glextlib.glAttachShader(self.program, self.fragment_shader);
    	
        glextlib.glLinkProgram(self.program);
        lStatus = 0x7FFFFFFF
        lStatus = glextlib.glGetProgramiv(self.program,
                                           glextlib.GL_LINK_STATUS,
                                           lStatus )
        if lStatus == 0:
            log = ""
            charsWritten  = 0
            progInfoLog = '\0' * 2048
            charsWritten, infoLog = glextlib.glGetProgramInfoLog(self.program, len(progInfoLog), 
                                                             charsWritten, progInfoLog)
            print "shader didn't link"
            print "log ",progInfoLog
            return 0
#        else:
#            glextlib.glValidateProgram(self.program)
#            lStatus = 0x7FFFFFFF
#            lStatus = glextlib.glGetProgramiv(self.program,
#                                               glextlib.GL_VALIDATE_STATUS,
#                                               lStatus )
#            if lStatus == 0:
#                log = ""
#                charsWritten  = 0
#                progInfoLog = '\0' * 2048
#                charsWritten, infoLog = glextlib.glGetProgramInfoLog(self.program, len(progInfoLog), 
#                                                                 charsWritten, progInfoLog)
#                print "ball shader did not validate, status:", lStatus
#                print "ball shader not validated",progInfoLog
#                return 0

        	#//
        	#// Setup Bond GLSL
        	#//
        if (self.bondVertexShaderProgramName==""):
            self.Bvertex_shader = self.LoadShaderFromString(glextlib.GL_VERTEX_SHADER, 
            									  self.bondVertexShaderProgramSource);
            print ":: Using internal bond vertex shaders" ;
        else:
            self.Bvertex_shader = LoadShader(glextlib.GL_VERTEX_SHADER, self.bondVertexShaderProgramName);
            print ":: Using external bond vertex shaders " + self.bondVertexShaderProgramName ;

        if(self.Bvertex_shader == 0):
            print( "Bonds vertex  shader error\n");
            return 0
        
        if (self.bondFragmentShaderProgramName==""):
            self.Bfragment_shader = self.LoadShaderFromString(glextlib.GL_FRAGMENT_SHADER, 
            self.bondFragmentShaderProgramSource);
            print ":: Using internal bond fragment shaders";
        else:
            self.Bfragment_shader = self.LoadShader(glextlib.GL_FRAGMENT_SHADER, self.bondFragmentShaderProgramName);
            print ":: Using external bond fragment shaders " + self.bondFragmentShaderProgramName;
        if(self.Bfragment_shader == 0): 
            print( "Bonds fragment  shader error\n");
            return 0
        
        self.Bprogram = glextlib.glCreateProgram();
        if (self.Bvertex_shader) :glextlib.glAttachShader(self.Bprogram, self.Bvertex_shader);  
        if (self.Bfragment_shader): glextlib.glAttachShader(self.Bprogram, self.Bfragment_shader);
        
        glextlib.glLinkProgram(self.Bprogram);
        lStatus = 0x7FFFFFFF
        lStatus = glextlib.glGetProgramiv(self.Bprogram,
                                           glextlib.GL_LINK_STATUS,
                                           lStatus )
        if lStatus == 0:
            log = ""
            charsWritten  = 0
            progInfoLog = '\0' * 2048
            charsWritten, infoLog = glextlib.glGetProgramInfoLog(self.Bprogram, len(progInfoLog), 
                                                             charsWritten, progInfoLog)
            print "shader didn't link"
            print "log ",progInfoLog
            return 0
#        else:
#            glextlib.glValidateProgram(self.Bprogram)
#            lStatus = 0x7FFFFFFF
#            lStatus = glextlib.glGetProgramiv(self.Bprogram,
#                                               glextlib.GL_VALIDATE_STATUS,
#                                               lStatus )
#            if lStatus == 0:
#                log = ""
#                charsWritten  = 0
#                progInfoLog = '\0' * 2048
#                charsWritten, infoLog = glextlib.glGetProgramInfoLog(self.Bprogram, len(progInfoLog), 
#                                                                 charsWritten, progInfoLog)
#                print "stick shader did not validate, status:", lStatus
#                print "stick shader not validated ",progInfoLog
#                return 0        
