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

"""
#// $Id: AtomAndBondShader.py,v 1.6.4.1 2017/07/13 22:27:37 annao Exp $
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
#// Adapted to Python by Ludovic Autin                                 #//
#// May 2011                                                           #//
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
"""
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
    from opengltk.extent import _gllib
    from opengltk.extent import glextlib
    GL.GL_STREAM_DRAW = glextlib.GL_STREAM_DRAW_ARB
    GL.GL_ELEMENT_ARRAY_BUFFER = glextlib.GL_ELEMENT_ARRAY_BUFFER # 0x8893
    GL.GL_ARRAY_BUFFER = glextlib.GL_ARRAY_BUFFER_ARB #0x8892
    glBindBuffer = glextlib.glBindBufferARB#GL.glBindBuffer
    glBufferData = glextlib.glBufferDataARB
    GL.glGenBuffers = glextlib.glGenBuffersARB
    GL.glTexImage2D = _gllib.glTexImage2D
    for k in glextlib.__dict__:
        if not hasattr(GL,k):
            setattr(GL,k, glextlib.__dict__[k])
else :
    from OpenGL import GL
    from OpenGL import GLU

GL_TEXTURE_RECTANGLE_NV = 0x84F5
GL_RGB32F_ARB = 0x8815
GL_RGBA32F_ARB = 0x8814#from OpenGL import GL as oGL   
import sys

#from OpenGL.arrays import vbo

sizeof = sys.getsizeof


#in case 
#35040
#// Identifier of the different buffers  used by the shaders
BUF_VERTICE=0
BUF_INDICE=1
BUF_TCOORD=2
BUF_TCOORD0=3
BUF_TCOORD1=4
BUF_TCOORD2=5
BUF_STICKV=2
BUF_STICKI=3

#// idenfifier of the different textures used by the shaders
TEX_POSITIONS=0
TEX_COLORS=1
TEX_SIZES=2
TEX_ASCALES=3
TEX_BSCALES=4
TEX_SHRINKS=5
#    	#// vertex, indices, and 4 for coordinates
NB_BUFFERS=4#7+2;
#    	#// atom position, atom color, atom size, atom scale, bond scale and bond shrink,framebuffer
NB_TEXTURES=7;

CUBE_VERTICES =  np.array([
         [-1,  -1, -1],
         [1,  -1, -1],
         [1,  -1, 1],
         [-1,  -1, 1],
         [-1,  1, -1],
         [1,  1, -1],
         [1,  1, 1],
         [-1,  1, 1],
    ],'f')

#6 quad face
CUBE_INDICES = np.array([
         [0,  1, 2],
         [0,  2, 3],
         [1,  5, 6],
         [1,  6, 2],
         [4,  6, 5],
         [4,  7, 6],
         [0,  7, 4],
         [0,  3, 7],
         [0,  5, 1],
         [0,  4, 5],
         [3,  2, 6],
         [3,  6, 7],    ],np.uint32)

import array

v=np.array([[ 1.,  1.,  0.],
       [ 0.,  1.,  0.],
       [ 0.,  0.,  0.],
       [ 1.,  0.,  0.],
       [ 1.,  1.,  1.],
       [ 0.,  1.,  1.],
       [ 0.,  0.,  1.],
       [ 1.,  0.,  1.]])

verts = np.array( [  
-0.285437,-0.744976,-0.471429,  
-0.285437,-0.744976,-2.471429,  
1.714563,-0.744976,-2.471429,  
1.714563,-0.744976,-0.471429,  
-0.285437,1.255024,-0.471429,  
-0.285437,1.255024,-2.471429,  
1.714563,1.255024,-2.471429,  
1.714563,1.255024,-0.471429] ) 
  
faces = np.array( [  
4,5,1,  
5,6,2,  
6,7,3,  
4,0,7,  
0,1,2,  
7,6,5,  
0,4,1,  
1,5,2,  
2,6,3,  
7,0,3,  
3,0,2,  
4,7,5] ) 


def hexahedron(radius):
        """
        Create the mesh data of a hexahedron of a given radius
        
        @type  radius: float
        @param radius: radius of the embeding sphere
        
        @rtype:   array
        @return:  vertex,face, face normal of the hexahedron        
        """        
        faceIndices = ( (3,2,0), (1,2,3), #bot
                        (4,6,7), (7,6,5), #top
                        (0,2,6), (6,4,0), #front
                        (2,1,5), (5,6,2), #right
                        (1,3,7), (7,5,1), #rear
                        (3,0,4), (4,7,3), #left
                        )
        a = 1.0
        faceNormals = None
        diameter = radius * 2
        width  = diameter
        height = diameter
        depth  = radius*sin(radians(45.))
        boundingRadius = radius
        surfaceArea = (sqrt(3) / 2) * (radius ** 2)
        volume = (4/3) * (radius ** 3)
        _corners = (
            (-radius, 0.0, -depth), (radius, 0.0, -depth),
            (0.0, -radius, -depth), (0.0, radius, -depth),
            (-radius, 0.0, depth), (radius, 0.0, depth),
            (0.0, -radius, depth), (0.0, radius, depth),)

        return np.array(_corners,'f'),np.array(faceIndices,np.uint32),np.array(faceNormals,float)

#verts,faces,fn= hexahedron(1.0) 
verts=CUBE_VERTICES
faces=CUBE_INDICES
#verts = np.array([[-1,0,0],[1,0,0],[1,1,0],[-1,1,0]],'f')*5.0#+np.array([9.2688, 9.7873, 7.9671],'f')
#faces = np.array([ [0,1,2],[0,2,3]],np.uint32)
#from OpenGLContext.scenegraph.basenodes import Sphere
#v_vbo,i_vbo,count = Sphere(  radius = 1  ).compile()

##//Return power of 2 superior than nb
##//Used for setting size of textures
def nextPower2(nb):
    power = 2;
    while(power <= nb):
        power = power << 1;
    return power;

class AtomAndBondShaders:
    def __init__(self,) :    	
#    	#// GL variables
        self.buffer_id =np.ones(NB_BUFFERS,int);
        self.texture_id= np.ones(NB_TEXTURES,int);#GLuint
        self.vertice=None;
        self.vertice_sticks=None
        self.texturecoord=None;
        self.texturecoord0=None;
        self.texturecoord1=None;
        self.texturecoord2=None;
        self.vertice_offset=None
        self.vertice_offset1=None
        self.vertice_offset2=None
        self.indice=None;
        self.indice_sticks=None
        self.colors=None;
        self.positions=None;
        self.atomScales=None;
        self.bondScales=None;
        self.sizes=None;
        self.shrinks=None;

        self.aScales=0.05;
        self.bScales=0.05;
        self.bshrinks=0.1;

        self.nbAtoms=nbBonds=0;
        self.textureSize=0
        self.textureSizeBonds=0;
        #// Parameter passed to the shaders??
        self.normMax=0.0;
        self.useVBO=True
        self.links = None
        
        self.atomVertexShaderProgramName=""
        self.bondVertexShaderProgramName=""
        self.atomFragmentShaderProgramName=""
        self.bondFragmentShaderProgramName=""

        self.atomVertexShaderProgramSource=""
        self.bondVertexShaderProgramSource=""
        self.atomFragmentShaderProgramSource=""
        self.bondFragmentShaderProgramSource=""

        self.ivert = np.reshape(np.mgrid[-1:2:2,-1:2:2,-1:2:2].T, (8,3))#one cube coordinate
        #one cube faces indices        
        self.ifaces= np.array([[0,1,2],[0,2,3],
        [1,5,6],[1,6,2],
        [4,6,5],[4,7,6],
        [0,7,4],[0,3,7],
        [0,5,1],[0,4,5],
        [3,2,6],[3,6,7],
        ],int )                                        
        
##//Destructor
#AtomAndBondShaders::~AtomAndBondShaders() {
#	if(buffer_id) 		delete [] buffer_id;
#	if(texture_id) 		delete [] texture_id;
#	if(vertice) 		delete [] vertice;
#	if(texturecoord) 	delete [] texturecoord;
#	if(texturecoord0)	delete [] texturecoord0;
#	if(texturecoord1) 	delete [] texturecoord1;
#	if(texturecoord2) 	delete [] texturecoord2;
#	if(indice)			delete [] indice;
#	if(colors) 			delete [] colors;
#	if(positions) 		delete [] positions;
#	if(atomScales) 		delete [] atomScales;
#	if(bondScales) 		delete [] bondScales;
#	if(sizes) 			delete [] sizes;
#	if(shrinks) 		delete [] shrinks;
#}

#pragma mark -
#pragma mark Setters
##//-----------------------------------------------------------------------------
    def setColors(self,colors):#GLfloat
        self.colors = colors#memcpy((GLfloat*)self.colors,(GLfloat*) colors,4*nbAtoms*sizeof(GLfloat));

##//-----------------------------------------------------------------------------
    def setPositions(self,positions):#GLfloat
        self.positions = positions;

##//-----------------------------------------------------------------------------
    def setAtomScales(self,scales):
        self.atomScales = scales #   memcpy((GLfloat*)self.atomScales,(GLfloat*) scales,nbAtoms*sizeof(GLfloat));

##//-----------------------------------------------------------------------------
    def setBondScales(self,scales):
        self.bondScales=scales#memcpy((GLfloat*)self.bondScales,(GLfloat*) scales,nbBonds*sizeof(GLfloat));

##//-----------------------------------------------------------------------------
    def setSizes(self,sizes):
        self.sizes = sizes#memcpy((GLfloat*)self.sizes,(GLfloat*) sizes,nbAtoms*sizeof(GLfloat));

##//-----------------------------------------------------------------------------
    def setShrinks(self, shrinks):
        self.shrinks=shrinks#memcpy((GLfloat*)self.shrinks,(GLfloat*) shrinks,nbAtoms*sizeof(GLfloat));

##//-----------------------------------------------------------------------------
    def setAtomFragmentShaderProgramName(self,fragmentShaderProgramName):
        self.atomFragmentShaderProgramName = fragmentShaderProgramName

##//-----------------------------------------------------------------------------
    def setBondFragmentShaderProgramName(self, fragmentShaderProgramName):
        self.bondFragmentShaderProgramName = fragmentShaderProgramName;

##//-----------------------------------------------------------------------------
    def setAtomVertexShaderProgramName(self, vertexShaderProgramName):
        self.atomVertexShaderProgramName = vertexShaderProgramName;

##//-----------------------------------------------------------------------------
    def setBondVertexShaderProgramName(self,vertexShaderProgramName):
        self.bondVertexShaderProgramName = vertexShaderProgramName;

    #//-----------------------------------------------------------------------------
    def setAtomFragmentShaderProgramSource(self,src):
    	self.atomFragmentShaderProgramSource = str(src);
    
    #//-----------------------------------------------------------------------------
    def setAtomVertexShaderProgramSource(self,src):
    	self.atomVertexShaderProgramSource = str(src);
    
    #//-----------------------------------------------------------------------------
    def setBondFragmentShaderProgramSource(self,src):
    	self.bondFragmentShaderProgramSource = str(src);
    
    #//-----------------------------------------------------------------------------
    def setBondVertexShaderProgramSource(self,src):
    	self.bondVertexShaderProgramSource = str(src);
    
    
    #pragma mark -
    #pragma mark Initialization Methods
    
    #//-----------------------------------------------------------------------------
    #//
    #// Init atoms parameters
    #//
    def initAtomTextures(self,colors,radius):
        #need to check shape of colors thatshould be (nbAtoms, 4)
        self.colors = np.array(colors,'f')#.flatten()
        if self.colors.shape[1] != 4 :
            one = np.ones( (self.colors.shape[0], 1), self.colors.dtype.char )
            self.colors = np.concatenate( (self.colors, one), 1 )
        self.sizes=np.array(radius,'f')
#        print self.sizes
        self.atomScales=np.ones(self.nbAtoms,'f')*self.aScales
#        print self.atomScales
#        print self.aScales
        
    def initAtomTextures_mol(self,dataAtoms):
        self.colors = np.array(map(lambda x: [x.colors['cpk'][0],x.colors['cpk'][1],x.colors['cpk'][2],1.0], dataAtoms),float)#.flatten()
        self.sizes=np.array(dataAtoms.vdwRadius,'f')
        self.atomScales=np.ones(self.nbAtoms,'f')*0.05    
    #//-----------------------------------------------------------------------------
    #//
    #// Init bond parameters
    #//
    def initBondTextures(self,):
        self.shrinks=np.ones(self.nbBonds,'f')*self.bshrinks
        self.bondScales=np.ones(self.nbBonds,'f')*self.bScales

    def initBondTextures_mol(self, dataBonds):
        self.shrinks=np.ones(self.nbBonds,'f')*0.1
        self.bondScales=np.ones(self.nbBonds,'f')*0.05
#        for i in range(self.nbBonds):
##    	for(int i = 0; i < nbBonds; ++i){
#    		self.shrinks[i] = dataBonds[i].shrink;
#    		self.bondScales[i] = dataBonds[i].scale;

    #//-----------------------------------------------------------------------------
    #//
    #// Init vertex positions and storage in textures
    #//
    def initVertTCoordIndice(self,):
        #cube vertices and faces
        self.ivert,self.ifaces,fn= hexahedron(1.0) 
        #Ball data
        nvert = self.nbAtoms# if  self.nbAtoms > self.nbBonds else self.nbBonds
        if nvert > 0:
            incr = np.repeat(np.arange(nvert)*len(self.ivert),len(self.ifaces))
            idata = np.tile(self.ifaces,(nvert,1))
            self.indice=idata+incr.reshape((len(incr),1))
            self.vertice = np.array(np.tile(self.ivert,(nvert,1)),'f')#added dirrectly the pshere position?
            self.vertice_offset = np.array(np.repeat(self.positions,len(self.ivert),axis=0),'f')
            self.texturecoord = np.vstack([np.repeat(np.arange(self.nbAtoms,dtype='f') % self.textureSize,len(self.ivert)),
                                           np.repeat(np.arange(self.nbAtoms,dtype='f') / self.textureSize,len(self.ivert))]).transpose().astype('f')#/self.textureSize
            self.texturecoord = np.floor(self.texturecoord)#/(self.textureSize-1)
            self.indice=np.array(self.indice,np.uint32)
        #add color, cant get it work with Texture !
#        self.vcolors = np.array(np.repeat(self.colors,len(self.ivert),axis=0),'f')
        #sticks data    
        nvert = self.nbBonds# if  self.nbAtoms > self.nbBonds else self.nbBonds
        if nvert > 0:
            incr = np.repeat(np.arange(nvert)*len(self.ivert),len(self.ifaces))
            idata = np.tile(self.ifaces,(nvert,1))
            self.indice_sticks=idata+incr.reshape((len(incr),1))
            self.vertice_sticks = np.array(np.tile(self.ivert,(nvert,1)),'f')#added dirrectly the pshere position?
    #        self.vertice_sticks_offset = np.array(np.repeat(self.positions,len(self.ivert),axis=0),'f')        
            vtext0 = np.array(self.links,int).transpose()
            #to fix 
            self.texturecoord0=np.vstack([np.repeat(vtext0[0] % self.textureSize,len(self.ivert)),
                                          np.repeat(vtext0[0] / self.textureSize,len(self.ivert))]).transpose().astype('f')#/self.textureSize
            self.texturecoord1=np.vstack([np.repeat(vtext0[1] % self.textureSize,len(self.ivert)),
                                          np.repeat(vtext0[1] / self.textureSize,len(self.ivert))]).transpose().astype('f')#/self.textureSize
            self.texturecoord2 = np.vstack([np.repeat(np.arange(self.nbBonds,dtype='f') % self.textureSizeBonds,len(self.ivert)),
                                            np.repeat(np.arange(self.nbBonds,dtype='f') / self.textureSizeBonds,len(self.ivert))]).transpose().astype('f')#/self.textureSizeBonds
            self.indice_sticks=np.array(self.indice_sticks,np.uint32)
            self.vertice_offset1 = np.array(np.repeat(np.take(self.positions,vtext0[0].astype(int),axis=0),len(self.ivert),axis=0),'f')
            self.vertice_offset2 = np.array(np.repeat(np.take(self.positions,vtext0[1].astype(int),axis=0),len(self.ivert),axis=0),'f')
            self.texturecoord0 = np.floor(self.texturecoord0)#/(self.textureSize-1)
            self.texturecoord1 = np.floor(self.texturecoord1)#/(self.textureSize-1)
            self.texturecoord2 = np.floor(self.texturecoord2)#/(self.textureSize-1)
            
#        self.vcolors1 = np.array(np.repeat(np.take(self.colors,vtext0[0].astype(int),axis=0),len(self.ivert),axis=0),'f')
#        self.vcolors2 = np.array(np.repeat(np.take(self.colors,vtext0[1].astype(int),axis=0),len(self.ivert),axis=0),'f')
    #//-----------------------------------------------------------------------------
    #//
    #// Init arrays used to fill the textures
    #//
    def setupBuffersAndTextures(self, nbAtoms, nbBonds, coords, bonds, colors, radii):
    	#//initialize nbAtoms 
    	self.nbAtoms = nbAtoms;
    	self.nbBonds = nbBonds;
    	self.positions = np.array(coords,'f')
    	#print self.positions
    	nbMax=self.nbAtoms if self.nbAtoms>self.nbBonds else self.nbBonds
    	
    	#// textureSize contains the side size of the  textures for atoms
    	self.textureSize=int(ceil(sqrt(self.nbAtoms)));
    	#// textureSizeBonds contains the side size of the  textures for bonds
    	self.textureSizeBonds=int(ceil(sqrt(self.nbBonds)));
    	
    	self.links=np.array(bonds,np.uint32)#new LinkAtom[nbBonds];
    	#for(int z=0; z<nbBonds; z++){
    	#	links[z].atom1 = dataBonds[z].idsrc;
    	#	links[z].atom2 = dataBonds[z].iddest;
    	#}
    	print "setupBuffersAndTextures"
    	self.initAtomTextures(colors,radii);
    	self.initBondTextures();
    	self.initVertTCoordIndice();
    	print "OK setupBuffersAndTextures"
        
    def setupBuffersAndTextures_mol(self, nbAtoms, nbBonds, dataAtoms, dataBonds):
    	#//initialize nbAtoms 
    	self.nbAtoms = nbAtoms;
    	self.nbBonds = nbBonds;
    	self.positions = dataAtoms.coords
    	print self.positions
    	nbMax=self.nbAtoms if self.nbAtoms>self.nbBonds else self.nbBonds
    	
    	#// textureSize contains the side size of the  textures for atoms
    	self.textureSize=int(ceil(sqrt(self.nbAtoms)));
    	#// textureSizeBonds contains the side size of the  textures for bonds
    	self.textureSizeBonds=int(ceil(sqrt(self.nbBonds)));
    	
    	self.links=dataBonds#new LinkAtom[nbBonds];
    	#for(int z=0; z<nbBonds; z++){
    	#	links[z].atom1 = dataBonds[z].idsrc;
    	#	links[z].atom2 = dataBonds[z].iddest;
    	#}
    	self.initAtomTextures_mol(dataAtoms);
    	self.initBondTextures_mol(dataBonds);
    	self.initVertTCoordIndice();

    #pragma mark -
    #pragma mark Update Parameters Methods
    
    def updateBuffersAndTextures(self, nbAtoms, nbBonds, coords, bonds, colors, radii):
        #//initialize nbAtoms 
        self.nbAtoms = nbAtoms;
        self.nbBonds = nbBonds;
        self.positions = np.array(coords,'f')
        #print self.positions
        nbMax=self.nbAtoms if self.nbAtoms>self.nbBonds else self.nbBonds
        #// textureSize contains the side size of the  textures for atoms
        self.textureSize=int(ceil(sqrt(self.nbAtoms)));
        #// textureSizeBonds contains the side size of the  textures for bonds
        self.textureSizeBonds=int(ceil(sqrt(self.nbBonds)));
        self.links=np.array(bonds,np.uint32)#new LinkAtom[nbBonds];
        self.initAtomTextures(colors,radii);
        self.initBondTextures();
        self.initVertTCoordIndice();
        self.activateAllTextures(update=True)
        self.updateVBO()

    def updateTextureSize(self):
        nbMax=self.nbAtoms if self.nbAtoms>self.nbBonds else self.nbBonds
        #// textureSize contains the side size of the  textures for atoms
        self.textureSize=int(ceil(sqrt(self.nbAtoms)));
        #// textureSizeBonds contains the side size of the  textures for bonds
        self.textureSizeBonds=int(ceil(sqrt(self.nbBonds)));
        
    #//-----------------------------------------------------------------------------
    def updatePositions(self, positions):
        #update textureSize
        self.activeTexture(GL.GL_TEXTURE0,self.texture_id[TEX_POSITIONS],np.array(positions,'f'),self.textureSize,size=3)    	
        #radius ?

    def updateSizes(self, sizes):
        self.sizes = np.array(sizes,'f')
        self.activeTexture(GL.GL_TEXTURE2,self.texture_id[TEX_SIZES],self.sizes,self.textureSize)   
        
    #//-----------------------------------------------------------------------------
    def updateColors(self, colors):
        self.colors = np.array(colors,'f')#.flatten()
        if self.colors.shape[1] != 4 :
            one = np.ones( (self.colors.shape[0], 1), self.colors.dtype.char )
            self.colors = np.concatenate( (self.colors, one), 1 )        
        self.activeTexture(GL.GL_TEXTURE1,self.texture_id[TEX_COLORS],self.colors,self.textureSize,size=4)
        print  self.texture_id,self.texture_id[TEX_COLORS]

    #//-----------------------------------------------------------------------------
    def updateAtomScales(self, newScale=1,scales=None):
        if scales is None :
            self.atomScales=np.ones(self.nbAtoms,'f')*newScale
        else :
            self.atomScales=np.array(scales[:],'f')
#    	self.activeTexture(GL.GL_TEXTURE2,self.texture_id[TEX_SIZES],self.sizes,self.textureSize) 
    	self.activeTexture(GL.GL_TEXTURE3,self.texture_id[TEX_ASCALES],self.atomScales,self.textureSize) 
    
    #//-----------------------------------------------------------------------------
    def updateBondScales(self, newScale=1,scales=None):
        if scales is None :
            self.bondScales=np.ones(self.nbBonds,'f')*newScale
        else :
            self.bondScales=np.array(scales[:],'f')
    	self.activeTexture(GL.GL_TEXTURE4,self.texture_id[TEX_BSCALES],self.bondScales,self.textureSizeBonds) 
    
    #//-----------------------------------------------------------------------------
    def updateShrinks(self, newShrink, shrinks=None):
        if shrinks is None :
            self.shrinks=np.ones(self.nbBonds,'f')*newShrink
        else :
            self.shrinks=np.array(shrinks[:],'f')
    	self.activeTexture(GL.GL_TEXTURE5,self.texture_id[TEX_SHRINKS],self.shrinks,self.textureSizeBonds) 
    
    
    #pragma mark -
    #pragma mark Opengl Drawing Methods
    
    #//-----------------------------------------------------------------------------
    def drawAtoms(self):
        pass
    
    #//-----------------------------------------------------------------------------
    def drawBonds(self):
        pass



    def generateVBOsimple(self):
        #Create the Vertex Array Object
#        v_vbo,i_vbo,count
        self.vertice_vbo=v_vbo#vbo.VBO( verts )#self.vertice
        self.indice_vbo=i_vbo#vbo.VBO( faces, target=GL.GL_ELEMENT_ARRAY_BUFFER )#self.indice
#        self.texturecoord_vbo=vbo.VBO( self.texturecoord,usage = GL.GL_STREAM_DRAW)

    def makeVBO(self,np_data,buff_indice,target=GL.GL_ARRAY_BUFFER,usage=GL.GL_STREAM_DRAW):
        buffer_id=int(self.buffer_id[buff_indice])#GL.glGenBuffers()
        glBindBuffer(target,int(buffer_id));
        glBufferData(target,np_data.nbytes,np_data.flatten(),usage)#,GL.GL_STREAM_DRAW);
        return buffer_id
        #4*len(vertices)*len(vertices[0])
    
    def bindVBO(self,buffer_id,np_data,target=GL.GL_ARRAY_BUFFER,usage=GL.GL_STREAM_DRAW):    
        glBindBuffer(target,int(buffer_id))
        glBufferData(target,np_data.nbytes,np_data.flatten(),usage)#,GL.GL_STREAM_DRAW);
        
    def generateVBO(self):
        #Create the Vertex Array Object thatcarry everthing coords, uv, normal if any
        #ball_vbo
        self.ball_data = np.hstack([self.vertice,self.texturecoord])#,self.vertice_offset])#,self.vcolors])
        self.vertice_vbo=self.makeVBO(self.ball_data,BUF_VERTICE) #vbo.VBO( self.ball_data, usage = GL.GL_STREAM_DRAW)#self.vertice
        self.indice_vbo=self.makeVBO( self.indice,BUF_INDICE,usage = GL.GL_STREAM_DRAW, target=GL.GL_ELEMENT_ARRAY_BUFFER )#self.indice
        #stick_vbo
        self.stick_data = np.hstack([self.vertice_sticks,self.texturecoord0,self.texturecoord1,self.texturecoord2])#,
#                                     self.vertice_offset1,self.vertice_offset2,])
#                                     self.vcolors1,self.vcolors2])        
        self.stick_vertice_vbo=self.makeVBO( self.stick_data,BUF_STICKV)# usage = GL.GL_STREAM_DRAW)#self.vertice
        self.stick_indice_vbo=self.makeVBO( self.indice_sticks,BUF_STICKI)#usage = GL.GL_STREAM_DRAW, target=GL.GL_ELEMENT_ARRAY_BUFFER )#self.indice

    def updateVBO(self):
        #Create the Vertex Array Object thatcarry everthing coords, uv, normal if any
        #ball_vbo
        self.ball_data = np.hstack([self.vertice,self.texturecoord])#,self.vertice_offset])#,self.vcolors])
#        self.vertice_vbo.set_array(self.ball_data)#  =vbo.VBO( self.ball_data, usage = GL.GL_STREAM_DRAW)#self.vertice
#        self.indice_vbo.set_array(self.indice)#=vbo.VBO( self.indice,usage = GL.GL_STREAM_DRAW, target=GL.GL_ELEMENT_ARRAY_BUFFER )#self.indice
        #stick_vbo
        self.stick_data = np.hstack([self.vertice_sticks,self.texturecoord0,self.texturecoord1,self.texturecoord2])#,
#                                     self.vertice_offset1,self.vertice_offset2,])
#                                     self.vcolors1,self.vcolors2])        
#        self.stick_vertice_vbo.set_array(self.stick_data)#=vbo.VBO( self.stick_data, usage = GL.GL_STREAM_DRAW)#self.vertice
#        self.stick_indice_vbo.set_array(self.indice_sticks)#=vbo.VBO( self.indice_sticks,usage = GL.GL_STREAM_DRAW, target=GL.GL_ELEMENT_ARRAY_BUFFER )#self.indice

    def generateVBO_VBO(self):
        #Create the Vertex Array Object thatcarry everthing coords, uv, normal if any
        #ball_vbo
        self.ball_data = np.hstack([self.vertice,self.texturecoord])#,self.vertice_offset])#,self.vcolors])
        self.vertice_vbo=vbo.VBO( self.ball_data, usage = GL.GL_STREAM_DRAW)#self.vertice
        self.indice_vbo=vbo.VBO( self.indice,usage = GL.GL_STREAM_DRAW, target=GL.GL_ELEMENT_ARRAY_BUFFER )#self.indice
        #stick_vbo
        self.stick_data = np.hstack([self.vertice_sticks,self.texturecoord0,self.texturecoord1,self.texturecoord2])#,
#                                     self.vertice_offset1,self.vertice_offset2,])
#                                     self.vcolors1,self.vcolors2])        
        self.stick_vertice_vbo=vbo.VBO( self.stick_data, usage = GL.GL_STREAM_DRAW)#self.vertice
        self.stick_indice_vbo=vbo.VBO( self.indice_sticks,usage = GL.GL_STREAM_DRAW, target=GL.GL_ELEMENT_ARRAY_BUFFER )#self.indice

    def updateVBO_VBO(self):
        #Create the Vertex Array Object thatcarry everthing coords, uv, normal if any
        #ball_vbo
        self.ball_data = np.hstack([self.vertice,self.texturecoord,self.vertice_offset])#,self.vcolors])
        self.vertice_vbo.set_array(self.ball_data)#  =vbo.VBO( self.ball_data, usage = GL.GL_STREAM_DRAW)#self.vertice
        self.indice_vbo.set_array(self.indice)#=vbo.VBO( self.indice,usage = GL.GL_STREAM_DRAW, target=GL.GL_ELEMENT_ARRAY_BUFFER )#self.indice
        #stick_vbo
        self.stick_data = np.hstack([self.vertice_sticks,self.texturecoord0,self.texturecoord1,self.texturecoord2,
                                     self.vertice_offset1,self.vertice_offset2,])
#                                     self.vcolors1,self.vcolors2])        
        self.stick_vertice_vbo.set_array(self.stick_data)#=vbo.VBO( self.stick_data, usage = GL.GL_STREAM_DRAW)#self.vertice
        self.stick_indice_vbo.set_array(self.indice_sticks)#=vbo.VBO( self.indice_sticks,usage = GL.GL_STREAM_DRAW, target=GL.GL_ELEMENT_ARRAY_BUFFER )#self.indice


    def useTexture(self,gl_id,tex_id):
        GL.glActiveTexture(gl_id);
#        GL.glBindTexture(GL_TEXTURE_RECTANGLE_NV,int(tex_id));#GL.GL_TEXTURE_2D
        GL.glBindTexture(GL.GL_TEXTURE_2D,int(tex_id));#GL.GL_TEXTURE_2D
        

    def activeTexture(self,gl_id,tex_id,data,texSize,size=1,update=False):    
    	GL.glActiveTexture(gl_id);
#    	ttype = GL_TEXTURE_RECTANGLE_NV
    	ttype = GL.GL_TEXTURE_2D#or GL_TEXTURE_RECTANGLE_NV
    	print "texSize ",texSize,ttype,int(tex_id),data.shape#,data
    	#flatten the data ?
    	if size == 1 :
            GL.glEnable(ttype);
            GL.glBindTexture(ttype,int(tex_id))#int(self.texture_id[tex_id]));
            GL.glTexImage2D (ttype, 0, GL.GL_INTENSITY , 
                                  texSize, texSize, 0, GL.GL_LUMINANCE, GL.GL_FLOAT, 
                      data.flatten());
    	elif size == 4 :
#            ttype = GL_TEXTURE_RECTANGLE_NV#or GL_TEXTURE_RECTANGLE_NV
            GL.glEnable(ttype);
            GL.glBindTexture(ttype,int(tex_id))#int(self.texture_id[tex_id]));
            GL.glTexImage2D (ttype, 0, GL_RGBA32F_ARB , 
                      texSize, texSize, 0, GL.GL_RGBA, GL.GL_FLOAT, 
                      data.flatten());
    	else :#3
            try :
                GL.glEnable(ttype);
                GL.glBindTexture(ttype,int(tex_id))#int(self.texture_id[tex_id]));
                #print "try glTexImage2D",GL_RGB32F_ARB,GL.GL_RGB,data.flatten()
                GL.glTexImage2D (ttype, 0, GL_RGB32F_ARB , 
                      texSize, texSize, 0, GL.GL_RGB, GL.GL_FLOAT, 
                      data.flatten());
            except :
                print "something rong"
    	print "ok teximage"
    	if not update :
         GL.glTexParameteri(ttype,GL.GL_TEXTURE_MIN_FILTER,GL.GL_NEAREST);
         GL.glTexParameteri(ttype,GL.GL_TEXTURE_MAG_FILTER,GL.GL_NEAREST);
    	GL.glDisable(ttype);

    def useAllTextures(self):
    	#// atom positions
    	self.useTexture(GL.GL_TEXTURE0,self.texture_id[TEX_POSITIONS])#,self.positions,size=3)    	
    	#// atom colors
    	self.useTexture(GL.GL_TEXTURE1,self.texture_id[TEX_COLORS])#,self.colors,size=4)	
    	#// atom sizes
    	self.useTexture(GL.GL_TEXTURE2,self.texture_id[TEX_SIZES])#,self.sizes) 
    	#// atom scales
    	self.useTexture(GL.GL_TEXTURE3,self.texture_id[TEX_ASCALES])#,self.atomScales) 
    	#//   bond scales
    	self.useTexture(GL.GL_TEXTURE4,self.texture_id[TEX_BSCALES])#,self.bondScales) 
    	#//   bond shrinks
    	self.useTexture(GL.GL_TEXTURE5,self.texture_id[TEX_SHRINKS])#,self.shrinks) 

    def activateAllTextures(self,update=False):
    	if not update :
            self.texture_id=GL.glGenTextures(NB_TEXTURES).tolist()
    	print "activateAllTextures",self.texture_id
    	print  self.texture_id,self.texture_id[TEX_COLORS]
    	print "#// atom positions"
    	self.activeTexture(GL.GL_TEXTURE0,self.texture_id[TEX_POSITIONS],self.positions,self.textureSize,size=3)    	
    	print "#// atom colors"
    	self.activeTexture(GL.GL_TEXTURE1,self.texture_id[TEX_COLORS],self.colors,self.textureSize,size=4)	
    	print "#// atom sizes"
    	self.activeTexture(GL.GL_TEXTURE2,self.texture_id[TEX_SIZES],self.sizes,self.textureSize) 
    	print "#// atom scales"
    	self.activeTexture(GL.GL_TEXTURE3,self.texture_id[TEX_ASCALES],self.atomScales,self.textureSize) 
    	print "#//   bond scales"
    	self.activeTexture(GL.GL_TEXTURE4,self.texture_id[TEX_BSCALES],self.bondScales,self.textureSizeBonds) 
    	print "#//   bond shrinks"
    	self.activeTexture(GL.GL_TEXTURE5,self.texture_id[TEX_SHRINKS],self.shrinks,self.textureSizeBonds) 
        
    def bindAllBuffer(self):
    	nbMax=self.nbAtoms if self.nbAtoms>self.nbBonds else self.nbBonds;
    	print "nmax is ",  nbMax

#    	print BUF_VERTICE,int(self.buffer_id[BUF_VERTICE]),nbMax*3*8*sizeof(float)   
#    	glBindBuffer(GL.GL_ARRAY_BUFFER,int(self.buffer_id[BUF_VERTICE]));
#    	glBufferData(GL.GL_ARRAY_BUFFER,self.vertice.nbytes,self.vertice,GL.GL_STREAM_DRAW)#,GL.GL_STREAM_DRAW);
#    	print "bin buffer   BUF_VERTICE",self.vertice.flatten(),self.buffer_id[BUF_INDICE],BUF_INDICE,self.indice.flatten()     
#
#    	glBindBuffer(GL.GL_ELEMENT_ARRAY_BUFFER,int(self.buffer_id[BUF_INDICE]));
#    	glBufferData(GL.GL_ELEMENT_ARRAY_BUFFER,self.indice.nbytes,self.indice,GL.GL_STREAM_DRAW);
#    	print "bin buffer   BUF_INDICE"     

    	glBindBuffer(GL.GL_ARRAY_BUFFER,int(self.buffer_id[BUF_TCOORD]));
    	glBufferData(GL.GL_ARRAY_BUFFER,self.texturecoord.nbytes,self.texturecoord,GL.GL_STREAM_DRAW);
    	print "bin buffer   BUF_TCOORD"     

    	glBindBuffer(GL.GL_ARRAY_BUFFER,int(self.buffer_id[BUF_TCOORD0]));
    	glBufferData(GL.GL_ARRAY_BUFFER,self.texturecoord0.nbytes,self.texturecoord0,GL.GL_STREAM_DRAW);

    	print "bin buffer   BUF_TCOORD0"     
    	glBindBuffer(GL.GL_ARRAY_BUFFER,int(self.buffer_id[BUF_TCOORD1]));
    	glBufferData(GL.GL_ARRAY_BUFFER,self.texturecoord1.nbytes,self.texturecoord1,GL.GL_STREAM_DRAW);
         #GL_TEXTURE_BUFFER??
    	print "bin buffer   BUF_TCOORD1"     
    	glBindBuffer(GL.GL_ARRAY_BUFFER,int(self.buffer_id[BUF_TCOORD2]));
    	glBufferData(GL.GL_ARRAY_BUFFER,self.texturecoord2.nbytes,self.texturecoord2,GL.GL_STREAM_DRAW);
    	print "bin buffer   BUF_TCOORD2"   


    def setupInstanceMatrice(self,):
        self.matrix_loc=pos=GL.glGetAttribLocation(self.program, "transformmatrix");
        pos1 = pos + 0; 
        pos2 = pos + 1; 
        pos3 = pos + 2; 
        pos4 = pos + 3; 
#        GL.glEnableVertexAttribArray(pos1);
#        GL.glEnableVertexAttribArray(pos2);
#        GL.glEnableVertexAttribArray(pos3);
#        GL.glEnableVertexAttribArray(pos4);
#        
#        glBindBuffer(GL_ARRAY_BUFFER, VBO_containing_matrices);
#        glVertexAttribPointer(pos1, 4, GL_FLOAT, GL_FALSE, sizeof(GLfloat) * 4 * 4, (void*)(0));
#        glVertexAttribPointer(pos2, 4, GL_FLOAT, GL_FALSE, sizeof(GLfloat) * 4 * 4, (void*)(sizeof(float) * 4));
#        glVertexAttribPointer(pos3, 4, GL_FLOAT, GL_FALSE, sizeof(GLfloat) * 4 * 4, (void*)(sizeof(float) * 8));
#        glVertexAttribPointer(pos4, 4, GL_FLOAT, GL_FALSE, sizeof(GLfloat) * 4 * 4, (void*)(sizeof(float) * 12));
#        glVertexAttribDivisor(pos1, 1);
#        glVertexAttribDivisor(pos2, 1);
#        glVertexAttribDivisor(pos3, 1);
#        glVertexAttribDivisor(pos4, 1);        
    #//-----------------------------------------------------------------------------
    def _initGL(self) :
    	print "#// Generation of the buffers and textures ID"
    	self.buffer_id=GL.glGenBuffers(NB_BUFFERS)#GL.glGenBuffers(NB_BUFFERS,self.buffer_id);
    	self.texture_id=GL.glGenTextures(NB_TEXTURES).tolist()
    	print self.buffer_id
    	print self.texture_id
    	if self.useVBO :
         self.generateVBO()#simple()
#         self.activateAllTextures()
         return
    	self.buffer_id=GL.glGenBuffers(NB_BUFFERS);
    
    	print "#// loading the buffers ",self.buffer_id      
    	nbMax=self.nbAtoms if self.nbAtoms>self.nbBonds else self.nbBonds;
    	self.bindAllBuffer()
    	print "#//Generation of the textures used by the shaders"
    	self.texture_id=GL.glGenTextures(NB_TEXTURES);
    	print self.texture_id
    	#// atom positions
    	self.activeTexture(GL.GL_TEXTURE0,TEX_POSITIONS,self.positions)
    	
    	#// atom colors
    	self.activeTexture(GL.GL_TEXTURE1,TEX_COLORS,self.colors)
    	
    	#// atom sizes
    	self.activeTexture(GL.GL_TEXTURE2,TEX_SIZES,self.sizes)
    	
    	#// atom scales
    	self.activeTexture(GL.GL_TEXTURE3,TEX_ASCALES,self.atomScales)
    	
    	#//   bond scales
    	self.activeTexture(GL.GL_TEXTURE4,TEX_BSCALES,self.bondScales)
    	
    	#//   bond shrinks
    	self.activeTexture(GL.GL_TEXTURE5,TEX_SHRINKS,self.shrinks)  
    	print "all textures activated "       
