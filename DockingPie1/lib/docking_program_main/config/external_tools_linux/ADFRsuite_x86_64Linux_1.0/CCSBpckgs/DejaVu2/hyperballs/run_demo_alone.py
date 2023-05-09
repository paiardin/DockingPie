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

import os, sys
sys.path.insert(0,".")
#sys.path.append("/Users/ludo/DEV/")
#sys.path.append("/Users/ludo/Library/Preferences/MAXON/CINEMA 4D R12_C333CB6C/plugins/ePMV/mgl64/MGLToolsPckgs")

import numpy
import numpy as np 
#import Image

from opengltk.OpenGL.GL import *
from opengltk.OpenGL.GLU import *
from OpenGL.GLUT import *
from OpenGL import GL as oGL
#import pyQutemol.Qutemol.glew_wrap as glew
from pyQutemol.Qutemol.trackball import glTrackball
from pyQutemol.Qutemol.quaternion import quaternion

from Pmv.moleculeViewer import MoleculeViewer

#from IPython.Shell import IPShellEmbed
#ipshell = IPShellEmbed()

ERRGL_OK = 0
ERRGL_NO_FS = 1
ERRGL_NO_VS = 2
ERRGL_NO_FBO_SHADOWMAP = 4
ERRGL_NO_FBO_HALO = 8
ERRGL_NO_FBO_AO = 16
ERRGL_NO_GLEW = 32

width = 1024
height = 1024
oldX, oldY = 0, 0
mustDoHQ = True
win_id=0

xCenter=0
yCenter=0
zCenter=0
xLength=0
yLength=0
zLength=0
bBox=numpy.zeros((2,3))
xTrans=0
yTrans=0
zTrans=0
xRot=0
yRot=0
zRot=0

def initDisplay():
#    glewInit()
#    res = 0
    glClearColor(1.0, 1.0, 1.0, 0.0);
    glClearDepth(1.0)
    glDepthFunc(GL_LESS)
    glEnable(GL_DEPTH_TEST)
#    glEnable(GL_LIGHTING)
#    glShadeModel(GL_SMOOTH)
    #//OpenGL modes
#    glEnable(GL_CULL_FACE);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glDepthFunc(GL_LESS);

    light_diffuse = [1.0, 0.0, 0.0, 1.0];  #/* Red diffuse light. */
    light_position = [1.0, 1.0, 1.0, 0.0];  #/* Infinite light location. */

    #/* Enable a single OpenGL light. */
    glLightfv(GL_LIGHT0, GL_DIFFUSE, light_diffuse);
    glLightfv(GL_LIGHT0, GL_POSITION, light_position);
    glEnable(GL_LIGHT0);
    glEnable(GL_LIGHTING);

#    glEnable(GL_CLIP_PLANE0)


#    glMatrixMode(GL_MODELVIEW); 
#    glLoadIdentity();
	#//for the frame rate
	#currentTime=previousTime=glutGet(GLUT_ELAPSED_TIME);
	#//	nbFrames=0;
	#//	nbFramesMax=2;
	

def setLightDir(d):
    f = (d[0], d[1], d[2], 0)
    glLightfv(GL_LIGHT0, GL_POSITION, f)

# This is not necessary any more - but I really have to figure out what I'm
# doing with the model-view matrices
#def getDirFromTrackball(mol):
#    # XXX this is complete wrong, but it shows that shadowing works
#    glPushMatrix()
#    gluLookAt(1,-3,-5,   0,0,0,   0,1,0)

#    #glMultMatrixd((-1*glTrackball.quat * mol.orien).asRotation())
#    d = glGetFloatv(GL_MODELVIEW_MATRIX)
#    glPopMatrix()
#    res = numpy.array([-d[2,0], -d[2,1], -d[2,2]])
#    res /= numpy.linalg.norm(res)
#    return res

def getGlLightPos():
    pos = glGetLightfv(GL_LIGHT0,GL_POSITION)
    x = glGetFloatv(GL_MODELVIEW_MATRIX).reshape(4,4)
    res = numpy.inner(x,pos.T)
    res /= numpy.linalg.norm(res)
    return -res[:3]


def setProjection(res):
    winx = width
    winy = height
    glMatrixMode(GL_PROJECTION)
    glLoadIdentity()
    nearPlane = 0.01#0.1
    farPlane = 10000#500
    size = 1.2
    ratio = size*winx/winy
#    if cgSettings.projmode == cgSettings.PERSPECTIVE:
#    gluPerspective(60.0, ratio, nearPlane, farPlane)
    gluPerspective(70.0,float(winx)/float(winy),0.1,500.0);
#    else:
#        glOrtho(-winx/winy, winx/winy,-1,1,40-2,40+200)
    glViewport(0,0,winx,winy)
    glMatrixMode(GL_MODELVIEW);  
    glLoadIdentity();


def drawScene(mol,shaders):
    drawFrame(mol,shaders)
#    glutSwapBuffers()
#    glutPostRedisplay();
    
def drawFrame(mol,shaders):
#    glPolygonMode( GL_FRONT_AND_BACK, GL_LINE );
#    shaders.useAllTextures()
    glClearColor(1.0, 1.0, 1.0, 0.0);
    
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)

    glMatrixMode(GL_MODELVIEW)
    glLoadIdentity()
#    // translation and rotation of the molecule
#    glTranslatef(xTrans+xCenter,yTrans+yCenter,zTrans+zCenter);
    #// parameters definition for the lighting
#    // camera position
    print xCenter,yCenter,bBox[0][2]-zLength,xCenter,yCenter,zCenter+1,0,1,0
    gluLookAt(xCenter,yCenter,bBox[0][2]-zLength,xCenter,yCenter,zCenter+1,0,1,0);
#    gluLookAt(0,0,40,   0,0,0,   0,1,0)
#    // displays shader information
#    displayInfo();
    glColor3f(1,1,1)
#    setProjection(600)	
		
#    // compute the rotations    
#    glPushMatrix();
#  3  glLoadIdentity();
#// translation and rotation of the molecule
#    glTranslatef(xTrans+xCenter,yTrans+yCenter,zTrans+zCenter);

    glMultMatrixd((glTrackball.quat * mol.orien).asRotation().reshape(16,))
   #// place the center of the molecule in (0,0,0) 
#    glTranslatef(-xCenter,-yCenter,-zCenter);
    
        
#    //Necessary initial draw call for Mac (bug ???)
#    oGL.glEnableClientState(GL_VERTEX_ARRAY);
#    oGL.glDrawElements(GL_TRIANGLES,0,
#				   GL_UNSIGNED_INT,0);
#    oGL.glDisableClientState(GL_VERTEX_ARRAY);
    print "Draw Atoms"	
#    glLoadIdentity();		
#    oGL.glUseProgram(shaders.program);
#    glTranslatef(mol.allAtoms.coords[0][0],mol.allAtoms.coords[0][1],mol.allAtoms.coords[0][2]);
#    glutSolidSphere(2,20,20)
#    oGL.glUseProgram(0);
#    glutSwapBuffers()
#    // draw atoms and bonds
#    shaders.drawAtoms();
    shaders.drawBonds_vbo();
    shaders.drawAtoms_vbo()
#    glPopMatrix();
#    glFlush()    
    glutSwapBuffers()
#    glutPostRedisplay()
    
def onMouseButton(button, state, x, y):
    global oldX, oldY
    global isRotating, isZooming, isClipping
    global mustDoHQ
    oldX, oldY = x, y
    if (button == GLUT_LEFT_BUTTON):
        if (state == GLUT_DOWN):
            mustDoHQ = False
            isRotating = True
        elif (state == GLUT_UP):
#            mol.orien = glTrackball.quat * mol.orien
            glTrackball.reset()
            mustDoHQ = True
            isRotating = not isRotating
            glutPostRedisplay()
    elif (button == GLUT_RIGHT_BUTTON):
        keys = glutGetModifiers()
        if (keys & GLUT_ACTIVE_SHIFT):
            if (state == GLUT_DOWN):
                isClipping = True
            elif (state == GLUT_UP):
                isClipping = not isClipping
        else:
            if (state == GLUT_DOWN):
                isZooming = True
            elif (state == GLUT_UP):
                isZooming = not isZooming

isRotating = False
isZooming = False
isClipping = False
def onMouseDrag(x, y):
    print "drag",x,y
    global clipplane, oldY, oldX
    if isRotating:
        glTrackball.update(oldX, oldY, x, y, width, height)
    if isZooming:
        ydiff = y - oldY
        oldX, oldY = x, y
#        mol.scaleFactor += -0.1*ydiff*mol.scaleFactor
#        if mol.scaleFactor < 0.1: mol.scaleFactor = 0.1
    elif isClipping:
        ydiff = oldY - y
        oldX, oldY = x, y
        clipplane[1] = -1
        #clipplane[0] = -1
        clipplane += [0, 0, 0, 0.1*ydiff]
    glutPostRedisplay()


def printHelp():
    print '''\
Welcomd to pyQutemol

Left mouse button to rotate the system
Right mousebutton to scale the system
Right mouse button + Shift key to move the clipping plane along y axis

Keyboard options:

q or Esc      - quit pyQutemol
h             - help
r             - run through trajectory
j             - jump to a certain step
n             - next step in trajectory (includes averaging)
p             - previous step in trajectory
+             - increase the number of frames averaged together (doesn't account for periodicity)
-             - decrease the number of frames averaged together (minimum of 1)
g             - next visualization
G             - previous visualization
k             - change color for selection (must be in hex format)
o             - redo ambient occlusion shading
a             - toggle axes
x, y, or z    - align viewpoint along respective axes
s             - save snapshot to shapshot.png in local directory
m             - make movie (requires ffmpeg) - edit makeMovie() to script the movie
c             - change primary selection
e             - change excluded selection (shown even with the clipping plane)

The selections only work if you have loaded a psf/dcd combination
'''

run_trj = False
draw_axes = True
def keyfunc_mol(mol):
    shadowmap = None
    shader_i = [0]
    def keyfunc(k, x, y):
        global run_trj, draw_axes
        if k == "q" or ord(k) == 27: # Escape
            sys.exit(0)
        elif k == "h":
            printHelp()
        elif k == "i":
            ipshell()
        elif k == "n":
            mol.read_next_frame()
            glutPostRedisplay()
        elif k == "p":
            mol.read_previous_frame()
            glutPostRedisplay()
        elif k == "j":
            print "Current frame %d, Total frames %d, select frame:"%(mol.universe.dcd.ts.frame, mol.universe.dcd.numframes)
            selection = raw_input("> ")
            try:
                frameno = int(selection)
                mol.universe.dcd[frameno]
                glutPostRedisplay()
            except:
                print "Invalid frame"
        elif k == "+":
            mol.averaging += 1
        elif k == "-":
            mol.averaging -= 1
            if mol.averaging < 1: mol.averaging = 1
        elif k == "r":
            run_trj = not run_trj
        elif k == "v":
            mol.PrepareAOSingleView(shadowmap)
            glutPostRedisplay()
        elif k == "a":
            draw_axes = not draw_axes
            glutPostRedisplay()
        elif k == "s":
            saveSnapshot(mainCanvas.GetHardRes()*2, mol)
        elif k == "g":
            shader_i[0] += 1
            if shader_i[0] == len(shaders): shader_i[0] = 0
            shaders[shader_i[0]].set(cgSettings)
            cgSettings.UpdateShaders()
            glutPostRedisplay()
        elif k == "G":
            shader_i[0] -= 1
            if shader_i[0] == -1: shader_i[0] = len(shaders)-1
            shaders[shader_i[0]].set(cgSettings)
            cgSettings.UpdateShaders()
            glutPostRedisplay()
        elif k == "o":
            mol.ResetAO()
            glutPostRedisplay()
        elif k == "m":
            makeMovie(mol)
        elif k == "k":
            print "Make selection for color change:"
            selection = raw_input("> ")
            try:
                sel = mol.universe.selectAtoms(selection)
                idx = sel.indices()
                print "input color:"
                color = raw_input("> ")
                mol.colors[idx] = convert_color(int(color,0))
                glutPostRedisplay()
            except:
                print "Invalid selection"
        elif k == "e":
            print "Make selection for exclusion:"
            selection = raw_input("> ")
            try:
                sel = mol.universe.selectAtoms(selection)
                mol.excl = sel.indices()
                mol.ResetAO()
                glutPostRedisplay()
            except:
                print "Invalid selection"
        elif k == "c":
            print "Make new selection:"
            selection = raw_input("> ")
            try:
                mol.sel = mol.universe.selectAtoms(selection)
                mol.pos = mol.sel.centerOfGeometry()
                coor = mol.sel.coordinates()
                min, max = numpy.minimum.reduce(coor), numpy.maximum.reduce(coor)
                mol.r = 0.5*numpy.sqrt(numpy.sum(numpy.power(max-min-4,2)))
                mol.min, mol.max = min, max
                mol.idx = mol.sel.indices()
                mol.ResetAO()
                glutPostRedisplay()
            except:
                print "Invalid selection"
        elif k == "x":
            mol.orien *= 0
            v = numpy.sin(numpy.pi/4.)
            mol.orien = quaternion([-.5,-.5,0.5,0.5])
            glTrackball.reset()
            glutPostRedisplay()
        elif k == "y":
            mol.orien *= 0
            v = numpy.sin(numpy.pi/4.)
            mol.orien.array += [-v,-v,0,0]
            glTrackball.reset()
            glutPostRedisplay()
        elif k == "z":
            mol.orien *= 0
            mol.orien.array[2] = -1
            glTrackball.reset()
            glutPostRedisplay()
    return keyfunc


def idlefunc_mol(mol):
    def idlefunc():
#        glutSetWindow(win_id);
        global run_trj
        if run_trj:
            mol.read_next_frame()
            glutPostRedisplay()
    return idlefunc

def ReshapeFunc( width,  height):
    print "ReshapeFunc",win_id
#    glutSetWindow(win_id);
    glutReshapeWindow(width, height);  
    glMatrixMode(GL_PROJECTION);     
    glLoadIdentity();                
    gluPerspective(70.0,float(width)/float(height),0.1,500.0);
    glMatrixMode(GL_MODELVIEW);     
    glViewport(0, 0, width, height);
    glMatrixMode(GL_MODELVIEW);  
    glLoadIdentity();

def OpenGlutWindow(posx,posy,width,height,mol,shaders):
    global xCenter,yCenter,zCenter
    glutInitDisplayMode(GLUT_DOUBLE |  GLUT_RGB | GLUT_DEPTH )

    glutInitWindowPosition(posx, posy);
    glutInitWindowSize( width, height)
    win_id = glutCreateWindow("Hyperball Demo");
    print win_id
	 #//    InitDisplay(width, height);
    initDisplay();

    print "after initGl"
    shaders.initGL();
    print "after shader init GL"
#    shaders.updatePositions(mol.allAtoms.coords[:2])#-np.array([xCenter,yCenter,zCenter]));
    xCenter,yCenter,zCenter = [0,0,0]
    print "ok update position"
    
    dispfunc = lambda: drawScene(mol,shaders)
    print "after drawScene"
    glutKeyboardFunc(keyfunc_mol(mol))
    print "dsFunc"
    glutDisplayFunc(dispfunc)
    glutMouseFunc(onMouseButton)
    glutMotionFunc(onMouseDrag)
    glutReshapeFunc(ReshapeFunc)
    print "dilFunc"
    glutIdleFunc(idlefunc_mol(mol))
    print "mainloop glutMainLoop()"    
    

if __name__=="__main__":
    from pyHyperballs import gl_config
    gl_config.importGL()
    
    print "OK"
    global pmv
    # Get around a bug in GLUT on OS X
    cwd = os.getcwd()
    os.chdir(cwd)
    
    if len(sys.argv) != 4:
        print "Usage: %s [prefix] [is_trj] [is_coarsegrain]"%sys.argv[0]
        print "If viewing a trajectory, prefix should be the name of the psf/dcd combo without the extension"
        print "otherwise the name of the pdb file"
        sys.exit(0)

    prefix = sys.argv[1]
    istrj = sys.argv[2]
    iscoarse = sys.argv[3]

    istrj = (int(istrj) == 1)
    iscoarse = (int(iscoarse) == 1)
    print "init"
    
    #    print "beforeMol"
#    mol = Molecule(prefix, istrj, iscoarse)    
    if not istrj:
        print "init PMV"
#        customizer = ePMV.__path__[0]+os.sep+"epmvrc.py"
        pmv = MoleculeViewer(logMode = 'overwrite', #customizer=customizer, 
                                master=None,title='pmv', withShell= 0,
                                verbose=False, gui = False)
    
        if len(prefix) == 4 :
            pmv.fetch(prefix)
            pmvmol = pmv.Mols[-1]
        else :
            pmvmol = pmv.readMolecule(prefix)
            pmvmol = pmv.Mols[-1]
        pmv.displayCPK(pmvmol)
        pmv.buildBondsByDistance(pmvmol)
        pmv.colorByAtomType(pmvmol,["cpk"])
#        pmv.colorByChains(pmvmol,["cpk"])
        mol = pmv.Mols[0]#Molecule(prefix, istrj, iscoarse,pmvmol=pmvmol)
        mol.orien = quaternion([-.5,-.5,0.5,0.5])
        #mol.updateColors()
    else :
        print "beforeMol"
        mol = pmv.Mols[0]#Molecule(prefix, istrj, iscoarse)
    print "afterMol"
    xCenter,yCenter,zCenter = mol.getCenter()
    # Set the clipping plane
    glutInit(sys.argv)    
    from pyHyperballs.hyperballs.AtomAndBondGLSL import AtomAndBondGLSL 
    from pyHyperballs.hyperballs.ballimproved_frag import ballimproved_frag, ball_frag
    from pyHyperballs.hyperballs.ballimproved_vert import ballimproved_vert, ball_vert
    from pyHyperballs.hyperballs.stickimproved_frag import stickimproved_frag 
    from pyHyperballs.hyperballs.stickimproved_vert import stickimproved_vert 
    shaders = AtomAndBondGLSL()
    print shaders
    shaders.setAtomFragmentShaderProgramSource(ballimproved_frag);
    shaders.setAtomVertexShaderProgramSource(ballimproved_vert);
    shaders.setBondFragmentShaderProgramSource(stickimproved_frag);
    shaders.setBondVertexShaderProgramSource(stickimproved_vert);

    nbAtoms=len(mol.allAtoms)
    bonds, atnobnd = mol.allAtoms.bonds
    indices = map(lambda x: [x.atom1._bndIndex_,
                             x.atom2._bndIndex_], bonds)
    nbBonds = len(indices)
    shaders.setupBuffersAndTextures_mol(nbAtoms,nbBonds,mol.allAtoms,indices);
    print "OpenGlutWindow"
    OpenGlutWindow(150,10,1024,1024,mol,shaders)
    print "glutMainLoop()"
#    glutMainLoop()
    bBox = mol.geomContainer.geoms['cpk'].ComputeBB("test")
    xLength=bBox[1][0]-bBox[0][0]
    yLength=bBox[1][1]-bBox[0][1]
    zLength=bBox[1][2]-bBox[0][2]
    print bBox,xCenter,yCenter,zCenter
    gl_config.printGLInfo()
#    glutMainLoop()
#pmv.colorByResidueType(pmvmol,["cpk"]);mol.updateColors();drawFrame(mol)
#pmv.colorAtomsUsingDG(pmvmol,["cpk"]);mol.updateColors();drawFrame(mol)
#pmv.colorByChains(pmvmol,["cpk"]);mol.updateColors();drawFrame(
