#test.py - an example program, (rewritten in python from the GLE distribution demo).
from opengltk.OpenGL import GL, GLU, GLUT
import gle
import numpy
import sys

lastx=0
lasty=0
NPTS = 6
points =  numpy.zeros((NPTS, 3), "f")
colors = numpy.zeros((NPTS,3), "f")
idx = 0

def PNT(x,y,z):
    global idx
    points[idx][0] = x
    points[idx][1] = y
    points[idx][2] = z
    idx = idx+1		


def COL(r,g,b):
    global idx
    colors[idx][0] = r
    colors[idx][1] = g
    colors[idx][2] = b


# * Initialize a bent shape with three segments. 
# * The data format is a polyline.
# *
# * NOTE that neither the first, nor the last segment are drawn.
# * The first & last segment serve only to determine that angle 
# * at which the endcaps are drawn.


def InitStuff ():

    #  initialize the join style here 
    gle.gleSetJoinStyle (gle.TUBE_NORM_EDGE | gle.TUBE_JN_ANGLE | gle.TUBE_JN_CAP)
    
    COL (0.0, 0.0, 0.0)
    PNT (-6.0, 6.0, 0.0)
    
    COL (0.0, 0.8, 0.3)
    PNT (6.0, 6.0, 0.0)
    
    COL (0.8, 0.3, 0.0)
    PNT (6.0, -6.0, 0.0)
    
    COL (0.2, 0.3, 0.9)
    PNT (-6.0, -6.0, 0.0)
    
    COL (0.2, 0.8, 0.5)
    PNT (-6.0, 6.0, 0.0)
    
    COL (0.0, 0.0, 0.0)
    PNT (6.0, 6.0, 0.0)
    

#draw the cylinder shape 
def DrawStuff ():

    GL.glClear (GL.GL_COLOR_BUFFER_BIT | GL.GL_DEPTH_BUFFER_BIT)
    
    # set up some matrices so that the object spins with the mouse 
    GL.glPushMatrix ()
    GL.glTranslatef (0.0, 0.0, -80.0)
    GL.glRotatef (lastx, 0.0, 1.0, 0.0)
    GL.glRotatef (lasty, 1.0, 0.0, 0.0)
    
    # Phew. FINALLY, Draw the polycylinder  -- 
    #gle.glePolyCylinder (NPTS, points, colors, 2.3)
    gle.glePolyCylinder (points, colors, 2.3)
    GL.glPopMatrix ()
    
    GLUT.glutSwapBuffers ()

# get notified of mouse motions 

def MouseMotion (x, y):
    
    lastx = x
    lasty = y
    GLUT.glutPostRedisplay ()

def keyboard( key, x, y):
    if 'e' == chr( key):
        raise RuntimeError
    elif 'q' == chr( key):
        sys.exit()

def JoinStyle (msg): 
   # get the current joint style 
   style = gle.gleGetJoinStyle ()

   # there are four different join styles, 
   # * and two different normal vector styles 
   normalStyle = gle.TUBE_NORM_PATH_EDGE
   joinStyle = gle.TUBE_JN_ANGLE
   gle.gleSetJoinStyle ( normalStyle | joinStyle )
   GLUT.glutPostRedisplay ()

# set up a light 
lightOnePosition = [40.0, 40, 100.0, 0.0]
lightOneColor = [0.99, 0.99, 0.99, 1.0] 

lightTwoPosition = [-40.0, 40, 100.0, 0.0]
lightTwoColor = [0.99, 0.99, 0.99, 1.0] 

###    initialize glut 
GLUT.glutInit (sys.argv)
GLUT.glutInitDisplayMode (GLUT.GLUT_DOUBLE | GLUT.GLUT_RGB | GLUT.GLUT_DEPTH)
GLUT.glutCreateWindow ("join styles")
GLUT.glutDisplayFunc (DrawStuff)
#GLUT.glutMotionFunc (MouseMotion)


#    initialize GL 
GL.glClearDepth (1.0)
GL.glEnable (GL.GL_DEPTH_TEST)
GL.glClearColor (0.0, 0.0, 0.0, 0.0)
GL.glShadeModel (GL.GL_SMOOTH)

GL.glMatrixMode (GL.GL_PROJECTION)
# roughly, measured in centimeters 
GL.glFrustum (-9.0, 9.0, -9.0, 9.0, 50.0, 150.0)
GL.glMatrixMode(GL.GL_MODELVIEW)

# initialize lighting 
GL.glLightfv (GL.GL_LIGHT0, GL.GL_POSITION, lightOnePosition)
GL.glLightfv (GL.GL_LIGHT0, GL.GL_DIFFUSE, lightOneColor)
GL.glEnable (GL.GL_LIGHT0)
GL.glLightfv (GL.GL_LIGHT1, GL.GL_POSITION, lightTwoPosition)
GL.glLightfv (GL.GL_LIGHT1, GL.GL_DIFFUSE, lightTwoColor)
GL.glEnable (GL.GL_LIGHT1)
GL.glEnable (GL.GL_LIGHTING)
GL.glEnable (GL.GL_NORMALIZE)
GL.glColorMaterial (GL.GL_FRONT_AND_BACK, GL.GL_DIFFUSE)
GL.glEnable (GL.GL_COLOR_MATERIAL)

InitStuff ()
GLUT.glutKeyboardFunc(keyboard)
print "Type 'q' in the demo window to quit"
GLUT.glutMainLoop ()


