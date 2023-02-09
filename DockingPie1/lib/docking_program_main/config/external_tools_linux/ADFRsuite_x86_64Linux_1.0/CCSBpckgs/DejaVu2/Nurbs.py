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

from opengltk.OpenGL import GL, GLU
from opengltk.extent import _glulib

import numpy

import DejaVu2
from DejaVu2.Geom import Geom

class GLUSurfaceNurbs(Geom):
    """Class for GLU Nurbs surfaces"""


    keywords = Geom.keywords + [
        'knots',
        'ctrlPts',
        'displayMode',
        ]

    def getState(self, full=False):
        state = Geom.getState(self, full)
        # add knots, ctrlPts, displayMode to state
##         state['quality'] = self.quality

##         if full:
##             rad = self.vertexSet.radii.array
##             if len(rad):
##                 state['radii'] =  rad
            
        return state


    def init_surface(self):
        self.ctlpoints = numpy.zeros( (4,4,3), 'f')
        for u in range(4):
            for v in range(4):
                self.ctlpoints[u][v][0] = 2.0*(u - 1.5)
                self.ctlpoints[u][v][1] = 2.0*(v - 1.5)

                if ( u == 1 or u == 2) and (v == 1 or v == 2):
                    self.ctlpoints[u][v][2] = 3.0
                else:
                    self.ctlpoints[u][v][2] = -3.0
        

    def __init__(self, name=None, check=1, **kw):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
        #print "Spheres.__init__"

        self.knots = numpy.array((0.0, 0.0, 0.0, 0.0, .5, .5, .5, .5), 'f')
        apply( Geom.__init__, (self, name, check), kw )
        self.theNurb = None
        self.init_surface()
        self.immediateRendering = True
        
    def nurbsError(self, errorCode):
        estring = GLU.gluErrorString(errorCode);
        print "Nurbs Error: %s"%estring
        sys.exit(0)
        
    def Draw(self):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()

         GL.glEnable(GL.GL_AUTO_NORMAL)
         GL.glEnable(GL.GL_NORMALIZE)
         
         if self.theNurb is None:         
             self.theNurb = GLU.gluNewNurbsRenderer()

             GLU.gluNurbsProperty(self.theNurb, GLU.GLU_SAMPLING_TOLERANCE,
                                  25.0)
             GLU.gluNurbsProperty(self.theNurb, GLU.GLU_DISPLAY_MODE,
                                  GLU.GLU_OUTLINE_POLYGON)
             GLU.gluNurbsCallback(self.theNurb, GLU.GLU_ERROR, self.nurbsError)

         GLU.gluBeginSurface(self.theNurb)
         _glulib.gluNurbsSurface( self.theNurb,
                                  8, self.knots,
                                  8, self.knots,
                                  4*3,
                                  3,
                                  self.ctlpoints,
                                  4, 4,
                                  GL.GL_MAP2_VERTEX_3)

         GLU.gluEndSurface(self.theNurb)

         
         GL.glPointSize(5.0)
         GL.glDisable(GL.GL_LIGHTING)
         GL.glColor3f(1.0, 1.0, 0.0)
         GL.glBegin(GL.GL_POINTS)
         for i in range(4):
             for j in range(4):
                 GL.glVertex3fv(self.ctlpoints[i][j]) 

         GL.glEnd()
         GL.glDisable(GL.GL_AUTO_NORMAL)
         GL.glDisable(GL.GL_NORMALIZE)

if __name__=='__main__':
    from DejaVu2.Nurbs import GLUSurfaceNurbs

    from DejaVu2 import Viewer
    vi = Viewer()

    nurbsg = GLUSurfaceNurbs('nurbs')
    vi.AddObject(nurbsg)

    def move():
        from random import uniform
        for n in range(200):
            for i in range(1,3):
                for j in range(1,3):
                    nurbsg.ctlpoints[i][j][2] += uniform(-0.3, 0.3)
                    vi.OneRedraw()
                    vi.master.update()
    import numpy

    def move1():
        for n in range(1, 11):
            v = n*0.1
            nurbsg.knots = numpy.array( (0.0, 0.0, 0.0, 0.0, v, v, v, v), 'f')
            print nurbsg.knots
            vi.OneRedraw()
            vi.master.update()

    def move2():
        for n in range(0, 10):
            v = n*0.1
            nurbsg.knots = numpy.array( (v, v, v, v, 1, 1, 1, 1), 'f')
            print nurbsg.knots
            vi.OneRedraw()
            vi.master.update()
