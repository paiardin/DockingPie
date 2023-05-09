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

import Tkinter, os
import numpy
from PIL import Image
from opengltk.OpenGL import GL
from opengltk.extent import _gllib
from DejaVu2.EventHandler import EventManager
from DejaVu2.MaterialEditor import OGLWidget

import warnings
from warnings import warn

## TODO
##
## try using glPixelMap to clip and change the outlines
## try using FL_FLOAT format for images to have more resolution in Z
## try glPicelTransfer to scale bias contour
## try a bluring filter to antialias outline
## detect imaging subset presence
##

class ImageViewer(OGLWidget):


    def __init__(self, master=None, image=None, name='ImageViewer',
                 cnf = {}, **kw):

        if not kw.has_key('double'):
            kw['double'] = 1
	if not kw.has_key('depth'):
            kw['depth'] = 0
	if not kw.has_key('accum'):
            kw['accum'] = 0

        if master is None:
            master = Tkinter.Toplevel()
            master.title(name)

        self.width = 100
        self.height = 100
        self.imarray = None
##         self.scale= 1.0
##         self.bias = 0.0
        
        apply( OGLWidget.__init__, (self, master, cnf), kw)

        self.initGL()

        self.setImage(image)


    def Exit(self):
        self.master.destroy()

        
    def Configure(self, event=None):
        """Cause the opengl widget to redraw itself."""
        #print 'Configure 0'
        if self.imarray is None:
            return
        if isinstance(self.master, Tkinter.Tk) or \
           isinstance(self.master, Tkinter.Toplevel):
               geom = '%dx%d' % (self.width, self.height)
               self.master.geometry(geom)
        self.unbind('<Configure>')
        self.configure(width=self.width, height=self.height)
        self.bind('<Configure>', self.Configure)
        self.initProjection()
        

    def initProjection(self):
        if self.imarray is None:
            return
        self.tk.call(self._w, 'makecurrent')
        GL.glViewport(0, 0, self.width, self.height)
        GL.glMatrixMode(GL.GL_PROJECTION)
        GL.glLoadIdentity()
        GL.glOrtho(0, float(self.width), 0, float(self.height), -1.0, 1.0)
        GL.glMatrixMode(GL.GL_MODELVIEW)


    def setImage(self, image, width=None, height=None, mode=None):
        self.width = None
        self.height = None
        self.imarray = None
        if image is None:
            return
        elif isinstance(image, Image.Image):
            self.width = image.size[0]
            self.height = image.size[1]
            im = image.transpose(Image.FLIP_TOP_BOTTOM)
            self.imarray = numpy.fromstring(im.tostring(), 'B')
            self.mode = image.mode
        elif isinstance(image, numpy.ndarray):
            if mode == 'RGB':
                lenmode = 3
            elif mode in ['L','P']:
                lenmode = 1
            assert image.dtype.char=='b'
            if len(image.shape)==3:
                self.width = image.shape[0]
                self.height = image.shape[1]
                self.imarray = numpy.reshape( image, (-1,))
            elif len(image.shape)==1:
                self.width = width
                self.height = height
                self.imarray = image
                self.imarray.shape = (self.width, self,height)
                self.numimarray = numpy.array(self.imarray)
                self.imarray = numpy.reshape( image, (-1,))
            else:
                raise RuntimeError, "bad shape for image array"
            assert len(self.imarray)==self.width*self.height*lenmode
            self.mode = mode
        else:
            print 'Not surported yet'
        
        self.Configure()
        self.tkRedraw()
        

    def initGL(self):
        GL.glPixelStorei(GL.GL_UNPACK_ALIGNMENT, 1)
        GL.glClearColor(0.4, 0.4, 0.4, 0.0)


    def redraw(self, event=None, filter=None):
        if self.imarray is None:
            return
        if filter:
            self.master.lift()
            GL.glConvolutionFilter2D(GL.GL_CONVOLUTION_2D, GL.GL_LUMINANCE,
                                     3, 3, GL.GL_LUMINANCE, GL.GL_FLOAT,
                                     filter)
##             GL.glConvolutionParameterfv(GL.GL_CONVOLUTION_2D,
##                                        GL.GL_CONVOLUTION_FILTER_SCALE,
##                                        (3., 3.,3.,3.))
##             s=  self.scale
##             GL.glPixelTransferf(GL.GL_POST_CONVOLUTION_RED_SCALE, s)
##             GL.glPixelTransferf(GL.GL_POST_CONVOLUTION_BLUE_SCALE, s)
##             GL.glPixelTransferf(GL.GL_POST_CONVOLUTION_GREEN_SCALE, s)
##             s = self.bias
##             GL.glPixelTransferf(GL.GL_POST_CONVOLUTION_RED_BIAS, s)
##             GL.glPixelTransferf(GL.GL_POST_CONVOLUTION_BLUE_BIAS, s)
##             GL.glPixelTransferf(GL.GL_POST_CONVOLUTION_GREEN_BIAS, s)
##                                        GL.GL_CONVOLUTION_FILTER_SCALE,
##                                        (3., 3.,3.,3.))
##             GL.glConvolutionParameteriv(GL.GL_CONVOLUTION_2D,
##                                        GL.GL_CONVOLUTION_FILTER_BIAS,
##                                        (1500, 1500, 1500, 1500))
            GL.glEnable(GL.GL_CONVOLUTION_2D)
            
        self.tk.call(self._w, 'makecurrent')
        GL.glClearColor(0.0, 0.0, 0.0, 0.0)
        GL.glClear(GL.GL_COLOR_BUFFER_BIT)
        GL.glRasterPos2i( 0, 0)
        GL.glPixelStorei(GL.GL_UNPACK_ALIGNMENT, 1)
        if self.mode=='RGB':
            _gllib.glDrawPixels( self.width, self.height,
                                 GL.GL_RGB, GL.GL_UNSIGNED_BYTE, self.imarray)
        elif self.mode in ['L','P']:
            _gllib.glDrawPixels( self.width, self.height,
                                 GL.GL_LUMINANCE, GL.GL_UNSIGNED_BYTE,
                                 self.imarray)

        GL.glDisable(GL.GL_CONVOLUTION_2D)


    def smooth(self):
        average = (numpy.ones( (3,3) , 'f')/9.).astype('f')
#        average[1][1] = .0
        print average
        self.redraw(filter=average)
        image = (self.imageAsArray()).astype('B')
        return image


    def secondDerivativeNum(self):
        sndDeriv = numpy.array([ [-0.125, -0.125, -0.125,],
                                    [-0.125,    1.0, -0.125,],
                                    [-0.125, -0.125, -0.125,] ], 'f')
        a = self.numimarray
        c = a.copy()
        numpy.convolve.Convolve2d(sndDeriv, a, c)
        return c


    def firstDerivative(self):
        if self.imarray is None:
            return
        fstDeriveV = numpy.array([ [-0.125,  -0.25, -0.125],
                                      [ 0.0  ,    0.0,  0.0  ],
                                      [ 0.125,   0.25,  0.125] ], 'f')

        a = numpy.array(self.imarray)
        c = a.copy()
        numpy.convolve.Convolve2d(fstDeriveV, a, c)
        
        fstDeriveV = numpy.array([ [ 0.125,   0.25,  0.125],
                                     [ 0.0  ,    0.0,  0.0  ],
                                     [-0.125,  -0.25, -0.125] ], 'f')
        d = a.copy()
        numpy.convolve.Convolve2d(fstDeriveV, a, d)

        fstDeriveH = numpy.array([ [-0.125,    0.0, 0.125],
                                     [-0.25 ,    0.0, 0.25  ],
                                     [-0.125,    0.0, 0.125] ], 'f')
        e = a.copy()
        numpy.convolve.Convolve2d(fstDeriveH, a, e)

        fstDeriveH = numpy.array([ [ 0.125,    0.0, -0.125],
                                     [ 0.25 ,    0.0, -0.25  ],
                                     [ 0.125,    0.0, -0.125] ], 'f')
        f = a.copy()
        numpy.convolve.Convolve2d(fstDeriveH, a, f)

        deriv = numpy.fabs(c+d*0.5)+numpy.fabs(e+f*0.5)

        return deriv.astype('B')


    def secondDerivative(self):
        sndDeriv = numpy.array([ [-0.125, -0.125, -0.125,],
                                   [-0.125,    1.0, -0.125,],
                                   [-0.125, -0.125, -0.125,] ], 'f')
        self.redraw(filter=sndDeriv)
        deriv = numpy.fabs(self.imageAsArray()).astype('B')
        return deriv
    

    def firstDerivative(self):
        fstDeriveV = numpy.array([ [-0.125,  -0.25, -0.125],
                                     [ 0.0  ,    0.0,  0.0  ],
                                     [ 0.125,   0.25,  0.125] ], 'f')

        self.redraw(filter=fstDeriveV)
        derivV = self.imageAsArray()
        if derivV is None:
            return None
        fstDeriveV = numpy.array([ [ 0.125,   0.25,  0.125],
                                     [ 0.0  ,    0.0,  0.0  ],
                                     [-0.125,  -0.25, -0.125] ], 'f')

        self.redraw(filter=fstDeriveV)
        derivV += self.imageAsArray()

        fstDeriveH = numpy.array([ [-0.125,    0.0, 0.125],
                                     [-0.25 ,    0.0, 0.25  ],
                                     [-0.125,    0.0, 0.125] ], 'f')
        self.redraw(filter=fstDeriveH)
        derivH = self.imageAsArray()

        fstDeriveH = numpy.array([ [ 0.125,    0.0, -0.125],
                                     [ 0.25 ,    0.0, -0.25  ],
                                     [ 0.125,    0.0, -0.125] ], 'f')
        self.redraw(filter=fstDeriveH)
        derivH += self.imageAsArray()

        deriv = numpy.fabs(derivH*0.5)+numpy.fabs(derivV*0.5)

        return deriv.astype('B')


    def getHistogram(self):
        GL.glHistogram(GL.GL_HISTOGRAM, 256, GL.GL_RGB, GL.GL_FALSE)
        GL.glEnable(GL.GL_HISTOGRAM)
        self.redraw()
        values = numpy.zeros( (256,3), numpy.uint16)
        #seg faults
        GL.glGetHistogram(GL.GL_HISTOGRAM, GL.GL_TRUE, GL.GL_RGB,
                          GL.GL_UNSIGNED_SHORT, values)
        return values

    
    def imageAsArray(self):
        width = self.width
        height = self.height
        if width is None or height is None:
            return None
        nar = numpy.zeros(3*width*height, numpy.uint8)

        self.tk.call(self._w, 'makecurrent')
        GL.glPixelStorei(GL.GL_PACK_ALIGNMENT, 1)

#        self.tk.call(self._w, 'swapbuffers')
        GL.glReadBuffer(GL.GL_BACK)
        _gllib.glReadPixels( 0, 0, width, height, GL.GL_RGB,
                             GL.GL_UNSIGNED_BYTE, nar)
#        self.tk.call(self._w, 'swapbuffers')
        return nar


if __name__=='__main__':
    im = Image.open('lena.jpg')
    vi = ImageViewer(image=im, name='lena.jpg')
    vi.redraw()
    im2 = Image.open('cAMP_mac.jpg')
    im3 = Image.open('spectrum.jpg')
    vi.setImage(im3)
    #vi.setImage(im2)
