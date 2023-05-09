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
# Date: Febuary 2006 Authors: Michel Sanner
#
#    vareille@scripps.edu
#    sanner@scripps.edu
#
#       The Scripps Research Institute (TSRI)
#       Molecular Graphics Lab
#       La Jolla, CA 92037, USA
#
# Copyright: Michel Sanner and TSRI
#
# Revision: Guillaume Vareille
#
#########################################################################
#
# $Header: /mnt/raid/services/cvs/DejaVu2/MaterialEditor.py,v 1.1.1.1.4.1 2017/07/13 22:28:32 annao Exp $
#
# $Id: MaterialEditor.py,v 1.1.1.1.4.1 2017/07/13 22:28:32 annao Exp $
#

import Tkinter, os, math
import numpy
import Pmw
import tkFileDialog

from opengltk.OpenGL import GL
from opengltk.extent.utillib import extractedGlutSolidSphere
from colorTool import ToHSV, ToRGB, OneColor, TkColor
from DejaVu2.extendedSlider import ExtendedSlider
from DejaVu2.Materials import Materials
import DejaVu2

from DejaVu2.Tk import loadTogl


class OGLWidget(Tkinter.Widget, Tkinter.Misc):
    """This class implements an OpenGL widget in a Tk Frame used to display 
color scales in the color editor.
"""

    def __init__(self, master, title, cnf={}, expand=1, **kw):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()

        if not kw.has_key('width'):
            kw['width']=150
        if not kw.has_key('height'):
            kw['height']=150
        if not kw.has_key('double'):
            kw['double']=1
        if not kw.has_key('depth'):
            kw['depth']=1

		# for unexplained reason anti-aliasing doesn't work
		# on mesa if we don't set the accumulation buffer
        if not kw.has_key('accum'):
            kw['accum'] = 1

        self.width = kw['width']
        self.height = kw['height']
        # load the TK-OpenGL extension (Togl)
        loadTogl(master)

        # create an Tk-OpenGL widget
        from opengltk.exception import GLerror
        from Tkinter import TclError
        try:
            Tkinter.Widget.__init__(self, master, 'togl', cnf, kw)
            try:
                GL.glAccum(GL.GL_LOAD, 1.0)
                self.accumBuffersError = False
            except GLerror:
                self.accumBuffersError = True
        except TclError, e:
            self.accumBuffersError = True
            kw.pop('accum')
            Tkinter.Widget.__init__(self, master, 'togl', cnf, kw)

        #currentcontext = self.tk.call(self._w, 'contexttag')
        #print "OGLWidget.__init__ currentcontext", currentcontext, title
        
        self.bind('<Expose>', self.tkExpose)
        self.bind('<Enter>', self.Enter_cb)
        self.bind('<Configure>', self.Configure)

        self.pack(side='left')


    def initProjection(self):
        """ guillaume: prepare the opengl matrices 
for drawing the widget 
(for instance the sliders in the case of the MaterialEditor) 
"""
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
        
        GL.glMatrixMode (GL.GL_PROJECTION)
        GL.glLoadIdentity ()
        GL.glOrtho(0., 1., 0., 1., -2.5, 2.5)
        GL.glMatrixMode(GL.GL_MODELVIEW)
        GL.glLoadIdentity()
        GL.glTranslatef(0, 0, 2.0)

        
    def tkExpose(self, *dummy):
        #guillaume: (commented out) redundant,this is allready in self.tkRedraw
        #self.tk.call(self._w, 'makecurrent') 
        #self.initProjection() 
        
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
        self.tkRedraw()
        

    def Activate(self):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()

        self.tk.call(self._w, 'makecurrent')
        

    def Enter_cb(self, event):
        """Call back function trigger when the mouse enters the camera"""
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()

        self.tk.call(self._w, 'makecurrent')


    def Configure(self, *dummy):
        """guillaume: set the opengl viewport"""
        # unbind <configure, because changing size of togl widget creates such
        # an event, hence causing an endless loop

        self.unbind('<Configure>')
        self.configure(width=self.width, height=self.height)
        self.bind('<Configure>', self.Configure)
        GL.glViewport(0, 0, self.width, self.height)
 
    def tkRedraw(self, *dummy):
        """Cause the opengl widget to redraw itself.
guillaume: call the preparation of the opengl matrices for drawing the widget
then draw the widget 
"""
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()

        #if not self.winfo_ismapped(): return
        self.update_idletasks()
        self.tk.call(self._w, 'makecurrent')
        self.initProjection()
        GL.glPushMatrix()
        self.redraw()
        GL.glFlush()
        GL.glPopMatrix()
        self.tk.call(self._w, 'swapbuffers')


    def setupLightModel(self):
        """
"""
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
        # this method has to be called explicitly by the derived classes if
        # a default lighting model is wanted

        GL.glLightfv(GL.GL_LIGHT0, GL.GL_AMBIENT,  [.5, .5, .5, 1.0])
        GL.glLightfv(GL.GL_LIGHT0, GL.GL_DIFFUSE,  [.5, .5, .5, 1.0])
        GL.glLightfv(GL.GL_LIGHT0, GL.GL_SPECULAR, [.5, .5, .5, 1.0])
        GL.glLightfv(GL.GL_LIGHT0, GL.GL_POSITION, [1.0, 1.0, 1.0, 0.0]);   

        GL.glLightfv(GL.GL_LIGHT1, GL.GL_AMBIENT,  [.5, .5, .5, 1.0])
        GL.glLightfv(GL.GL_LIGHT1, GL.GL_DIFFUSE,  [.5, .5, .5, 1.0])
        GL.glLightfv(GL.GL_LIGHT1, GL.GL_SPECULAR, [.5, .5, .5, 1.0])
        GL.glLightfv(GL.GL_LIGHT1, GL.GL_POSITION, [-1.0, 1.0, 1.0, 0.0]);   

        GL.glLightModelfv(GL.GL_LIGHT_MODEL_AMBIENT, [0.2, 0.2, 0.2, 1.0])
        #GL.glEnable(GL.GL_LIGHTING)
        if self.viewer is not None:
            self.viewer.enableOpenglLighting()
        GL.glEnable(GL.GL_LIGHT0)
        GL.glEnable(GL.GL_LIGHT1)



class MaterialEditor(OGLWidget):
    """The class implements the material editor widget allowing to modify the ambient, ....
The editor maintains a list of colors for each property in the .materials attribute.  
This list can be set using the defineMaterial(material) method. 
(FIXME passing a material object here ties the Material Editor to DejaVu2)
A material editor instance can be bound to an object to edit either its front 
or back facing polygons material. This is done using the setObject(object, face) method.  
The object bound to the matrial editor is stored in  the .object attribute.
To use a Material editor on objects other than a DejaVu2 Geometry 
the setObject() method should be overwritten.
the sphere drawn in the upper left corner shows the effect of changing 
the various components of  a material. Chekcbuttons allow to ...
etc ...
"""
    
    ambi = 0
    diff = 1
    emis = 2
    spec = 3
    shini = 4
    opac = 5


    def _destroy(self):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
        self.peFrame.destroy()
        self.destroy()
            
    
    def Configure(self, *dummy):
        """guillaume: set the OpenGL Viewport for the material editor"""
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
        
        #print 'Configure 1'
        width = self.winfo_width()
        height = self.winfo_height()
        if width > height:
          xo = (width-height)/2
          yo = 0
          size = height
        else:
          xo = 0
          yo = (height-width)/2
          size = width
        GL.glViewport(xo, yo, size, size)


    def initProjection(self):
        """ guillaume: prepare the opengl matrices to display the sphere"""
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
        
        GL.glMatrixMode (GL.GL_PROJECTION)
        GL.glLoadIdentity ()
        halfSize = self.halfSize
        GL.glOrtho(float(-halfSize), float(halfSize),
                   float(-halfSize), float(halfSize),
                   float(-halfSize),
                   float(halfSize))
        GL.glMatrixMode(GL.GL_MODELVIEW)
        GL.glLoadIdentity()
        GL.glTranslatef(0, 0, -2.0)


    def tkExpose(self, *dummy):
        """ guillaume: override the one in OGLWidget
the only addition is a call to OGLWidget.configure()
"""
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()

        #print 'Expose 1'
        self.tk.call(self._w, 'makecurrent')
        self.initProjection() 
        self.Configure()
        self.tkRedraw()


    def redraw(self):
        """ redraw the Material editor opengl sphere that shows the effect
of the modifications
"""
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()

        self.tk.call(self._w, 'makecurrent')
        GL.glClearColor(0,0,0,0)
        self.initProjection()
        self.Configure()
        GL.glClear(GL.GL_COLOR_BUFFER_BIT | GL.GL_DEPTH_BUFFER_BIT)

        """ draw the black background squares """
        GL.glDisable( GL.GL_LIGHTING )
        GL.glColor3f( 0.1, 0.1, 0.1 )
        GL.glBegin(GL.GL_QUADS)
        GL.glVertex3f(-2, 0, 2); GL.glVertex3f(-2, 2, 2)
        GL.glVertex3f( 0, 2, 2); GL.glVertex3f( 0, 0, 2)
        GL.glEnd()
        GL.glBegin(GL.GL_QUADS)
        GL.glVertex3f( 0,-2, 2); GL.glVertex3f( 0, 0, 2)
        GL.glVertex3f( 2, 0, 2); GL.glVertex3f( 2,-2, 2)
        GL.glEnd()

        """ draw the grey background squares """
        GL.glColor3f( 0.3, 0.3, 0.3 )
        GL.glBegin(GL.GL_QUADS)
        GL.glVertex3f(-2,-2, 2); GL.glVertex3f(-2, 0, 2)
        GL.glVertex3f( 0, 0, 2); GL.glVertex3f( 0,-2, 2)
        GL.glEnd()
        GL.glBegin(GL.GL_QUADS)
        GL.glVertex3f( 0, 0, 2); GL.glVertex3f( 0, 2, 2)
        GL.glVertex3f( 2, 2, 2); GL.glVertex3f( 2, 0, 2)
        GL.glEnd()

        """ enable the sphere transparancy """
        GL.glEnable(GL.GL_BLEND)
        GL.glDepthMask(GL.GL_FALSE)
        GL.glBlendFunc(GL.GL_SRC_ALPHA, GL.GL_ONE_MINUS_SRC_ALPHA)

        """ draw the sphere """
        #GL.glEnable(GL.GL_LIGHTING)
        if self.viewer is not None:
            self.viewer.enableOpenglLighting()
        self.setMaterial()
        extractedGlutSolidSphere(1.6, 30, 30, 0)


    def __init__(self, master, viewer, title='MaterialEditor',cnf={}, **kw):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()

        self.ownsMaster = False
        if master is None:
            master = Tkinter.Toplevel(viewer.master)#screenName=screenName)
            master.title(title)
            master.protocol('WM_DELETE_WINDOW', self.dismiss )
            self.ownsMaster = True

        # DO NOT USE master because OGLWidget.__init__ will define
        # self.master as a Frame
        self.root = master
        
        if not kw.has_key('width'):
            kw['width'] = 200

        if not kw.has_key('height'):
            kw['height'] = 200
            
        self.frame = Tkinter.Frame(self.root, borderwidth=3, relief='ridge')

        self.halfSize = 2
        self.material = [
            numpy.array([[0.4, 0.4, 0.4, 1.0]],'f'),
            numpy.array([[0.4, 0.4, 0.4, 1.0]],'f'),
            numpy.array([[1., 1., 1., 1.0]],'f'),
            numpy.array([[0., 0., 0., 1.0]],'f'),
            numpy.array([50],'f'),
            numpy.array([1.0],'f'),
        ]

        apply( OGLWidget.__init__, (self, self.frame, title, cnf), kw )
#        self.frame.pack(pady=5, padx=10)
        self.frame.grid(row=0, column=0, padx=5, pady=5)

        self.object = None
        self.viewer = viewer

        #self.redraw = self.redraw

        f = self.peFrame = Tkinter.Frame(self.root)
        self.ambipe = MaterialPropertyEditor(f, self, 'Ambi: ')
        self.ambipe.setRGB(self.material[self.ambi][0])
        self.ambipe.setCallback(self.setAmbient)
        
        self.diffpe = MaterialPropertyEditor(f, self, 'Diff: ')
        self.diffpe.setRGB(self.material[self.diff][0])
        self.diffpe.setCallback(self.setDiffuse)

        self.emispe = MaterialPropertyEditor(f, self, 'Emis: ')
        self.emispe.setRGB(self.material[self.emis][0])
        self.emispe.setCallback(self.setEmission)

        self.specpe = MaterialPropertyEditor(f, self, 'Spec: ')
        self.specpe.setRGB(self.material[self.spec][0])
        self.specpe.setCallback(self.setSpecular)

        self.opacitySl = ExtendedSlider(f, label='  Opacity', minval=0.0,
                   maxval=1.0, init=self.material[self.diff][0][3],
                   onButton=1, sd='left', withValue=0, immediate=1)
        self.opacitySl.AddCallback(self.opacity_cb)

        self.shiniSl = ExtendedSlider(f, label='Shininess', minval=0.0,
                   maxval=128.0, init=self.material[self.shini][0],
                   onButton=1, sd='left', withValue=0, immediate=1)
        self.shiniSl.AddCallback(self.shininess_cb)
        self.shiniSl.frame.pack(side='top', expand=1, fill='y')
        self.opacitySl.frame.pack(side='top', expand=1, fill='y')
        f.grid(row=1, column=0, padx=5, pady=5)
#        f.pack(pady=5)

        self.currentPalette = None
        self.paletteNum = 0
        self.lastDir = './'
        
        f = Tkinter.Frame(self.root)




#        # guillaume's work in progress
#        # target
#        f0 = Tkinter.Frame(f, borderwidth=3, relief = 'ridge')
#        self.targetOption = Pmw.OptionMenu(f0, 
#                                          label_text = 'target:',
#                                          labelpos = 'w',
#                                          initialitem = 'Front',
#                                          #command = self.setTarget,
#                                          items = ('Front', 'Back'),
#                                          menubutton_pady=0,
#                                         )
#        self.targetOption.grid(row=0, column=0, padx=2, pady=2)
#        f0.pack(padx=2, pady=2)



        f1 = Tkinter.Frame(f, borderwidth=3, relief = 'ridge')
        
        l = ['artdeco','autumn','glass','metal','neon','rococo','santafe',
             'sheen','silky','spring','summer','tropical','winter']

        self.paletteBt = Pmw.ComboBox(f1, label_text='Palettes:',
                                      labelpos='nw',
                                      entryfield_entry_width=10,
                                      scrolledlist_items=l,
                                      selectioncommand=self.selectPalette)
        self.paletteBt.grid(row=0, column=0, padx=2, pady=2)

        self.paletteNumBt = Pmw.Counter(f1, labelpos = 'nw',
                    label_text = 'Prop number:', entry_width = 2,
                    entryfield_value = self.paletteNum,
                    entryfield_validate = {'validator' : 'integer',
                                           'min' : 0, 'max' : 34} )
        self.paletteNumBt.grid(row=0, column=1, padx=2, pady=2)

        l = Tkinter.Button(f1, text="Load ",command=self.loadPaletteMaterial)
        l.grid(row=1, column=0, padx=2, pady=2, sticky='ew')

        l = Tkinter.Button(f1, text="Apply ", command=self.broadcast)
        l.grid(row=1, column=1, padx=2, pady=2, sticky='ew')

        f1.pack(padx=2, pady=2)

        f1 = Tkinter.Frame(f, borderwidth=3, relief = 'ridge')

        l = Tkinter.Button(f1, text="Read ...",command=self.readMaterial)
        l.grid(row=0, column=0, padx=2, pady=2, sticky='ew')

        l = Tkinter.Button(f1, text="Write ..", command=self.writeMaterial)
        l.grid(row=0, column=1, padx=2, pady=2, sticky='ew')
        
        f1.pack(padx=2, pady=2)

        f.grid(row=0, column=1)
        
        f = self.ccFrame = Tkinter.Frame(self.root)
        self.cc = ColorChooser(f, title=None)
        f.grid(row=1, column=1)
#        f.pack(pady=5)
        
        self.dismissTk = Tkinter.Button(self.root, text='DISMISS',
                                        command=self.dismiss)
        self.dismissTk.grid(row=2, column=0, columnspan=2, padx=5, pady=5)

        self.tk.call(self._w, 'makecurrent')
        
        self.initLight(self.viewer)


#    def setTarget(self, face):
#        obj = self.viewer.currentObject
#        self.setObject(obj, face)
#        delf.defineMaterial(obj.materials[face].prop, face)


    def readMaterial(self):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
        types = [('materials', '*_prop.py'), ('all files', '*.*')]
        title = "Read material definition"
        file = tkFileDialog.askopenfilename(
            filetypes=types, initialdir=self.lastDir, initialfile=None,
            title=title)
        
        if len(file)==0: return
        d = {}
        execfile(file, d)
        for k,v in d.items():
            if k=='__builtins__': continue
            break
        self.defineMaterial(v)
            

    def writeMaterial(self):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
        types = [('materials', '*_prop.py'), ('all files', '*.*')]
        title = "Save material definition"
        file = tkFileDialog.asksaveasfilename(
            filetypes=types, initialdir=self.lastDir,
            initialfile=None, title=title)
        if len(file)==0: return
        f = open(file, 'w')
        f.write("material = [\n")
        f.write(repr([self.material[self.ambi][0].tolist()])+', # ambient\n')
        f.write(repr([self.material[self.diff][0].tolist()])+', # diffuse\n')
        f.write(repr([self.material[self.emis][0].tolist()])+', # emission\n')
        f.write(repr([self.material[self.spec][0].tolist()])+', # specular\n')
        f.write(repr([self.material[self.shini][0].tolist()])+', # shininess\n')
        f.write(repr([self.material[self.opac][0].tolist()])+' # opacity\n')
        f.write("]\n")
        f.close()
        
        
    def loadPaletteMaterial(self, event=None):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
        if self.currentPalette is not None:
            num = int(self.paletteNumBt.get())
            mat = Materials(self.currentPalette, num)
            self.defineMaterial(mat.prop)
            
            
    def selectPalette(self, palette):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
        self.currentPalette = palette
        
        
    def initLight(self, viewer):
        """
"""        
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
        if viewer:
            viewer.lightModel.applyTo.append(self)
            viewer.lightModel.apply()
            for l in viewer.lights:
                if l.enabled:
                    self.tk.call(self._w, 'makecurrent')
                    l.apply()
        else:
            self.setupLightModel()
            

    def setObject(self, object, face):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
        self.object = object
        self.defineMaterial(object.materials[face].prop, face)
        
    def dismiss(self, event=None):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
        self.frame.master.withdraw()


    def show(self, event=None):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
        self.frame.master.deiconify()
      

    def broadcast(self):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
        if self.object is None: return
        mask = [ self.ambipe.onTk.get(), self.diffpe.onTk.get(),
                 self.emispe.onTk.get(), self.specpe.onTk.get(),
                 self.shiniSl.onTk.get(), self.opacitySl.onTk.get() ]
        material = {'ambient':self.material[self.ambi],
                    'diffuse':self.material[self.diff],
                    'emission':self.material[self.emis],
                    'specular': self.material[self.spec],
                    'shininess':self.material[self.shini],
                    'opacity':self.material[self.opac],
                    'binding':None}
        
        if self.objectFace == GL.GL_BACK:
            self.object.Set(rawMaterialB=material, matMask=mask,
                            transparent='implicit', redo=1)
        else:
            self.object.Set(rawMaterialF=material, matMask=mask,
                            transparent='implicit', redo=1)

        if self.viewer:
            self.tk.call(self.viewer.currentCamera._w, 'makecurrent')
            self.viewer.Redraw()

        
    def setAmbient(self, rgb):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
        #self.material[self.ambi][0][:3] = list(rgb[:3])
        self.material[self.ambi][0][:3] = numpy.array(rgb[:3],'f',copy=1)
        self.tkRedraw(self)
        self.broadcast()
        
    def setDiffuse(self, rgb):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
        #self.material[self.diff][0][:3] = list(rgb[:3])
        self.material[self.diff][0][:3] = numpy.array(rgb[:3],'f',copy=1)
        self.tkRedraw(self)
        self.broadcast()
        
    def setEmission(self, rgb):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
        #self.material[self.emis][0][:3] = list(rgb[:3])
        self.material[self.emis][0][:3] = numpy.array(rgb[:3],'f',copy=1)
        self.tkRedraw(self)
        self.broadcast()
        
    def setSpecular(self, rgb):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
        #self.material[self.spec][0][:3] = list(rgb[:3])
        self.material[self.spec][0][:3] = numpy.array(rgb[:3],'f',copy=1)
        self.tkRedraw(self)
        self.broadcast()


    def defineMaterial(self, material, face=GL.GL_FRONT):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
        # This method is used to set the material from an object to which
        # the material editor is bound
        self.objectFace = face
        #self.material[self.ambi] = material[self.ambi]
        self.material[self.ambi] = [material[self.ambi][0]]
        self.ambipe.setRGB(self.material[self.ambi][0])
        if len(material[self.ambi]) > 1:
            self.ambipe.disable()  # disable editing button
        self.ambipe.tkRedraw(self)
            
        #self.material[self.diff] = material[self.diff]
        self.material[self.diff] = [material[self.diff][0]]
        self.diffpe.setRGB(self.material[self.diff][0])
        if len(material[self.diff]) > 1:
            self.diffpe.disable()
        self.diffpe.tkRedraw(self)

        #self.material[self.emis] = material[self.emis]
        self.material[self.emis] = [material[self.emis][0]]
        self.emispe.setRGB(self.material[self.emis][0])
        if len(material[self.emis]) > 1:
            self.emispe.disable()
        self.emispe.tkRedraw(self)

        #self.material[self.spec] = material[self.spec]
        self.material[self.spec] = [material[self.spec][0]]
        self.specpe.setRGB(self.material[self.spec][0])
        if len(material[self.spec]) > 1:
            self.specpe.disable()
        self.specpe.tkRedraw(self)

        #self.material[self.shini] = material[self.shini]
        self.material[self.shini] = [material[self.shini][0]]
        self.shiniSl.set(self.material[self.shini][0], update=0)
        if len(material[self.shini]) > 1:
            self.shiniSl.disable()

        #self.material[self.opac] = material[self.opac]
        self.material[self.opac] = [material[self.opac][0]]
        self.opacitySl.set(self.material[self.opac][0], update=0)
        if len(material[self.opac]) > 1:
            self.opacitySl.disable()

        self.tkRedraw(self)

    
    def opacity_cb(self, val):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
        self.material[self.opac][0] = val
        self.tkRedraw(self)
        self.broadcast()


    def shininess_cb(self, val):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
        self.material[self.shini][0] = val
        self.tkRedraw(self)
        self.broadcast()


    def setMaterial(self):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
        mat = self.material
        GL.glMaterialfv(GL.GL_FRONT, GL.GL_AMBIENT, mat[self.ambi][0])
        diff = mat[self.diff][0]
        diff[3] = mat[self.opac][0]
        GL.glMaterialfv(GL.GL_FRONT, GL.GL_DIFFUSE, mat[self.diff][0])
        GL.glMaterialfv(GL.GL_FRONT, GL.GL_SPECULAR, mat[self.spec][0])
        GL.glMaterialfv(GL.GL_FRONT, GL.GL_EMISSION, mat[self.emis][0])
        GL.glMaterialf(GL.GL_FRONT, GL.GL_SHININESS, float(mat[self.shini][0]))

        GL.glColor4fv(mat[self.diff][0])



class MaterialPropertyEditor(OGLWidget):
    """This object provides a slider over a hue/saturation WISIWYG ramp"""


    def disable(self):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
        self.onTk.set(0)
        self.setVal.configure(state='disabled', bg='#AAAAAA')
        self.setcw.configure(state='disabled')
        self.edit.configure(state='disabled')
        self.unbind('<ButtonPress-1>')


    def enable(self):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
        self.onTk.set(1)
        self.setVal.configure(state='normal', bg='#CC9999')
        self.setcw.configure(state='normal')
        self.edit.configure(state='normal')
        self.bind('<ButtonPress-1>', self.mouse1Down)


    def on_cb(self, event=None):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
        if self.onTk.get(): self.enable()
        else: self.disable()

        
    def __init__(self, master, editor, label, title='MaterialPropertyEditor',
                 width=None, height=None, cnf={}, **kw):

        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
        if width is None: kw['width'] = 180
        if height is None: kw['height'] = 19
        self.width = kw['width']
        self.callback = None
        self.editor = editor

        tf = self.topFrame = Tkinter.Frame(master, borderwidth=3,
                                           relief='ridge')
        self.frame = Tkinter.Frame(tf, borderwidth=3, relief='sunken')
        
        apply( OGLWidget.__init__, (self, self.frame, title, cnf, 0), kw )

        self.forget()
        self.pack(side='left')
        self.valueTk = Tkinter.StringVar()
        self.setVal = Tkinter.Entry(self.frame, textvariable=self.valueTk,
                                   width=5, bg='#CC9999')
        self.setVal.bind('<Return>', self.setVal_cb)
        self.setVal.pack(side='left')

        f = Tkinter.Frame(tf)
        self.onTk = Tkinter.IntVar()
        cb = Tkinter.Checkbutton(f, text=label, variable=self.onTk,
                                 command=self.on_cb)
        self.onTk.set(1)
        self.editTk = Tkinter.IntVar()
        self.edit = Tkinter.Checkbutton(f, variable=self.editTk,
                                        text='Edit', command=self.edit_cb)
        self.setcwTk = Tkinter.IntVar()
        self.setcw = Tkinter.Button(f, text='SetCW', pady=0,
                                    command=self.setcw_cb)

        cb.pack(side='left', padx=3)
        self.edit.pack(side='left', padx=3)
        self.setcw.pack(side='left', padx=3)
        f.pack(side='top')
        self.frame.pack(side='bottom', pady=2)        
        tf.pack(padx=2, pady=2)
        self.bind('<ButtonPress-1>', self.mouse1Down)

        self.setRGB((1,1,1))


    def setVal_cb(self, event=None):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
        self.setValue(float(self.valueTk.get()))
        

    def setcw_cb(self, event=None):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
        self.editor.cc.Set(self.rgb)
        

    def edit_cb(self, event=None):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
        val = self.editTk.get()
        if val:
            self.editor.cc.callbacks.append(self.update)
        else:
            self.editor.cc.callbacks.remove(self.update)

            
    def setCallback(self, func):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
        assert callable(func)
        self.callback = func


    def callCallback(self):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
        # called when value changes
        if self.callback:
            self.callback(self.rgb)


    def update(self, rgb):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
        # call back called by colorwheel
        self.h, self.s, v = ToHSV(rgb[:3])
        self.rgbMax = ToRGB( (self.h, self.s, 1.0) )
        self.rgb = ToRGB( (self.h, self.s, self.v) )
        self.tkRedraw()
        self.callCallback()

       
    def setRGB(self, rgb):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
        # used to set the color for a given rgb tuple
        self.h, self.s, self.v = ToHSV(rgb[:3])
        self.valueTk.set(float('%.4f'%self.v))
        self.rgbMax = ToRGB( (self.h, self.s, 1.0) )
        self.rgb = rgb
        

    def setValue(self, value):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
        # call back for value (as v in hsv) slider
        if value < 0.0: value = 0.0
        if value > 1.0: value = 1.0
        self.v = value
        self.valueTk.set(float('%.4f'%self.v))
        self.rgb = ToRGB( (self.h, self.s, self.v) )
        self.tkRedraw()
        self.callCallback()


    def mouse1Down(self, event):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
        self.bind('<B1-Motion>', self.mouse1Move)
        self.bind('<ButtonRelease-1>', self.mouse1Up)
        self.setValue(float(event.x)/self.width)

        
    def mouse1Move(self, event):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
        self.setValue(float(event.x)/self.width)


    def mouse1Up(self, event):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
        self.setValue(float(event.x)/self.width)
        self.unbind('<B1-Motion>')
        self.unbind('<ButtonRelease-1>')

        
    def redraw(self):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
        self.tk.call(self._w, 'makecurrent')

        GL.glDisable( GL.GL_DEPTH_TEST )
        GL.glDisable( GL.GL_LIGHTING )
        #GL.glPolygonMode(GL.GL_FRONT, GL.GL_FILL)
        GL.glBegin(GL.GL_QUADS)
        GL.glColor3f(0.,0.,0.)
        GL.glVertex2f(0., 1.); GL.glVertex2f(0., 0.)
        GL.glColor3f(float(self.rgbMax[0]),float(self.rgbMax[1]),float(self.rgbMax[2]))
        GL.glVertex2f( 1., 0.); GL.glVertex2f( 1., 1.)
        GL.glEnd()

        GL.glEnable(GL.GL_COLOR_LOGIC_OP)
        GL.glLogicOp(GL.GL_XOR)
        GL.glLineWidth(2)
        GL.glColor3f(.5,.5,.5)
        GL.glBegin(GL.GL_LINES)
        x1 = self.v-0.01
        x2 = self.v+0.01
        GL.glVertex2f(float(x1), 1.); GL.glVertex2f(float(x1), 0.)
        GL.glVertex2f(float(x2), 0.); GL.glVertex2f(float(x2), 1.)
        GL.glEnd()
        GL.glDisable(GL.GL_COLOR_LOGIC_OP)



class ColorWheel:        

    def __init__(self, master, title=None):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()

        if not master:
            master = Tkinter.Toplevel()#screenName=screenName)

        if title is not None:
            master.title(title)

        f = self.frame = Tkinter.Frame(master)

        path = __import__('DejaVu2').__path__
        iconfile = os.path.join(path[0],'cw.ppm')
        self.cwim = Tkinter.PhotoImage(master=master, file=iconfile)
        self.width = self.cwim.width()
        self.height = self.cwim.height()
        self.cwcanvas = Tkinter.Canvas(f, width=self.width,
                       height=self.height, relief='sunken', borderwidth=3 )
        self.cwcanvas.create_image(3, 3, anchor=Tkinter.NW, image=self.cwim)
        self.cwcanvas.pack()
        self.frame.pack()

        self.callback = None
        self.immediate = 1
        self.x = 0
        self.y = 0
        self.radius = 55
        cx = self.cx = self.width/2 + 3
        cy = self.cy = self.height/2 + 3

        self.cursor = self.cwcanvas.create_line(
            cx-3, cy-3, cx-3, cy+3, cx+3,cy+3, cx+3, cy-3, cx-3, cy-3 )

        self.hsvColor = [1.,1.,1.]
        
        self.cwcanvas.bind('<ButtonPress-1>', self.mouse1Down)
#        self.cursor = self.cwcanvas.create_line(cx, cy, cx-55, cy)


    def _MoveCursor(self, x, y):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
	# find the saturation based on distance
	s = math.sqrt(x*x + y*y) / self.radius
	if s > 1.0:
	    x = x / s
	    y = y / s
	    s = 1.0

	# now find the hue based on the angle 
	if x or y:
	    angle = math.atan2(y, x)
	    if angle < 0.0:
                angle = angle + (2.*math.pi)
	    h = 1. - angle / (2.0 * math.pi)
	else:
	    h = 0

	# check if redraw and callback are needed
	if self.hsvColor[0] != h or self.hsvColor[1] != s:
	    self.hsvColor[0] = h
	    self.hsvColor[1] = s
            cx = self.cx+x
            cy = self.cy+y
            self.cwcanvas.coords( self.cursor, cx-3, cy-3, cx-3, cy+3,
                                  cx+3,cy+3, cx+3, cy-3, cx-3, cy-3 )

	    if self.immediate and self.callback:
                self.callback(self.Get('RGB'))


    def mouse1Down(self, event=None):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
        self.cwcanvas.bind('<B1-Motion>', self.mouse1Move)
        self.cwcanvas.bind('<ButtonRelease-1>', self.mouse1Up)
        self._MoveCursor(event.x - self.cx, event.y - self.cy)

        
    def mouse1Move(self, event=None):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
        self._MoveCursor(event.x - self.cx, event.y - self.cy)

        
    def mouse1Up(self, event=None):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
        self._MoveCursor(event.x - self.cx, event.y - self.cy)
        self.cwcanvas.unbind('<B1-Motion>')
        self.cwcanvas.unbind('<ButtonRelease-1>')


    def Get(self, mode='HSV'):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
	"""Get the current color"""
	if mode == 'RGB':
	    rgb = ToRGB(self.hsvColor)
	    return OneColor(rgb)
	elif mode == 'HSV':
	    return OneColor(self.hsvColor)
	elif mode == 'TkRGB':
	    col = numpy.array(ToRGB(self.hsvColor[:3]), 'f') * 255
	    return TkColor(col)


    def Set(self, color, mode='HSV'):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
	"""Set the current color"""
	assert len(color) in (3,4)
	color = OneColor(color)
	if mode=='RGB': color[:3] = ToHSV(color[:3])
	self.hsvColor = color[:3]

        # update cursor
	rad = self.hsvColor[1] * self.radius
	angle = 2.0 * math.pi * (1. - self.hsvColor[0])
	cx = self.cx + int(rad * math.cos(angle))
	cy = self.cy + int(rad * math.sin(angle))
        self.cwcanvas.coords( self.cursor, cx-3, cy-3, cx-3, cy+3,
                              cx+3,cy+3, cx+3, cy-3, cx-3, cy-3 )
        
	if self.immediate and self.callback:
            self.callback(self.Get('RGB'))


class ColorChooser:

    def __init__(self, master, title='Color Chooser'):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()

        if not master:
            master = Tkinter.Toplevel()#screenName=screenName)

        if title is not None:
            master.title(title)

        self.frame = Tkinter.Frame(master, borderwidth=3, relief='ridge')

        self.wheelFrame = Tkinter.Frame(self.frame)
        self.cw = ColorWheel(self.wheelFrame)
        self.wheelFrame.pack(pady=5, padx=10)

        f = self.chipFrame = Tkinter.Frame(self.frame, relief='ridge')
        self.chip1 = Tkinter.Frame(f, relief='sunken', borderwidth=3,
                                   width=60, height=30, bg='white')
        self.chip2 = Tkinter.Frame(f, relief='sunken', borderwidth=3,
                                   width=60, height=30, bg='white')
        self.chip1.pack(side='left', anchor='e')
        self.chip2.pack(side='right', anchor='w')
        f.pack(pady=2)

        f = self.btFrame = Tkinter.Frame(self.frame, relief='ridge')
        self.saveBt = Tkinter.Button(f, text='->', command=self.save)
        self.swapBt = Tkinter.Button(f, text='<->', command=self.swap)
        self.restoreBt = Tkinter.Button(f, text='<-', command=self.restore)

        self.saveBt.pack(side='left')
        self.swapBt.pack(side='left')
        self.restoreBt.pack(side='left')
        f.pack(pady=2)

        self.frame.pack(padx=5, pady=5)

        self.callbacks = []
        self.cw.callback = self.updateCurrent
        self.currentColor = [1., 1., 1.]
        self.savedColor = [1., 1., 1.]


    def Set(self, color):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
        self.updateCurrent(color)
        self.cw.Set(color, 'RGB')

        
    def updateCurrent(self, color):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
        self.chip1.configure( bg = TkColor( numpy.array(color, 'f')*255 ) )
        self.currentColor[:3] = list(color[:3])
        for f in self.callbacks:
            f(self.currentColor)
            

    def updateSaved(self, color):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
        self.chip2.configure( bg = TkColor( numpy.array(color, 'f')*255 ) )
        self.savedColor[:3] = list(color[:3])

        
    def save(self, event=None):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
        self.updateSaved(self.currentColor)


    def swap(self, event=None):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
        saved = self.savedColor[:]
        self.updateSaved(self.currentColor)
        self.cw.Set(saved, 'RGB')


    def restore(self, event=None):
        if __debug__:
         if hasattr(DejaVu2, 'functionName'): DejaVu2.functionName()
        self.cw.Set(self.savedColor, 'RGB')


if __name__=='__main__':
    import Tkinter
    master = Tkinter.Toplevel()

    from DejaVu2 import Viewer
    vi = Viewer()
    from DejaVu2.colorMap import ColorMap
    
    test = MaterialEditor(master, vi)
    
    from DejaVu2.Spheres import Spheres
    s = Spheres('test', centers = [[0.,0,0], [5,0,0.], [0,5,0]], quality=30,
                matName='tropical', matInd=0)
    vi.AddObject(s)

    from DejaVu2.materialsDef.tropical import tropical
    test.defineMaterial(tropical[0])
    test.setObject(s)

##      root1 = Tkinter.Tk()
##      cc = ColorChooser(root1, title='Color chooser')

