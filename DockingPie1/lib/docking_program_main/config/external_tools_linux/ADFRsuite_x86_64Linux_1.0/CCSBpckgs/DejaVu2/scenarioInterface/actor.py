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

from Scenario2.actor import  Actor, CustomActor
from Scenario2.datatypes import DataType
from Scenario2.interpolators import Interpolator, FloatScalarInterpolator,\
     BooleanInterpolator
from SimPy.Simulation import now
import numpy

## class ActionWithRedraw(Action):
##     # subclass Action to have these actions set a redraw flag in the director
##     # so that the redraw process only triggers a redraw if something changed

##     def __init__(self, *args, **kw):
##         Action.__init__(self, *args, **kw)
##         self.postStep_cb = (self.setGlobalRedrawFlag, (), {})


##     def setGlobalRedrawFlag(self):
##         #print 'setting needs redraw at', now(), self._actor().name
##         self._actor()._maa().needsRedraw = True


class RedrawActor(Actor):

    def __init__(self, vi, *args, **kw):
        Actor.__init__(self, "redraw", vi)
        #self.addAction( Action() )
        self.visible = False
        self.recording = False
        self.scenarioname =  "DejaVu2Scenario"
        self.guiVar = None


    def setValueAt(self, frame=None, off=0):
        self.setValue()


    def setValue(self, value=None):

        needsRedraw = True
        if self._maa:
            if self._maa().needsRedraw == False:
                needsRedraw = False
                if self._maa()._director and self._maa()._director().needsRedraw:
                    needsRedraw = True
        if not needsRedraw: return
        #print "oneRedraw from redraw actor"
        self.object.OneRedraw()
        # handle recording a frame if need be
        camera = self.object.currentCamera
        if hasattr(camera, 'videoRecordingStatus'):
            if camera.videoRecordingStatus == 'recording':
                camera.recordFrame()


    def onAddToDirector(self):
        gui = self._maa().gui
        if gui:
            try:
                from DejaVu2.Camera import RecordableCamera
                isrecordable = True
            except:
                isrecordable = False
            if isrecordable:
                camera = self.object.currentCamera
                if isinstance(camera, RecordableCamera):
                    gui.createRecordMovieButton()

    def startRecording(self, file):
        camera = self.object.currentCamera
        camera.setVideoOutputFile(file)
        camera.setVideoParameters()
        camera.videoRecordingStatus = 'recording'


    def stopRecording(self):
        self.object.currentCamera.stop()


    def clone(self):
        return RedrawActor( self.object )


def getActorName(object, propname):
    # create a name for object's actor
    objname = object.name
    import string
    if hasattr(object, "fullName"):
        names = object.fullName.split("|")
        if len(names)> 1:
            objname = string.join(names[1:], ".")
    if objname.count(" "):
        # remove whitespaces from the object's name
        objname = string.join(objname.split(), "")
    return "%s.%s"%(objname, propname)


class DejaVu2Actor(CustomActor):
    """
    class for DejaVu2 actors.
        initialValue= None  - initial value of the object's attribute (object.name), 
        interp = None       - interpolator class,
        setFunction = None  - function to set the value on the object,
                              if None, object.Set(name=value) will be used
        setMethod:         method of the object to be called at each time step.
                           The function will be called using  obj.method(value)
        getFunction = None  - function that can be called to get the
                              current value of the attribute (object.name)
                              The function and its arguments have to be specified as a
                              3-tuple (func, *args, **kw). It will be called using
                              func(*(object,)+args), **kw) if it is a function
                              or func(*args, **kw) if it is a method;
                              if None, getattr(object, name) will be used to get the value
     to set the value and getattr(geom, name) to get the value
    """
    def __init__(self, name, object, initialValue=None, datatype=None,
                 interp=None, setFunction=None, setMethod=None,
                 getFunction=None):

        assert isinstance(name, str)
        self.varname = name
        if not getFunction:
            if hasattr(object, name):
                getFunction = (getattr, (name,), {})
        actorname = getActorName(object, name)
        ## objname = object.name
##         if objname.count(" "):
##                 # remove whitespaces from the object's name
##                 import string
##                 objname = string.join(objname.split(), "")
        CustomActor.__init__(self, actorname, object,
                             initialValue, datatype, interp, setFunction=setFunction, setMethod=setMethod, getFunction=getFunction)
        if self.hasGetFunction:
            self.recording = True
        self.guiVar = None
        self.scenarioname =  "DejaVu2Scenario" 


    def setValue(self, value):
        if self.setFunction:
            self.setFunction( *(self, value) )
        elif self.setMethod:
            self.setMethod(value)
        else:
            self.object.Set( **{self.varname:value} )


    def setValueAt(self, frame, off=0):
        #print 'GGGGGGGGGGGGGGG', self.name
        value = self.actions.getValue(frame-off)
        if value != 'Nothing There':
            self.setValue(value)
            maa = self._maa()
            maa.needsRedraw = True
            #print 'HHHHHHHHHHHHHHHHHH', maa.name
##             if maa._director is not None:
##                 print maa._director(), 'LLLLLLLLLLLLLLLL' 
##                 maa._director().needsRedraw = True

    def clone(self):
        if self.setMethod is not None:
            setMethod = self.setMethod.__name__
        else:
            setMethod = None
        datatype = None
        if self.datatype is not None:
            datatype = self.datatype.__class__
        newActor = self.__class__(
            self.varname, self.object, initialValue=self.initialValue,
            datatype=datatype, interp=self.interp,
            setFunction=self.setFunction, setMethod=setMethod)

        newActor.getFuncTuple = self.getFuncTuple
        if self.getFuncTuple:
            newActor.hasGetFunction = True
        newActor.behaviorList.easeIn = self.behaviorList.easeIn
        newActor.behaviorList.easeOut = self.behaviorList.easeOut        
        newActor.behaviorList.constant = self.behaviorList.constant            

        return newActor
 
from interpolators import MaterialInterpolator
from Scenario2.datatypes import FloatType, IntType, BoolType,IntVectorType,\
     FloatVectorType, IntVarType, FloatVarType, VarVectorType


class DejaVu2CameraActor(DejaVu2Actor):
    
    def __init__(self, name, object, initialValue=None, datatype=DataType,
                 interp=BooleanInterpolator, setFunction=None,
                 setMethod=None, getFunction=None):
        
        assert  hasattr(object, 'viewer')
        self.cameraind = 0
        for i, c in enumerate(object.viewer.cameras):
            if c == object:
                self.cameraind=i
                break
        DejaVu2Actor.__init__(self, name, object, initialValue, datatype,
                             interp, setFunction=setFunction, setMethod=setMethod,
                             getFunction=getFunction)
        


    def getCameraObject(self):
        vi = self.object.viewer
        camera = vi.cameras[self.cameraind]
        if camera != self.object:
            self.object = camera
            if self.setMethod is not None:
                method = getattr(self.object, self.setMethod.__name__)
                self.setMethod = method
        return camera


    def setValue(self, value):
        self.getCameraObject()
        DejaVu2Actor.setValue(self, value)


    def getValueFromObject(self):
        self.getCameraObject()
        return DejaVu2Actor.getValueFromObject(self)
        
        
    

class DejaVu2GeomVisibilityActor(DejaVu2Actor):
    """
    Actor to set geometry visibility flag
    when set to 1 we need to make sure each parent's visibility flag is 1
    else the object will not appear
    """
    def __init__(self, name, object, initialValue=None, datatype=DataType,
                 interp=BooleanInterpolator, setFunction=None,
                 setMethod=None, getFunction=None):
        
        DejaVu2Actor.__init__(self, name, object, initialValue, datatype,
                             interp, getFunction=self.getValue)


    def setValue(self, value):
        obj = self.object
        self.object.Set( visible = value)
        if value: # set all parents visibility to 1
            while obj.parent:
                obj = obj.parent
                obj.Set( visible=value )
        #else:
        #    for child in obj.AllObjects():
        #       if child!=self:
        #           child.Set( visible = value)
                

    def getValue(self):
        # MS: maybe we should return 0 if 1 parent is not visible
        return self.object.visible
    


class DejaVu2GeomEnableClipPlaneActor(DejaVu2Actor):
    """
    Actor to enable clipping plane for a geometry
    """
    def __init__(self, name, object, initialValue=None, datatype=DataType,
                 interp=None, setFunction=None,
                 setMethod=None, getFunction=None):
        
        DejaVu2Actor.__init__(self, name, object, initialValue, datatype,
                             interp, getFunction=None)


    def setValue(self, value):
        if not value:
            return
        geom = self.object
        for clip in value:
            assert len(clip) == 4
            if clip[3] == True:
                #print "addclipplane:", clip[0], clip[0].num
                geom.AddClipPlane(clip[0], clip[1], clip[2])
            else:
                #print "removeClipPlane:", clip[0], clip[0].num
                geom.RemoveClipPlane(clip[0])


    def getValue(self):
        geom = self.object
        val = []
        clip = []
        inh = False
        if len(geom.clipP): clip = geom.clipP
        elif len(geom.clipPI):
            clip = geom.clipPI
            inh = True
        for cp in clip:
            val.append([cp, geom.clipSide[cp.num], inh, True])
        return val



class DejaVu2MaterialActor(DejaVu2Actor):

    def __init__(self, name, object, initialValue=None, datatype=DataType,
                 interp=MaterialInterpolator, setFunction=None, setMethod=None, getFunction=None, mode="front"):
        self.mode = mode
        DejaVu2Actor.__init__(self, name, object, initialValue, datatype,
                             interp, getFunction=self.getValue)
        self.hasGetFunction = True


    def setValue(self, value):
        object = self.object
       ##  i=0
##         for v,name in zip(value, ('ambi', 'emis', 'spec', 'shini')):
##             #if self.activeVars[i]:
##             object.Set(propName=name, materials=v, transparent='implicit', redo=1)
##             i +=1
        #mat = [ value[0], None, value[1], value[2], value[3], None]
        mat = {'ambient':value[0], 'emission': value[1], 'specular':value[2], 'shininess':value[3]}
        mask = [1, 0, 1, 1, 1, 0]
        if self.mode == "front":
            object.Set(rawMaterialF=mat, matMask=mask,transparent='implicit', redo=1)
        else:
            object.Set(rawMaterialB=mat, matMask=mask,transparent='implicit', redo=1)
 
    def getValue(self):
        if self.mode == "front":
            mat = self.object.materials[1028].prop
        else:
            mat = self.object.materials[1029].prop
        return [mat[0], mat[2], mat[3], mat[4]]


    def compareValues(self, oldval, newval):
        vvt = VarVectorType()
        fvt = FloatVectorType()
        for i in range(3):
            if not vvt.isEqual(oldval[i], newval[i]):
                return False
        if not fvt.isEqual(oldval[3], newval[3]):
            return False
        return True


    def clone(self):
        newactor = DejaVu2Actor.clone(self)
        newactor.mode = self.mode
        newactor.hasGetFunction  = self.hasGetFunction 
        return newactor
       

from Scenario2.interpolators import FloatVectorInterpolator, VarVectorInterpolator

class DejaVu2ScissorActor(DejaVu2Actor):
    """ This actor manages resizing of DejaVu2 object's scissor"""

    def __init__(self, name, object, initialValue=None, datatype=FloatVectorType,
                 interp=FloatVectorInterpolator, setFunction=None, setMethod=None, getFunction=None):

        DejaVu2Actor.__init__(self, name, object, initialValue, datatype, interp, getFunction=self.getValue)
        self.hasGetFunction = True
        self.varnames = ['scissorH', 'scissorW', 'scissorX', 'scissorY']
        
        self.activeVars = ['scissorH', 'scissorW', 'scissorX', 'scissorY']

    def setValue(self, value):
        object = self.object
        kw = {}
        for i, var in enumerate (self.varnames):
            if var in self.activeVars:
                kw[var] = value[i]
        object.Set(**kw)

    def getValue(self):
        obj = self.object
        return [float(obj.scissorH), float(obj.scissorW),
                float(obj.scissorX), float(obj.scissorY)] 


from Scenario2.interpolators import RotationInterpolator
from mglutil.math.transformation import UnitQuaternion, Transformation
from mglutil.math.rotax import mat_to_quat
from DejaVu2.Camera import Camera
from math import cos, acos, sin, pi
from Scenario2.interpolators import matToQuaternion, quatToMatrix

class DejaVu2RotationActor(DejaVu2Actor):
    """
    This actor manages rotation of DejaVu2 object by setting the rotation
    """

    def __init__(self, name, object, initialValue=None, datatype=FloatVectorType,
                 interp=RotationInterpolator, setFunction=None, setMethod=None, getFunction=None):
        
        DejaVu2Actor.__init__(self, name, object, initialValue, datatype, interp, getFunction=self.getValue)
        self.hasGetFunction = True

##     def setValue(self, value):
##         object = self.object
##         q = (value[0], numpy.array(value[1:], 'f'))  
##         t = Transformation(quaternion=q)
##         object.Set(rotation = t.getRotMatrix(shape=(16,), transpose=1))
##         if isinstance(object, Camera):
##             object.Redraw()

##     def getValue(self):
##         obj = self.object
##         value = self.object.rotation
##         if len(value)==16:
##             q = UnitQuaternion(mat_to_quat(value))
            
##             value = [q.real, q.pure[0], q.pure[1], q.pure[2]]
##             #print "in DejaVu2RotationActor.getValue: ", value
##             return value


    def setValue(self, value):
        object = self.object
        mat = quatToMatrix(value)
        #object._setRotation(mat)
        object.Set(rotation=mat)
        if isinstance(object, Camera):
            object.Redraw()


    def getValue(self):
        mat = self.object.rotation
        quat = matToQuaternion(mat)
        return quat
        

        
class DejaVu2ConcatRotationActor(DejaVu2RotationActor):
    """
    This actor manages rotation of DejaVu2 object by concatenating the rotation
    """

    def __init__(self, name, object, initialValue=None,
                 datatype=FloatVectorType,
                 interp=RotationInterpolator, setFunction=None,
                 setMethod=None, getFunction=None):
        
        DejaVu2RotationActor.__init__(
            self, name, object, initialValue, datatype, interp,
            getFunction=self.getValue)
        self.hasGetFunction = True


    def setValue(self, value):
        object = self.object
        mat = quatToMatrix(value)
        object.ConcatRotation(mat)
        if isinstance(object, Camera):
            object.Redraw()
        


from mglutil.math.rotax import rotax
import math
from DejaVu2.scenarioInterface.interpolators import  RotationConcatInterpolator

class DejaVu2AxisConcatRotationActor(DejaVu2Actor):
    """
    This actor manages rotation of DejaVu2 object by concatenating the rotation
    """

    def __init__(self, name, object, initialValue=None,
                 datatype=FloatVectorType, axis=(0,1,0),
                 interp=RotationConcatInterpolator, setFunction=None,
                 setMethod=None, getFunction=None):
        
        DejaVu2Actor.__init__(self, name, object, initialValue, datatype,
                             interp)
        self.axis = axis
        self.hasGetFunction = False


    def setValue(self, value):
        #print 'Rotate by', value
        mat = rotax( (0,0,0), self.axis, value*math.pi/180.)
        object = self.object
        object.ConcatRotation(mat.flatten())
        if isinstance(object, Camera):
            object.Redraw()


    def clone(self):
        newactor = DejaVu2Actor.clone(self)
        newactor.axis = self.axis
        newactor.hasGetFunction  = self.hasGetFunction 
        return newactor



class DejaVu2ClipZActor(DejaVu2CameraActor):
    """ This actor manages the near and far camera's atributes"""

    def __init__(self, name, object, initialValue=None, datatype=FloatVectorType,
                 interp=FloatVectorInterpolator, setFunction=None, setMethod=None, getFunction=None):

        DejaVu2CameraActor.__init__(self, name, object, initialValue, datatype,
                             interp, getFunction=self.getValue)
        self.hasGetFunction = True
        self.varnames = ['near', 'far']
        
        self.activeVars = ['near', 'far']


    def setValue(self, value):
        
        camera = self.getCameraObject()
        kw = {}
        for i, var in enumerate (self.varnames):
            if var in self.activeVars:
                kw[var] = value[i]
        camera.Set(**kw)
        camera.Redraw()


    def getValue(self):
        c = self.getCameraObject()
        return numpy.array([c.near, c.far,], "f")



class DejaVu2FogActor(DejaVu2CameraActor):
    """ This actor manages the start and end atributes of camera's fog"""

    def __init__(self, name, object, attr='start', initialValue=None, datatype=FloatType,
                 interp=FloatScalarInterpolator, setFunction=None, setMethod=None, getFunction=None):
        
        from DejaVu2.Camera import Fog
        if isinstance(object, Fog):
           assert object.camera is not None
           object = object.camera()
        self.attr = attr
        DejaVu2CameraActor.__init__(self, name, object, initialValue, datatype,
                             interp, getFunction=self.getValue)
        self.hasGetFunction = True


    def setValue(self, value):
        camera = self.getCameraObject()
        kw = {self.attr: value}
        camera.fog.Set(**kw)
        #camera.Redraw()


    def getValue(self):
        c = self.getCameraObject()
        return getattr(c.fog, self.attr)

    def clone(self):
        newactor = DejaVu2CameraActor.clone(self)
        newactor.attr = self.attr
        return newactor

    
from Scenario2.interpolators import FloatVarScalarInterpolator
from interpolators import LightColorInterpolator

class DejaVu2LightColorActor(DejaVu2Actor):

    def __init__(self, name, object, initialValue=None, datatype = DataType, interp = LightColorInterpolator, setFunction=None, setMethod=None, getFunction=None):
        
        DejaVu2Actor.__init__(self, name, object, initialValue, datatype, interp,
                             getFunction=self.getValue)
        self.varnames = ['ambient', 'diffuse', 'specular']
        self.activeVars = ['ambient', 'diffuse', 'specular']
        self.hasGetFunction = True

    def setValue(self, value):
        object = self.object
        kw = {}
        for v,name in zip(value, ('ambient', 'diffuse', 'specular')):
            if name in self.activeVars:
               kw[name] = v 
        object.Set(**kw)
 
    def getValue(self):
        obj = self.object
        return [obj.ambient, obj.diffuse, obj.specular]

    def compareValues(self, oldval, newval):
        fvt = FloatVectorType()
        for i in range(3):
            if not fvt.isEqual(oldval[i], newval[i]):
                return False
        return True

    

class DejaVu2SpheresRadiiActor(DejaVu2Actor):
    """ This actor manages the raduis attribute of spheres"""

    def __init__(self, name, object, initialValue=None, datatype=FloatVarType,
                 interp=FloatVarScalarInterpolator, setFunction=None, setMethod=None, getFunction=None):
            
        DejaVu2Actor.__init__(self, name, object, initialValue,
                             datatype=datatype, interp=interp,
                             getFunction=self.getValue)
        self.hasGetFunction = True
        self.varnames = ['radii']
        self.activeVars = ['radii']

    def setValue(self, value):
        object = self.object
        object.Set(radii=value)

    def getValue(self):
        object = self.object
        if object.oneRadius:
            return object.radius
        else:
            return object.vertexSet.radii.array.ravel()

        

from interpolators import TransformInterpolator

class DejaVu2TransformationActor(DejaVu2Actor):

    def __init__(self, name, object, initialValue=None, datatype=DataType,
                 interp=TransformInterpolator, setFunction=None, setMethod=None, getFunction=None):
        DejaVu2Actor.__init__(self, name, object, initialValue, datatype,
                             interp, getFunction=self.getValue)
        self.hasGetFunction = True


    def setValue(self, value):
        object = self.object
        rotation = quatToMatrix(value[0])
        redo = False
        if object != object.viewer.rootObject:
            # need to rebuild display list if the object is not root 
             redo = True
        object.Set(rotation=rotation,translation = value[1],scale = value[2],
                   pivot = value[3], redo=redo)

 
    def getValue(self):
        obj = self.object
        rotation = matToQuaternion(obj.rotation)
        return [rotation, obj.translation[:], obj.scale[:], obj.pivot[:] ]


    def compareValues(self, oldval, newval):
        fvt = FloatVectorType()
        for i in range(len(oldval)):
            if not fvt.isEqual(oldval[i], newval[i]):
                return False
        return True


    def clone(self):
        newactor = DejaVu2Actor.clone(self)
        newactor.hasGetFunction  = self.hasGetFunction 
        return newactor

        


def getAllSubclasses(klass):
    # recursively build a list of all sub-classes
    bases = klass.__bases__
    klassList = list(bases)
    for k in bases:
        klassList.extend( getAllSubclasses(k) )
    return klassList


import inspect
from actorsDescr import actorsDescr
from mglutil.gui.BasicWidgets.Tk.thumbwheel import ThumbWheel
from Scenario2.interpolators import FloatScalarInterpolator, IntScalarInterpolator
def getAnimatableAttributes(object):
    # return two dicts that contain a dict for each attribute of object
    # that can be animated. The first dict is for attribute foudn in
    # actorsDescr, the second is for attributes picked up on the fly

    # merge the dict of attribute for all base classes
    d1 = {}
    for klass, d2 in actorsDescr.items():
        if isinstance(object, klass):
            d1.update(d2)

    d2 = {}
    # find all attribute that are float
    attrs = inspect.getmembers(object, lambda x: isinstance(x, float))
    for name,value in attrs:
        if d1.has_key(name): continue
        d2[name] = {
            'interp': FloatScalarInterpolator,
            'interpKw': {'values':[value, value]},
            'valueWidget': ThumbWheel,
            'valueWidgetKw': {'type': 'float', 'initialValue':value},
            }

    # find all attribute that are bool or int
    attrs = inspect.getmembers(object, lambda x: isinstance(x, int))
    for name,value in attrs:
        if d1.has_key(name): continue
        d2[name] = {
            'interp': IntScalarInterpolator,
            'interpKw': {'values':[value, value]},
            'valueWidget': ThumbWheel,
            'valueWidgetKw': {'type':'int', 'initialValue':value},
            }
    return d1, d2

    
def getDejaVu2Actor(object, attribute):
    # return a DejaVu2 Actor given a DejaVu2 object and attribute name
    baseClasses = [object.__class__]
    baseClasses.extend( getAllSubclasses(object.__class__) )
    #print 'getDejaVu2Actor', object,attribute
    for klass in baseClasses:
        d1 = actorsDescr.get(klass, None)
        if d1:
            d2 = d1.get(attribute, None)
            if d2:
                actorCalss, args, kw = d2['actor']
                args = (attribute,object) + args
                actor = actorCalss( *args, **kw )

                return actor
##             else: # attribute not found in dictionary
##                 if attribute in object.keywords: # it is setable
##                     if hasattr(object, attribute):
##                         actor = DejaVu2ActorSetGetattr(
##                             object, name=actorName, setName=attribute,
##                             getName=attribute)
##                     else:
##                         return None # FIXME
##                 else:
                    
                
    return None

