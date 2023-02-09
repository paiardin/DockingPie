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

from DejaVu2.Transformable import Transformable
from DejaVu2.Displayable import Displayable
from DejaVu2.Geom import Geom
from DejaVu2.IndexedGeom import IndexedGeom
from DejaVu2.Cylinders import Cylinders
from DejaVu2.Light import Light
from DejaVu2.Spheres import Spheres
from DejaVu2.Clip import ClippingPlane
from DejaVu2.Camera import Camera, Fog
from DejaVu2.Materials import propertyNum

from opengltk.OpenGL import GL

from Scenario2.interpolators import VarVectorInterpolator, FloatVectorInterpolator,\
     IntScalarInterpolator, RotationInterpolator,\
     BooleanInterpolator, FloatVarScalarInterpolator, FloatScalarInterpolator
##     ReadDataInterpolator
from Scenario2.datatypes import FloatType, IntType, BoolType,IntVectorType,\
     FloatVectorType, IntVarType, FloatVarType, VarVectorType

from actor import DejaVu2Actor, DejaVu2MaterialActor,  DejaVu2ScissorActor, \
     DejaVu2ClipZActor, DejaVu2FogActor, DejaVu2LightColorActor, \
     DejaVu2SpheresRadiiActor, DejaVu2RotationActor, DejaVu2ConcatRotationActor, \
     DejaVu2GeomVisibilityActor, DejaVu2TransformationActor, DejaVu2GeomEnableClipPlaneActor,\
     DejaVu2CameraActor

from mglutil.gui.BasicWidgets.Tk.thumbwheel import ThumbWheel
from DejaVu2.states import getRendering, setRendering, setObjectRendering
import numpy

id4x4 = numpy.identity(4).astype('f')
actorsDescr = {

    Transformable: {
        'rotation': {
            'actor': (DejaVu2RotationActor, (), {})
             } ,

        'concatRotation': {
            'actor': (DejaVu2ConcatRotationActor, (), {})
             } ,

    
        'translation': {
            'actor': (DejaVu2Actor, (), {'interp':FloatVectorInterpolator,
                                        #'setMethod': '_setTranslation',
                                        'datatype': FloatVectorType} )
            },
        'scale': {
            'actor': (DejaVu2Actor, (), {'interp':FloatVectorInterpolator,
                                        #'setMethod': '_setScale',
                                        'datatype': FloatVectorType} )
            },
        'pivot': {
            'actor': (DejaVu2Actor, (), {'interp': FloatVectorInterpolator,
                                        #'setMethod': '_setPivot',
                                        'datatype': FloatVectorType} )
            },
        'transformation':{
            'actor': (DejaVu2TransformationActor, (), {})
             } ,
        },
    

   Displayable: {
        'colors':{
            'actor': (DejaVu2Actor, (),
              {
                'setFunction':\
                #lambda actor, value: actor.object.Set(materials=value, inheritMaterial=0),
                lambda actor, value: actor.object.Set(materials=value, redo=1, tagModified=False,
                                                      transparent='implicit', inheritMaterial=0),
                'getFunction':(lambda x,y: numpy.copy(x.materials[y].prop[1][:, :3]), (GL.GL_FRONT,), {}),
                'interp': VarVectorInterpolator,
                'datatype':VarVectorType
               },
                      ),
                 },
        'colorsB':{
            'actor': (DejaVu2Actor, (),
              {
                'setFunction':\
                #lambda actor, value: actor.object.Set(materials=value, inheritMaterial=0),
                lambda actor, value: actor.object.Set(materials=value, redo=1,
                                                      tagModified=False, polyFace = GL.GL_BACK,
                                                      transparent='implicit', inheritMaterial=0),
                'getFunction':(lambda x,y: numpy.copy(x.materials[y].prop[1][:, :3]), (GL.GL_BACK,), {}),
                'interp': VarVectorInterpolator,
                'datatype':VarVectorType
               },
                      ),
                 },
        'lineWidth': {
        'actor': (DejaVu2Actor, (),
                  {'setFunction': \
               lambda actor, value: actor.object.Set(lineWidth=value, inheritLineWidth=0),
                   'interp': IntScalarInterpolator, 'datatype':IntType})
                     },
        
        'pointWidth': {
              'actor': (DejaVu2Actor, (),
                        {'setFunction': \
               lambda actor, value: actor.object.Set(pointWidth=value, inheritPointWidth=0),
                         'interp': IntScalarInterpolator, 'datatype':IntType})
              },
        
        'rendering': {
              'actor': (DejaVu2Actor, (),
                        {'getFunction': lambda object: getRendering(object.viewer,checkAnimatable=True).get(object.fullName, None),
                         'setFunction': lambda actor, value: setObjectRendering(actor.object.viewer, actor.object, value)})
              },
    
        },
    
    Geom:  {
       'material': {
            'actor': (DejaVu2MaterialActor, (), {} ),

            },
        'opacity': {
            'actor': (DejaVu2Actor, (),
              {'setFunction': \
               lambda actor, value: actor.object.Set(opacity=value, transparent='implicit', polyFace=GL.GL_FRONT_AND_BACK,
                                                     inheritMaterial=0),
               'getFunction':(
                   lambda x,y: x.materials[y].prop[propertyNum['opacity']], (GL.GL_FRONT,), {}
                             ),
               'interp': FloatVarScalarInterpolator,
               'datatype': FloatVarType
               }
                      ),
                  },
       'opacityF': {
            'actor': (DejaVu2Actor, (),
              {'setFunction': \
               lambda actor, value: actor.object.Set(opacity=value, transparent='implicit',
                                                     polyFace=GL.GL_FRONT,
                                                     inheritMaterial=0),
               'getFunction':(
                   lambda x,y: x.materials[y].prop[propertyNum['opacity']], (GL.GL_FRONT,), {}
                             ),
               'interp': FloatVarScalarInterpolator,
               'datatype': FloatVarType
               }
                      ),
                  },
        'opacityB': {
            'actor': (DejaVu2Actor, (),
              {'setFunction': \
               lambda actor, value: actor.object.Set(opacity=value, transparent='implicit',
                                                     polyFace=GL.GL_BACK,
                                                     inheritMaterial=0),
               'getFunction':(
                   lambda x,y: x.materials[y].prop[propertyNum['opacity']], (GL.GL_BACK,), {}
                             ),
               'interp': FloatVarScalarInterpolator,
               'datatype': FloatVarType
               }
                      ),
                  },
       
        'visible': {
            'actor': (DejaVu2GeomVisibilityActor, (),
                      {'interp': BooleanInterpolator, 'datatype': BoolType})
            },
        
        'depthMask': {
            'actor': (DejaVu2Actor, (),
                      {'interp': BooleanInterpolator, 'datatype': BoolType})
            },
        
        'scissor': {
              'actor': (DejaVu2Actor, (), {'interp': BooleanInterpolator, 'datatype': BoolType})
                     }, # turns the scissor on/off.
        
        'scissorResize': {
              'actor': ( DejaVu2ScissorActor, (), {}),             
              },
       
        'xyz':{
            'actor': (DejaVu2Actor, (),
              {
                'setFunction': lambda actor, value: actor.object.Set(vertices=value),
                'getFunction':(lambda obj: obj.vertexSet.vertices.array, (), {}),
                'interp': VarVectorInterpolator,
                'datatype':VarVectorType
               },
                      ),
           },

        'clipEnable': {
              'actor': (DejaVu2GeomEnableClipPlaneActor, (), {}),
              },

##        'function':{
##     'actor':(DejaVu2Actor, (), {'setFunction': None, 'getFunction': None,
##                                'interp':ReadDataInterpolator, 'datatype': None})  
##                   }

        },

    IndexedGeom:  {
     'faces':{
            'actor': (DejaVu2Actor, (),
              {
                'setFunction': lambda actor, value: actor.object.Set(faces=value),
                'getFunction':(lambda obj: obj.faceSet.faces.array, (), {}),
                'interp': VarVectorInterpolator,
                'datatype':VarVectorType
               },),
           },

     'vertface':{
            'actor': (DejaVu2Actor, (),
              {
                'setFunction': lambda actor, value: actor.object.Set(vertices=value[0], faces=value[1]),
                'getFunction':(lambda obj: [obj.vertexSet.vertices.array, obj.faceSet.faces.array], (), {}),
               },),
           },
     
        },


        Spheres: {
        'radii':   {
             'actor': (DejaVu2SpheresRadiiActor, (), {}),
                    },
        'quality':   {
             'actor': (DejaVu2Actor, (),
                       {
                        'setFunction': lambda actor, value: actor.object.Set(quality=value),
                        'interp': IntScalarInterpolator, 'datatype':IntType
                       }),
                    },
        'vertices':{
                     'actor': (DejaVu2Actor, (), 
                {'setFunction': lambda actor, value: actor.object.Set(vertices=value),
                 'getFunction':(lambda obj: obj.vertexSet.vertices.array, (), {}) } ),
                    }
         },

        Cylinders: {
         'radii':   {
             'actor': (DejaVu2Actor, (),
                       {'setFunction': lambda actor, value: actor.object.Set(radii=list(value)),
                        'getFunction':(lambda obj: obj.vertexSet.radii.array, (), {}),
                        'interp': FloatVectorInterpolator,
                        'datatype': FloatVectorType} 
                       ),
                    },
         'quality':   {
             'actor': (DejaVu2Actor, (),
                       {
                        'setFunction': lambda actor, value: actor.object.Set(quality=value),
                        'interp': IntScalarInterpolator, 'datatype':IntType
                       }),
                    },
        'vertices':{
                     'actor': (DejaVu2Actor, (), 
                {'setFunction': lambda actor, value: actor.object.Set(vertices=value),
                 'getFunction':(lambda obj: obj.vertexSet.vertices.array, (), {}) } ),
                    }         
          },

        Fog:{
##           'start': {
##             'actor': (DejaVu2Actor, (),
##                       {
##                           'getFunction':(getattr, ('start',), {}),
##                           'interp': FloatScalarInterpolator,
##                           'datatype': FloatType
##                           },
##                       ),
##             },
##           'end': {
##             'actor': (DejaVu2Actor, (),
##                       {
##                           'getFunction':(getattr, ('end',), {}),
##                           'interp': FloatScalarInterpolator,
##                           'datatype': FloatType
##                           },
##                       ),
##             },
             'start': {'actor': (DejaVu2FogActor, (), {'attr':'start'},)},
             'end': {'actor': (DejaVu2FogActor, (), {'attr':'end'},) },

        },
    
        Camera:{
         'clipZ':{
             'actor': (DejaVu2ClipZActor, (), {} ),
                       },
         
         'backgroundColor': {
              'actor': (DejaVu2CameraActor, (),
                         #(lambda c, value: c.Set(color=value),  ),
               {'setFunction': lambda actor, value: (actor.object.Set(color=value), actor.object.fog.Set(color=value)),  #,actor.object.Redraw()),
                        
                'getFunction':(getattr, ('backgroundColor',), {}),
                'interp': FloatVectorInterpolator,
                'datatype': FloatVectorType
                },
                        ),
              
              },

         'fieldOfView': {
              'actor': (DejaVu2CameraActor, (),
                        {'setMethod': '_setFov',
                        'getFunction':(getattr, ('fovy',), {}),
                        'interp': FloatScalarInterpolator,
                        'datatype': FloatType
                         },
                        ),
              
              },
          'near': {
              'actor': (DejaVu2CameraActor, (),
                        {
                        'getFunction':(getattr, ('near',), {}),
                        'interp': FloatScalarInterpolator,
                        'datatype': FloatType
                         },
                        ),
              
              },
         'far': {
              'actor': (DejaVu2CameraActor, (),
                        {
                        'getFunction':(getattr, ('far',), {}),
                        'interp': FloatScalarInterpolator,
                        'datatype': FloatType
                         },
                        ),
              
              },
          'width': {
              'actor': (DejaVu2CameraActor, (),
                        {
                        'getFunction':(getattr, ('width',), {}),
                        'interp': IntScalarInterpolator,
                        'datatype': IntType
                         },
                        ),
              
              },
         'heigt': {
              'actor': (DejaVu2CameraActor, (),
                        {
                        'getFunction':(getattr, ('heigt',), {}),
                        'interp': IntScalarInterpolator,
                        'datatype': IntType
                         },
                        ),
              
              },
         'antialiased': {
              'actor': (DejaVu2CameraActor, (),
                        {
                        'getFunction':(getattr, ('antiAliased',), {}),
                        'interp': IntScalarInterpolator,
                        'datatype': IntType
                         },
                        ),
              
              },
         'boundingbox':   {
             'actor': (DejaVu2CameraActor, (),
                       {
                        'interp': IntScalarInterpolator,
                        'getFunction':(getattr, ('drawBB',), {}),
                        'datatype':IntType
                       }),
                    },
         'projectionType': {
            'actor': (DejaVu2CameraActor, (),
                      {'interp': BooleanInterpolator, 'datatype': BoolType})
            },
         'drawThumbnail': {
            'actor': (DejaVu2CameraActor, (),
                      {'interp': BooleanInterpolator, 'datatype': BoolType})
            },
         'contours': {
            'actor': (DejaVu2CameraActor, (),
                      {'interp': BooleanInterpolator,
                       'getFunction':(getattr, ('drawThumbnailFlag',), {}),
                       'datatype': BoolType})
            },
         
         'lookFrom': {
              'actor': (DejaVu2CameraActor, (),
                        {'setMethod': '_setLookFrom',
                         'getFunction':(getattr, ('lookFrom',), {}),
                         'interp':FloatVectorInterpolator,
                         'datatype': FloatVectorType} )
              },
         'lookAt': {
              'actor': (DejaVu2CameraActor, (),
                        {'getFunction':(getattr, ('lookAt',), {}),
                         'interp':FloatVectorInterpolator,
                         'datatype': FloatVectorType} )
              },
         
         }, 
    
        ClippingPlane: {
            'visible': {
              'actor': (DejaVu2Actor, (), {'interp': BooleanInterpolator, 'datatype': BoolType
                                          })
                     },
            'color': {
              'actor': (DejaVu2Actor, (), {'interp': FloatVectorInterpolator,
                                          'datatype': FloatVectorType 
                                          })
    
                            },
            },

    Light: {
        'visible': {
              'actor': (DejaVu2Actor, (), {'interp': BooleanInterpolator, 'datatype': BoolType
                                          })
                     },

        'color': {
               'actor':(DejaVu2LightColorActor, (), {}),
               },
          
    'direction': { 'actor': (DejaVu2Actor, (), {'interp': FloatVectorInterpolator,
                                               'datatype': FloatVectorType} )
                   }
        }
    }

