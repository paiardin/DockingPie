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

def getRendering(viewer, checkAnimatable=False):
    """
    build and return a dictionary storing the state of all geometries

    renderingDict <- getRendering(viewer)

    """
    statesMem = {}
    for g in viewer.rootObject.AllObjects():
        if checkAnimatable:
            if not g.animatable:
                continue
        state = g.getState(full=1)
        if state.has_key('rotation'):
            del state['rotation']
        if state.has_key('translation'):
            del state['translation']
        if state.has_key('scale'):
            del state['scale']
        if state.has_key('pivot'):
            del state['pivot']
        del state['name']
        if g.parent is not None and g.parent.visible == 0:
            state['visible'] = 0
        #add clipping planes info
        clip = []
        if hasattr(g, 'clipP'):
            for cp in g.clipP:
                clip.append([cp.num, g.clipSide[cp.num], False])
            for cp in g.clipPI:
                clip.append([cp.num, g.clipSide[cp.num], True])
        state['clipPlanes'] = clip
        statesMem[g.fullName] = state
        
    # get Camera state:
    cam = viewer.currentCamera
    state = {#'height':cam.height, 'width':cam.width,
             'near':cam.near, 'far':cam.far,
             'color':cam.backgroundColor,
             'antialiased':cam.antiAliased,
             'boundingbox':cam.drawBB,
             'rotation':list(cam.rotation),
             'translation':list(cam.translation),
             'scale':list(cam.scale),
             'pivot':list(cam.pivot),
             'direction':list(cam.direction),
             'lookAt':list(cam.lookAt),
             'projectionType': cam.projectionType,
             'drawThumbnail':cam.drawThumbnailFlag,
             'contours': cam.contours}
    statesMem['camera'] = state
    # get light state
    for i,l in enumerate(viewer.lights):
            if not l._modified: continue
            statesMem['light%d'%i] = l.getState()
    for i,c in enumerate(viewer.clipP):
            if not c._modified: continue
            statesMem['clipplane%d'%i] = c.getState()
    return statesMem

from opengltk.OpenGL import GL

def setRendering( viewer, rendering, redraw=False):
    """
    restore a given rendering for a viewer

    None <- setRendering( viewer, rendering)
    """
    old = viewer.suspendRedraw 
    viewer.suspendRedraw = True
    for g in viewer.rootObject.AllObjects():
        if not rendering.has_key(g.fullName):
            g.Set(visible = 0)
        else:
            vis = rendering[g.fullName].get('visible', None)
            if vis == 0:
                g.Set(visible = 0)
            else:
                tmp = rendering[g.fullName].copy()
                mat = tmp.pop('rawMaterialF', None)
                if mat: g.materials[GL.GL_FRONT].Set(**mat)
                mat = tmp.pop('rawMaterialB', None)
                if mat: g.materials[GL.GL_BACK].Set(**mat)
                if tmp.has_key('clipPlanes'):
                    tmp.pop('clipPlanes')
                g.Set( **tmp )

    cam = viewer.currentCamera
    if rendering.has_key('camera'):
        cam.Set(**rendering['camera'])
    viewer.suspendRedraw = old
    if redraw:
        viewer.OneRedraw()

    
def setObjectRendering(viewer, object, rendering):
    if not rendering: return
    tmp = rendering.copy()
    mat = tmp.pop('rawMaterialF', None)
    if mat: object.materials[GL.GL_FRONT].Set(**mat)
    mat = tmp.pop('rawMaterialB', None)
    if mat: object.materials[GL.GL_BACK].Set(**mat)
    if tmp.has_key('clipPlanes'):
        tmp.pop('clipPlanes')
    object.Set( **tmp )



def getOrientation( geom , camera=True):
    """
    build and return a dictionary storing the orientation for geom

    orientDict <- getOrientation( geom )

    """
    orientMem = {}
    orientMem['rotation'] = geom.rotation[:]
    orientMem['translation'] = geom.translation.copy()
    orientMem['scale'] = geom.scale[:]
    orientMem['pivot'] = geom.pivot[:]
    if camera:
        cam =  geom.viewer.currentCamera
        orientMem['fieldOfView'] = cam.fovy
        orientMem['lookFrom'] = cam.lookFrom.copy()
    
    return orientMem


def setOrientation(geom, orient):
    """
    restore a given orientation for a geometry and viewer

    None <- setOrientation( geom, orient)
    """
    geom.Set(rotation=orient['rotation'],
             translation=orient['translation'],
             scale=orient['scale'],
             pivot=orient['pivot'])
    if orient.has_key('fieldOfView'):
        geom.viewer.currentCamera._setFov(orient['fieldOfView'])
    if orient.has_key('lookFrom'):
        geom.viewer.currentCamera._setLookFrom(orient['lookFrom'])
