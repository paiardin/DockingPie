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

ball_vert_simple="""
//attribute vec3 position;
attribute mat4 transformmatrix;
varying vec3 vertex_color;

void main (void)
{
    mat4 mvp = gl_ModelViewProjectionMatrix * transformmatrix;
    //vec4 p = vec4(position,1.0);
    gl_Position = mvp * gl_Vertex;//gl_ModelViewProjectionMatrix * p;
    vertex_color = vec3(1.0,0.0,0.0);//gl_Normal.xyz;
}
//#version 120
//#extension GL_EXT_gpu_shader4 : enable

"""
ball_vert="""
attribute vec3 position;
attribute vec2 texco;
attribute vec3 posoffset;

//attribute vec3 offset;
//uniform sampler2DRect texturePosition;\n
uniform sampler2DRect textureColors;\n//or in the vbo ?
uniform sampler2DRect textureSizes;\n
uniform sampler2DRect textureScale;\n
//uniform int textureSize;
varying vec3 vertex_color;

void main (void)
{
    //vec2 coord = texco;
    vec4 color = texture2DRect(textureColors, texco.xy);\n
    vec4 p = vec4(position+posoffset,1.0);
    /*
    vec4 spaceposition;
    vec4 sphereposition;
    vec2 coord = texco;//vec2 ( float(gl_InstanceID % textureSize),float(gl_InstanceID)/float(textureSize));
    
    float radius = texture2DRect(textureSizes, coord.xy).x *\n
    	 texture2DRect(textureScale, texco.xy).x * 10.0 ;\n

    spaceposition = texture2DRect(texturePosition, coord.xy);\n
    spaceposition.w = 1.0;\n
    sphereposition = texture2DRect(texturePosition, coord.xy);\n
    sphereposition.w = 1.0;\n
    if (radius < 1.0)\n
        spaceposition.xyz += p.xyz;\n
    else\n
        spaceposition.xyz += p.xyz*radius*radius;\n
    */
    gl_Position = gl_ModelViewProjectionMatrix*p;//spaceposition;\n
    vertex_color = color.xyz;//vec3(texco.xy,0.0);//vec3(tex_co.xy,0.0);//gl_Normal.xyz;color.xyz;//
}
"""
ballimproved_vert = """
\n
#extension GL_ARB_texture_rectangle : enable\n
\n
attribute vec3 position;
attribute vec2 texco;
//attribute vec3 posoffset;

varying vec4 i_near;\n
varying vec4 i_far;\n
varying vec4 sphereposition;\n
varying vec4 color;\n
varying float radius;\n
uniform sampler2D texturePosition;\n
uniform sampler2D textureColors;\n
uniform sampler2D textureSizes;\n
uniform sampler2D textureScale;\n
uniform int textureSize;
\n
\n

void main()  {\n
    //(2i + 1)/(2N)
    vec4 spaceposition;
    vec2 c = vec2(texco.x/float(textureSize),texco.y/float(textureSize));
    vec2 mult=(2.0*texco.xy + 1.0)/(2.0*float(textureSize));
    color =  texture2D(textureColors, mult.xy);\n////vec4(texco.xy/float(textureSize),0.0,1.0);//
    //sphereposition = vec4(posoffset,1.0);
    //spaceposition = vec4(posoffset,1.0);
    \n
    //why not one map ? is faster here 
    radius = texture2D(textureSizes, mult.xy).x *\n
    		 texture2D(textureScale, mult.xy).x * 10.0 ;\n//
    spaceposition = texture2D(texturePosition, mult.xy);\n
    spaceposition.w = 1.0;\n
    sphereposition = texture2D(texturePosition, mult.xy);\n
    sphereposition.w = 1.0;\n
    if (radius < 1.0)\n
        spaceposition.xyz += position;\n
    else\n
        spaceposition.xyz += position*radius*radius;\n
    gl_Position = (gl_ModelViewProjectionMatrix*spaceposition);\n
      // Calcul near from position\n
      vec4 near = gl_Position ;\n
      near.z = 0.0 ;\n
      i_near = (gl_ModelViewProjectionMatrixInverse*near) ;\n
      //i_near = near;\n
      // Calcul far from position\n
      vec4 far = gl_Position ;\n
      far.z = far.w ;\n
      i_far = (gl_ModelViewProjectionMatrixInverse*far) ;\n}\n
"""
