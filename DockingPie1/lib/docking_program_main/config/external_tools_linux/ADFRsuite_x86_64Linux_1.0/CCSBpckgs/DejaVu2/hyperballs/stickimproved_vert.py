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

#ifndef __STICKIMPROVED_VERT_H__
#define __STICKIMPROVED_VERT_H__
stickimproved_vert = """
\n
#extension GL_ARB_texture_rectangle : enable\n
\n
attribute vec3 position;
attribute vec2 texco1;
attribute vec2 texco2;
attribute vec2 texco3;
//attribute vec3 position1;
//attribute vec3 position2;

varying mat4  	matrix_near;\n
varying vec4  	color_atom1;\n
varying vec4  	color_atom2;\n
varying float   	shrink;\n
\n
varying vec4 		prime1;\n
varying vec4		prime2;\n
\n
uniform sampler2D texturePosition;\n//atom pos
uniform sampler2D textureColors;\n//atom color
uniform sampler2D textureSizes;\n//atom size ie vdw
uniform sampler2D textureShrink;\n//bonds shrink
uniform sampler2D textureScale;\n//bonds scale

uniform int textureSize;
uniform int textureBondsSize;

\n
void main()\n
{\n
vec2 mult1=(2.0*texco1.xy + 1.0)/(2.0*float(textureSize));
vec2 mult2=(2.0*texco2.xy + 1.0)/(2.0*float(textureSize));
vec2 mult3=(2.0*texco3.xy + 1.0)/(2.0*float(textureBondsSize));
\n
vec4 spaceposition;\n
vec3 position_atom1;\n
vec3 position_atom2;\n
vec4 vertex_position;\n
\n
color_atom1 = texture2D(textureColors, mult1.xy);\n
color_atom2 = texture2D(textureColors, mult2.xy);\n
\n
shrink = texture2D(textureShrink, mult3.xy).x;\n
\n
float radius1, radius2;\n
float size, scale;\n
size = texture2D(textureSizes, mult1.xy).x;\n
scale = texture2D(textureScale, mult3.xy).x;\n
radius1 = size * scale * 10.0;\n
size = texture2D(textureSizes, mult2.xy).x;\n
radius2 = size * scale * 10.0;\n
\n
position_atom1 = texture2D(texturePosition, mult1.xy).xyz;//position1.xyz;//\n
position_atom2 = texture2D(texturePosition, mult2.xy).xyz;//position2.xyz;//\n
\n
\n
float distance = sqrt( (position_atom1.x - position_atom2.x)*(position_atom1.x - position_atom2.x) + (position_atom1.y - position_atom2.y)*(position_atom1.y - position_atom2.y) + (position_atom1.z - position_atom2.z)*(position_atom1.z - position_atom2.z) );\n
\n
spaceposition.z = position.z * distance;\n
\n
if (radius1 > radius2) {\n
  spaceposition.y = position.y * 1.5 * radius1;\n
  spaceposition.x = position.x * 1.5 * radius1;\n
} else {\n
  spaceposition.y = position.y * 1.5 * radius2;\n
  spaceposition.x = position.x * 1.5 * radius2;\n
}\n
  spaceposition.w = 1.0;\n
\n
\n
\n
\n
  vec4 e3;\n
  vec3 e1, e1_temp, e2, e2_temp;\n
\n
  // Calculation of bond direction: e3\n
  e3.xyz = normalize(position_atom1-position_atom2);\n
\n
  // little hack to avoid some problems of precision due to graphic card limitation using float: To improve soon\n
  if (e3.z == 0.0) { e3.z = 0.0000000000001;}\n
  if ( (position_atom1.x - position_atom2.x) == 0.0) { position_atom1.x += 0.001;}\n
  if ( (position_atom1.y - position_atom2.y) == 0.0) { position_atom1.y += 0.001;}\n
  if ( (position_atom1.z - position_atom2.z) == 0.0) { position_atom1.z += 0.001;}\n
\n
  // Focus calculation\n
  vec4 focus;\n
  focus.x = ( position_atom1.x*position_atom1.x - position_atom2.x*position_atom2.x + ( radius2*radius2 - radius1*radius1 )*e3.x*e3.x/shrink )/(2.0*(position_atom1.x - position_atom2.x));\n
  focus.y = ( position_atom1.y*position_atom1.y - position_atom2.y*position_atom2.y + ( radius2*radius2 - radius1*radius1 )*e3.y*e3.y/shrink )/(2.0*(position_atom1.y - position_atom2.y));\n
  focus.z = ( position_atom1.z*position_atom1.z - position_atom2.z*position_atom2.z + ( radius2*radius2 - radius1*radius1 )*e3.z*e3.z/shrink )/(2.0*(position_atom1.z - position_atom2.z));\n
\n
 // e1 calculation\n
 e1.x = 1.0;\n
 e1.y = 1.0;\n
 e1.z = ( (e3.x*focus.x + e3.y*focus.y + e3.z*focus.z) - e1.x*e3.x - e1.y*e3.y)/e3.z;\n
 e1_temp = e1 - focus.xyz;\n
 e1 = normalize(e1_temp);\n
\n
 // e2 calculation\n
 e2_temp = e1.yzx * e3.zxy - e1.zxy * e3.yzx;\n
 e2 = normalize(e2_temp);\n
\n
 //ROTATION:\n
 // final form of change of basis matrix:\n
 mat3 R= mat3(e1.xyz, e2.xyz, e3.xyz);\n
 // Apply rotation and translation to the bond primitive\n
 vertex_position.xyz = R*spaceposition.xyz;\n
 vertex_position.w = 1.0;\n
\n
  // TRANSLATION:\n
  vertex_position.x +=  (position_atom1.x+position_atom2.x)/2.0;\n
  vertex_position.y +=  (position_atom1.y+position_atom2.y)/2.0;\n
  vertex_position.z +=  (position_atom1.z+position_atom2.z)/2.0;\n
\n
  // New position\n
  gl_Position = (gl_ModelViewProjectionMatrix*vertex_position);\n
\n
\n
\n
\n
\n
vec4 i_near, i_far;\n
\n
  // Calcul near from position\n
  vec4 near = gl_Position ;\n
  near.z = 0.0 ;\n
  near = (gl_ModelViewProjectionMatrixInverse*near) ;\n
  i_near = near;\n
  //i_near = vec4(1.0,1.0,1.0,1.0);\n
\n
  // Calcul far from position\n
  vec4 far = gl_Position ;\n
  far.z = far.w ;\n
  i_far = (gl_ModelViewProjectionMatrixInverse*far) ;\n
  //i_far = vec4(1.0,1.0,1.0,1.0);\n
\n
\n
prime1.xyz = position_atom1 - (position_atom1 - focus.xyz)*shrink;\n
prime2.xyz = position_atom2 - (position_atom2 - focus.xyz)*shrink;\n
\n
float Rsquare  = (radius1*radius1/shrink) - ( (position_atom1.x - focus.x)*(position_atom1.x - focus.x) + (position_atom1.y - focus.y)*(position_atom1.y - focus.y) + (position_atom1.z - focus.z)*(position_atom1.z - focus.z) );\n
\n
focus.w = Rsquare;\n
\n
matrix_near = mat4(i_near,i_far, focus, e3);\n
}\n
\n
void main2(){\n
\n
}\n
"""

