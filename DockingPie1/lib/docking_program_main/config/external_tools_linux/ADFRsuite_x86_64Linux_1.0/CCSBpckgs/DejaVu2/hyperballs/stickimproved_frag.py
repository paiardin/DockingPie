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

stickimproved_frag = """

#extension GL_ARB_texture_rectangle : enable


struct Ray {
   vec3 origin ;
   vec3 direction ;
};


varying mat4  	matrix_near;
varying vec4  		color_atom1;
varying vec4  		color_atom2;
varying float   	shrink;

varying vec4 		prime1;
varying vec4		prime2;

vec4 lit(float NdotL, float NdotH, float m) {

  float ambient = 1.0;
  float diffuse = max(NdotL, 0.0);
  float specular = pow(NdotH,m);
  if(NdotL < 0.0 || NdotH < 0.0)
  	specular = 0.0;

  return vec4(ambient, diffuse, specular, 1.0);
}


bool cutoff_plane (vec3 M, vec3 cutoff, vec3 x3){
  float a = x3.x;
  float b = x3.y;
  float c = x3.z;
  float d = -x3.x*cutoff.x-x3.y*cutoff.y-x3.z*cutoff.z;
  float l = a*M.x+b*M.y+c*M.z+d;
  if (l<0.0) {return true;}
  else{return false;}
}


vec3 isect_surf(Ray r, mat4 matrix_coef){
   vec4 direction = vec4(r.direction, 0.0);
   vec4 origin = vec4(r.origin, 1.0);
   float a = dot(direction,(matrix_coef*direction));
   float b = dot(origin,(matrix_coef*direction));
   float c = dot(origin,(matrix_coef*origin));
   float delta =b*b-a*c;
   if (delta<0.0) discard;
   float t1 =(-b-sqrt(delta))/a;

   // Second solution not necessary if you don't want
   // to see inside spheres and cylinders, save some fps
   //float t2 = (-b+sqrt(delta)) / a  ;
   //float t =(t1<t2) ? t1 : t2;

   return r.origin+t1*r.direction;
}


Ray primary_ray(vec4 near1, vec4 far1){
    vec3 near=near1.xyz/near1.w;
    vec3 far=far1.xyz/far1.w;
    return Ray(near,far-near);
}


float update_z_buffer(vec3 M, mat4 ModelViewP){
    float  depth1;
    vec4 Ms=(ModelViewP*vec4(M,1.0));
    return depth1=(1.0+Ms.z/Ms.w)/2.0;
}


void main()
{

vec4 i_near, i_far, focus;
vec3 e3, e1, e1_temp, e2;

i_near = vec4(matrix_near[0][0],matrix_near[0][1],matrix_near[0][2],matrix_near[0][3]);
i_far  = vec4(matrix_near[1][0],matrix_near[1][1],matrix_near[1][2],matrix_near[1][3]);
focus = vec4(matrix_near[2][0],matrix_near[2][1],matrix_near[2][2],matrix_near[2][3]);
e3 = vec3(matrix_near[3][0],matrix_near[3][1],matrix_near[3][2]);

e1.x = 1.0;
e1.y = 1.0;
e1.z = ( (e3.x*focus.x + e3.y*focus.y + e3.z*focus.z) - e1.x*e3.x - e1.y*e3.y)/e3.z;
e1_temp = e1 - focus.xyz;
e1 = normalize(e1_temp);

e2 = normalize(cross(e1,e3));


vec4 equation = focus;

float shrinkfactor = shrink;
float t1 = -1.0/(1.0-shrinkfactor);
float t2 = 1.0/(shrinkfactor);

vec4 colonne1, colonne2, colonne3, colonne4;
mat4 mat;

vec3 equation1 = vec3(t2,t2,t1);


  float A1 = - e1.x*equation.x - e1.y*equation.y - e1.z*equation.z;
  float A2 = - e2.x*equation.x - e2.y*equation.y - e2.z*equation.z;
  float A3 = - e3.x*equation.x - e3.y*equation.y - e3.z*equation.z;

  float A11 = equation1.x*e1.x*e1.x +  equation1.y*e2.x*e2.x + equation1.z*e3.x*e3.x;
  float A21 = equation1.x*e1.x*e1.y +  equation1.y*e2.x*e2.y + equation1.z*e3.x*e3.y;
  float A31 = equation1.x*e1.x*e1.z +  equation1.y*e2.x*e2.z + equation1.z*e3.x*e3.z;
  float A41 = equation1.x*e1.x*A1   +  equation1.y*e2.x*A2   + equation1.z*e3.x*A3;

  float A22 = equation1.x*e1.y*e1.y +  equation1.y*e2.y*e2.y + equation1.z*e3.y*e3.y;
  float A32 = equation1.x*e1.y*e1.z +  equation1.y*e2.y*e2.z + equation1.z*e3.y*e3.z;
  float A42 = equation1.x*e1.y*A1   +  equation1.y*e2.y*A2   + equation1.z*e3.y*A3;

  float A33 = equation1.x*e1.z*e1.z +  equation1.y*e2.z*e2.z + equation1.z*e3.z*e3.z;
  float A43 = equation1.x*e1.z*A1   +  equation1.y*e2.z*A2   + equation1.z*e3.z*A3;

  float A44 = equation1.x*A1*A1 +  equation1.y*A2*A2 + equation1.z*A3*A3 - equation.w;

  colonne1 = vec4(A11,A21,A31,A41);
  colonne2 = vec4(A21,A22,A32,A42);
  colonne3 = vec4(A31,A32,A33,A43);
  colonne4 = vec4(A41,A42,A43,A44);

  mat = mat4(colonne1,colonne2,colonne3,colonne4);



 // Ray calculation using near and far
 Ray ray = primary_ray(i_near,i_far) ;

 // Intersection between ray and surface for each pixel
 vec3 M;
 M = isect_surf(ray, mat);

 // Recalculate the depth in function of the new pixel position
 gl_FragDepth = update_z_buffer(M, gl_ModelViewProjectionMatrix) ;

 // cut the extremities of bonds to superimpose bond and spheres surfaces
 if (cutoff_plane(M, prime1.xyz, -e3) || cutoff_plane(M, prime2.xyz, e3)){discard;}


  // Transform normal to model space to view-space
  vec4 M1 = vec4(M,1.0);
  vec4 M2 =  mat*M1;
  vec3 normal = normalize((gl_ModelViewMatrixInverseTranspose*M2).xyz);

  // Give light vector position perpendicular to the screen
  vec3 lightvec = normalize(vec3(0.0,0.0,1.2));
  vec3 eyepos = vec3(0.0,0.0,1.0);

  // calculate half-angle vector
  vec3 halfvec = normalize(lightvec + eyepos);

  // Parameters used to calculate per pixel lighting
  // see http://http.developer.nvidia.com/CgTutorial/cg_tutorial_chapter05.html
  float diffuse = dot(normal,lightvec);
  float specular = dot(normal, halfvec);
  vec4 lighting = lit(diffuse, specular, 32.0);

  // Mix the color bond in function of the two atom colors
  float distance_ratio = ((M.x-prime2.x)*e3.x + (M.y-prime2.y)*e3.y +(M.z-prime2.z)*e3.z)/distance(prime2.xyz,prime1.xyz);
  // lerp function not in GLSL. Find something else ...
  vec3 diffusecolor = mix( color_atom2.xyz, color_atom1.xyz, distance_ratio );



  vec3 specularcolor = vec3(1.0,1.0,1.0);

  // Give color parameters to the Graphic card
  gl_FragColor.rgb = lighting.y * diffusecolor + lighting.z * specularcolor;
  gl_FragColor.a = 1.0;

  	 // ############## Fog effect #####################################################
	 // and uncomment the next lines.
	 // Color of the fog: white
	 //float fogDistance  = update_z_buffer(M, gl_ModelViewMatrix) ;
  	 //float fogExponent  = fogDistance * fogDistance * 0.007;
	 //vec3 fogColor   = vec3(1.0, 1.0, 1.0);
	 //float fogFactor   = exp2(-abs(fogExponent));
	 //fogFactor = clamp(fogFactor, 0.0, 1.0);

	 //vec3 final_color = lighting.y * diffusecolor + lighting.z * specularcolor;
	 //gl_FragColor.rgb = mix(fogColor,final_color,fogFactor);
  	 //gl_FragColor.a = 1.0;
	 // ##################################################################################

}
"""

