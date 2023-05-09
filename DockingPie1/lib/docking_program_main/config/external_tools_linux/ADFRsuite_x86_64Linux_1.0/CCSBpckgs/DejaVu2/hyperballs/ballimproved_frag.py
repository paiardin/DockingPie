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

#ifndef __BALLIMPROVED_FRAG_H__
#define __BALLIMPROVED_FRAG_H__
ball_frag_simple="""
varying vec3 vertex_color;
void main( ) {
    gl_FragColor = vec4(vertex_color,1.0);
}
"""
ball_frag="""
varying vec3 vertex_color;
void main( ) {
    gl_FragColor = vec4(vertex_color,1.0);
}
"""

ballimproved_frag = """

#extension GL_ARB_texture_rectangle : enable


varying vec4 i_near;
varying vec4 i_far;
varying vec4 sphereposition;
varying vec4 color;
varying float radius;


struct Ray{
   vec3 origin;
   vec3 direction;
};

vec3 isect_surf(Ray r, mat4 matrix_coef){
   vec4 direction = vec4(r.direction, 0.0);
   vec4 origin = vec4(r.origin, 1.0);
   float a = dot(direction,(matrix_coef*direction));
   float b = dot(origin,(matrix_coef*direction));
   float c = dot(origin,(matrix_coef*origin));
   float delta =b*b-a*c;
   if (delta<0.0) discard;//discard ?
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

vec4 lit(float NdotL, float NdotH, float m) {
  float ambient = 1.0;
  float diffuse = max(NdotL, 0.0);
  float specular = pow(NdotH,m);
  if(NdotL < 0.0 || NdotH < 0.0)
  	specular = 0.0;
  return vec4(ambient, diffuse, specular, 1.0);
}


void main()
{
 // Create matrix for the quadric equation of the sphere
 vec4 colonne1, colonne2, colonne3, colonne4;
 mat4 mat;
 vec4 equation = vec4(1.0,1.0,1.0,radius*radius);


 colonne1 = vec4(equation.x,0.0,0.0,-equation.x*sphereposition.x);
 colonne2 = vec4(0.0,equation.y,0.0,-equation.y*sphereposition.y);
 colonne3 = vec4(0.0,0.0,equation.z,-equation.z*sphereposition.z);
 colonne4 = vec4(-equation.x*sphereposition.x,
				 -equation.y*sphereposition.y,
				 -equation.z*sphereposition.z,
				 -equation.w +  equation.x*sphereposition.x*sphereposition.x + equation.y*sphereposition.y*sphereposition.y +equation.z*sphereposition.z*sphereposition.z);

 mat = mat4(colonne1,colonne2,colonne3,colonne4);

 // Ray calculation using near and far
 Ray ray = primary_ray(i_near,i_far) ;
 // Intersection between ray and surface for each pixel
 vec3 M;
 M = isect_surf(ray, mat);

 // Recalculate the depth in function of the new pixel position
  gl_FragDepth = update_z_buffer(M, gl_ModelViewProjectionMatrix) ;


  // Transform normal to model space to view-space
  vec4 M1 = vec4(M,1.0);
  vec4 M2 =  (mat*M1);
  vec3 normal = normalize((gl_ModelViewMatrixInverseTranspose*M2).xyz);

  // Give light vector position perpendicular to the screen
  vec3 lightvec = normalize(vec3(0.0,0.0,1.2));
  vec3 eyepos = vec3(0.0,0.0,1.0);

  // calculate half-angle vector
  vec3 halfvec = normalize(lightvec + eyepos);

  // Parameters used to calculate per pixel lighting
  // see http://http.developer.nvidia.com/CgTutorial/cg_tutorial_chapter05.html

  float shininess = 0.5;
  float diffuse = dot(normal,lightvec);
  float specular = dot(halfvec, normal);
  vec4 lighting = lit(diffuse, specular, 32.0) ;

  vec3 diffusecolor = color.xyz;
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
void main2()
{
 // Create matrix for the quadric equation of the sphere
 vec4 colonne1, colonne2, colonne3, colonne4;
 mat4 mat;
 vec4 equation = vec4(1.0,1.0,1.0,radius*radius);


 colonne1 = vec4(equation.x,0.0,0.0,-equation.x*sphereposition.x);
 colonne2 = vec4(0.0,equation.y,0.0,-equation.y*sphereposition.y);
 colonne3 = vec4(0.0,0.0,equation.z,-equation.z*sphereposition.z);
 colonne4 = vec4(-equation.x*sphereposition.x,
				 -equation.y*sphereposition.y,
				 -equation.z*sphereposition.z,
				 -equation.w +  equation.x*sphereposition.x*sphereposition.x + equation.y*sphereposition.y*sphereposition.y +equation.z*sphereposition.z*sphereposition.z);

     mat = mat4(colonne1,colonne2,colonne3,colonne4);
    // Give color parameters to the Graphic card
 // Ray calculation using near and far
 Ray ray = primary_ray(i_near,i_far) ;
 // Intersection between ray and surface for each pixel
 vec3 M;
 M = isect_surf(ray, mat);

 // Recalculate the depth in function of the new pixel position
 gl_FragDepth = update_z_buffer(M, gl_ModelViewProjectionMatrix) ;
    gl_FragColor = color;
}
"""

