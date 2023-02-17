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

import numpy

def blurSpheres(coords, elist, totmass, box=0, apix=1.0,res=1.5, solv=False):
  totchg = 0.0
  xr  = [[0.,0.,0.],[0.,0.,0.],[0.,0.,0.]]
  rmax = 0.0
  g = [0., 0., 0.]
  for i, c in enumerate(coords):
    e = elist[i]

    x,y,z = c
    g[0] += x
    g[1] += y
    g[2] += z
    elist.append(e)
    tr=x*x + y*y + z*z
    if tr > rmax:
      rmax=tr

    if x<xr[0][0]: xr[0][0]=x
    if x>xr[1][0]: xr[1][0]=x
    if x<xr[0][1]: xr[0][1]=y
    if x>xr[1][1]: xr[1][1]=y
    if x<xr[0][2]: xr[0][2]=z
    if x>xr[1][2]: xr[1][2]=z

    xr[2][0]+=x*e
    xr[2][1]+=y*e
    xr[2][2]+=z*e
    totchg+=e

  xr[2][0]/=totchg
  xr[2][1]/=totchg
  xr[2][2]/=totchg

  from math import sqrt, ceil, pi, exp, fabs
  rmax = max( fabs(xr[1][0]-xr[0][0]),
              fabs(xr[1][1]-xr[0][1]),
              fabs(xr[1][2]-xr[0][2]) )

  g[0] /= len(coords)
  g[1] /= len(coords)
  g[2] /= len(coords)

  if box<5:
    box=ceil(rmax)
    box/=apix
    box=4*((box-1)/4+2) # round up to nearest 4 then add 4
    box*=2
  box = int(box)
  
  # the constant a for a Gaussian y=exp(-a*x^2) is
  #	a=ln(2)/res^2	for res defined at y=0.5
  #	a=1/res^2	for res defined at y=1/e
  #
  #  here we use the res defined in Fourier space: eg. Fourier transform
  #  of the above function: Y=exp(-pi^2*k^2/a)
  # a=ln(2)*pi^2/res^2	for res defined at Y=0.5
  # a=pi^2/res^2	for res defined at Y=1/e is used in this program
  
  rp=res/apix
  rp=sqrt(pi/rp)  # constant for gaussian falloff
  kn=(rp/pi)**1.5
  w=round(res*3.0/apix)
##   if w<3:
##     print "insufficient sampling for this resolution. Decrease apix."
##     return

  d = numpy.zeros( (box, box, box), 'f')
  if solv:
    rd = (numpy.ones( (box, box, box) )*100.).astype('f')

  b2 = box/2.
  for c,e in zip(coords, elist):
    if e== -1:
      continue
    x,y,z = c
    x-=g[0]
    y-=g[1]
    z-=g[2]
    xx=(x/apix)+b2
    yy=(y/apix)+b2
    zz=(z/apix)+b2
    xmin=int(round(xx-w))
    xmax=int(round(xx+w))
    ymin=int(round(yy-w))
    ymax=int(round(yy+w))
    zmin=int(round(zz-w))
    zmax=int(round(zz+w))

    if (xmin>=box or ymin>=box or zmin>=box or xmax<0 or ymax<0 or zmax<0):
      continue

    if xmin<0: xmin=0
    if xmax>=box: xmax=box-1
    if ymin<0: ymin=0
    if ymax>=box: ymax=box-1
    if zmin<0: zmin=0
    if zmax>=box: zmax=box-1

    cnst = kn*e
    for k in range(zmin, zmax):
      for j in range(ymin, ymax):
        for i in range(xmin, xmax):
	  r= (i-xx)*(i-xx) + (j-yy)*(j-yy) + (k-zz)*(k-zz)
          if solv:
            if rd[i][j][k]>r:
              rd[i][j][k] = r
          d[i][j][k] += cnst*exp(-r*rp)

  origin = ( g[0]-b2*apix, g[1]-b2*apix, g[2]-b2*apix)

  return d, origin, apix


if __name__=='__main__':
  from MolKit import Read
  mol = Read('cv.pdb')
  atoms = mol.chains.residues.atoms
  from MolKit.chargeMass import getChargeMass
  elist, totmass = getChargeMass(atoms)
  volarr, origin, boxsize = synthMap(atoms.coords, elist, totmass)
