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
from string import split

class ReadGamessOrbitals:
    """Read a rawiv binary file"""
    
    def mapArray(self, values, mini, maxi):

        # allocate numeric array that is **2
        dim = int(round(len(values)**(1./3)))
        print dim, dim**3, len(values)
        if dim**3!=len(values):
            raise ValueError, "Length of array is not **3"
        dim2 = 2
        while dim2<dim:
            dim2 = dim2*2
        print 'dim: padding from:', dim, 'to:', dim2
        intval = numpy.zeros( (dim2,dim2,dim2), numpy.uint8 )

        valuesN = numpy.array(values)
        
        minim = min(values)
        maxim = max(values)
        scale = (maxi - mini)/float(maxim - minim)
        valuesN.shape = (dim, dim, dim)

        for i in range(dim):
            face = valuesN[i]
            for j in range(dim):
                line = face[j]
                for k in range(dim):
                    intval[i][j][k] = ((line[k]-minim)*scale) + mini
        return intval, minim, maxim


    def read(self, filename, normalize=True):
        
        myfile = open(filename)
        data = myfile.readlines()
        myfile.close()
        self.header = {}

        dataspl = map( split, data )
        values = []
        for l in dataspl:
            values.extend( map(float, l) )

        # Fixme this should be done in separate nodes
        self.data, mini, maxi = self.mapArray(values, 1, 255)
        w, h, d = self.data.shape
        self.header['width'] = w
        self.header['height'] = h
        self.header['depth'] = d
        self.header['minif'] = mini
        self.header['maxif'] = maxi

        return self.header, self.data
