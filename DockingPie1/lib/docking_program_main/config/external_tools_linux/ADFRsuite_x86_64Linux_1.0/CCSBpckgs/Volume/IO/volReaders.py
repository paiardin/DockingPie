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

import math
from struct import pack, unpack, calcsize
import numpy
from Volume.Grid3D import Grid3DUC, Grid3DSI, Grid3DF, Grid3DD
from mglutil.math.crystal import Crystal
from mglutil.math.stats import stats
from time import time
import warnings

class VolumeReaderBase:

    def getDataElemSize(self, mode):
        """For a given mode, returns the size in number of bytes of a data
        elements, the character to be used in the unpack format and the
        Numeric array type"""
        
        # FIXME should use sizeof
        
        if mode==0:
            size = 1 #sizeof(char) *nDataVar;
            unpackType = 'c'
            arraytype = numpy.int8
            gridType = Grid3DUC
        elif mode==1:
            size = 2 #sizeof(short)*nDataVar;
            unpackType = 'h'
            arraytype = numpy.int16
            gridType = Grid3DSI
        elif mode==2:
            size = 4 #sizeof(float)*nDataVar;
            unpackType = 'f'
            arraytype = numpy.float32
            gridType = Grid3DF
        elif mode==3:
            size = 2 #sizeof(short)*nDataVar;
            unpackType = 'h'
            arraytype = numpy.int16
            gridType = Grid3DSI            
        elif mode==4:
            size = 4 #sizeof(float)*nDataVar;
            unpackType = 'f'
            arraytype = numpy.float32
            gridType = Grid3DF
        return size, unpackType, arraytype, gridType

    def getGRDElemSize(self,mode):
        """For GRD maps, returns the size in number of bytes of a data
        elements, the character to be used in the unpack format and the
        Numeric array type"""

        if mode==0:
            size = 1
            unpackType = 'c'
            arraytype = numpy.int8
            gridType = Grid3DUC
        elif mode==1:  # unsigned char
            size = 1
            unpackType = 'B'
            arraytype = numpy.uint8
            gridType = Grid3DUC
        elif mode==2:
            size = 1
            unpackType = 'b'
            arraytype = numpy.int8
            gridType = Grid3DUC
        elif mode==3:
            size = 2
            unpackType = 'H'
            arraytype = numpy.uint16
            gridType = Grid3DSI
        elif mode==4:
            size = 2
            unpackType = 'h'
            arraytype = numpy.int16
            gridType = Grid3DSI
        elif mode==5:
            size = 2
            unpackType = 'I'
            arraytype = numpy.uint16
            gridType = Grid3DSI
        elif mode==6:
            size = 2
            unpackType = 'i'
            arraytype = numpy.int16
            gridType = Grid3DSI
        elif mode==7:
            size = 4
            unpackType = 'L'
            arraytype = numpy.uint32
            gridType = Grid3DSI
        elif mode==8:
            size = 4
            unpackType = 'l'
            arraytype = numpy.int32
            gridType = Grid3DSI
        elif mode==9:
            size = 4
            unpackType = 'f'
            arraytype = numpy.float32
            gridType = Grid3DF
        elif mode==10:
            size = 8
            unpackType = 'd'
            arraytype = numpy.float64
            gridType = Grid3DD
        return size, unpackType, arraytype, gridType

    def describe(self, hd):
        strg = "%s file: %s\n"%(hd['mapType'], hd['filename'])
        if hd.has_key('mode'): strg += "mode: %d\n"%hd['mode']
        if hd.has_key('mapc'):
            strg += "indices (1:fastest, 3:slowest) %d %d %d\n"%(
                hd['mapc'], hd['mapr'], hd['maps'] )
        strg += "nc: %d, nr: %d, ns: %d\n"%(
            hd['nc'], hd['nr'], hd['ns'])
        strg += "dimensions: %s\n"%hd['dimensions']
        strg += "origin: %s, %s, %s\n"%(
            hd['origin'][0], hd['origin'][1], hd['origin'][2])
        if hd.has_key('ccell'):
            strg += "acell: %10.6f, bcell: %10.6f, ccell: %10.6f\n"%(
                hd['acell'], hd['bcell'], hd['ccell'])
            strg += "alpha: %10.6f, beta: %10.6f, gamma: %10.6f\n"%(
                hd['alpha'], hd['beta'], hd['gamma'])
        if hd.has_key('nsymbt'):
            strg += "symmetry ops: %d bytes in length\n"%hd['nsymbt']
        if hd.has_key('lskflg'):
            strg += "Skew transformation Flag = %d\n"%hd['lskflg']
        strg += "min:  %10.6f\nmax:   %10.6f\nmean: %10.6f\nstdev: %10.6f\n"%(
            hd['amin'], hd['amax'], hd['amean'], hd['arms'])
        return strg
    
class ReadCCP4(VolumeReaderBase):
    """Read a CCP4 file and return a Grid3D object
This parser handles mapc, mapr and maps. There is a performance difference
depending on on the values of mapc mapr and maps.

1 - C-style (i.e. mapc=3 mapr=2 maps=1) is the fastest and is also the ordering
    used by write CCP4
2 - Fortran style (i.e. mapc=1 mapr=2 maps=3) is also efficient although
    not as fast as C-style and will cause slow writting to file
3 - other ordering such as results or cutting maps with mapin are the slowest
    (factor of 3 to 4 compared to C-style)

A crystal attribute is added to the grid object to enable conversion between
cartesian and fractional space.  The origin ans stepSize are kept in fractional
space since isocontouring and volume rendering cannot handle embedded data.

By default the data is normalized.
"""
    
    def read(self, filename, normalize=True, disp_out=True):
        t1 = time()
        sw = "" # used in format to swap bytes if necessary
        self.filename = filename
        
        # open file to read in binary mode
        f = open(filename, 'rb')
        # read the 1024 bytes of the header
        data = f.read(1024)
        h = self.header = {}
        h['mapType']='CCP4'
        h['filename'] = filename
        
        nc, nr, ns, mode = unpack("4i", data[:16])
        #print "nc=", nc, "nr=", nr, "ns=", ns, "mode=", mode
        
        # find out about byteswapping (mode has to be 0-5 or 10)
        # if mode is not in 0-5 or 10 then we have to swap
        # first determine if byte order is big-endian
        if mode not in (0,1,2,3,4,5,10):
            sw=">"
            nc, nr, ns, mode = unpack(sw+"4i", data[:16])

        # if not big-endian, check if little-endian
        if mode not in (0,1,2,3,4,5,10):
            sw="<"
            nc, nr, ns, mode = unpack(sw+"4i", data[:16])

        # else this is an unreadable ccp4 file
        if mode not in (0,1,2,3,4,5,10):
            f.close()
            print "ReadCCP4: ERROR: %s is not a valid CCP4 file"%filename
            return
            
        # Collect header data:
        h.update( {'nc':nc, 'nr':nr, 'ns':ns, 'mode':mode } )
        
        ncstart,nrstart,nsstart,nx,ny,nz = unpack(sw+"6i", data[16:40])
        # Note: nx,ny,nz represent grid sampling

        #print "ncstart=",ncstart, "nrstart=", nrstart ,"nsstart=", nsstart , "nx:", nx, "ny=" ,ny,"nz=", nz
        h.update( {'ncstart':ncstart,'nrstart':nrstart,'nsstart':nsstart })
        h.update( {'nx':nx,'ny':ny,'nz':nz })

        acell,bcell,ccell,alpha, beta, gamma = unpack(sw+"6f", data[40:64])
        h.update( {'acell':acell,'bcell':bcell,'ccell':ccell })
        h.update( {'alpha':alpha,'beta':beta,'gamma':gamma })
        
        mapc, mapr, maps = unpack(sw+"3i", data[64:76])
        h.update( {'mapc':mapc,'mapr':mapr,'maps':maps })

        amin, amax, amean = unpack(sw+"3f", data[76:88])
        h.update( {'amin':amin,'amax':amax,'amean':amean })

        ispg, nsymbt = unpack(sw+"2i", data[88:96])
        h.update( {'ispg':ispg,'nsymbt':nsymbt})

        lskflg = unpack(sw+"1i", data[96:100])[0]
        h.update( {'lskflg':lskflg})

        skwmat = unpack(sw+"9f", data[100:136])
        h.update( {'skwmat':skwmat})
        
        skwtrn = unpack(sw+"3f", data[136:148])
        h.update( {'skwtrn':skwtrn})
        
        future_words = unpack(sw+"15i", data[148:208])
        h.update( {'future_words':future_words })

        # we don't use this, but put it here for completion
        mapstring = unpack(sw+"4c", data[208:212])

        # Machine stamp
        machst = unpack(sw+"4c", data[212:216])

        arms = unpack(sw+"1f", data[216:220])[0]
        h.update( {'arms':arms })
        
        nlabl = unpack(sw+"i", data[220:224])[0]
        labl = unpack(sw+"800c", data[224:1025])

        size, unpackType, arraytype, gtype = self.getDataElemSize(mode)

        # read the data
        # sometimes the symmetry info is not included in the file,
        # even though the header dictates the # of bytes allocated for it,
        # so it is safer to read the data from the end of the file
        data = f.read(  ) # read to end
        f.close()
        #print 'size difference', nc*nr*ns*size, len(data)
        
        ndata = numpy.array( unpack(sw+"%d%c"%(nc*nr*ns,unpackType),
                                      data[-nc*nr*ns*size:]), arraytype )

        self.data = ndata
        
        # compute min, max, mean rms and replace values if needed
        mymin, mymax, mymean, mystdev = stats(ndata)

        if mymean != amean:
            amean = mymean
            h['amean'] = mymean
        if mystdev != arms:
            arms = mystdev
            h['arms'] = mystdev
        if mymin != amin:
            amin = mymin
            h['amin'] = mymin
        if mymax != amax:
            amax = mymax
            h['amax'] = mymax

        #print 'indices speeds: mapc mapr maps', mapc, mapr, maps, nc, nr, ns
        cmapc = mapc-1
        cmapr = mapr-1
        cmaps = maps-1
        
        if (mapc==1 and mapr==2 and maps==3):  #fortran style x-fastest
            #print 'FORTRAN style'
            ndata.shape = (ns,nr,nc)
            #transpose the scalar data due to FORTRAN style
            self.data = numpy.ascontiguousarray(numpy.transpose(ndata), ndata.dtype.char)
            if normalize:
                tcode = self.data.dtype.char
                self.data = ((self.data-amean)/arms).astype(tcode)
        elif (mapc==3 and mapr==2 and maps==1):  #C style z-fastest
            #print 'C style'
            ndata.shape = (ns,nr,nc)
            if normalize:
                tcode = self.data.dtype.char
                self.data = ((self.data-amean)/arms).astype(tcode)
        else:
            #print 'Generic style'
            dims = [0,0,0]
            dims[cmapc]=nc
            dims[cmapr]=nr
            dims[cmaps]=ns
            nndata = numpy.zeros( dims, 'f')
            cc = [0,0,0]
            l = 0
            for i in range(ns):
                cc[cmaps] = i
                for j in range(nr):
                    cc[cmapr] = j
                    for k in range(nc):
                        cc[cmapc] = k
                        val = ndata.flat[l]
                        if normalize:
                            val = (val-amean)/arms
                        nndata[cc[0], cc[1], cc[2]] = val
                        l+=1
            self.data = nndata

        if normalize:
            h['amean'] = 0.0
            h['arms'] = 1.0
            h['amin'] = (h['amin']-amean)/arms
            h['amax'] = (h['amax']-amean)/arms

        # we compute the origin and stepSize in Crystal space coordinates
        # the isocontour operates in this space and the result of
        # iso-contouring has then to be converted to real space
        
        #nx,ny,nz represent grid sampling
        assert nx != 0
        assert ny != 0
        assert nz != 0
        dims = (nx, ny, nz)

##         dims = [0,0,0]
##         dims[cmapc]=nc
##         dims[cmapr]=nr
##         dims[cmaps]=ns
        #print 'dims', dims
        start = [0.,0.,0.]
        start[cmapc] = ncstart
        start[cmapr] = nrstart
        start[cmaps] = nsstart
        #print 'start', start
        origin = [0.,0.,0.]
        origin[cmapc] = start[cmapc]/float(dims[cmapc])
        origin[cmapr] = start[cmapr]/float(dims[cmapr])
        origin[cmaps] = start[cmaps]/float(dims[cmaps])
        h['origin'] = origin 

        stepSize = [0.,0.,0.]
        stepSize[cmapc] = 1./float(dims[cmapc])
        stepSize[cmapr] = 1./float(dims[cmapr])
        stepSize[cmaps] = 1./float(dims[cmaps])
        #print 'stepSize', stepSize

        crystal = Crystal( (acell, bcell, ccell), (alpha, beta, gamma))
        g = gtype( self.data, origin, stepSize, h, crystal)
        h['dimensions'] = str(self.data.shape)

        if disp_out is True:
            print self.describe(h)
            print 'time to read file: ', time()-t1
        
        return g
    

    def info(self):
        txt = """http://www.ccp4.ac.uk/dist/html/maplib.html#description
1) INTRODUCTION

The standard map file format used by the CCP4 programs is the map/image file
format devised at the MRC LMB Cambridge, originally by David Agard. The
advantages of the format include the following:

a) The information in the header describes the relationship of the map to the
crystal cell, and other
information useful in crystallographic calculations (eg crystal symmetry) 

b) The file may be written and read in different ways e.g. written section by
section and read line by line. 

c) The format is suitable for both crystallographic work and for image
processing so that Fourier and plotting programs can be used for both
purposes. 

2) DETAILED DESCRIPTION OF THE MAP FORMAT
The overall layout of the file is as follows: 

a) File header (256 longwords)
b) Symmetry information
c) Map, stored as a 3-dimensional array 

The files are read & written using the diskio package, which allows the file
to be treated as a direct-access byte stream, essentially as in C fread &
fwrite calls. 

The header is organised as 56 words followed by space for ten 80 character
text labels as follows: 

 
 1      NC              # of Columns    (fastest changing in map)
 2      NR              # of Rows
 3      NS              # of Sections   (slowest changing in map)
 4      MODE            Data type
                          0 = envelope stored as signed bytes (from
                              -128 lowest to 127 highest)
                          1 = Image     stored as Integer*2
                          2 = Image     stored as Reals
                          3 = Transform stored as Complex Integer*2
                          4 = Transform stored as Complex Reals
                          5 == 0        
 
                          Note: Mode 2 is the normal mode used in
                                the CCP4 programs. Other modes than 2 and 0
                                may NOT WORK
 
 5      NCSTART         Number of first COLUMN  in map
 6      NRSTART         Number of first ROW     in map
 7      NSSTART         Number of first SECTION in map
 8      NX              Number of intervals along X
 9      NY              Number of intervals along Y
10      NZ              Number of intervals along Z
11      X length        Cell Dimensions (Angstroms)
12      Y length                     ''
13      Z length                     ''
14      Alpha           Cell Angles     (Degrees)
15      Beta                         ''
16      Gamma                        ''
17      MAPC            Which axis corresponds to Cols.  (1,2,3 for X,Y,Z)
18      MAPR            Which axis corresponds to Rows   (1,2,3 for X,Y,Z)
19      MAPS            Which axis corresponds to Sects. (1,2,3 for X,Y,Z)
20      AMIN            Minimum density value
21      AMAX            Maximum density value
22      AMEAN           Mean    density value    (Average)
23      ISPG            Space group number
24      NSYMBT          Number of bytes used for storing symmetry operators
25      LSKFLG          Flag for skew transformation, =0 none, =1 if foll
26-34   SKWMAT          Skew matrix S (in order S11, S12, S13, S21 etc) if
                        LSKFLG .ne. 0.
35-37   SKWTRN          Skew translation t if LSKFLG .ne. 0.
                        Skew transformation is from standard orthogonal
                        coordinate frame (as used for atoms) to orthogonal
                        map frame, as
 
                                Xo(map) = S * (Xo(atoms) - t)
 
38      future use       (some of these are used by the MSUBSX routines
 .          ''            in MAPBRICK, MAPCONT and FRODO)
 .          ''   (all set to zero by default)
 .          ''
52          ''

53      MAP             Character string 'MAP ' to identify file type
54      MACHST          Machine stamp indicating the machine type
                        which wrote file
55      ARMS            Rms deviation of map from mean density
56      NLABL           Number of labels being used
57-256  LABEL(20,10)    10  80 character text labels (ie. A4 format)

Symmetry records follow - if any - stored as text as in International Tables,
operators separated by * and grouped into 'lines' of 80 characters (ie.
symmetry operators do not cross the ends of the 80-character 'lines' and the
'lines' do not terminate in a *). 

Map data array follows.\n"""
        print txt


class ReadMRC(VolumeReaderBase):
    """Read an MRC file."""

    def read(self, filename, disp_out=True, normalize=True):

        t1 = time()
        sw = "" # used in format to swap bytes if necessary
        self.filename = filename

        # open file to read in binary mode
        f = open(filename, 'rb')
        # read the 1024 bytes of the header
        data = f.read(1024)
        h = self.header = {}
        h['mapType'] = 'MRC'
        h['filename'] = filename

        nc, nr, ns, mode = unpack("4i", data[:16])
        
        # find out about byteswapping (mode has to be 0-4)
        # first determine if byte order is big-endian
        if mode<0 or mode>4:
            sw=">"
            nc, nr, ns, mode = unpack(sw+"4i", data[:16])

        # if not big-endian, check if little-endian
        if mode<0 or mode>4:
            sw="<"
            nc, nr, ns, mode = unpack(sw+"4i", data[:16])

        # else this is an unreadable ccp4 file
        if mode<0 or mode>4:
            f.close()
            print "ReadMRC: ERROR: %s is not a valid MRC file"%filename
            return

        #Collect header data
        h.update( {'nc':nc, 'nr':nr, 'ns':ns, 'mode':mode } )

        ncstart,nrstart,nsstart,mx,my,mz = unpack(sw+"6i", data[16:40])
        h.update( {'ncstart':ncstart,'nrstart':nrstart,'nsstart':nsstart })
        h.update( {'mx':mx,'my':my,'mz':mz })

        acell,bcell,ccell,alpha, beta, gamma = unpack(sw+"6f", data[40:64])
        h.update( {'acell':acell,'bcell':bcell,'ccell':ccell })
        h.update( {'alpha':alpha,'beta':beta,'gamma':gamma })

        mapc, mapr, maps = unpack(sw+"3i", data[64:76])
        h.update( {'mapc':mapc,'mapr':mapr,'maps':maps })

        amin, amax, amean = unpack(sw+"3f", data[76:88])
        h.update( {'amin':amin,'amax':amax,'amean':amean })

        ispg, nsymbt = unpack(sw+"2i", data[88:96])
        h.update( {'ispg':ispg,'nsymbt':nsymbt})

        extra = unpack(sw+"25i", data[96:196])
        h.update( {'extra':extra })

        xorigin, yorigin, zorigin = unpack(sw+"3f", data[196:208])
        h.update( {'xorigin':xorigin,'yorigin':yorigin, 'zorigin':zorigin})

        map = unpack(sw+"4c", data[208:212])
        machinestamp = unpack(sw+"i", data[212:216])[0]
        arms = unpack(sw+"f", data[216:220])[0]
        nlabl = unpack(sw+"i", data[220:224])[0]
        label = unpack(sw+"800c", data[224:1024])
        h.update( {'map':map,'machinestamp':machinestamp, 'arms':arms,
                   'nlabl':nlabl,'label':label})

        size, unpackType, arraytype, gtype = self.getDataElemSize(mode)
        # read the data
        #data = f.read( nx*ny*nz*size )
        data = f.read(  ) # read to end
        f.close()
        ndata = numpy.array( unpack(sw+"%d%c"%(nc*nr*ns,unpackType),
                                          data), arraytype )

        self.data = ndata
        
        mymin, mymax, mymean, mystdev = self.data.min(), self.data.max(), self.data.mean(), self.data.std()
#        print mymin, mymax,mymean,mystdev
        if mymean != amean:
#            print "mean changed from %f to %f"%(amean,mymean)
            amean = mymean
            h['amean'] = mymean
        if mystdev != arms:
#            print "rms changed from %f to %f"%(arms,mystdev)
            arms = mystdev
            h['arms'] = mystdev
        if mymin != amin:
#            print "min changed from %f to %f"%(amin,mymin)
            amin = mymin
            h['amin'] = mymin
        if mymax != amax:
#            print "max changed from %f to %f"%(amax,mymax)
            amax = mymax
            h['amax'] = mymax

#        print 'indices speeds: mapc mapr maps', mapc, mapr, maps, nc, nr, ns
        cmapc = mapc-1
        cmapr = mapr-1
        cmaps = maps-1

        if (mapc==1 and mapr==2 and maps==3):  #fortran style x-fastest
            #print 'FORTRAN style'
            ndata.shape = (ns,nr,nc)
            self.data = numpy.ascontiguousarray(numpy.transpose(ndata), ndata.dtype.char)
            if normalize:
                tcode = self.data.dtype.char
                self.data = ((self.data-amean)/arms).astype(tcode)
        elif (mapc==3 and mapr==2 and maps==1):  #C style z-fastest
            #print 'C style'
            ndata.shape = (ns,nr,nc)
            if normalize:
                tcode = self.data.dtype.char
                self.data = ((self.data-amean)/arms).astype(tcode)
        else:
            #print 'Generic style'
            dims = [0,0,0]
            dims[cmapc]=nc
            dims[cmapr]=nr
            dims[cmaps]=ns
            nndata = numpy.zeros( dims, 'f')
            cc = [0,0,0]
            l = 0
            for i in range(ns):
                cc[cmaps] = i
                for j in range(nr):
                    cc[cmapr] = j
                    for k in range(nc):
                        cc[cmapc] = k
                        val = ndata.flat[l]
                        if normalize:
                            val = (val-amean)/arms
                        nndata[cc[0], cc[1], cc[2]] = val
                        l+=1
            self.data = nndata

        if normalize:
            h['amean'] = 0.0
            h['arms'] = 1.0
            h['amin'] = (h['amin']-amean)/arms
            h['amax'] = (h['amax']-amean)/arms


        # we compute the origin and stepSize in Crystal space coordinates
        # the isocontour operates in this space and the result of
        # iso-contouring has then to be converted to real space

        dims = (nc, nr, ns)
        #print 'dims', dims
        center = [0,0,0]
        center[cmapc] = xorigin
        center[cmapr] = yorigin
        center[cmaps] = zorigin
        
        start = [0.,0.,0.]
        start[cmapc] = ncstart
        start[cmapr] = nrstart
        start[cmaps] = nsstart
        #print 'start', start
        origin = [0.,0.,0.]
        origin[cmapc] = (start[cmapc]-center[cmapc])/float(dims[cmapc])
        origin[cmapr] = (start[cmapr]-center[cmapr])/float(dims[cmapr])
        origin[cmaps] = (start[cmaps]-center[cmaps])/float(dims[cmaps])
        h['origin'] = origin
        #print 'origin', origin

        stepSize = [0.,0.,0.]
        stepSize[cmapc] = 1./float(dims[cmapc])
        stepSize[cmapr] = 1./float(dims[cmapr])
        stepSize[cmaps] = 1./float(dims[cmaps])
        #print 'stepSize', stepSize

        crystal = Crystal( (acell, bcell, ccell), (alpha, beta, gamma))
        g = gtype( self.data, origin, stepSize, h, crystal)
        h['dimensions'] = str(self.data.shape)
        
        if disp_out is True:
            print self.describe(h)
            print 'time to read file: ', time()-t1
        
        return g
    
import re
import array
class ReadCNS(VolumeReaderBase):
    """Read a CNS/XPLOR file and return a Grid3D object"""

    def read(self, filename, disp_out=True, normalize=True):
        
        t1 = time()
        sw = "" # used in format to swap bytes if necessary
        self.filename = filename
        h = self.header = {}
        h['mapType'] = 'CNS/XPLOR'
        h['filename'] = filename

        # open file to read in ascii mode
        f = open(filename, 'r')
        numOfLines=len(f.readlines())
        f.close
        # read the header lines
        #  content = f.readlines()
        count = 0
        f = open(filename, 'r')
        cns_file_error = False
        # determine location of cell edges in header
        i=0
        while i < numOfLines:
            line=f.readline()
            # gather grid information
            if re.search("^(\s+-?\d+){9}", line):
                edges=re.split("\s+", line)
                na=int(edges[1])
                amin=int(edges[2])
                amax=int(edges[3])
                nb=int(edges[4])
                bmin=int(edges[5])
                bmax=int(edges[6])
                nc=int(edges[7])
                cmin=int(edges[8])
                cmax=int(edges[9])

                # calculate edge lengths
                alen, blen, clen = (amax-amin+1, bmax-bmin+1, cmax-cmin+1)
                break
            i+=1
        if i==numOfLines: cns_file_error=True
        
        # gather crystal cell dimensions
        cryst_cell=f.readline()
        if re.search("^([\s-]\d.+){6}",cryst_cell):
            acell=float(cryst_cell[:12])
            bcell=float(cryst_cell[12:24])
            ccell=float(cryst_cell[24:36])
            alpha=float(cryst_cell[36:48])
            beta=float(cryst_cell[48:60])
            gamma=float(cryst_cell[60:72])
            #compute interval step size
            sx=1./(na-1)
            sy=1./(nb-1)
            sz=1./(nc-1)
        else: cns_file_error=True

        # the following line should contain "ZYX"
        if not re.search("ZYX",f.readline()): cns_file_error=True
        if cns_file_error is True:
            print "ReadCNS: ERROR: %s is not a valid CNS/XPLOR file"%filename
            return

        h.update({'nc':alen, 'nr':blen, 'ns':clen, 'mode':2})
        h.update({'mapc':1, 'mapr':2,'maps':3})
        h.update({'ncstart':amin, 'nrstart':bmin, 'nsstart':cmin })
        h.update({'acell':acell, 'bcell':bcell, 'ccell':ccell})
        h.update({'alpha':alpha, 'beta':beta, 'gamma':gamma})

        # allocate array for density
        ndata = numpy.zeros((alen, blen, clen), 'f')

        # loop over sections, then lines, then columns and read
        # a new line everytime we read 6 float (each stored in 12 charaters)
        # after each section (i.e. XY plane) there is a 1 line section header
        for z in range(clen):
            # skip section header line
            line = f.readline()

            # read first line
            line = f.readline()
            ci = 0 # initialize character counter on this line
            for y in range(blen):
                if ci==72:  # if we read all the float
                    line = f.readline() # read next line
                    ci = 0              # and reset character pointer
                for x in range(alen):
                    if ci==72:  # if we read all the float
                        line = f.readline() #  read next line
                        ci = 0              #  and reset character pointer
                    # put the denity intothe array
                    ndata[x,y,z] = float(line[ci:ci+12])
                    ci+=12 # increment the character pointer

        self.data = ndata
        
        line = f.readline() # footer - int value must be -9999
        assert float(line)==-9999, "penultimate line should hold -9999', got %d"%line
        line = f.readline() # density average and standard dev
        f.close()
        mean, stddev = map(float, line.split())
        h['amean'] = mean
        h['arms'] = stddev

        dims = (alen, blen, clen)

        #print 'dims', dims
        start = [0.,0.,0.]
        start[0] = amin
        start[1] = bmin
        start[2] = cmin
        #print 'start', start
        origin = [amin*sx, bmin*sy, cmin*sz]
        h['origin'] = origin

        #print 'origin', origin
        stepSize = [sx,sy,sz]
        
        size, unpackType, arraytype, gtype = self.getDataElemSize(2)
        crystal = Crystal( (acell, bcell, ccell), (alpha, beta, gamma))
        g = gtype(self.data, origin, stepSize, h, crystal)       

        # compute min, max, mean rms
        amin, amax, amean, arms = g.stats()
        h.update( {'amin':amin, 'amax':amax, 'amean':amean, 'arms':arms })

        h['dimensions'] = str(self.data.shape)

        if disp_out is True:
            print self.describe(h)
            print 'time: ', time()-t1

        return g

    def info(self):
        txt = """http://www.sinica.edu.tw/~scimath/msi/xplor981/formats.html
        The X-PLOR map file begins with an eight-line header.
1. Line 1

An empty line written by the `/ ` FORTRAN format
descriptor in the formatted map file.


2. Lines 2- 5

Title information written as character strings.
These lines are written as 80-character strings
in the formatted file map.


3. Line 6

A series of nine integers:
NA, AMIN, AMAX, NB, BMIN, BMAX, NC, CMIN, CMAX
The values NA, NB and NC indicate the total number
of grid points along the a,b, and c cell edges.
The items AMIN, AMAX, BMIN, BMAX, CMIN, CMAX
indicate the starting and stopping grid points
along each cell edge in the portion of the map that
is written. In the formatted map file this line is
written using the FORTRAN format statement (9I8).


4. Line 7

A series of six double-precision items corresponding to
the crystal cell dimensions a, b, c, alpha, beta, gamma.
In the formatted map file these items are written using
the FORTRAN format statement (6E12.5).


5. Line 8

A three-letter character string which always reads `ZXY'.
"""
        print txt

class ReadGRD(VolumeReaderBase):
    """Read a GRD file and return a Grid3D object"""
    
    def read(self, filename, normalize=True, disp_out=True):
        
        t1 = time()
        sw = "" # used in format to swap bytes if necessary
        self.filename = filename
        
        # open file to read in binary mode
        f = open(filename, 'rb')
        # read the 512 bytes of the header
        data = f.read(512)
        h = self.header = {}
        h['mapType'] = 'GRD'
        h['filename'] = filename
        
        file_num, processor, mode = unpack("3i", data[:12])
        
        # find out about byteswapping (mode has to be 1-3)
        # if mode is not in 0-3 then we have to swap
        # first determine if byte order is big-endian
        if processor not in (1,2,3):
            sw=">"
            file_num, processor, mode = unpack(sw+"3i", data[:12])

        # if not big-endian, check if little-endian
        if processor not in (1,2,3):
            sw="<"
            file_num, processor, mode = unpack(sw+"3i", data[:12])

        # else this is an unreadable GRD file
        if processor not in (1,2,3):
            f.close()
            print "ReadGRD: ERROR: %s is not a valid GRD file"%filename
            return
            
        offset, nc, nr, ns  = unpack(sw+"4i", data[12:28])
        cmapc, cmapr, cmaps = unpack(sw+"3i", data[28:40])
        mapc = cmapc+1
        mapr = cmapr+1
        maps = cmaps+1
        xtra = unpack(sw+"472c", data[40:512])
        
        h.update( {'file_num':file_num,'processor':processor })
        h.update( {'mode':mode })
        h.update( {'offset':offset,'nc':nc,'nr':nr,'ns':ns })
        h.update( {'mapc':mapc,'mapr':mapr,'maps':maps })
        h.update( {'cmapc':cmapc,'cmapr':cmapr,'cmaps':cmaps })

        size, unpackType, arraytype, gtype = self.getGRDElemSize(mode)
        
        data = f.read(  ) # read to end
        f.close()

        ndata = numpy.array( unpack(sw+"%d%c"%(nc*nr*ns,unpackType),
                                      data[:nc*nr*ns*size]), arraytype )

        self.data = ndata

        if (mapc==1 and mapr==2 and maps==3):  #fortran style x-fastest
            #print 'FORTRAN style'
            ndata.shape = (ns,nr,nc)
            #transpose the scalar data due to FORTRAN style
            self.data = numpy.ascontiguousarray(numpy.transpose(ndata), ndata.dtype.char)
        elif (mapc==3 and mapr==2 and maps==1):  #C style z-fastest
            #print 'C style'
            ndata.shape = (ns,nr,nc)
        else:
            #print 'Generic style'
            dims = [0,0,0]
            dims[cmapc]=nc
            dims[cmapr]=nr
            dims[cmaps]=ns
            nndata = numpy.zeros( dims, 'f')
            cc = [0,0,0]
            l = 0
            for i in range(ns):
                cc[cmaps] = i
                for j in range(nr):
                    cc[cmapr] = j
                    for k in range(nc):
                        cc[cmapc] = k
                        val = ndata.flat[l]
                        nndata[cc[0], cc[1], cc[2]] = val
                        l+=1
            self.data = nndata

        dims = (nc, nr, ns)
        origin = [0.,0.,0.]
        stepSize = [1.,1.,1.]
        h['origin'] = origin

        g = gtype( self.data, origin, stepSize, h)

        # compute min, max, mean rms
        amin, amax, amean, arms = g.stats()
        h.update( {'amin':amin, 'amax':amax, 'amean':amean, 'arms':arms })

        h['dimensions'] = str(self.data.shape)

        if disp_out is True:
            print self.describe(h)
            print 'time: ', time()-t1
        
        return g

    def info(self):
        txt = """
        GRD is a general purpose grid file format developed by C. Henn and
        R. Buerki at the M.E. Mueller Institute, Basel. It contains a 
        default 512 byte header containing 64 four byte records. Only
        header records 0-10 are occupied so far.

                [0]     magic file number
                [1]     processor where the file was generated:
                                1 = Data_from_VAX
                                2 = Data_from_MIPS
                                3 = Data_from_Convex

                [2]     data record format:
                                0  = Data_in_free_fromat
                                1  = Data_in_u_char
                                2  = Data_in_char
                                3  = Data_in_u_short
                                4  = Data_in_short
                                5  = Data_in_u_int
                                6  = Data_in_int
                                7  = Data_in_u_long
                                8  = Data_in_long
                                9  = Data_in_float
                                10 = Data_in_double

                [3]     offset of data section after default header

                [4]     number of columns
                [5]     number of rows
                [6]     number of sections
                
                [8]     fastest coordinate
                [9]     medium coordinate
                [10]    slowest coordinate

        The rest of the header section is unused.
"""
        print txt

class ReadBRIX(VolumeReaderBase):
    """Read a BRIX or MAPPAGE (DSN6) crystallographic density
       map used by O, and return a Grid3D object"""

    def read(self, filename, normalize=True, disp_out=True):
        
        t1 = time()
        sw = "" # used in format to swap bytes if necessary
        self.filename = filename
        
        # open file to read in binary mode
        f = open(filename, 'rb')

        # determine if the file is of type BRIX or DSN6 by looking
        # at the beginning of the header
        O_format = f.read(3)
        if O_format == ':-)':
            O_format = 'BRIX'
        else: O_format = 'DSN6'
        
        h = self.header = {}
        h['mapType'] = O_format
        h['filename'] = filename
        
        #READ DSN6 FORMATTED HEADER
        if O_format == 'DSN6':
            #go back to beginning of file
            f.seek(0,0)        

            # store header data
            data = f.read(512)

            prod_norm = unpack('h', data[36:38])[0]

            # find out about byteswapping (19th element should be 100)
            # first determine if byte order is big-endian
            if prod_norm != 100:
                sw=">"
                prod_norm = unpack(sw+'h', data[34:36])

            # if not big-endian, check if little-endian
            if prod_norm != 100:
                sw="<"
                prod_norm = unpack(sw+'h', data[34:36])

            # if neither, then this is an invalid file
            if prod_norm != 100:                
                print "ReadBRIX: ERROR: %s is not a valid DSN6 file"%filename
                return

            xorigin,yorigin,zorigin=unpack(sw+'3h', data[:6])
            h.update( {'xorigin':xorigin,'yorigin':yorigin, 'zorigin':zorigin})
            nc,nr,ns=unpack(sw+'3h', data[6:12])
            h.update( {'nc':nc, 'nr':nr, 'ns':ns })
            nx,ny,nz=unpack(sw+'3h', data[12:18])
            h.update( {'nx':nx,'ny':ny,'nz':nz })

            acell,bcell,ccell,alpha,beta,gamma=unpack(sw+'6h', data[18:30])
            prod,plus=unpack(sw+'2h', data[30:34])
            cell_norm=unpack(sw+'h', data[34:36])[0]

            #normalize cell parameters
            acell=acell/float(cell_norm)
            bcell=bcell/float(cell_norm)
            ccell=ccell/float(cell_norm)
            alpha=alpha/float(cell_norm)
            beta=beta/float(cell_norm)
            gamma=gamma/float(cell_norm)
            h.update( {'acell':acell,'bcell':bcell,'ccell':ccell })
            h.update( {'alpha':alpha,'beta':beta,'gamma':gamma })
            h.update( {'xorigin':xorigin,'yorigin':yorigin, 'zorigin':zorigin})

            #normalize prod
            prod=float(prod)/prod_norm
            h.update( {'prod':prod,'plus':plus })

        # READ BRIX FORMATTED HEADER
        if O_format == 'BRIX':
            headerData=f.read(509)
            info=re.search("^\s*origin\s+(\S+)\s+(\S+)\s+(\S+)\s+extent\s+(\d+)\s+(\d+)\s+(\d+)\s+grid\s+(\d+)\s+(\d+)\s+(\d+)\s+cell\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+prod\s+(\S+)\s+plus\s+(\d+)\s+sigma\s+(\S+)",headerData,re.I)
            if info:
                xorigin=int(info.group(1))
                yorigin=int(info.group(2))
                zorigin=int(info.group(3))
                nc=int(info.group(4))
                nr=int(info.group(5))
                ns=int(info.group(6))
                nx=int(info.group(7))
                ny=int(info.group(8))
                nz=int(info.group(9))
                acell=float(info.group(10))
                bcell=float(info.group(11))
                ccell=float(info.group(12))
                alpha=float(info.group(13))
                beta=float(info.group(14))
                gamma=float(info.group(15))
                prod=float(info.group(16))
                plus=int(info.group(17))
                sigma=float(info.group(18))
                h.update( {'xorigin':xorigin,'yorigin':yorigin, 'zorigin':zorigin})
                h.update( {'nc':nc, 'nr':nr, 'ns':ns })
                h.update( {'nx':nx,'ny':ny,'nz':nz })
                h.update( {'acell':acell,'bcell':bcell,'ccell':ccell })
                h.update( {'alpha':alpha,'beta':beta,'gamma':gamma })
                h.update( {'prod':prod,'plus':plus })
            else:
                print "ReadBRIX: ERROR: %s is not a valid BRIX file"%filename
                return
            
        size, unpackType, arraytype, gtype = self.getGRDElemSize(1)

        # allocate array for density
        ndata = numpy.zeros((ns, nr, nc), arraytype)
            
        #data is stored in 8x8x8 cubes, within
        #an INTEGER*2 array with 256 elements
        for zcube in range (int((math.ceil(float(ns)/8)))):
            zstart = zcube*8
            if zstart+8 > ns:
                zfwd = ns-zstart
            else: zfwd = 8
            for ycube in range (int((math.ceil(float(nr)/8)))):  
                ystart = ycube*8
                if ystart+8 > nr:
                    yfwd = nr-ystart
                else: yfwd = 8
                for xcube in range (int((math.ceil(float(nc)/8)))):
                    xstart = xcube*8
                    if xstart+8 > nc:
                        xfwd = nc-xstart
                    else: xfwd = 8
                    cubeData = f.read(512) #read this 8x8x8 cube of data
                    # if BRIX, no byteswapping necessary
                    if O_format == 'BRIX':
                        fixedData=(numpy.fromstring(cubeData, arraytype))
                    else: #do byteswapping
                        rawData = numpy.fromstring(cubeData, numpy.Int16)
                        swapData = rawData.byteswap().tostring()
                        fixedData = numpy.fromstring(swapData, arraytype)
                    matrix = numpy.reshape(fixedData, (8,8,8))
                    ndata[zstart:zstart+zfwd, ystart:ystart+yfwd, xstart:xstart+xfwd] = \
                          matrix[0:zfwd, 0:yfwd, 0:xfwd]
        f.close()                

        self.data=numpy.ascontiguousarray(numpy.transpose(ndata), ndata.dtype.char)
        
        sx = 1./(nx-1)
        sy = 1./(ny-1)
        sz = 1./(nz-1)
        stepSize = [sx,sy,sz]

        origin=[xorigin*sx,yorigin*sy,zorigin*sz]
        h['origin']=origin
        h['dimensions'] = str(self.data.shape)

        crystal = Crystal( (acell, bcell, ccell), (alpha, beta, gamma))
        g = gtype(self.data, origin, stepSize, h, crystal)       

        #compute stats
        amin, amax, amean, arms = g.stats()
        h.update( {'amin':amin,'amax':amax,'amean':amean,'arms':arms} )
    
        if disp_out is True:
            print self.describe(h)
            print 'time: ', time()-t1
        
        return g

    def info(self):
        txt = """http://www.uoxray.uoregon.edu/tnt/manual/node104.html
The format of the bricked files made by BRIX is almost the same
as the old MAPPAGE (DNS6) format, except that no byte swapping is required.

The header record

The first 512 bytes of the file is the header information, containing the following information

1)  x start
2)  y start
3)  z start
4)  x extent
5)  y extent
6)  z extent
7)  x sampling rate
8)  y sampling rate
9)  z sampling rate
10) Header(18) * A Cell Edge
11) Header(18) * B Cell Edge
12) Header(18) * C Cell Edge
13) Header(18) * alpha
14) Header(18) * beta
15) Header(18) * gamma
16) Header(19) * (253-3)/(rhomax - rhomin)
17) (3rhomax - 253rhomin)/(rhomax - rhomin)
18) Cell Constant Scaling Factore
19) 100

In the original DSN6 format, this information was stored in
the first elements of an INTEGER*2 array with 256 elements.
The problem with that is that such an array is stored differently
in big endian (for example SGI) and small endian (for example
ALPHA) machines. BRIX overcomes this problem by storing the header
as a character string. The first few bytes of the character string
is a signature (a "smiley") that enables O and other programs to
recognize the format of the file. Then comes a keyworded sequence
giving the origin, extent, unit cell parameters etc.
"""
        print txt
    
class ReadSPIDER(VolumeReaderBase):
    """Read a SPIDER formatted file"""

    def read(self, filename, normalize=True, disp_out=True):
        
        t1 = time()
        sw = "" # used in format to swap bytes if necessary
        self.filename = filename
        
        # open file to read in binary mode
        f = open(filename, 'rb')

        h = self.header = {}
        h['mapType'] = 'SPIDER'
        h['filename'] = filename

        size, unpackType, arraytype, gtype = self.getDataElemSize(4)

        data = f.read(1052)
        iform = unpack('f', data[16:20])[0]

        # find out about byteswapping (iform should be 3)
        # first determine if byte order is big-endian
        if iform != 3:
            sw=">"
            iform = unpack('f', data[16:20])[0]

        # if not big-endian, check if little-endian
        if iform != 3:
            sw="<"
            prod_norm = unpack('f', data[16:20])[0]

        # if neither, then this is an invalid file
        if iform != 3:                
            print "ReadSPIDER: ERROR: %s is not a valid SPIDER file"%filename
            return
        
        ns, nr, irec, nhrec, iform = unpack(sw+'5f', data[:20])
        h.update({'ns':ns, 'nr':nr, 'irec':irec, 'nhrec':nhrec, 'iform':iform})
        imami, amax, amin, amean, arms = unpack(sw+'5f', data[20:40])
        h.update({'imami':imami, 'amax':amax, 'amin':amin,
                  'amean':amean, 'arms':arms})
        nc, labrec, iangle = unpack(sw+'3f', data[44:56])
        h.update({'nc':nc, 'labrec':labrec, 'iangle':iangle})
        if iangle == 0: alpha, beta, gamma = (90.,90.,90.)
        else: alpha, beta, gamma = unpack(sw+'3f', data[56:68])
        h.update({ 'alpha':alpha, 'beta':beta, 'gamma':gamma })
        xoff, yoff, zoff = unpack(sw+'3f', data[68:80])
        h.update({'xoff':xoff, 'yoff':yoff, 'zoff':zoff })
        scale, labbyt, lenbyt, istack = unpack(sw+'4f', data[80:96])
        h.update({'scale':scale, 'labbyt':labbyt,
                  'lenbyt':lenbyt, 'istack':istack })
        maxim, imgnum, lastindx = unpack(sw+'3f', data[100:112])
        h.update({'maxim':maxim, 'imgnum':imgnum, 'lastindx':lastindx })
        kangle, phi1, theta1, psi1 = unpack(sw+'4f', data[120:136])
        h.update({ 'kangle':kangle, 'phi1':phi1, 'theta1':theta1, 'psi1':psi1 })
        phi2, theta2, psi2 = unpack(sw+'3f', data[136:148])
        h.update({ 'phi2':phi2, 'theta2':theta2, 'psi2':psi2 })

        # go to the where the data starts (specified in header)
        f.seek(labbyt)
        
        data_size=int(nc*nr*ns)
        data=f.read( ) # read to end of file
        f.close()

        ndata = numpy.array( unpack(sw+"%d%c"%(data_size,unpackType),
                                      data[:data_size*size]), arraytype )

        self.data = numpy.reshape(ndata, [ns,nr,nc])

        origin = [0.,0.,0.]
        stepSize = [1.,1.,1.]
        h['origin'] = origin

        g = gtype( self.data, origin, stepSize, h)

        h['dimensions'] = str(self.data.shape)

        if disp_out is True:
            print self.describe(h)
            print 'time: ', time()-t1
        
        return g

    def info(self):
        txt = """http://www.wadsworth.org/spider_doc/spider/docs/image_doc.html
Layout of the SPIDER header is as follows:

Word No.
1.  nslice = number of slices (planes) in volume (=1 for an image)
2.  nrow = number of rows per slice.
3.  irec = total number of records in the file (unused)
4.  nhistrec = (obsolete, unused)
5.  iform = file type specifier.
    Obsolete file types d, 8, 11, 12, 16, -1, -3, -7, and -9
    are no longer supported in SPIDER.
    iform (type)  data type
     1    (r)     2D image.
     3    (r)     3D volume.
    -11   (fo)    2D Fourier, mixed radix odd.
    -12   (fe)    2D Fourier, mixed radix even.
    -21   (fo)    3D Fourier, mixed radix odd.
    -22   (fe)    3D Fourier, mixed radix even.
6.  imami = maximum/minimum flag.
    Is set at 0 when the file is created, and at 1 when the
    maximum, minimum, average, and standard deviation have been
    computed and stored into this header record (see following locations).
7.  fmax = maximum value.
8.  fmin = minimum value.
9.  av = average value.
10. sig = standard deviation. A value of -1.0 indicates
    that sig has not been computed previously.
11. ihist = (obsolete, no longer used).
12. nsam = number of pixels per line.
13. labrec = number of records in file header (label).
14. iangle = flag that tilt angles are present.
15. phi = tilt angle (See note #2 below).
16. theta = tilt angle.
17. gamma = tilt angle (also called psi).
18. xoff = x translation.
19. yoff = y translation.
20. zoff = z translation.
21. scale = scale factor.
22. labbyt = total number of bytes in header.
23. lenbyt = record length in bytes.
24. istack = This position has a value of 0 in simple 2D or 3D (non-stack)
    files. In an "image stack" there is one overall stack header followed
    by a stack of images in which each image has its own image header.
    (An image stack differs from a simple 3D image in that each stacked
    image has its own header.) A value of >0 in this position in the overall
    stack header indicates a stack of images. A value of <0 in this position
    in the overall stack header indicates an indexed stack of images and
    gives the maximum image number allowed in the index.
25. NOTUSED = This position is unused now! Prior to release 9.0,
    a -1 at this location in an overall stack indicated a valid stack
    and in the stacked images, a value of 1 indicated that this image
    was in use (existed).
26. maxim = This position is only used in the overall header for a stacked
    image file. There, this position contains the number of the highest
    image currently used in the stack. This number is updated, if necessary,
    when an image is added or deleted from the stack.
27. imgnum = This position is only used in a stacked image header.
    There, this position contains the number of the current image or
    zero if the image is unused.
28. lastindx = This position is only used in the overall header of indexed
    stacks. There, this position is the highest index currently in use.
29. unused
30. unused
31. Kangle = flag that additional angles are present in header.
    1 = one additional rotation is present,
    2 = additional rotation that preceeds the rotation that was
    stored in words 15..20.
32. phi1
33. theta1
34. psi1
35. phi2
36. theta2
37. psi2
50-76 == reserved for Jose Maria's transforms 
212-214 == cdat = character * 11 containing creation date e.g. 27-MAY-1999 
215-216 -- ctim = character * 8 containing creation time e.g. 09:43:19 
217-256 -- ctit = character * 160 containing title

Note#1 :
All character arrays are retrieved from the floating point buffer array
containing the header record(s) by equivalence assignments. Thus character
arrays are stored in the header without any conversion.

Note#2 :
The angle, offset & scale factor locations contained in the SPIDER header
are available to communicate between different SPIDER operations. Currently
they are NOT used in the code distributed with SPIDER, but some outside labs
make extensive use of these positions. The angles are usually in Euler format
and are given in degrees.

Note#3 :
SGI, most IBM, Macintosh, and Sun Unix machines use a different byte
ordering from GNU/Linux on Intel, or HP Alpha machines. SPIDER contains
the "CP TO OPEND" operation to interconvert these files. However SPIDER
can read/write either byte ordering now.
"""
        print txt


class ReadRawiv(VolumeReaderBase):
    """Read a rawiv binary file"""
    
    def read(self, filename, normalize=0, verbose=True):
        
        myfile = open(filename, 'rb')
        self.header = {'title': 'rawiv from %s'%filename}

        #the header
        header=myfile.read(68)
        # unpack header following rawiv format,big endian
        h = unpack('>6f5I6f',header)
        width=int(h[8])
        self.header['width'] = width
        height=int(h[9])
        self.header['height'] = height
        depth =int(h[10])
        self.header['depth'] = depth
        nverts = int(h[6])
        ncells = int(h[7])
        origin = h[11:14]
        step = h[14:17]
        self.header['nverts'] = nverts
        self.header['ncells'] = ncells
        size = width*height*depth
        if verbose:
            print "header: ", h
            print "nverts: ", nverts, "ncells: ", ncells
        # load the data
        #l = myfile.read(width*height*depth)
        l = myfile.read() # read the rest of the file
        myfile.close()
        nbytes = len(l)
        elsize = nbytes/size
        if elsize == 1:
            mode = 1
        elif elsize == 4:
            mode = 9
        elif elsize == 8:
            mode = 10
        else:
            print "Error: in ReadRawiv - unsupported data type, size= %d"%elsize
            return None
        
        size, unpackType, arraytype, gtype = self.getGRDElemSize(mode)
        #self.data = numpy.fromstring(l, arraytype, (width*height*depth))
        #self.data.shape = (width, height, depth)
        frm = ">%i%s"%(width*height*depth, unpackType)
        self.data = numpy.reshape(numpy.array(unpack(frm, l), arraytype), (width,height,depth))
        grid = gtype( self.data, origin, step, self.header)
        return grid


class ReadEM(VolumeReaderBase):
    # (c) Daniel Stoffler, Biozentrum, University of Basel, Switzerland,
    # July 2005
    
    """Read an EM volume file and return a Grid3D object. """

    def read(self, filename, disp_out=True, normalize=True):

        t1 = time()
        sw = "" # used in format to swap bytes if necessary
        self.filename = filename
        
        # open file to read in binary mode
        f = open(filename, 'rb')
        # read the 512 bytes of the header
        header = f.read(512)

        self.header = {}
        self.header['mapType'] = 'EM'
        self.header['filename'] = filename

        machine, dummy1, dummy2, coding = unpack("4b", header[0:4])
        
        if machine in [5,6]:
            swap = "<" # PC or Mac: Little Endian
        else:
            swap = ">" # VAX, Convex, SGI, Sun: Big Endian

        self.header['machine'] = machine
        self.header['coding'] = coding

        xDim, yDim, zDim = unpack(swap+"3L", header[4:16])
        self.header["mapc"] = xDim
        self.header["mapr"] = yDim
        self.header["maps"] = zDim

        comment = unpack(swap+"80c", header[16:96])
        self.header["comment"] = comment

        # important map info follows:
        userData1 = []
        for i in range(96,256,4):
            userData1.append( unpack(swap+"1l", header[i:i+4])[0] )
        self.header["userData1"] = userData1

        userData2 = ""
        for i in range(256,512,1):
            userData2 +=  unpack(swap+"1c", header[i:i+1])[0]
        self.header["userData2"] = userData2


        data = f.read(  ) # read to end
        f.close()

        assert coding in [1,2,4,5,8,9]

        # set the correct data coding mode
        if coding == 1:
            mode = 2 # byte
        elif coding == 2:
            mode = 4 # short
        elif coding == 4:
            mode = 8 # long int
        elif coding == 5:
            mode = 9 # float
        elif coding == 8:  # COMPLEX IS NOT SUPPORTED CURRENTLY
            warnings.warn(
                "READ EM: unsupported data type COMPLEX!\nMAP IMPORT FAILED!")
            return None
        elif coding == 9:
            mode = 10 # double

        size, unpackType, arraytype, gtype = self.getGRDElemSize(mode)

        nc = xDim
        nr = yDim
        ns = zDim
        
        ndata = numpy.array( unpack(swap+"%d%c"%(nc*nr*ns,unpackType),
                                      data[:nc*nr*ns*size]), arraytype )

        # 'FORTRAN style' (x fastest, then y, then z)
        ndata.shape = (ns,nr,nc)
        #transpose the scalar data due to FORTRAN style
        self.data = numpy.ascontiguousarray(numpy.transpose(ndata), ndata.dtype.char)

        origin = [0.,0.,0.]   # origin is fixed in this data set
        stepSize = [1.,1.,1.] # so is the stepsize 
        h = {}
        h['origin'] = origin
        h.update( {'nc':nc, 'nr':nr, 'ns':ns })

        g = gtype( self.data, origin, stepSize, h)
        # compute min, max, mean rms
        amin, amax, amean, arms = g.stats()
        h.update( {'amin':amin, 'amax':amax, 'amean':amean, 'arms':arms })

        h['dimensions'] = str(self.data.shape)
        self.header.update(h)

        if disp_out is True:
            print self.describe(self.header)
            print 'time: ', time()-t1
        
        return g
       
    def info(self):
        txt = """
This map format was developed at the MPI for Biochemistry, Munich, in
Prof. Wolfgang Baumeister's group. 

For problems with this reader, please contact Daniel Stoffler.

BYTE 1: Machine Coding:
        OS-9        0
        VAX         1
        Convex      2
        SGI         3
        SUN         4
        MAC         5
        PC          6

BYTE 2: General purpose. 0 is old version, 1 is new version

BYTE 3: Not used. If 1, the header is abandoned

BYTE 4: Data Type Coding:
         byte     (1 byte)   1
         short    (2 byte)   2
         long int (4 byte)   4
         float    (4 byte)   5
         complex  (8 byte)   8
         double   (8 byte)   9

BYTE 5-16: Three long integers (3x4byte) are image size in x,y,z dimensions

BYTE 17-96: 80 Characters comment

BYTE 97-256 40 long integers (4x40byte) with user-defined values:
 
Organization of userData1:
0: accelerating voltage (Volt)
1: Cs of objective lens (mm)
2: apperture (mrad)
3: end magnification
4: postmagnification of CCD (fixed 1000)
5: exposure time in seconds
6: pixel size in object plane (nm)
7: EM Code: EM420=1; CM12=2; CM200=3; CM120/Biofilter=4; CM300=5; Polara=6;
            extern=0
8: pyhsical pixelsize on CCD (um)  u=micro
9: phys_pixel_siz * nr_pixels (um)
10: defocus, underfocus is negative (Angstr)
11: astigmatism (Angstr)
12: angle of astigmatism (deg)
13: focus increment for focus-series (Angstr)
14: counts per primary electron, sensitivity of CCD
15-39: many more, less important values...

BYTE 257-512: 256 Bytes with user data, first 20 char username, 8 chars date

BYTE 513-...: raw data: x fastest dimension, then y, then z

""" 
        return txt

class ReadFLDBinary(VolumeReaderBase):
    """Reads AVS/Express Field binary data
    http://help.avs.com/Express/doc/help/books/dv/dvfdtype.html"""
    
    def read(self, filename, normalize=True):
        myfile = open(filename, 'r')
        head = myfile.readlines(15) #limits the number of header lines
        self.header = ''
        self.comment = ''
        self.gridType = None
        self.min_ext = []
        self.max_ext = []
        for line in head:
            if line[0] == '#':
                self.comment += line
                continue
            if line[:4] == 'ndim':
                str = line.split('=')
                str= str[1].split('#')
                self.ndim = int(str[0])
                self.comment += line
                continue
            if line[:4] == 'dim1':
                str = line.split('=')
                str= str[1].split('#')
                self.dim1 = int(str[0])                
                self.comment += line
                continue
            if line[:4] == 'dim2':
                str = line.split('=')
                str= str[1].split('#')
                self.dim2 = int(str[0])                
                self.comment += line
                continue
            if line[:4] == 'dim3':
                str = line.split('=')
                str= str[1].split('#')
                self.dim3 = int(str[0])                
                self.comment += line
                continue
            if line[:6] == 'nspace':
                str = line.split('=')
                str= str[1].split('#')
                self.nspace = int(str[0])                
                self.comment += line
                continue
            if line[:6] == 'veclen':
                str = line.split('=')
                str= str[1].split('#')
                self.veclen = int(str[0])                
                self.comment += line
                continue
            if line[:4] == 'data':
                if line[5:10] == 'float':
                    dtypec='f'
                    ntype=numpy.float32
                    self.gridType = Grid3DF
                    self.comment += line
                    continue
                if line[5:12] == 'double':
                    dtypec='d'
                    ntype=numpy.float64
                    self.gridType = Grid3DD
                    self.comment += line
                    continue
                if line[5:13] == 'integer':
                    dtypec='i'
                    ntype=numpy.int32
                    self.gridType = Grid3DUI
                    self.comment += line
                    continue
                if line[5:9] == 'byte':
                    dtypec='B'
                    ntype=numpy.uint8
                    self.gridType = Grid3DUC
                    self.comment += line
                    continue
            if line[:7] == 'min_ext':            
                str = line.split('=')
                str= str[1].split('#')
                str= str[0].split()
                self.min_ext.append(float(str[0]))
                self.min_ext.append(float(str[1]))
                self.min_ext.append(float(str[2]))
                self.comment += line
                continue
            if line[:7] == 'max_ext':            
                str = line.split('=')
                str= str[1].split('#')
                str= str[0].split()
                self.max_ext.append(float(str[0]))
                self.max_ext.append(float(str[1]))
                self.max_ext.append(float(str[2]))
                self.comment += line
                continue
            if line[:5] == 'field':            
                str = line.split('=')
                str= str[1].split('#')
                self.field = str[0]
                self.comment += line
                continue

        flag_mead = False #this flag used to identify fld generated by MEAD

        myfile.close()
        myfile = open(filename, 'rb')
        data = myfile.read()
        start = data.find(pack('bb',12,12))
        data = data[start+2:]
        size = self.dim1*self.dim2*self.dim3
        #size += len(self.max_ext)+len(self.min_ext)

        data = unpack('%d%1s'%(size,dtypec), data[:size*calcsize(ntype)])
        data = numpy.array(data,ntype)

        data = numpy.reshape(data,(self.dim1,self.dim2,self.dim3))
        origin = []
        origin.append(self.min_ext[0]) 
        origin.append(self.min_ext[1]) 
        origin.append(self.min_ext[2]) 

        stepSize = []
        stepSize.append((self.max_ext[0]-self.min_ext[0])/(self.dim1-1)) 
        stepSize.append((self.max_ext[1]-self.min_ext[1])/(self.dim2-1)) 
        stepSize.append((self.max_ext[2]-self.min_ext[2])/(self.dim3-1)) 
        g = self.gridType(data, origin, stepSize, self.comment)
        return g
if __name__=='__main__':
#    reader = ReadCNS()
#    header = reader.read("Tests/Data/agrin.cns")
    reader = ReadSPIDER()
    header = reader.read("Tests/Data/ccmv.DAT")
    
##      #import isocontour
##      from UTpackages.UTisocontour import isocontour
##      vol.shape = (1, 1) + vol.shape # add time steps and variables dimensions
##      orig = [-xlen/2., -ylen/2., -zlen/2.]
##      span = (xlen/nx, ylen/ny, zlen/nx)

##      the_data = isocontour.newDatasetRegFloat3D(vol, orig, span )
##      isovar   = 0
##      timestep = 0
##      isovalue = 5.0

##      isoc = isocontour.getContour3d(the_data, isovar, timestep, isovalue,
##                                     isocontour.NO_COLOR_VARIABLE)
##      vert = numpy.zeros((isoc.nvert,3)).astype('f')
##      norm = numpy.zeros((isoc.nvert,3)).astype('f')
##      col = numpy.zeros((isoc.nvert)).astype('f')
##      tri = numpy.zeros((isoc.ntri,3)).astype('i')
##      print "nvert:", isoc.nvert
##      print "ntri:", isoc.ntri
##      isocontour.getContour3dData(isoc, vert, norm, col, tri, 1)

##      from DejaVu import Viewer
##      vi = Viewer()
##      from DejaVu.IndexedPolygons import IndexedPolygons
##      pol = IndexedPolygons('iso', vertices=vert, vnormals=norm, faces=tri)
##      #pol.Set(vertices=vert, vnormals=norm, faces=tri)
##      vi.AddObject(pol)
