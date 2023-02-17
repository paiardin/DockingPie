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

#############################################################################
#
# Author: Anna Omelchenko
#
# Copyright: M. Sanner TSRI
#
#############################################################################

import struct, os
from SimpleDialog import SimpleDialog
import Tkinter

class WriteVolumeToFile:
    """class that has a method to write data into a file file"""

    def __init__(self, nx, ny, nz, data = None, file = None, voxsize = 16,
                 endian = "L"):
        self.file = file
        self.size = (nx, ny, nz)
        self.data = data
        if self.data:
           assert len(self.data.shape)==1
        assert voxsize == 16 or voxsize == 8
        self.voxsize = voxsize
        assert endian == "L" or endian == "B"
        self.endian = endian


    def checkIfExists(self):
        """Check if self.file exists. If yes - creates a dialog,
        asking to ovwrwrite the file.
        Function returns: 1 if file exists - no overwritting;
                          0 if no file or file exists - overwrite. """
        if os.path.exists(self.file):
            root = Tkinter.Toplevel()
            root.withdraw()
            d = SimpleDialog(root,
                             title = "Overwrite Existing File Question",
                             text = "File %s exists. Overwrite it?"% self.file,
                             buttons=["Yes", "Cancel"],
                             default = 1,
                             cancel = 1).go()
            root.destroy()
            return d
        else: return 0
        
class WriteVolumeToVox(WriteVolumeToFile):
    
    def isVox(self):
        """Checks if self.file is .vox file"""
        suff = os.path.splitext(self.file)[-1]
        if suff != '.vox':
            print "Error: expected .vox file"
            return 0
        else : return 1


    def write_file(self):
        """writes data in .vox file (VOX file format)"""
        #assert len(self.data.shape)==1
        assert self.data.flags.contiguous
        size = self.size[0]*self.size[1]*self.size[2]
        of = open(self.file,"wb")
        of.write("Vox1999a\n")
        of.write("VolumeCount 1\n")
        of.write("##\f\n")
        of.write("##\n")
        of.write("VolumeSize %d %d %d\n" % self.size)
        of.write("VoxelSize %d \n" % self.voxsize)
        #of.write("Endian L\n")
        of.write("Endian %s\n" % self.endian)
        if self.voxsize == 16:
            # lower 12 bits of data - range of values: 0 .....4095
            of.write("Field0 (Position 0 Size 12 Name \"Density\")\n")
            fmt = "<%dH" % size  # H for unsigned short
        else:
            of.write("Field0 (Position 0 Size 8 Name \"Density\")\n")
            fmt = ">%dB" % size  # B for unsigned char
        of.write("##\f\n")
        of.write( apply(struct.pack, (fmt,)+tuple(self.data.ravel())) )
        of.close()

class WriteVolumeToRawiv(WriteVolumeToFile):
    
    def write_file(self, packType = 'B'):
        """Writes data in .rawiv file."""
        nx,ny,nz = self.size
        print "writing to file: %s, data size: %d, %d, %d" % (self.file, nx,ny,nz)
        size = nx*ny*nz
        size1 = (nx-1)*(ny-1)*(nz-1)
        of = open(self.file,"wb")
        #Header: FIXME - not sure if this is right
        st=(0.0, 0.0, 0.0, float(nx-1), float(ny-1),float(nz-1),
            size, size1, nx, ny, nz, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0)
        of.write(apply(struct.pack, ('>6f5I6f',)+st))
        #fmt = ">%dB"%size
        fmt = ">%d%s"%(size, packType)
        of.write( apply(struct.pack, (fmt,)+tuple(self.data.ravel())))
        of.close()
