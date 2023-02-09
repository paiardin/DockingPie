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

import struct, string
def readVox(filename):
    f = open(filename)
    data = f.readlines()
    f.close()
    volsize = 0
    voxsize = 0
    endian = None
    i = 0
    while(1):
        i = i+1
        if data[i]=='##\014\012': break
    while(1):
        i = i+1
        if data[i]=='##\014\012': break
        else:
            s = string.split(data[i])
            if s[0] == "VolumeSize":
                volsize = int(s[1])*int(s[2])*int(s[3])
            elif s[0] == "VoxelSize":
                voxsize = int(s[1])
            elif s[0] == "Endian":
                endian = s[1]
    s = reduce( lambda x,y: x+y, data[i+1:] )
    print "type:", type(s)
    print "len(s):", len(s)
    print "volsize:", volsize, "voxsize:", voxsize, "endian:", endian
    if voxsize == 8:
        fmt = "%s%dB"%(endian, volsize)
    elif voxsize == 16:
        fmt = "%s%dH"%(endian, volsize)
    # 'H' - unsigned short
    #fmt = "%dH"%len(s)
    # 'B' - unsigned char
    #fmt = "%dB"%volsize
    #fmt = "%dH"%volsize
    values = struct.unpack( fmt, s)
    return values
