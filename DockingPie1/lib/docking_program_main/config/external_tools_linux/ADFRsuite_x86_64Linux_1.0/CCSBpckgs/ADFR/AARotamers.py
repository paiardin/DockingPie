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

########################################################################
#
# Date: Jan 2015 Author: Michel Sanner
#
#   sanner@scripps.edu
#       
#   The Scripps Research Institute (TSRI)
#   Molecular Graphics Lab
#   La Jolla, CA 92037, USA
#
# Copyright: Michel Sanner and TSRI
#
#########################################################################

#
# $Header: /mnt/raid/services/cvs/ADFR/AARotamers.py,v 1.2 2016/12/07 00:38:32 sanner Exp $
#
# $Id: AARotamers.py,v 1.2 2016/12/07 00:38:32 sanner Exp $
#
import os

import ADFR
ADFRPATH = ADFR.__path__[0]

class RotamerLib:

    def __init__(self, libFile=None):
        """
libFile: rotamer library file name
defFile: rotamer angle definition file name
confFile: conformation of amino acids. 
"""
        
        self.angleList= None
        self.angleDef = None
        if libFile is None:
            libFile = os.path.join(ADFRPATH, 'Data', 'bbind02.May.lib')
        defFile = os.path.join(ADFRPATH, 'Data', 'rotamer.def')
        self.loadRotamerLib(libFile, defFile)

    def loadRotamerLib(self, libFile, defFile):
        """ load backbone-independent rotamer library May 15, 2002
http://dunbrack.fccc.edu/bbdep/bbind02.May.lib.gz
libFile: bbind02.May.lib, all descriptive lines begin with '#'
defFile: define the chi angles and the atoms that are moved by the chi rotation

Roland L. Dunbrack, Jr., Ph. D.
Institute for Cancer Research, Fox Chase Cancer Center                  
7701 Burholme Avenue, Philadelphia PA 19111
"""
        try:
            input_file = file(libFile, 'r')
            lines=input_file.readlines()
        except:
            print "Error in opening ",libFile
            return 
        try:
            input_file = file(defFile, 'r')
            defLines =input_file.readlines()
        except:
            print "Error in opening ",defFile
            return 

        self.angleList = {}
        self.angleDef = {}
        self.angleStddev = {}
        rotlib = self.angleList
        angdef = self.angleDef
        angledev = self.angleStddev
        # parsing angle definition
        for line in defLines:
            data= line.split()
            if len(data) ==0:
                continue
            if data[0][0] == '#':
                continue                
            name = data[0]
                        
            X = data[1:5]
            moving = data[5:]
            if angdef.has_key(name):
                angdef[name].append([X, moving])
            else:
                angdef[name]= [[X, moving]]

        # parsing angleList
        name=''
        for line in lines:
            data= line.split()
            if len(data) ==0:
                continue
            if data[0][0] == '#':
                continue                
            name = data[0]            
            # hardwired code for parsing Dunbrack's lib
            number = (len(data) - 11)/2 # number of CHI angles
            chi_angles=[]
            stdev = []
            for i in range(number):
                chi_angles.append(float(data[11 + i*2]))
                stdev.append(float(data[12+ i*2]))
            if rotlib.has_key(name):
                rotlib[name].append(chi_angles)
                angledev[name].append(stdev)                
            else:
                rotlib[name]= [chi_angles]
                angledev[name]= [stdev]

    def get(self, residueName ):
        """ returns the angle definition and the angle lists of the residue (residueName) """
        
        if self.angleDef.has_key(residueName) and \
               self.angleList.has_key(residueName):
            return self.angleDef[residueName], self.angleList[residueName], \
                   self.angleStddev[residueName]
        
        else:
            print residueName, "is not found in the library"
            print self.angleDef.keys(), 'are defined'
            print self.angleDef.keys(), 'angles are loaded '
            resListStr = ""
            for key in self.angleDef.keys():
                resListStr += key + ", "
            if len(resListStr): resListStr = resListStr[:-2] # to remove last ", "
            raise KeyError ("residue %s cannot be made flexible. Only residues of type: %s can be made flexible"%(
                residueName, ''.join(resListStr)))
