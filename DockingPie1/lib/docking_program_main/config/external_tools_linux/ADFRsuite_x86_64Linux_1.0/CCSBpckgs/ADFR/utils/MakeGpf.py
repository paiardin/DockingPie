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
# Date: 2015 Author: Michel F. SANNER
#
#    sanner@scripps.edu
#
#       The Scripps Research Institute (TSRI)
#       Molecular Graphics Lab
#       La Jolla, CA 92037, USA
#
# Copyright: Michel Sanner and TSRI
#
#
#########################################################################
#
# $Header: /mnt/raid/services/cvs/ADFR/utils/MakeGpf.py,v 1.13 2017/04/20 23:29:10 sanner Exp $
#
# $Id: MakeGpf.py,v 1.13 2017/04/20 23:29:10 sanner Exp $
#
import os

class ADGPF:
    """
    Generate AutoDock Grid Parameter File (GPF)
    """

    KEYWORDS = [
        ('outlev', '%d'),
        ('parameter_file', '%s'),
        ('npts', ''),
        ('gridfld', ""),
        ('spacing', "%f"),
        ('receptor_atypes', "%s"),
        ('ligand_atypes', "%s"),
        ('receptor', "%s"),
        ('gridcenter', ''),
        ('smooth', "%f"),
        ('map', ''),
        ('dielectric', "%f"),
        ]

    def __init__(self, npts=None, spacing=0.375,
                 receptor_atypes=None, ligand_atypes=None, receptor=None,
                 gridcenter=None, smooth=0.5, dielectric=-0.1456,
                 atypesOnly=False, paramFile=None, outlev=1):

        self.npts = npts
        self.spacing = spacing
        self.receptor_atypes = receptor_atypes
        self.ligand_atypes = ligand_atypes
        self.receptor = receptor
        self.gridcenter = gridcenter
        self.smooth = smooth
        self.dielectric = dielectric
        self.map = None
        self.gridfld = None
        self.atypesOnly = atypesOnly # when True e and d maps are not computed
        self.parameter_file = paramFile
        self.outlev = outlev
        
        if npts:
            self.setNpts(npts)
        if spacing:
            self.setSpacing(spacing)
        if receptor_atypes:
            self.setReceptor_atypes(receptor_atypes)
        if ligand_atypes:
            self.setLigand_atypes(ligand_atypes)
        if receptor:
            self.setReceptor(receptor)
        if gridcenter:
            self.setGridcenter(gridcenter)
        if smooth:
            self.setSmooth(smooth)
        if dielectric:
            self.setDielectric(dielectric)

    def setNpts(self, npts):
        assert len(npts)==3
        assert isinstance(npts[0], int), npts
        assert npts[0] > 1
        assert npts[0] < 2048
        assert isinstance(npts[1], int), npts
        assert npts[1] > 1
        assert npts[1] < 2048
        assert isinstance(npts[2], int), npts
        assert npts[2] > 1
        assert npts[2] < 2048
        self.npts = npts
        
    def setSpacing(self, spacing):
        self.spacing = spacing
        
    def setReceptor_atypes(self, receptor_atypes):
        assert len(receptor_atypes) > 0
        for t in receptor_atypes:
            assert isinstance(t, (unicode,str))
        self.receptor_atypes = receptor_atypes
        
    def setLigand_atypes(self, ligand_atypes):
        #assert len(ligand_atypes) >0
        for t in ligand_atypes:
            assert isinstance(t, (unicode,str))
        self.ligand_atypes = ligand_atypes
        
    def setReceptor(self, receptor):
        assert isinstance(receptor, (unicode,str))
        assert os.path.exists(receptor)
        self.receptor = os.path.abspath(receptor)
        
    def setGridcenter(self, gridcenter):
        assert len(gridcenter)==3
        float(gridcenter[0])
        float(gridcenter[1])
        float(gridcenter[2])
        self.gridcenter = gridcenter
        
    def setSmooth(self, smooth):
        self.smooth = smooth
        
    def setDielectric(self, dielectric):
        self.dielectric = dielectric   

    def isDataValid(self):
        if self.npts is None: return False, 'missing npts[0]'
        if self.receptor_atypes is []: return False,'missing receptor_atypes'
        #if self.ligand_atypes is []: return False,'missing ligand_atypes'
        if self.receptor is '': return False,'missing receptor'
        if self.gridcenter is None: return False,'missing gridcenter[0]'
        return True, ''
    
    def getGPFlines(self):
        valid, msg = self.isDataValid()
        if not valid:
            raise ValueError("invalid data %s"%msg)
        lines = []

        rec_name = os.path.basename(os.path.splitext(self.receptor)[0])

        for k, _format in self.KEYWORDS:
            value = getattr(self, k)            
            if k == 'gridfld':
                lines.append("gridfld %s.maps.fld # avs field file" % rec_name)
            elif k == 'gridcenter':
                lines.append('%s %s'%(k,"%2.3f %2.3f %2.3f" % tuple(value)))
            elif k == 'npts':
                lines.append('%s %s'%(k,"%d %d %d" % tuple(value)))
            elif k == 'map':
                valueStr = ''
                for a in self.ligand_atypes:
                    lines.append("map %s.%s.map" % (rec_name, a))
                if not self.atypesOnly:
                    lines.append("elecmap %s.e.map" % rec_name)
                    lines.append("dsolvmap %s.d.map" % rec_name)
            elif k == 'receptor_atypes':
                line = "receptor_types "
                for t in self.receptor_atypes:
                    line += t+' '
                lines.append(line)
            elif k == 'receptor':
                lines.append("receptor %s"%os.path.basename(self.receptor))
            elif k == 'ligand_atypes':
                line = "ligand_types "
                for t in self.ligand_atypes:
                    line += t+' '
                lines.append(line)
            #elif k == 'parameter_file':
            #    if self.paramFile:
            #      lines.append('%s %s'%(k, self.paramFile))  
            else:
                lines.append('%s %s'%(k,_format%getattr(self, k)))
        return lines

if __name__=='__main__':
    obj=ADGPF()
    assert obj.isDataValid()[0] == False
    obj.setNpts([128, 128, 128])
    assert obj.isDataValid()[0] == False
    obj.setSpacing(0.375)
    assert obj.isDataValid()[0] == False
    obj.setReceptor_atypes(['C', 'N', 'O'])
    assert obj.isDataValid()[0] == False
    obj.setLigand_atypes(['C', 'N', 'O'])
    assert obj.isDataValid()[0] == False
    obj.setReceptor('foo')
    assert obj.isDataValid()[0] == False
    obj.setSmooth(0.5)
    assert obj.isDataValid()[0] == False
    obj.setDielectric(-0.1456)
    assert obj.isDataValid()[0] == False
    obj.setGridcenter([0, 0.,0])
    assert obj.isDataValid()[0] == True
    lines = obj.getGPFlines()
    
    for l in lines:
        print l
        
