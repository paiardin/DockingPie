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

#
# $Header: /mnt/raid/services/cvs/PmvApp/pmvPalettes.py,v 1.1.4.1 2017/07/13 20:55:28 annao Exp $
#
# $Id: pmvPalettes.py,v 1.1.4.1 2017/07/13 20:55:28 annao Exp $
#
import numpy

RasmolAmino = { 'ALA':(.781, 0.781, 0.781), 'ARG':(0.078, 0.352, 0.999),
                'ASN':(0.,.898, 0.898), 'ASP':(0.898, 0.039, 0.039),
                'CYS':(0.898, 0.898, 0.0), 'GLN':(0.,0.898, 0.898), 
                'GLU':(0.898, 0.039, 0.039),
                'GLY':(0.918, 0.918, 0.918),
                'HIS':(0.508, 0.508, 0.820),
                'ILE':(0.059, 0.508, 0.059),
                'LEU':(0.059, 0.508, 0.059),
                'LYS':(0.078, 0.352, 0.999),
                'MET':(0.898, 0.898, 0.0), 'PHE':(0.195, 0.195, 0.664),
                'PRO':(0.859, 0.586, 0.508), 'SER':(0.977,0.586,0.),
                'THR':(0.977,0.586,0.0), 'TRP':(0.703,0.352,0.703),
                'TYR':(0.195,0.195,0.664), 'VAL':(0.059,0.508,0.059),
                '  A': (1,0,0), '  G': (0,0,1), '  C': (1,1,0),  '  T': (0,1,0),  '  U': (1,0.5,0)}

RasmolAminoSortedKeys = [ 'ALA', 'GLY', 'PRO', 'LEU', 'ILE', 'VAL',
                          'MET', 'CYS', 'THR', 'SER', 'GLU', 'ASP', 'TRP', 
                         'HIS','GLN', 'ASN', 'ARG', 'LYS', 'TYR', 'PHE', '  A', '  C', '  G', '  T', '  U' ]

#AtomElements = { 'N':(0.,0.,1.), 'C':(.54,.54,.54), 'O':(1.,0.,0.),
#                 'H':(1.,1.,1.), 'S':(1.,1.,0.), 'CA':(1.,0.,1.),
#                 'P':(1.,0.,1.), 'A':(0.,1.,0.) }

# New Graham-Stefano
AtomElements = { 'H':[0.85, 0.85, 0.85], 'C':[.54,.54,.54], 'O':[0.90, 0.10, 0.10], 'N':[0.20, 0.33, 1.00], 
                'P':[1.00, 0.59, 0.00], 'S':[0.95, 0.85, 0.00], 'B':[0.95, 0.73, 0.56],

                # halogens
                'F':[0.51, 1.00, 0.00], 'Cl':[0.51, 1.00, 0.00], 'Br':[0.60, 0.23, 0.08],
                'CL':[0.51, 1.00, 0.00], 'BR':[0.60, 0.23, 0.08], 'I':[0.42, 0.22, 0.76],

                # transition metals
                'Mn':[0.95, 0.73, 0.56], 'Co':[0.95, 0.73, 0.56], 'Ni':[0.95, 0.73, 0.56],
                'Cu':[0.95, 0.73, 0.56], 'Ru':[0.95, 0.73, 0.56], 'Pd':[0.95, 0.73, 0.56],
                'Ag':[0.95, 0.73, 0.56], 'Cd':[0.95, 0.73, 0.56], 'Pt':[0.95, 0.73, 0.56],
                'Au':[0.95, 0.73, 0.56], 'Hg':[0.95, 0.73, 0.56], 

                # non-metals
                'Si':[0.95, 0.73, 0.56], 'Al':[0.95, 0.73, 0.56], 

                # alkaline metals
                'Li':[0.63, 0.50, 1.00], 'Na':[0.63, 0.50, 1.00], 'K':[0.63, 0.50, 1.00],
                'Rb':[0.63, 0.50, 1.00], 'Cs':[0.63, 0.50, 1.00], 
                
                # alkaline earth metals
                'Be':[0.27, 0.52, 0.11], 'Mg':[0.27, 0.52, 0.11], 'Ca':[0.27, 0.52, 0.11],
                'Sr':[0.27, 0.52, 0.11], 'Ba':[0.27, 0.52, 0.11], 'Ra':[0.27, 0.52, 0.11],
                'Zn':[0.34, 0.79, 1.00],

                # iron, selenium
                'Fe':[0.83, 0.40, 0.00],  'Se':[0.70, 0.90, 0.20],
                
                }

DavidGoodsell = { 'N':(0.7,0.7,1.), 'C':(1.,1.,1.), 'O':(1.0,.7,.7),
                  'H':(1.,1.,1.), 'HN':(0.7,0.7,1.0), 
                  'S':(0.9,0.85,0.1), 'ASPOD1':(1.,0.2,0.2),
                  'ASPOD2':(1.,0.2,0.2), 'GLUOE1':(1.,0.2,0.2),
                  'GLUOE2':(1.,0.2,0.2), 'SERHG':(1.0,.7,.7),
                  'THRHG1':(1.0,.7,.7), 'TYROH':(1.0,.7,.7),
                  'TYRHH':(1.0,.7,.7), 'LYSNZ':(0.2,0.3,1.0), 
                  'LYSHZ1':(0.2,0.3,1.0), 'LYSHZ2':(0.2,0.3,1.0), 
                  'LYSHZ3':(0.2,0.3,1.0), 'ARGNE':(0.2,0.3,1.0),
                  'ARGNH1':(0.2,0.3,1.0), 'ARGNH2':(0.2,0.3,1.0),
                  'ARGHH11':(0.2,0.3,1.0), 'ARGHH12':(0.2,0.3,1.0), 
                  'ARGHH21':(0.2,0.3,1.0), 'ARGHH22':(0.2,0.3,1.0), 
                  'ARGHE':(0.2,0.3,1.0), 
                  'GLNHE21':(0.7,0.7,1.0), 'GLNHE22':(0.7,0.7,1.0),
                  'ASNHD21':(0.7,0.7,1.0), 'ASNHD22':(0.7,0.7,1.0),
                  'HISHD1':(0.7,0.7,1.0), 'HISHE2':(0.7,0.7,1.0),
                  'GLNHE2':(0.7,0.7,1.0),'ASNHD2':(0.7,0.7,1.0),
                  'CYSHG':(0.9,0.85,0.1), 'HH':(1.,1.,1.) }
DavidGoodsellSortedKeys = [
    'C', 'HN', 'HH', 'H', 'O', 'SERHG', 'TYROH', 'TYRHH', 'THRHG1',
    'N', 'ASNHD21', 'ASNHD22', 'GLNHE21', 'GLNHE22', 'GLNHE2', 'ASNHD2',
    'HISHE2', 'HISHD1', 'S', 'CYSHG', 'GLUOE1', 'GLUOE2', 'ASPOD1',
    'ASPOD2', 'LYSNZ', 'LYSHZ1', 'LYSHZ2','LYSHZ3','ARGNE', 'ARGNH1', 'ARGNH2',
    'ARGHH11','ARGHH12','ARGHH21','ARGHH22','ARGHE', ]

SecondaryStructureSides = { 'up': (1., 0., 0.), 'down': (0., 1., 0.),
                            'left': (0., 0., 1.),
                            'right': (0., 0., 1.) }

SecondaryStructureType = {'Helix': (0.937, 0.0, 0.5),
                          'Strand': (1.0, 1.0, 0.0),
                          'Turn': (0.375, .5, 1.0),
                          'Coil': (1.0, 1.0, 1.0) }

Shapely = { 'ALA':(0.117, 0.703, 0.117), 'ARG':(0.222, 0.222, 0.406),
            'ASN':(0.734, 0.558, 0.558), 'ASP':(0.996, 0.117, 0.117),
            'CYS':(0.996, 0.996, 0.000), 'GLN':(0.996, 0.558, 0.558), 
            'GLU':(0.703, 0.078, 0.039), 'GLY':(0.937, 0.937, 0.937),
            'HIS':(0.390, 0.757, 0.996), 'ILE':(0.183, 0.308, 0.183),
            'LEU':(0.417, 0.554, 0.136), 'LYS':(0.195, 0.195, 0.996),
            'MET':(0.578, 0.578, 0.066), 'PHE':(0.722, 0.281, 0.242),
            'PRO':(0.750, 0.750, 0.750), 'SER':(0.996, 0.484, 0.054),
            'THR':(0.996, 0.328, 0.054), 'TRP':(0.429, 0.105, 0.097),
            'TYR':(0.855, 0.574, 0.437), 'VAL':(0.000, 0.996, 0.000) }

DnaElements = { 'DA':(0.898, 0.039, 0.039),
                'DT':(0.508, 0.508, 0.820),
                'DU':(0.078, 0.352, 0.999),
                'DG':(0.898, 0.898, 0.0),
                'DC':(0.722, 0.281, 0.242),
                'DH':(1.,1.,1.)  }

                 

#Rainbow= {0: (0.,0.,1.),  1: (0.,1.,0.), 2: (1.,0.,0.),
#          3: (0.,1.,1.),  4: (1.,1.,0.), 5: (1.,0.,1.),
#          6: (0.,.75,1.), 7: (0.,1.,.5), 8: (.6,1.,0.),
#          9: (1.,.5,0.), 10: (1.,0.,.5), 11:(.5,0.,1.),
#          12:(.5,0.,.2), 13: (0.,.5,.2), 14:(0.75,0,1),
#          15:(1.,0.75,0),16:(1.,0.,0.75,),17:(0., 1., .75),
#          18:(1.,0.,0.25,), 19: (0.,.5,1.)}
#                    
#RainbowSortedKey = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19]

Rainbow= {'0': (0.,0.,1.),  '1': (0.,1.,0.), '2': (1.,0.,0.),
          '3': (0.,1.,1.),  '4': (1.,1.,0.), '5': (1.,0.,1.),
          '6': (0.,.75,1.), '7': (0.,1.,.5), '8': (.6,1.,0.),
          '9': (1.,.5,0.), '10': (1.,0.,.5), '11':(.5,0.,1.),
          '12':(.5,0.,.2), '13': (0.,.5,.2), '14':(0.75,0,1),
          '15':(1.,0.75,0),'16':(1.,0.,0.75,),'17':(0., 1., .75),
          '18':(1.,0.,0.25,), '19': (0.,.5,1.)}
                    
RainbowSortedKey = ['0','1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19']

##
## build color arrays
##

from MolKit2.babelElements import babel_elements
import numpy

## elementColorsArray provides a colot for an element number i.e. 6 is Carbon
elementColors = numpy.array( ([0.,0.,0.],)*len(babel_elements), 'f')
for k,v in babel_elements.items():
    try:
        col = AtomElements[k]
    except KeyError:
        col = (0.,1.,0.)
    elementColors[v['num']] = col

# SequentialColors provides an array of 20 colors where consecutive colors are
# different from each other
sequentialColors = numpy.array([
    (1.,0.5,1.), (0., 0.509,0.), (1.,0.5,0.5),  (0.25,0.75,1.), (0.517,0.2549,0.2588),
    (1.,0.,1.),  (0.,.75,1.),    (0.,1.,.5),    (.6,1.,0.),     (1.,.5,0.),
    (1.,0.,.5),  (.5,0.,1.),     (.5,0.,.2),    (0.,.5,.2),     (0.75,0,1),
    (1.,0.75,0), (1.,0.,0.75),  (0., 1., .75), (1.,0.,0.25),  (0.,.5,1.)] , 'f')

from DejaVu2.colorTool import RGBRamp
rainbow256Colors = RGBRamp(size=256)

rainbowColors = numpy.array(
    [(0.,0.,1.),   (0.,1.,0.),   (1.,0.,0.),
     (0.,1.,1.),   (1.,1.,0.),   (1.,0.,1.),
     (0.,.75,1.),  (0.,1.,.5),   (.6,1.,0.),
     (1.,.5,0.),   (1.,0.,.5),   (.5,0.,1.),
     (.5,0.,.2),   (0.,.5,.2),   (0.75,0,1),
     (1.,0.75,0),  (1.,0.,0.75,),(0., 1., .75),
     (1.,0.,0.25,),(0.,.5,1.)], 'f')
 
