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

from DejaVu.colorMap import ColorMap
from numpy import array
cm = ColorMap('rgb256')
cfg = {'name': 'rgb256', 'ramp': array([[ 0.      ,  0.      ,  1.      ,  1.      ],
       [ 0.      ,  0.015625,  1.      ,  1.      ],
       [ 0.      ,  0.03125 ,  1.      ,  1.      ],
       [ 0.      ,  0.046875,  1.      ,  1.      ],
       [ 0.      ,  0.0625  ,  1.      ,  1.      ],
       [ 0.      ,  0.078125,  1.      ,  1.      ],
       [ 0.      ,  0.09375 ,  1.      ,  1.      ],
       [ 0.      ,  0.109375,  1.      ,  1.      ],
       [ 0.      ,  0.125   ,  1.      ,  1.      ],
       [ 0.      ,  0.140625,  1.      ,  1.      ],
       [ 0.      ,  0.15625 ,  1.      ,  1.      ],
       [ 0.      ,  0.171875,  1.      ,  1.      ],
       [ 0.      ,  0.1875  ,  1.      ,  1.      ],
       [ 0.      ,  0.203125,  1.      ,  1.      ],
       [ 0.      ,  0.21875 ,  1.      ,  1.      ],
       [ 0.      ,  0.234375,  1.      ,  1.      ],
       [ 0.      ,  0.25    ,  1.      ,  1.      ],
       [ 0.      ,  0.265625,  1.      ,  1.      ],
       [ 0.      ,  0.28125 ,  1.      ,  1.      ],
       [ 0.      ,  0.296875,  1.      ,  1.      ],
       [ 0.      ,  0.3125  ,  1.      ,  1.      ],
       [ 0.      ,  0.328125,  1.      ,  1.      ],
       [ 0.      ,  0.34375 ,  1.      ,  1.      ],
       [ 0.      ,  0.359375,  1.      ,  1.      ],
       [ 0.      ,  0.375   ,  1.      ,  1.      ],
       [ 0.      ,  0.390625,  1.      ,  1.      ],
       [ 0.      ,  0.40625 ,  1.      ,  1.      ],
       [ 0.      ,  0.421875,  1.      ,  1.      ],
       [ 0.      ,  0.4375  ,  1.      ,  1.      ],
       [ 0.      ,  0.453125,  1.      ,  1.      ],
       [ 0.      ,  0.46875 ,  1.      ,  1.      ],
       [ 0.      ,  0.484375,  1.      ,  1.      ],
       [ 0.      ,  0.5     ,  1.      ,  1.      ],
       [ 0.      ,  0.515625,  1.      ,  1.      ],
       [ 0.      ,  0.53125 ,  1.      ,  1.      ],
       [ 0.      ,  0.546875,  1.      ,  1.      ],
       [ 0.      ,  0.5625  ,  1.      ,  1.      ],
       [ 0.      ,  0.578125,  1.      ,  1.      ],
       [ 0.      ,  0.59375 ,  1.      ,  1.      ],
       [ 0.      ,  0.609375,  1.      ,  1.      ],
       [ 0.      ,  0.625   ,  1.      ,  1.      ],
       [ 0.      ,  0.640625,  1.      ,  1.      ],
       [ 0.      ,  0.65625 ,  1.      ,  1.      ],
       [ 0.      ,  0.671875,  1.      ,  1.      ],
       [ 0.      ,  0.6875  ,  1.      ,  1.      ],
       [ 0.      ,  0.703125,  1.      ,  1.      ],
       [ 0.      ,  0.71875 ,  1.      ,  1.      ],
       [ 0.      ,  0.734375,  1.      ,  1.      ],
       [ 0.      ,  0.75    ,  1.      ,  1.      ],
       [ 0.      ,  0.765625,  1.      ,  1.      ],
       [ 0.      ,  0.78125 ,  1.      ,  1.      ],
       [ 0.      ,  0.796875,  1.      ,  1.      ],
       [ 0.      ,  0.8125  ,  1.      ,  1.      ],
       [ 0.      ,  0.828125,  1.      ,  1.      ],
       [ 0.      ,  0.84375 ,  1.      ,  1.      ],
       [ 0.      ,  0.859375,  1.      ,  1.      ],
       [ 0.      ,  0.875   ,  1.      ,  1.      ],
       [ 0.      ,  0.890625,  1.      ,  1.      ],
       [ 0.      ,  0.90625 ,  1.      ,  1.      ],
       [ 0.      ,  0.921875,  1.      ,  1.      ],
       [ 0.      ,  0.9375  ,  1.      ,  1.      ],
       [ 0.      ,  0.953125,  1.      ,  1.      ],
       [ 0.      ,  0.96875 ,  1.      ,  1.      ],
       [ 0.      ,  0.984375,  1.      ,  1.      ],
       [ 0.      ,  1.      ,  1.      ,  1.      ],
       [ 0.      ,  1.      ,  0.984375,  1.      ],
       [ 0.      ,  1.      ,  0.96875 ,  1.      ],
       [ 0.      ,  1.      ,  0.953125,  1.      ],
       [ 0.      ,  1.      ,  0.9375  ,  1.      ],
       [ 0.      ,  1.      ,  0.921875,  1.      ],
       [ 0.      ,  1.      ,  0.90625 ,  1.      ],
       [ 0.      ,  1.      ,  0.890625,  1.      ],
       [ 0.      ,  1.      ,  0.875   ,  1.      ],
       [ 0.      ,  1.      ,  0.859375,  1.      ],
       [ 0.      ,  1.      ,  0.84375 ,  1.      ],
       [ 0.      ,  1.      ,  0.828125,  1.      ],
       [ 0.      ,  1.      ,  0.8125  ,  1.      ],
       [ 0.      ,  1.      ,  0.796875,  1.      ],
       [ 0.      ,  1.      ,  0.78125 ,  1.      ],
       [ 0.      ,  1.      ,  0.765625,  1.      ],
       [ 0.      ,  1.      ,  0.75    ,  1.      ],
       [ 0.      ,  1.      ,  0.734375,  1.      ],
       [ 0.      ,  1.      ,  0.71875 ,  1.      ],
       [ 0.      ,  1.      ,  0.703125,  1.      ],
       [ 0.      ,  1.      ,  0.6875  ,  1.      ],
       [ 0.      ,  1.      ,  0.671875,  1.      ],
       [ 0.      ,  1.      ,  0.65625 ,  1.      ],
       [ 0.      ,  1.      ,  0.640625,  1.      ],
       [ 0.      ,  1.      ,  0.625   ,  1.      ],
       [ 0.      ,  1.      ,  0.609375,  1.      ],
       [ 0.      ,  1.      ,  0.59375 ,  1.      ],
       [ 0.      ,  1.      ,  0.578125,  1.      ],
       [ 0.      ,  1.      ,  0.5625  ,  1.      ],
       [ 0.      ,  1.      ,  0.546875,  1.      ],
       [ 0.      ,  1.      ,  0.53125 ,  1.      ],
       [ 0.      ,  1.      ,  0.515625,  1.      ],
       [ 0.      ,  1.      ,  0.5     ,  1.      ],
       [ 0.      ,  1.      ,  0.484375,  1.      ],
       [ 0.      ,  1.      ,  0.46875 ,  1.      ],
       [ 0.      ,  1.      ,  0.453125,  1.      ],
       [ 0.      ,  1.      ,  0.4375  ,  1.      ],
       [ 0.      ,  1.      ,  0.421875,  1.      ],
       [ 0.      ,  1.      ,  0.40625 ,  1.      ],
       [ 0.      ,  1.      ,  0.390625,  1.      ],
       [ 0.      ,  1.      ,  0.375   ,  1.      ],
       [ 0.      ,  1.      ,  0.359375,  1.      ],
       [ 0.      ,  1.      ,  0.34375 ,  1.      ],
       [ 0.      ,  1.      ,  0.328125,  1.      ],
       [ 0.      ,  1.      ,  0.3125  ,  1.      ],
       [ 0.      ,  1.      ,  0.296875,  1.      ],
       [ 0.      ,  1.      ,  0.28125 ,  1.      ],
       [ 0.      ,  1.      ,  0.265625,  1.      ],
       [ 0.      ,  1.      ,  0.25    ,  1.      ],
       [ 0.      ,  1.      ,  0.234375,  1.      ],
       [ 0.      ,  1.      ,  0.21875 ,  1.      ],
       [ 0.      ,  1.      ,  0.203125,  1.      ],
       [ 0.      ,  1.      ,  0.1875  ,  1.      ],
       [ 0.      ,  1.      ,  0.171875,  1.      ],
       [ 0.      ,  1.      ,  0.15625 ,  1.      ],
       [ 0.      ,  1.      ,  0.140625,  1.      ],
       [ 0.      ,  1.      ,  0.125   ,  1.      ],
       [ 0.      ,  1.      ,  0.109375,  1.      ],
       [ 0.      ,  1.      ,  0.09375 ,  1.      ],
       [ 0.      ,  1.      ,  0.078125,  1.      ],
       [ 0.      ,  1.      ,  0.0625  ,  1.      ],
       [ 0.      ,  1.      ,  0.046875,  1.      ],
       [ 0.      ,  1.      ,  0.03125 ,  1.      ],
       [ 0.      ,  1.      ,  0.015625,  1.      ],
       [ 0.      ,  1.      ,  0.      ,  1.      ],
       [ 0.015625,  1.      ,  0.      ,  1.      ],
       [ 0.03125 ,  1.      ,  0.      ,  1.      ],
       [ 0.046875,  1.      ,  0.      ,  1.      ],
       [ 0.0625  ,  1.      ,  0.      ,  1.      ],
       [ 0.078125,  1.      ,  0.      ,  1.      ],
       [ 0.09375 ,  1.      ,  0.      ,  1.      ],
       [ 0.109375,  1.      ,  0.      ,  1.      ],
       [ 0.125   ,  1.      ,  0.      ,  1.      ],
       [ 0.140625,  1.      ,  0.      ,  1.      ],
       [ 0.15625 ,  1.      ,  0.      ,  1.      ],
       [ 0.171875,  1.      ,  0.      ,  1.      ],
       [ 0.1875  ,  1.      ,  0.      ,  1.      ],
       [ 0.203125,  1.      ,  0.      ,  1.      ],
       [ 0.21875 ,  1.      ,  0.      ,  1.      ],
       [ 0.234375,  1.      ,  0.      ,  1.      ],
       [ 0.25    ,  1.      ,  0.      ,  1.      ],
       [ 0.265625,  1.      ,  0.      ,  1.      ],
       [ 0.28125 ,  1.      ,  0.      ,  1.      ],
       [ 0.296875,  1.      ,  0.      ,  1.      ],
       [ 0.3125  ,  1.      ,  0.      ,  1.      ],
       [ 0.328125,  1.      ,  0.      ,  1.      ],
       [ 0.34375 ,  1.      ,  0.      ,  1.      ],
       [ 0.359375,  1.      ,  0.      ,  1.      ],
       [ 0.375   ,  1.      ,  0.      ,  1.      ],
       [ 0.390625,  1.      ,  0.      ,  1.      ],
       [ 0.40625 ,  1.      ,  0.      ,  1.      ],
       [ 0.421875,  1.      ,  0.      ,  1.      ],
       [ 0.4375  ,  1.      ,  0.      ,  1.      ],
       [ 0.453125,  1.      ,  0.      ,  1.      ],
       [ 0.46875 ,  1.      ,  0.      ,  1.      ],
       [ 0.484375,  1.      ,  0.      ,  1.      ],
       [ 0.5     ,  1.      ,  0.      ,  1.      ],
       [ 0.515625,  1.      ,  0.      ,  1.      ],
       [ 0.53125 ,  1.      ,  0.      ,  1.      ],
       [ 0.546875,  1.      ,  0.      ,  1.      ],
       [ 0.5625  ,  1.      ,  0.      ,  1.      ],
       [ 0.578125,  1.      ,  0.      ,  1.      ],
       [ 0.59375 ,  1.      ,  0.      ,  1.      ],
       [ 0.609375,  1.      ,  0.      ,  1.      ],
       [ 0.625   ,  1.      ,  0.      ,  1.      ],
       [ 0.640625,  1.      ,  0.      ,  1.      ],
       [ 0.65625 ,  1.      ,  0.      ,  1.      ],
       [ 0.671875,  1.      ,  0.      ,  1.      ],
       [ 0.6875  ,  1.      ,  0.      ,  1.      ],
       [ 0.703125,  1.      ,  0.      ,  1.      ],
       [ 0.71875 ,  1.      ,  0.      ,  1.      ],
       [ 0.734375,  1.      ,  0.      ,  1.      ],
       [ 0.75    ,  1.      ,  0.      ,  1.      ],
       [ 0.765625,  1.      ,  0.      ,  1.      ],
       [ 0.78125 ,  1.      ,  0.      ,  1.      ],
       [ 0.796875,  1.      ,  0.      ,  1.      ],
       [ 0.8125  ,  1.      ,  0.      ,  1.      ],
       [ 0.828125,  1.      ,  0.      ,  1.      ],
       [ 0.84375 ,  1.      ,  0.      ,  1.      ],
       [ 0.859375,  1.      ,  0.      ,  1.      ],
       [ 0.875   ,  1.      ,  0.      ,  1.      ],
       [ 0.890625,  1.      ,  0.      ,  1.      ],
       [ 0.90625 ,  1.      ,  0.      ,  1.      ],
       [ 0.921875,  1.      ,  0.      ,  1.      ],
       [ 0.9375  ,  1.      ,  0.      ,  1.      ],
       [ 0.953125,  1.      ,  0.      ,  1.      ],
       [ 0.96875 ,  1.      ,  0.      ,  1.      ],
       [ 0.984375,  1.      ,  0.      ,  1.      ],
       [ 1.      ,  1.      ,  0.      ,  1.      ],
       [ 1.      ,  0.984375,  0.      ,  1.      ],
       [ 1.      ,  0.96875 ,  0.      ,  1.      ],
       [ 1.      ,  0.953125,  0.      ,  1.      ],
       [ 1.      ,  0.9375  ,  0.      ,  1.      ],
       [ 1.      ,  0.921875,  0.      ,  1.      ],
       [ 1.      ,  0.90625 ,  0.      ,  1.      ],
       [ 1.      ,  0.890625,  0.      ,  1.      ],
       [ 1.      ,  0.875   ,  0.      ,  1.      ],
       [ 1.      ,  0.859375,  0.      ,  1.      ],
       [ 1.      ,  0.84375 ,  0.      ,  1.      ],
       [ 1.      ,  0.828125,  0.      ,  1.      ],
       [ 1.      ,  0.8125  ,  0.      ,  1.      ],
       [ 1.      ,  0.796875,  0.      ,  1.      ],
       [ 1.      ,  0.78125 ,  0.      ,  1.      ],
       [ 1.      ,  0.765625,  0.      ,  1.      ],
       [ 1.      ,  0.75    ,  0.      ,  1.      ],
       [ 1.      ,  0.734375,  0.      ,  1.      ],
       [ 1.      ,  0.71875 ,  0.      ,  1.      ],
       [ 1.      ,  0.703125,  0.      ,  1.      ],
       [ 1.      ,  0.6875  ,  0.      ,  1.      ],
       [ 1.      ,  0.671875,  0.      ,  1.      ],
       [ 1.      ,  0.65625 ,  0.      ,  1.      ],
       [ 1.      ,  0.640625,  0.      ,  1.      ],
       [ 1.      ,  0.625   ,  0.      ,  1.      ],
       [ 1.      ,  0.609375,  0.      ,  1.      ],
       [ 1.      ,  0.59375 ,  0.      ,  1.      ],
       [ 1.      ,  0.578125,  0.      ,  1.      ],
       [ 1.      ,  0.5625  ,  0.      ,  1.      ],
       [ 1.      ,  0.546875,  0.      ,  1.      ],
       [ 1.      ,  0.53125 ,  0.      ,  1.      ],
       [ 1.      ,  0.515625,  0.      ,  1.      ],
       [ 1.      ,  0.5     ,  0.      ,  1.      ],
       [ 1.      ,  0.484375,  0.      ,  1.      ],
       [ 1.      ,  0.46875 ,  0.      ,  1.      ],
       [ 1.      ,  0.453125,  0.      ,  1.      ],
       [ 1.      ,  0.4375  ,  0.      ,  1.      ],
       [ 1.      ,  0.421875,  0.      ,  1.      ],
       [ 1.      ,  0.40625 ,  0.      ,  1.      ],
       [ 1.      ,  0.390625,  0.      ,  1.      ],
       [ 1.      ,  0.375   ,  0.      ,  1.      ],
       [ 1.      ,  0.359375,  0.      ,  1.      ],
       [ 1.      ,  0.34375 ,  0.      ,  1.      ],
       [ 1.      ,  0.328125,  0.      ,  1.      ],
       [ 1.      ,  0.3125  ,  0.      ,  1.      ],
       [ 1.      ,  0.296875,  0.      ,  1.      ],
       [ 1.      ,  0.28125 ,  0.      ,  1.      ],
       [ 1.      ,  0.265625,  0.      ,  1.      ],
       [ 1.      ,  0.25    ,  0.      ,  1.      ],
       [ 1.      ,  0.234375,  0.      ,  1.      ],
       [ 1.      ,  0.21875 ,  0.      ,  1.      ],
       [ 1.      ,  0.203125,  0.      ,  1.      ],
       [ 1.      ,  0.1875  ,  0.      ,  1.      ],
       [ 1.      ,  0.171875,  0.      ,  1.      ],
       [ 1.      ,  0.15625 ,  0.      ,  1.      ],
       [ 1.      ,  0.140625,  0.      ,  1.      ],
       [ 1.      ,  0.125   ,  0.      ,  1.      ],
       [ 1.      ,  0.109375,  0.      ,  1.      ],
       [ 1.      ,  0.09375 ,  0.      ,  1.      ],
       [ 1.      ,  0.078125,  0.      ,  1.      ],
       [ 1.      ,  0.0625  ,  0.      ,  1.      ],
       [ 1.      ,  0.046875,  0.      ,  1.      ],
       [ 1.      ,  0.03125 ,  0.      ,  1.      ],
       [ 1.      ,  0.015625,  0.      ,  1.      ]],'f'), 'maxi': 10.0, 'mini': 0.0}
cm.configure(**cfg)
