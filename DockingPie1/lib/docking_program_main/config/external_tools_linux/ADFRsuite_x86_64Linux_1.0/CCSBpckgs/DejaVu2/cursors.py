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
# Author: Michel F. SANNER
#
# Copyright: M. Sanner TSRI 2010
#
#############################################################################

#
# $Header: /mnt/raid/services/cvs/DejaVu2/cursors.py,v 1.1.1.1.4.1 2017/07/13 22:28:33 annao Exp $
# 
# $Id: cursors.py,v 1.1.1.1.4.1 2017/07/13 22:28:33 annao Exp $
#

import sys, os

## define cursors
cursorsDict = {
    'default':'',
    'busy':'watch',
    'None':'pirate',
    'default':'',
    'addToSelection':'',
    'removeFromSelection':'',
    'rotation':'exchange',
    'XYtranslation':'fleur',
    'Ztranslation':'double_arrow',
    'scale':'',
    'zoom':'sizing',
    'pivotOnPixel':'cross_reverse',
    'picking':'dotbox'
#    '':'',
    }

import DejaVu2
cpath = os.path.join(DejaVu2.__path__[0], 'Cursors')

if sys.platform=='linux2':
    cursorsDict['addToSelection'] = (
        '@'+os.path.join(cpath, 'AddSel.xbm'),
        os.path.join(cpath, 'AddSelMask.xbm'),
        'black', 'white')
    cursorsDict['removeFromSelection'] = (
        '@'+os.path.join(cpath, 'MinusSel.xbm'),
        os.path.join(cpath, 'MinusSelMask.xbm'),
        'black', 'white')
    cursorsDict['Ztranslation'] = (
        '@'+os.path.join(cpath, 'Ztrans.xbm'),
        os.path.join(cpath, 'ZtransMask.xbm'),
        'black', 'white')
elif os.name == 'nt':
    pass
    ## cursorsDict['None']='@'+os.path.join(cpath, 'NoAction.cur')
    ## cursorsDict['addToSelection']='@'+os.path.join(cpath, 'AddSel.cur')
    ## cursorsDict['removeFromSelection']='@'+os.path.join(cpath,
    ##                                                      'MinusSel.cur')
    ## cursorsDict['XYtranslation']='@'+os.path.join(cpath, 'XYtrans.cur')
    ## cursorsDict['Ztranslation']='@'+os.path.join(cpath, 'Ztrans.cur')

elif sys.platform=='darwin':
    cursorsDict['addToSelection'] = (
        '@'+os.path.join(cpath, 'AddSel.xbm'),
        os.path.join(cpath, 'AddSelMask.xbm'),
        'black', 'white')
    cursorsDict['removeFromSelection'] = (
        '@'+os.path.join(cpath, 'MinusSel.xbm'),
        os.path.join(cpath, 'MinusSelMask.xbm'),
        'black', 'white')
    cursorsDict['Ztranslation'] = (
        '@'+os.path.join(cpath, 'Ztrans.xbm'),
        os.path.join(cpath, 'ZtransMask.xbm'),
        'black', 'white')
    ## cursorsDict['None']='@'+os.path.join(cpath, 'NoAction.cur')
    ## cursorsDict['addToSelection']='@'+os.path.join(cpath, 'AddSel.cur')
    ## cursorsDict['removeFromSelection']='@'+os.path.join(cpath,
    ##                                                      'MinusSel.cur')
    ## cursorsDict['XYtranslation']='@'+os.path.join(cpath, 'XYtrans.cur')
    ## cursorsDict['Ztranslation']='@'+os.path.join(cpath, 'Ztrans.cur')
