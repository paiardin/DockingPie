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
## Copyright (c) MGL TSRI 2016
##
################################################################################

#
# The contents of this file are subject to the Mozilla Public
# License Version 1.1 (the "License"); you may not use this file
# except in compliance with the License. You may obtain a copy of
# the License at http://www.mozilla.org/MPL/
# 
# Software distributed under the License is distributed on an "AS
# IS" basis, WITHOUT WARRANTY OF ANY KIND, either express or
# implied. See the License for the specific language governing
# rights and limitations under the License.
# The Original Code is "Java-Python Extension libplus (JPE-libplus)".
# 
# The Initial Developer of the Original Code is Frederic Bruno Giacometti.
# Portions created by Frederic Bruno Giacometti are
# Copyright (C) 2001-2002 Frederic Bruno Giacometti. All Rights Reserved.
# 
# Contributor(s): frederic.giacometti@arakne.com
# 
# Acknowledgments:
# Particular gratitude is expressed to the following parties for their
# contributing support to the development of JPE-libplus:
#     - The Molecular Graphics Laboratory (MGL)
#       at The Scripps Research Institute (TSRI), in La Jolla, CA, USA;
#       and in particular to Michel Sanner and Arthur Olson.
#

"""helper for generating C callback files
"""

import string, re, types,sys
platformname = sys.platform
winfpref = platformname == 'win32' and 'WSTDCALL ' or ''

def cb_sigargs( args):
    if args == ['void']:
        return 'void'
    else:
        result = ' '
        for arg in args:
            if '...' == arg:
                raise NotImplementedError( arg, args)
            if result != ' ':
                result += ', '
            result += '%s %s' % arg
        return result

def cb_auxsigargs( args):
    result = cb_sigargs( args)
    if result == 'void':
        result = ''
    return result

def cb_callargs( args):
    if args == ['void']:
        return ''
    else:
        result = ' '
        for arg in args:
            if '...' == arg:
                raise NotImplementedError( arg, args)
            if result != ' ':
                result += ', '
            result += '%s' % arg[ 1]
        return result
        

def cb_name( rettype, args):
    def subfun( name):
        return re.sub( r'\*', 'star',
                       re.sub( r' ', '', name))
    return string.join( [subfun( rettype)]
                         + [subfun( isinstance( x, types.StringType)
                                    and x or x[ 0])
                            for x in args], '_'
                        )

def cb_interface( headername, nslots, pyargs, rettype, args, pythread):

    sigargs = cb_sigargs( args)
    auxsigargs = cb_auxsigargs( args)
    callargs = cb_callargs( args)
    cbname = cb_name( rettype, args)

    result = '/* %s %s callback */\n\n' % (rettype, args)

    result += '%%{\n%s %s_callback( int idx%s%s);\n\n'\
             % (rettype, cbname, auxsigargs and ',' or '',
                auxsigargs)

    for idx in range( nslots):
        result += 'static %s %s%s_%i( %s)\n' % (rettype, winfpref, cbname,
                                                idx, sigargs)\
                  + '{\n  %s_callback( %i%s%s);\n}\n\n'\
                  % (cbname, idx, callargs and ',' or '', callargs)
        
    result += 'static %s_f %s_array[] = {\n' % ((cbname,) * 2)
    
    for idx in range( nslots):
        result += '  %s_%i,\n' % (cbname, idx)
        
    result += '};\n\n'

    result += '#define %s_DIM (sizeof %s_array / sizeof *%s_array)\n\n'\
              % ((cbname,) * 3)

    result += 'static %s_f %s_array_get( int idx)\n' % ((cbname,) * 2)
    result += '{\n  assert( 0 <= idx);\n  assert( idx < %s_DIM);\n' % cbname
    result += '  return %s_array[ idx];\n}\n\n' % cbname
    
    result += 'static %s %s_callback( int idx%s%s)\n{\n'\
          % (rettype, cbname, auxsigargs and ',' or '', auxsigargs)

    result += '  static PyObject_t s_runcallback = NULL;\n'
    result += '  PyObject_t result;\n\n'

    result += '  if (%s) PyEval_AcquireThread( %s);\n'\
              % (pythread, pythread)

    result += '  if (NOT s_runcallback)\n'\
              '    s_runcallback = PypImport_ModuleAttr('\
              ' "%s", "swigcallback");\n' % __name__
    result += '  if (NOT s_runcallback) goto FUNEXIT;\n'

    result += '  result = PyObject_CallFunction( s_runcallback,'\
              ' "si%s", "%s", idx%s%s);\n'\
                % (isinstance( pyargs, types.StringType)
                   and (pyargs,
                        cbname,
                        callargs and ',' or '',
                        callargs)
                   or (pyargs[ 0], cbname, ',', pyargs[ 1]))
    result += '  if (NOT result) goto FUNEXIT;\n'
    result += '  Py_DECREF( result);\n'
  
    result += 'FUNEXIT:\n'
    result += '  if (PyErr_Occurred()) PypCallback_ProcessErr( "%s");\n'\
              % cbname
    result += '  if (%s) PyEval_ReleaseThread( %s);\n' % (pythread, pythread)
    result += '  return;\n}\n%}\n\n'

    result += '%%include gltypedef.i\n%%include %s\n' % headername

    result += '%s_f %s_array_get( int idx);\n\n' % ((cbname,) * 2)

    result += '%%constant int %s_dim = %s_DIM;\n' % ((cbname,) * 2)

    result += '%%constant %s_f %s_NULL = NULL;\n\n' % ((cbname,) * 2)

    return result

def cb_header( rettype, args):
    return 'typedef %s (%s*%s_f)(%s);\n\n'\
          % (rettype, winfpref, cb_name(rettype, args), cb_sigargs( args))

######################

swigfuns = {}

def swigfun( moduledef, sig, cbfun):
    assert cbfun is None or callable( cbfun)
    ftype = apply( cb_name, sig)
    if cbfun is None:
        return getattr( moduledef, '%s_NULL' % ftype)
    else:
        try:
            flist = swigfuns[ ftype]
        except KeyError:
            flist = swigfuns[ ftype] = []
        try:
            idx = flist.index( cbfun)
        except ValueError:
            idx = len( flist)
            nmax = getattr( moduledef, '%s_dim' % ftype)
            if nmax <= idx:
                raise ValueError( 'internal limit on registered python'
                                  'callback reached', nmax, idx)
            flist.append( cbfun)
        assert idx < getattr( moduledef, '%s_dim' % ftype)
        return getattr( moduledef, '%s_array_get' % ftype)( idx)


def swigcallback( ftype, idx, *args):
    '''callback from C (callback.i)
    '''
    try:
        return apply( swigfuns[ ftype][ idx], args)
    except SystemExit:
        raise
    except:
        import sys, traceback
        traceback.print_stack()
        traceback.print_exc()
        sys.exit( 1)

def isswigfun( obj):
    return not (obj is None or callable( obj)  )  

def swigptrcallback( moduledef, sig, callback):
    if not isswigfun( callback):
        callback = swigfun( moduledef, sig, callback)
    return callback

