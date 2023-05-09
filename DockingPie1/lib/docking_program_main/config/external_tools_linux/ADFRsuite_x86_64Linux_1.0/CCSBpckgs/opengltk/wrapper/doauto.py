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
# copyright_notice
#

from __future__ import nested_scopes

import string

typmap = {
    'b': 'GLbyte',
    's': 'GLshort',
    'i': 'GLint',
    'f': 'GLfloat',
    'd': 'GLdouble',
    'ub': 'GLubyte',
    'us': 'GLushort',
    'ui': 'GLuint',
    
    'Boolean': 'GLboolean',
    'Integer': 'GLint',
    'Float': 'GLfloat',
    'Double': 'GLdouble',
    }

pixeltypes = ['f', 'ui', 'us']

def gl_auto():

    print '"""automatically generated"""\n'
    print 'import numpy\n'
    print 'from opengltk.extent import _gllib\n'
    print 'from opengltk.util import glGetXXDim, %s\n'\
          % string.join( typmap.values(), ', ')
    
    allnames = []

    glGets = [
        ('glGet', 'Boolean', []),
        ('glGet', 'Integer', []),
        ('glGet', 'Float', []),
        ('glGet', 'Double', []),
        ]

    glGets += [('glGetLight', x, [('light', 'GLenum')]) for x in 'fi']
    glGets += [('glGetMaterial', x, [('face', 'GLenum')]) for x in 'fi']
    glGets += [('glGetTexEnv', x, [('target', 'GLenum')]) for x in 'fi']
    glGets += [('glGetTexGen', x, [('coord', 'GLenum')]) for x in 'dfi']
    glGets += [('glGetTexLevelParameter', x,
                [('target', 'GLenum'), ('level', 'GLint')])
               for x in 'fi']
    glGets += [('glGetTexParameter', x, [('target', 'GLenum')])
               for x in 'fi']

    def glGetNames( glgetlist):
        from operator import add
        return reduce( add, [['%s%s' % x[ :2], '%s%sv' % x[ :2]]
                             for x in glgetlist])

    def argsdocpr( args):
        return string.join( ['%s,' % x[ 0] for x in args]),\
               string.join( ['%s - %s' % x for x in args], '\n    ')

    def glGetFunv( name, typdesc, args):
        typ = typmap[ typdesc]
        argprint, argsdoc = argsdocpr( args)
        fname = '%s%s' % (name, typdesc)
        vcall = '%sv(%s pname' % (fname, argprint)
        return string.join( [
            '\ndef %s):' % vcall,
            '    """%s' % argsdoc,
            '\nreturn numpy array( %s)\n"""' % typ,
            '    result = numpy.zeros( glGetXXDim[ pname], %s)' % typ,
            '    _gllib.%s, result)' % vcall,
            '    return result',
            '\ndef %s(%s pname):' % (fname, argprint),
            '    """\n    %s' % argsdoc,
            'return numpy array( %s) or %s if singleton\n"""' %(typ,typ),
            '    result = %s)' % vcall,
            '    if 1 == len( result):',
            '        return result[ 0]',
            '    else:',
            '        return result',
            ] , '\n')
    
    allnames += glGetNames( glGets)
    print string.join( [ apply( glGetFunv, x) for x in glGets], '\n')

    coltypes = ['b', 'd', 'f', 'i', 's', 'ub', 'ui', 'us']

    sztypv = []
    sztypv += [('glColor', dim, typ) for dim in [3, 4] for typ in coltypes]
    sztypv += [('glEvalCoord', dim, typ) for dim in [1, 2] for typ in 'df']
    sztypv += [('glRasterPos', dim, typ)
               for dim in [2, 3, 4] for typ in 'dfis']
    sztypv += [('glTexCoord', dim, typ)
               for dim in [1, 2, 3, 4] for typ in 'dfis']
    sztypv += [('glVertex', dim, typ)
               for dim in [2, 3, 4] for typ in 'dfis']

    def sztypvname( name, dim, typ):
        return '%s%i%sv' % (name, dim, typ)

    allnames += [apply( sztypvname, x) for x in sztypv]

    for name, dim, typ in sztypv:
        fname = sztypvname( name, dim, typ)
        ctype = typmap[ typ]
        print '\ndef %s( v):' % fname
        print '    """\n    v - seq( %s, %i)\n    """' % (ctype, dim)
        print '    if %i != len( v):' % dim
        print '        raise TypeError( len( v), "%i-array expected")' % dim
        print '    _gllib.%s( v)' % fname

    indimarrays = [
        ('glDeleteTextures', 'GLuint'),
        ]
    allnames += [x[ 0] for x in indimarrays]
    for fun, typ in indimarrays:
        print '\ndef %s( seq):' % fun
        print '    """\n    vseq - sequence( %s)\n    """' % typ
        print '    return _gllib.%s( len( seq), seq)' % fun


    outdimarrays = [
        ('glGenTextures', 'GLuint'),
        ]
    allnames += [x[ 0] for x in outdimarrays]
    for fun, typ in outdimarrays:
        print '\ndef %s( n):' % fun
        print '    """\nreturn - sequence( %s, n)\n    """' % typ
        print '    result = numpy.zeros( n, %s)' % typ
        print '    _gllib.%s( n, result)' % fun
        print '    return result'


    ptrninputs = [('glIndex', 1, styp) for styp in ['d', 'f', 'i', 's', 'ub']]
    def ptrninputname( name, dim, styp):
        return '%s%sv' % (name, styp)
    def ptrninputfun( name, dim, styp):
        fname = ptrninputname( name, dim, styp)
        ctype = typmap[ styp]
        strs = [
            '\ndef %s( v):' % fname,
            '    """\n    v - seq( %s, %i)\n"""' % (ctype, dim),
            '    if %i != len( v):' % dim,
            '        raise TypeError( len( v), "%i-array expected")' % dim,
            '    _gllib.%s( v)' % fname,
            ]
        return string.join( strs, '\n')
    allnames += [apply( ptrninputname, x) for x in ptrninputs]
    for name, dim, styp in ptrninputs:
            print ptrninputfun( name, dim, styp)


    rectvs = [('glRect', 2, styp) for styp in 'dfis']
    def rectvname( name, dim, styp):
        return '%s%sv' % (name, styp)
    def rectvfun( name, dim, styp):
        fname = rectvname( name, dim, styp)
        ctype = typmap[ styp]
        strs = [
            '\ndef %s( v1, v2):' % fname,
            '    """\n    v1, v2 - seq( %s, %i)\n"""' % (ctype, dim),
            '    if %i != len( v1):' % dim,
            '        raise TypeError( len( v1), "%i-array expected for v1")'
                % dim,
            '    if %i != len( v2):' % dim,
            '        raise TypeError( len( v2), "%i-array expected for v2")'
                % dim,
            '    _gllib.%s( v1, v2)' % fname,
            ]
        return string.join( strs, '\n')
    allnames += [apply( rectvname, x) for x in rectvs]
    for name, dim, styp in rectvs:
            print rectvfun( name, dim, styp)


    inparms = []
    inparms += [('glLight', x, [('light', 'GLenum')]) for x in 'fi']
    inparms += [('glMaterial', x, [('face', 'GLenum')]) for x in 'fi']
    inparms += [('glLightModel', x, []) for x in 'fi']
    inparms += [('glFog', x, []) for x in 'fi']
    inparms += [('glTexEnv', x, [('target', 'GLenum')]) for x in 'fi']
    inparms += [('glTexGen', x, [('coord', 'GLenum')]) for x in 'fi']
    inparms += [('glTexParameter', x, [('target', 'GLenum')]) for x in 'fi']

    def inparmvname( name, typ):
        return '%s%sv' % (name, typ)

    def inparmfun( name, typ, args):
        ctype = typmap[ typ]
        argprint, argsdoc = argsdocpr( args)
        fname = inparmvname(name, typ)
        fcall = '%s(%s pname, parms)' % (fname, argprint)
        return string.join( [
            '\ndef %s:' % fcall,
            '    """%s' % argsdoc,
            '    parms - sequence\n"""',
            '    if len( parms) != glGetXXDim[ pname]:',
            '        raise TypeError( len( parms), glGetXXDim[ pname],',
            '                         "wrong size of parms")',
            '    _gllib.%s' % fcall,
            ] , '\n')
    allnames += [apply( inparmvname, x[ :2]) for x in inparms]
    for inparm in inparms:
        print apply( inparmfun, inparm)

        
##      invectfuns = [('glNormal', 3, x) for x in 'bdfis']
    invectfuns = []
    def invectname( name, dim, typ):
        return '%s%i%sv' % (name, dim, typ)
    def invectf( name, dim, typ):
        fname = invectname( name, dim, typ)
        ctype = typmap[ typ]
        return string.join( [
            '\ndef %s( v):' % fname,
            '    """\n    v - seq( %i, %s)\n    """' % (dim, ctype),
            '    if %i != len( v):' % dim,
            '        raise TypeError( len( v), "v must be a 3-array")',
            '    _gllib.%s( v)' % fname,
            ] , '\n')
    allnames += [apply( invectname, x) for x in invectfuns]
    for x in invectfuns:
        print apply( invectf, x)
        

    invtyp = [('glPixelMap', x, [('map', 'GLenum')]) for x in pixeltypes]
    def invtypname( name, typ, args):
        return '%s%sv' % (name, typ)
    allnames += [apply( invtypname, x) for x in invtyp]
    for name, typdesc, args in invtyp:
        typ = typmap[ typdesc]
        argprint, argsdoc = argsdocpr( args)
        fname = invtypname( name, typdesc, args)
        print '\ndef %s(%s values):' % (fname, argprint)
        print '    """\n    %s' % argsdoc
        print '    values - seq( %s)\n    """' % typ
        print '    _gllib.%s(%s values)' % (fname, argprint)

    print '\n__all__ = ["%s"]' % string.join( allnames, '",\n    "')

def gl_deprec():
    deprfuns = []
    
    deprfuns += [('glGetMap%sv' % x, '( target, query, v)', '') for x in 'dfi']
    
    deprfuns += [('glGetPixelMap%sv' % x, '( map, values)',
                  'use glGet to get value sz to return')
                 for x in pixeltypes]
    
    deprfuns += [('glGetPointerv', '( pname, params)',
                  'convert cobject to appropriate numpy array (need size...)')]
    
##      deprfuns += [('glGetPolygonStipple', '( mask)',
##                    'what is the expact output format ??')]

    deprfuns += [('glGetTexImage', '( target, level, format, type, pixels)',
                  '')]

##      deprfuns += [('glPolygonStipple', '( mask)',
##                    'what is the stipple format ??')]

    print '"""automatically generated"""\n'
    print '\n__all__ = ["%s"]' % string.join( [x[ 0] for x in deprfuns],
                                              '",\n    "')
    print '\nfrom warnings import warn\n'
    print 'from opengltk.extent import _gllib\n'
    print '\nclass ToWrap( RuntimeWarning):\n    pass\n'

    for fun, args, desc in deprfuns:
        print '\ndef %s%s:' % (fun, args)
        print '    """To Wrap: %s"""' % desc
        print '    warn( desc, ToWrap)'
        print '    return _gllib.%s%s' % (fun, args)
