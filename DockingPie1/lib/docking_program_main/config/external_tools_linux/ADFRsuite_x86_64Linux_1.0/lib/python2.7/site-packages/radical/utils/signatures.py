
__author__    = "Radical.Utils Development Team (Andre Merzky)"
__copyright__ = "Copyright 2013, RADICAL@Rutgers"
__license__   = "MIT"

import sys
import os


# ------------------------------------------------------------------------------
#
# Method call parameters/return value type checking decorators.
# (c) 2006-2007, Dmitry Dvoinikov <dmitry@targeted.org>
# Distributed under BSD license.
#
# code adjustements by 'The SAGA Project', 2013
# See also 
# http://code.activestate.com/recipes/577065-type-checking-function-overloading-decorator/
#
#
# Samples:
#
# from typecheck import *
#
# @takes(int, str) # takes   int, str, upon a problem throws TypeError
# @returns(int)    # returns int,      upon a problem throws TypeError
# def foo(i, s): 
#     return i + len(s)
#
# @takes((int, long), by_regex("^[0-9]+$")) # int or long, numerical string
# def foo(i, s, anything):                  # and the third parameter is not checked
#     ...
#
# @takes(int, int, foo = int, bar = optional(int)) # keyword argument foo must be int
# def foo(a, b, **kwargs):                         # bar may be int or missing
#     ...
#
# Note: @takes for positional arguments, @takes for keyword arguments and @returns
# all support the same checker syntax, for example for the following declaration
# 
# @takes(C)
# def foo(x):
#     ...
#
# then C may be one of the simple checkers:
#
# --------- C ---------     ------------- semantics -------------
# typename              ==> ok if x is is an instance of typename
# "typename"            ==> ok if x is is an instance of typename
# with_attr("a", "b")   ==> ok if x has specific attributes
# some_callable         ==> ok if some_callable(x) is True
# one_of(1, "2")        ==> ok if x is one of the literal values
# by_regex("^foo$")     ==> ok if x is a matching basestring
# nothing               ==> ok if x is None
# anything              ==> always ok
#
# simple checkers can further be combined with OR semantics using tuples:
#
# --------- C ---------     ------------- semantics -------------
# (checker1, checker2)      ==> ok if x conforms with either checker
#
# be optional:
#
# --------- C ---------     ------------- semantics -------------
# optional (checker)        ==> ok if x is checker-conformant or None
#
# or nested recursively into one of the following checkers
#
# --------- C ---------     ------------- semantics -------------
# list_of  (checker)        ==> ok if x is a list  of checker-conformant values
# tuple_of (checker)        ==> ok if x is a tuple of checker-conformant values
# set_of   (checker)        ==> ok if x is a set   of checker-conformant values
# dict_of  (key_checker, 
#           value_checker)  ==> ok if x is a dict  mapping key_checker-
#                           conformant keys to value_checker-conformant values
#
# More samples:
#
# class foo (object):
#     @takes ("foo", optional (int)) # foo, maybe int, but foo is yet incomplete
#     def __init__ (self, i = None): # and is thus specified by name
#         ...
#     @takes ("foo", int)            # foo, and int if presents in args,
#     def bar (self, *args):         # if args is empty, the check passes ok
#         ...
#     @takes ("foo")                 
#     @returns (object)              # returns foo which is fine, because
#     def biz (self):                # foo is an object
#         return self
#     @classmethod                   # classmethod's and staticmethod's
#     @takes (type)                  # go same way
#     def baz (cls):
#         ...
#    
# @takes   (int)
# @returns (optional ("int", foo))   # returns either int, foo or NoneType
# def bar (i):                       # "int" (rather than just int) is for fun
#     if i > 0: 
#         return i
#     elif i == 0:
#         return foo()               # otherwise returns NoneType
#
# @takes   (callable)                # built-in functions are treated as predicates
# @returns (lambda x: x == 123)      # and so do user-defined functions or lambdas
# def exec  (f, *args, **kwargs):
#     return f (*args, **kwargs)
#
# assert execute (execute, execute, execute, lambda x: x, 123) == 123
#
# def readable (x):                  # user-defined type-checking predicate
#     return hasattr (x, "read")
#
# anything is an alias for predicate lambda: True,
# nothing is an alias for NoneType, as in:
#
# @takes   (callable, readable, 
#           optional (anything), 
#           optional (int))
# @returns (nothing)
# def foo (f, r, x = None, i = None):
#     ...
#
#                                    # another way of protocol checking
# @takes (with_attr ("read", "write")) 
# def foo (pipe):
#     ...
#
# @takes (list_of (int))             # list of ints
# def foo (x):
#     print x[0]
#
# @takes (tuple_of (callable))       # tuple of callables
# def foo(x):
#     print x[0]()
#
#                                    # dict mapping strs to lists of int
# @takes (dict_of (str, list_of (int))) 
# def foo (x):
#     print sum (x["foo"])
#
# @takes (by_regex ("^[0-9]{1,8}$")) # integer-as-a-string regex
# def foo (x):
#     i = int (x)
#
# @takes (one_of (1, 2))             # must be equal to either one
# def set_version (version):
#     ...
#
# The (3 times longer) source code with self-tests is available from:
# http://www.targeted.org/python/recipes/typecheck.py
#
# ------------------------------------------------------------------------------

__all__ = [ "takes",    "returns",   "optional", "nothing",   
            "anything", "list_of",   "tuple_of", "dict_of",
            "by_regex", "with_attr", "one_of",   "set_of" ]

no_check        = True  # set this to True to turn all checks off
no_return_check = True   # set this to True to turn return value cchecks off

if 'RADICAL_DEBUG_SIG' in os.environ:
    no_check    = False  # set this to True to turn all checks off

# ------------------------------------------------------------------------------

from traceback import extract_stack
from inspect   import getargspec, isclass
from types     import NoneType
from re        import compile as regex
from functools import wraps

# ------------------------------------------------------------------------------

def base_names (C):
    "Returns list of base class names for a given class"
    return [ x.__name__ for x in C.__mro__ ]
    
# ------------------------------------------------------------------------------

def type_name (v):
    "Returns the name of the passed value's type"
    return type (v).__name__

# ------------------------------------------------------------------------------

class Checker (object):

    def __init__ (self, reference):
        self.reference = reference
        self.spectype  = reference

    def check (self, value): # abstract
        pass

    _registered = [] # a list of registered descendant class factories

    @staticmethod
    def create (value): # static factory method
        for f, t in Checker._registered:
            if f (value):
                return t (value)
        else:
            return None

# ------------------------------------------------------------------------------

class TypeChecker (Checker):

    def check (self, value):
        return isinstance (value, self.reference)

Checker._registered.append ((isclass, TypeChecker))
nothing = NoneType

# ------------------------------------------------------------------------------

class StrChecker (Checker):

    def check (self, value):
        value_base_names = base_names (type (value))
        return self.reference in value_base_names or \
               "instance"     in value_base_names
   
Checker._registered.append ((lambda x: isinstance (x, str), StrChecker))

# ------------------------------------------------------------------------------

class TupleChecker (Checker):

    def __init__ (self, reference):
        self.reference = map (Checker.create, reference)
        self.spectype  = reference

    def check (self, value):
        return reduce (lambda r, c: r or c.check (value), self.reference, False)

Checker._registered.append ((lambda x: isinstance (x, tuple) and not
                     filter (lambda y: Checker.create (y) is None, x), 
                     TupleChecker))

optional = lambda *args: args + (NoneType, )

# ------------------------------------------------------------------------------

class CallableChecker (Checker):

    def check (self, value):
        return self.reference (value)

# note that the callable check is the most relaxed of all, therefore it should
# be registered last, after all the more specific cases have been registered

Checker._registered.append ((callable, CallableChecker))

anything = lambda *args: True

# ------------------------------------------------------------------------------

class ListOfChecker (Checker):

    def __init__ (self, reference):
        self.reference = Checker.create (reference)

    def check (self, value):
        return isinstance (value, list) and \
               not filter (lambda e: not self.reference.check (e), value)

list_of = lambda *args: ListOfChecker (*args).check

# ------------------------------------------------------------------------------

class TupleOfChecker (Checker):

    def __init__ (self, reference):
        self.reference = Checker.create (reference)

    def check (self, value):
        return isinstance (value, tuple) and \
               not filter (lambda e: not self.reference.check (e), value)

tuple_of = lambda *args: TupleOfChecker (*args).check

# ------------------------------------------------------------------------------

class SetOfChecker (Checker):

    def __init__ (self, reference):
        self.reference = Checker.create (reference)

    def check (self, value):
        return isinstance (value, set) and \
               not filter (lambda e: not self.reference.check (e), value)

set_of = lambda *args: SetOfChecker (*args).check

# ------------------------------------------------------------------------------

class DictOfChecker (Checker):

    def __init__ (self, key_reference, value_reference):
        self.key_reference   = Checker.create (key_reference)
        self.value_reference = Checker.create (value_reference)

    def check (self, value):
        return isinstance (value, dict) and \
               not filter (lambda e: not self.key_reference.check   (e), value.iterkeys  ()) and \
               not filter (lambda e: not self.value_reference.check (e), value.itervalues())

dict_of = lambda *args: DictOfChecker (*args).check

# ------------------------------------------------------------------------------

class RegexChecker (Checker):

    def __init__(self, reference):
        self.reference = regex (reference)

    def check (self, value):
        return isinstance (value, basestring) and self.reference.match (value)

by_regex = lambda *args: RegexChecker (*args).check

# ------------------------------------------------------------------------------

class AttrChecker (Checker):

    def __init__ (self, *attrs):
        self.attrs = attrs

    def check (self, value):
        return reduce (lambda r, c: r and c, 
                  map (lambda a: hasattr (value, a), self.attrs), True)

with_attr = lambda *args: AttrChecker (*args).check

# ------------------------------------------------------------------------------

class OneOfChecker (Checker):

    def __init__ (self, *values):
        self.values = values

    def check (self, value):
        return value in self.values

one_of = lambda *args: OneOfChecker (*args).check


def create_return_exception (method, spectype, result) :

    stack = extract_stack ()
    for f in stack : 
        if  'utils/signatures.py' in f[0] :
            break
        frame = f

    msg   = "\nSignature Mismatch\n"
    msg  += "  in function   : %s\n"       % (frame[2])
    msg  += "  in file       : %s +%s\n"   % (frame[0], frame[1])
    msg  += "  on line       : %s\n"       % (frame[3])
    msg  += "  method        : %s\n"       % (method.__name__)
    msg  += "  returned type : %s\n"       % (type(result).__name__)
    msg  += "  instead  of   : %s\n"       % (spectype.__name__)
    msg  += "  This is an internal error!"

    return TypeError (msg)


def create_type_exception (method, arg0, i, arg, spectype, kwname="") :

    narg = 0
    
    if  arg0 and isinstance (arg0, object) :
        narg = i
    else :
        narg = i+1
    
    stack = extract_stack ()
    for f in stack : 
        if  'utils/signatures.py' in f[0] :
            break
        frame = f

    if isinstance(spectype,tuple) or \
       isinstance(spectype,list )  :
        expected = list()
        for t in spectype :
            if isinstance(t, type):
                expected.append (t.__name__)
            else:
                expected.append (str(t))

    elif isinstance(spectype,type):
        expected = spectype.__name__

    else:
        expected = str(spectype)

    msg   = "\nSignature Mismatch\n"
    msg  += "  in function        : %s\n"       % (frame[2])
    msg  += "  in file            : %s +%s\n"   % (frame[0], frame[1])
    msg  += "  on line            : %s\n"       % (frame[3])
    msg  += "  method             : %s\n"       % (method.__name__)
    if  not kwname :
        msg  += "  argument           : #%s\n"  % (narg)
    else :
        msg  += "  parameter          : %s"     % (kwname)
    msg  += "  has incorrect type : %s\n"       % (type (arg).__name__)
    msg  += "  instead of         : %s"         % (str(expected))

    return TypeError (msg)


# ------------------------------------------------------------------------------

def takes (*args, **kwargs):
    "Method signature checking decorator"

    if no_check: # no type checking is performed, return decorated method itself

        def takes_proxy (method):
            return method

    else:

        # convert decorator arguments into a list of checkers

        checkers = []
        for i, arg in enumerate (args):
            checker = Checker.create (arg)
            if checker is None:
                raise TypeError ("@takes decorator got parameter %d of unsupported "
                                 "type %s" % (i + 1, type_name(arg)))
            checkers.append (checker)

        kwcheckers = {}
        for kwname, kwarg in kwargs.iteritems ():
            checker = Checker.create (kwarg)
            if checker is None:
                raise TypeError ("@takes decorator got parameter %s of unsupported "
                                 "type %s" % (kwname, type_name (kwarg)))
            kwcheckers[kwname] = checker


        def takes_proxy (method):
            
            method_args, method_defaults = getargspec (method)[0::3]

            @wraps(method)
            def signature_check (*pargs, **pkwargs):
    
                # append the default parameters

                if  method_defaults is not None and len (method_defaults) > 0 \
                    and len (method_args) - len (method_defaults) <= len (pargs) < len (method_args):
                    pargs += method_defaults[len (pargs) - len (method_args):]

                # check the types of the actual call parameters

                for i, (arg, checker) in enumerate (zip (pargs, checkers)):
                    if  not checker.check (arg):
                      # print 'Checker.spectype %s' % checker.spectype
                        raise (create_type_exception (method, pargs[0], i,
                                                      arg, checker.spectype))

                for kwname, checker in kwcheckers.iteritems ():
                    if  not checker.check (pkwargs.get (kwname, None)):
                      # print 'checker.spectype %s' % checker.spectype
                        raise (create_type_exception (method, pargs[0], i,
                                                      arg, checker.spectype, kwname))

                try :
                    return method(*pargs, **pkwargs)
                except :
                    # remove signature decorator from exception call stack
                    et, ei, tb = sys.exc_info()
                    raise et, ei, tb.tb_next

            signature_check.__name__ = method.__name__
            return signature_check
    
    return takes_proxy

# ------------------------------------------------------------------------------

def returns (sometype):
    "Return type checking decorator"

    if no_check: # no type checking is performed, return decorated method itself

        def returns_proxy (method):
            return method

    else:

        # convert decorator argument into a checker
        
        checker = Checker.create(sometype)
        if checker is None:
            raise TypeError ("@returns decorator got parameter of unsupported "
                             "type %s" % type_name (sometype))

        def returns_proxy (method):
            
            @wraps(method)
            def signature_check (*args, **kwargs):
                
                try :
                    result = method (*args, **kwargs)
                except :
                    # remove signature decorator from exception call stack
                    et, ei, tb = sys.exc_info()
                    raise et, ei, tb.tb_next
                
                if  not checker.check (result):
                    
                    if not  no_return_check:
                      # print 'Checker.spectype %s' % checker.spectype
                        raise (create_return_exception (method, checker.spectype, result))

                return result
    
            signature_check.__name__ = method.__name__
            return signature_check
        
    return returns_proxy


# ------------------------------------------------------------------------------


