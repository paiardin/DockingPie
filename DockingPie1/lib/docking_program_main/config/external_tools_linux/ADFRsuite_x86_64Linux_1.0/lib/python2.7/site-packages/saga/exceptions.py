
__author__    = "Andre Merzky, Ole Weidner"
__copyright__ = "Copyright 2012-2013, The SAGA Project"
__license__   = "MIT"


""" Exception classes
"""

import sys
import weakref
import operator
import traceback

# We have the choice of doing signature checks in exceptions, or to raise saga
# exceptions on signature checks -- we cannot do both.  At this point, we use
# the saga.exceptions in signatures, thus can *not* have signature checks
# here...
#
# import saga.base             as sb


# ------------------------------------------------------------------------------
#
class SagaException (Exception) :
    """
    The Exception class encapsulates information about error conditions
    encountered in SAGA.

    Additionally to the error message (e.message), the exception also provides
    a trace to the code location where the error condition got raised
    (e.traceback).

    B{Example}::

      try :
          file = saga.filesystem.File ("sftp://alamo.futuregrid.org/tmp/data1.dat")

      except saga.Timeout as to :
          # maybe the network is down?
          print "connection timed out"
      except saga.Exception as e :
          # something else went wrong
          print "Exception occurred: %s %s" % (e, e.traceback)

    There are cases where multiple backends can report errors at the same time.
    In that case, the saga-python implementation will collect the exceptions,
    sort them by their 'rank', and return the highest ranked one.  All other
    catched exceptions are available via :func:`get_all_exceptions`, or via the
    `exceptions` property.

    The rank of an exception defines its explicity: in general terms: the higher
    the rank, the better defined / known is the cause of the problem.
    """

    # --------------------------------------------------------------------------
    #
    def __init__ (self, msg, parent=None, api_object=None, from_log=False) :
        """ 
        Create a new exception object.

        :param msg:         The exception message.
        :param parent:      Original exception
        :param api_object:  The object that has caused the exception, default is
                            None.
        :param from_log:    Exception c'tor originates from the static log_
                            member method (ignore in exception stack!)
        """
        Exception.__init__(self, msg)

        self._plain_message = msg
        self._exceptions    = [self]
        self._top_exception = self
        self._ptype         = type(parent).__name__   # parent exception type
        self._stype         = type(self  ).__name__   # own exception    type 

        ignore_stack = 2
        if  from_log : 
            ignore_stack += 1


        if api_object : 
            self._object    = weakref.ref (api_object)
        else :
            self._object    = None


        # did we get a parent exception?
        if  parent :

            # if so, then this exception is likely created in some 'except'
            # clause, as a reaction on a previously catched exception (the
            # parent).  Thus we append the message of the parent to our own
            # message, but keep the parent's traceback (after all, the original
            # exception location is what we are interested in).
            #
            if  isinstance (parent, SagaException) :
                # that all works nicely when parent is our own exception type...
                self._traceback = parent.traceback

                frame           = traceback.extract_stack ()[- ignore_stack]
                line            = "%s +%s (%s)  :  %s" % frame 
                self._message   = "  %-20s: %s (%s)\n%s" \
                                % (self._stype, msg, line, parent.msg)

            else :
                if self._stype != "NoneType" :
                    # ... but if parent is a native (or any other) exception type,
                    # we don't have a traceback really -- so we dig it out of
                    # sys.exc_info. 
                    trace           = sys.exc_info ()[2]
                    stack           = traceback.extract_tb  (trace)
                    traceback_list  = traceback.format_list (stack)
                    self._traceback = "".join (traceback_list)

                    # the message composition is very similar -- we just inject the
                    # parent exception type inconspicuously somewhere (above that
                    # was part of 'parent.message' already).
                    frame           = traceback.extract_stack ()[- ignore_stack]
                    line            = "%s +%s (%s)  :  %s" % frame 
                    self._message   = "  %-20s: %s (%s)\n  %-20s: %s" \
                                    % (self._stype, msg, line, self._ptype, parent)

        else :

            # if we don't have a parent, we are a 1st principle exception,
            # i.e. a reaction to some genuine code error.  Thus we extract the
            # traceback from exactly where we are in the code (the last stack
            # frame will be the call to this exception constructor), and we
            # create the original exception message from 'stype' and 'message'.
            stack           = traceback.extract_stack ()
            traceback_list  = traceback.format_list (stack)
            self._traceback = "".join (traceback_list[:-1])
            frame           = traceback.extract_stack ()[- ignore_stack -1]
            line            = "%s +%s (%s)  :  %s" % frame 
            self._message   = "%s (%s)" % (msg, line)

        # we can't do that earlier as _msg was not set up before
        self._messages = [self._message]


    # --------------------------------------------------------------------------
    #
    def __str__ (self) :
        return self.get_message ()


    # --------------------------------------------------------------------------
    #
    def __repr__ (self) :
        return "%s\n%s" % (self._message, self._traceback)


    # --------------------------------------------------------------------------
    #
    def _clone (self) :
        """ This method is used internally -- see :func:`_get_exception_stack`."""

        clone = self.__class__ ("")

        clone._parent    = self._parent
        clone._object    = self._object
        clone._message   = self._message
        clone._messages  = self._messages
        clone._exception = self._exceptions
        clone._traceback = self._traceback
        clone._stype     = self._stype
        clone._ptype     = self._ptype

        return clone

    # --------------------------------------------------------------------------
    #
    @classmethod
    def _log (cls, logger, msg, parent=None, api_object=None, level='error'):
        """ this class method allows to log the exception message while
            constructing a SAGA exception, like::

              # raise an exception, no logging
              raise saga.IncorrectState ("File is not open")

              # raise an exception, log as error event (error level is default)
              raise saga.IncorrectState._log (self._logger, "File is not open")

              # raise an exception, log as warning event
              raise saga.IncorrectState._log (self._logger, "File is not open", level=warning)
              raise saga.IncorrectState._log (self._logger, "File is not open", warning) # same

            This way, the 'raise' remains clearly in the code, as that is the
            dominating semantics of the call.
        """

        log_method = logger.error

        try :
            log_method = getattr (logger, level.lower())
        except :
            sys.stderr.write ("unknown log level '%s'"  %  level)

        log_method ("%s: %s" % (cls.__name__, msg))

        return cls (msg, parent=parent, api_object=api_object, from_log=True)



    # --------------------------------------------------------------------------
    #
    def get_message (self) :
        """ Return the exception message as a string.  That message is also
        available via the 'message' property."""
        return self._message


    # --------------------------------------------------------------------------
    #
    def _get_plain_message (self) :
        """ Return the plain error message as a string. """
        return self._message


    # --------------------------------------------------------------------------
    #
    def get_type (self):
        """ Return the type of the exception as string.
        """
        return self._stype


    # --------------------------------------------------------------------------
    #
    def get_object (self) :
        """ Return the object that raised this exception. An object may not
        always be available -- for example, exceptions raised during object
        creation may not have the option to keep an incomplete object instance
        around.  In those cases, this method will return 'None'.  Either way,
        the object is also accessible via the 'object' property.
        """

        # object is a weak_ref, and may have been garbage collected - we simply
        # return 'None' then
        return self._object ()


    # --------------------------------------------------------------------------
    #
    def _add_exception (self, e) :
        """
        Some sub-operation raised a SAGA exception, but other exceptions may
        be catched later on.  In that case the later exceptions can be added to
        the original one with :func:`_add_exception`\(e).  Once all exceptions are
        collected, a call to :func:`_get_exception_stack`\() will return a new
        exception which is selected from the stack by rank and order.  All other
        exceptions can be accessed from the returned exception by
        :func:`get_all_exceptions`\() -- those exceptions are then also ordered 
        by rank.
        """

        self._exceptions.append (e)
        self._messages.append   (e.message)

        if e._rank > self._top_exception._rank :
            self._top_exception = e


    # --------------------------------------------------------------------------
    #
    def _get_exception_stack (self) :
        """ 
        This method is internally used by the saga-python engine, and is only
        relevant for operations which (potentially) bind to more than one
        adaptor.
        """

        if self._top_exception == self :
            return self
        
        # damned, we can't simply recast ourself to the top exception type -- so
        # we have to create a new exception with that type, and copy all state
        # over...

        # create a new exception with same type as top_exception
        clone = self._top_exception._clone ()
        clone._exceptions = []
        clone._messages   = []

        # copy all state over
        for e in sorted (self._exceptions, key=operator.attrgetter ('_rank'), reverse=True) :
            clone._exceptions.append (e)
            clone._messages.append (e._message)

        return clone

                
    # --------------------------------------------------------------------------
    #
    def get_all_exceptions (self) :
        return self._exceptions


    # --------------------------------------------------------------------------
    #
    def get_all_messages (self) :
        return self._messages


    # --------------------------------------------------------------------------
    #
    def get_traceback (self) :
        return self._traceback


    # --------------------------------------------------------------------------
    #
    message        = property (get_message)         # string
    object         = property (get_object)          # object type
    type           = property (get_type)            # exception type
    exceptions     = property (get_all_exceptions)  # list [Exception]
    messages       = property (get_all_messages)    # list [string]
    traceback      = property (get_traceback)       # string


# ------------------------------------------------------------------------------
#
class NotImplemented(SagaException):
    """ SAGA-Python does not implement this method or class. (rank: 11)"""

    _rank = 11

    def __init__ (self, msg, parent=None, api_object=None, from_log=False) :
        SagaException.__init__ (self, msg, parent, api_object, from_log)


# ------------------------------------------------------------------------------
#
class IncorrectURL(SagaException): 
    """ The given URL could not be interpreted, for example due to an incorrect
        / unknown schema. (rank: 10)"""

    _rank = 10
    
    def __init__ (self, msg, parent=None, api_object=None, from_log=False) :
        SagaException.__init__ (self, msg, parent, api_object, from_log)


# ------------------------------------------------------------------------------
#
class BadParameter(SagaException): 
    """ A given parameter is out of bound or ill formatted. (rank: 9)"""

    _rank = 9
    
    def __init__ (self, msg, parent=None, api_object=None, from_log=False) :
        SagaException.__init__ (self, msg, parent, api_object, from_log)


# ------------------------------------------------------------------------------
#
class AlreadyExists(SagaException): 
    """ The entity to be created already exists. (rank: 8)"""

    _rank = 8
    
    def __init__ (self, msg, parent=None, api_object=None, from_log=False) :
        SagaException.__init__ (self, msg, parent, api_object, from_log)


# ------------------------------------------------------------------------------
#
class DoesNotExist(SagaException):
    """ An operation tried to access a non-existing entity. (rank: 7)"""

    _rank = 7
    
    def __init__ (self, msg, parent=None, api_object=None, from_log=False) :
        SagaException.__init__ (self, msg, parent, api_object, from_log)


# ------------------------------------------------------------------------------
#
class IncorrectState(SagaException): 
    """ The operation is not allowed on the entity in its current state. (rank: 6)"""

    _rank = 6
    
    def __init__ (self, msg, parent=None, api_object=None, from_log=False) :
        SagaException.__init__ (self, msg, parent, api_object, from_log)


# ------------------------------------------------------------------------------
#
class PermissionDenied(SagaException): 
    """ The used identity is not permitted to perform the requested operation. (rank: 5)"""

    _rank = 5
    
    def __init__ (self, msg, parent=None, api_object=None, from_log=False) :
        SagaException.__init__ (self, msg, parent, api_object, from_log)


# ------------------------------------------------------------------------------
#
class AuthorizationFailed(SagaException): 
    """ The backend could not establish a valid identity. (rank: 4)"""

    _rank = 4
    
    def __init__ (self, msg, parent=None, api_object=None, from_log=False) :
        SagaException.__init__ (self, msg, parent, api_object, from_log)


# ------------------------------------------------------------------------------
#
class AuthenticationFailed(SagaException): 
    """ The backend could not establish a valid identity. (rank: 3)"""

    _rank = 3
    
    def __init__ (self, msg, parent=None, api_object=None, from_log=False) :
        SagaException.__init__ (self, msg, parent, api_object, from_log)


# ------------------------------------------------------------------------------
#
class Timeout(SagaException): 
    """ The interaction with the backend times out. (rank: 2)"""

    _rank = 2
    
    def __init__ (self, msg, parent=None, api_object=None, from_log=False) :
        SagaException.__init__ (self, msg, parent, api_object, from_log)


# ------------------------------------------------------------------------------
#
class NoSuccess(SagaException): 
    """ Some other error occurred. (rank: 1)"""

    _rank = 1
    
    def __init__ (self, msg, parent=None, api_object=None, from_log=False) :
        SagaException.__init__ (self, msg, parent, api_object, from_log)




