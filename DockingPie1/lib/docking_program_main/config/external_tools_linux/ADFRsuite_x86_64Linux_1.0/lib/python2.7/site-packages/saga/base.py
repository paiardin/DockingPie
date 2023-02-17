
__author__    = "Andre Merzky, Ole Weidner"
__copyright__ = "Copyright 2012-2013, The SAGA Project"
__license__   = "MIT"


import sys
import string
import inspect

import radical.utils              as ru
import radical.utils.signatures   as rus
import radical.utils.logger       as rul

import saga.engine.engine
import saga.adaptors.base         as sab

# ------------------------------------------------------------------------------
#
class SimpleBase (object) :
    """ This is a very simple API base class which just initializes
    the self._logger and self._engine members, but does not perform any further
    initialization, nor any adaptor binding.  This base is used for API classes
    which are not backed by multiple adaptors (no session, tasks, etc).
    """

    # --------------------------------------------------------------------------
    #
    @rus.takes   ('SimpleBase')
    @rus.returns (rus.nothing)
    def __init__  (self) :

        if  not hasattr (self, '_apitype') :
            self._apitype = self._get_apitype ()

        self._logger = ru.get_logger  ('radical.saga')
        self._id     = ru.generate_id (self._get_apitype (), mode=ru.ID_SIMPLE)

      # self._logger.debug ("[saga.Base] %s.__init__()" % self._apitype)


    # --------------------------------------------------------------------------
    #
    @rus.takes   ('SimpleBase')
    @rus.returns (basestring)
    def _get_apitype (self) :

        # apitype for saga.job.service.Service should be saga.job.Service --
        # but we need to make sure that this actually exists and is equivalent.

        mname_1 = self.__module__           # saga.job.service
        cname   = self.__class__.__name__   # Service

        mname_1_elems = mname_1.split ('.') # ['saga', 'job', 'service']
        mname_1_elems.pop ()                # ['saga', 'job']
        mname_2 = '.'.join (mname_1_elems)  #  'saga.job'

        for mod_cname, mod_class in inspect.getmembers (sys.modules[mname_2]) :
            if  mod_cname == cname           and \
                inspect.isclass (mod_class)  and \
                isinstance (self, mod_class)     :

                apitype = "%s.%s" % (mname_2, cname) # saga.job.Service
                return apitype

        apitype = "%s.%s" % (mname_1, cname) # saga.job.service.Service
        return apitype

    
# ------------------------------------------------------------------------------
#
class Base (SimpleBase) :

    # --------------------------------------------------------------------------
    #
    @rus.takes   ('Base',
                  basestring, 
                  rus.optional (sab.Base), 
                  rus.optional (dict), 
                  rus.optional (rus.anything),
                  rus.optional (rus.anything))
    @rus.returns (rus.nothing)
    def __init__  (self, schema, adaptor, adaptor_state, *args, **kwargs) :

        SimpleBase.__init__ (self)

        _engine       = saga.engine.engine.Engine ()

        self._adaptor = adaptor
        self._adaptor = _engine.bind_adaptor (self, self._apitype, schema, adaptor)

        # Sync creation (normal __init__) will simply call the adaptor's
        # init_instance at this point.  _init_task should *not* be evaluated,
        # ever, for __init__ based construction! (it is private, see?)
        #
        # For any async creation (create()), we need to return a task which
        # performs initialization.  We rely on the sync/async method decorators
        # on CPI level to provide the task instance itself, and point the task's
        # workload to the adaptor level init_instance method.

        self._init_task = self._adaptor.init_instance (adaptor_state, *args, **kwargs)

        if 'ttype' in kwargs and kwargs['ttype'] :
            # in this case we in in fact need the init_task later on, to return
            # on the create() call
            pass
        else :
            # in the sync case, we can get rid of the init_task reference, to
            # simplify garbage collection
            self._init_task = None


    # --------------------------------------------------------------------------
    #
    @rus.takes   ('Base')
    @rus.returns ('saga.Session')
    def get_session (self) :
        """ 
        Returns the session which is managing the object instance.  For objects
        which do not accept a session handle on construction, this call returns
        None.

        The object's session is also available via the `session` property.
        """
        return self._adaptor.get_session ()

    session = property (get_session)

#
# ------------------------------------------------------------------------------




