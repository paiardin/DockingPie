
__author__    = "Andre Merzky, Ole Weidner, Thomas Schatz"
__copyright__ = "Copyright 2012-2013, The SAGA Project"
__license__   = "MIT"


""" SAGA job description interface """

import radical.utils.signatures as rus

import saga


#-------------------------------------------------------------------------------
#
class Description (saga.Attributes) :
    """ The job description class. """

    # --------------------------------------------------------------------------
    #
    @rus.takes   ('Description')
    @rus.returns (rus.nothing)
    def __init__(self):

        # set attribute interface properties

        import saga.attributes as sa

        self._attributes_extensible  (True)
        self._attributes_camelcasing (True)

        # register properties with the attribute interface

        self._attributes_register  (saga.job.EXECUTABLE           , None, sa.STRING, sa.SCALAR, sa.WRITEABLE)
        self._attributes_register  (saga.job.PRE_EXEC             , None, sa.STRING, sa.VECTOR, sa.WRITEABLE)
        self._attributes_register  (saga.job.POST_EXEC            , None, sa.STRING, sa.VECTOR, sa.WRITEABLE)
        self._attributes_register  (saga.job.ARGUMENTS            , None, sa.STRING, sa.VECTOR, sa.WRITEABLE)
        self._attributes_register  (saga.job.ENVIRONMENT          , None, sa.STRING, sa.DICT,   sa.WRITEABLE)
        self._attributes_register  (saga.job.TOTAL_CPU_COUNT      , None, sa.INT,    sa.SCALAR, sa.WRITEABLE)
        self._attributes_register  (saga.job.NUMBER_OF_PROCESSES  , None, sa.INT,    sa.SCALAR, sa.WRITEABLE)
        self._attributes_register  (saga.job.PROCESSES_PER_HOST   , None, sa.INT,    sa.SCALAR, sa.WRITEABLE)
        self._attributes_register  (saga.job.THREADS_PER_PROCESS  , None, sa.INT,    sa.SCALAR, sa.WRITEABLE)
        self._attributes_register  (saga.job.WORKING_DIRECTORY    , None, sa.STRING, sa.SCALAR, sa.WRITEABLE)
        self._attributes_register  (saga.job.INTERACTIVE          , None, sa.BOOL,   sa.SCALAR, sa.WRITEABLE)
        self._attributes_register  (saga.job.INPUT                , None, sa.STRING, sa.SCALAR, sa.WRITEABLE)
        self._attributes_register  (saga.job.OUTPUT               , None, sa.STRING, sa.SCALAR, sa.WRITEABLE)
        self._attributes_register  (saga.job.ERROR                , None, sa.STRING, sa.SCALAR, sa.WRITEABLE)
        self._attributes_register  (saga.job.FILE_TRANSFER        , None, sa.STRING, sa.VECTOR, sa.WRITEABLE)
        self._attributes_register  (saga.job.CLEANUP              , None, sa.BOOL,   sa.SCALAR, sa.WRITEABLE)
        self._attributes_register  (saga.job.JOB_START_TIME       , None, sa.TIME,   sa.SCALAR, sa.WRITEABLE)
        self._attributes_register  (saga.job.WALL_TIME_LIMIT      , None, sa.INT,    sa.SCALAR, sa.WRITEABLE)
        self._attributes_register  (saga.job.TOTAL_PHYSICAL_MEMORY, None, sa.INT,    sa.SCALAR, sa.WRITEABLE)
        self._attributes_register  (saga.job.CPU_ARCHITECTURE     , None, sa.ENUM,   sa.SCALAR, sa.WRITEABLE)
        self._attributes_register  (saga.job.OPERATING_SYSTEM_TYPE, None, sa.ENUM,   sa.SCALAR, sa.WRITEABLE)
        self._attributes_register  (saga.job.CANDIDATE_HOSTS      , None, sa.STRING, sa.VECTOR, sa.WRITEABLE)
        self._attributes_register  (saga.job.QUEUE                , None, sa.STRING, sa.SCALAR, sa.WRITEABLE)
        self._attributes_register  (saga.job.NAME                 , None, sa.STRING, sa.SCALAR, sa.WRITEABLE)
        self._attributes_register  (saga.job.PROJECT              , None, sa.STRING, sa.SCALAR, sa.WRITEABLE)
        self._attributes_register  (saga.job.JOB_CONTACT          , None, sa.STRING, sa.VECTOR, sa.WRITEABLE)
        self._attributes_register  (saga.job.SPMD_VARIATION       , None, sa.ENUM,   sa.SCALAR, sa.WRITEABLE)
      # self._attributes_set_enums (saga.job.SPMD_VARIATION,      ['MPI', 'OpenMP', 'MPICH-G'])

        self._env_is_list = False

        self._attributes_set_getter (saga.job.ENVIRONMENT, self._get_env)
        self._attributes_set_setter (saga.job.ENVIRONMENT, self._set_env)


    # --------------------------------------------------------------------------
    #
    def _set_env (self, val) :
        if  isinstance (val, list) :
            self._env_is_list = True


    # --------------------------------------------------------------------------
    #
    def _get_env (self) :
        env = self.get_attribute (saga.job.ENVIRONMENT)
        if  self._env_is_list :
            self._env_is_list = False
            return ["%s=%s" % (key, val) for (key, val) in env.items ()]
        return env


    # --------------------------------------------------------------------------
    #
    @rus.takes   ('Description',
                  ('Description', dict))
    @rus.returns ('Description')
    def __deepcopy__ (self, memo) :
        other = saga.job.Description ()
        return self.clone (other)


    # --------------------------------------------------------------------------
    #
    @rus.takes   ('Description',
                  'Description')
    @rus.returns ('Description')
    def clone (self, other=None) :
        """
        clone()

        Implements deep copy. u

        Unlike the default python assignment (copy object reference),
        a deep copy will create a new object instance with the same state --
        after a deep copy, a change on one instance will not affect the other.
        """

        # a job description only has attributes - so create a new instance,
        # clone the attribs, and done.
        if not other :
            other = saga.job.Description ()

        return self._attributes_deep_copy (other)


# ------------------------------------------------------------------------------

