
__author__    = "Andre Merzky, Ole Weidner"
__copyright__ = "Copyright 2012-2013, The SAGA Project"
__license__   = "MIT"


""" SAGA job interface """

import radical.utils.signatures as rus

from   saga.constants        import SYNC, ASYNC, TASK
from   saga.job.constants    import *

import saga.adaptors.base    as sab
import saga.attributes       as sa
import saga.exceptions       as se
import saga.async            as sasync
import saga.task             as st
import saga.base             as sb
import saga.url              as surl

import description           as descr

# ------------------------------------------------------------------------------
#
class Job (sb.Base, st.Task, sasync.Async) :
    '''Represents a SAGA job as defined in GFD.90
    
    A 'Job' represents a running application instance, which may consist of one
    or more processes.  Jobs are created by submitting a Job description to
    a Job submission system -- usually a queuing system, or some other service
    which spawns jobs on the user's behalf.

    Jobs have a unique ID (see get_job_id()), and are stateful entities -- their
    'state' attribute changes according to a well defined state model:

    A job as returned by job.Service.create(jd) is in 'New' state -- it is not
    yet submitted to the job submission backend.  Once it was submitted, via
    run(), it will enter the 'Pending' state, where it waits to get actually
    executed by the backend (e.g. waiting in a queue etc).  Once the job is
    actually executed, it enters the 'Running' state -- only in that state is
    the job actually consuming resources (CPU, memory, ...).

    Jobs can leave the 'Running' state in three different ways: they finish
    successfully on their own ('Done'), they finish unsuccessfully on their own,
    or get canceled by the job management backend ('Failed'), or they get
    actively canceled by the user or the application ('Canceled').

    The methods defined on the Job object serve two purposes: inspecting the
    job's state, and initiating job state transitions.

    '''
    
    # --------------------------------------------------------------------------
    #
    @rus.takes   ('Job',
                  rus.optional (basestring),
                  rus.optional (sab.Base),
                  rus.optional (dict),
                  rus.optional (rus.one_of (SYNC, ASYNC, TASK)))
    @rus.returns (rus.nothing)
    def __init__ (self, _method_type='run', _adaptor=None, _adaptor_state={}, _ttype=None) : 
        '''
        _adaptor`` references the adaptor class instance which created this task
        instance.

        The ``_method_type`` parameter is flattened into the job constructor to
        satisfy the bulk optimization properties of the saga.Task class, whose
        interface is implemented by saga.job.Job.
        ``_method_type`` specifies the SAGA API method which task is
        representing.  For jobs, that is the 'run' method.

        We don't have a create classmethod -- jobs are never constructed by the user
        '''

        if not _adaptor :
            raise se.IncorrectState ("saga.job.Job constructor is private")

        self._valid = False


        # we need to keep _method_type around, for the task interface (see
        # :class:`saga.Task`)
        self._method_type = _method_type

        # We need to specify a schema for adaptor selection -- and
        # simply choose the first one the adaptor offers.
        schema = _adaptor.get_schemas()[0]
        if  'job_schema' in _adaptor_state :
            schema = _adaptor_state['job_schema']

        self._base = super  (Job, self)
        self._base.__init__ (schema, _adaptor, _adaptor_state, ttype=None)


        # set attribute interface properties
        self._attributes_allow_private (True)
        self._attributes_extensible    (False)
        self._attributes_camelcasing   (True)

        # register properties with the attribute interface 
        self._attributes_register   (STATE,            UNKNOWN, sa.ENUM,   sa.SCALAR, sa.READONLY)
        self._attributes_register   (EXIT_CODE,        None,    sa.INT,    sa.SCALAR, sa.READONLY)
        self._attributes_register   (CREATED,          None,    sa.INT,    sa.SCALAR, sa.READONLY)
        self._attributes_register   (STARTED,          None,    sa.INT,    sa.SCALAR, sa.READONLY)
        self._attributes_register   (FINISHED,         None,    sa.INT,    sa.SCALAR, sa.READONLY)
        self._attributes_register   (EXECUTION_HOSTS,  None,    sa.STRING, sa.VECTOR, sa.READONLY)
        self._attributes_register   (ID,               None,    sa.STRING, sa.SCALAR, sa.READONLY)
        self._attributes_register   (NAME,             None,    sa.STRING, sa.SCALAR, sa.READONLY)
        self._attributes_register   (SERVICE_URL,      None,    sa.URL,    sa.SCALAR, sa.READONLY)

        self._attributes_set_enums  (STATE, [UNKNOWN, NEW,     PENDING,  RUNNING,
                                             DONE,    FAILED,  CANCELED, SUSPENDED])

        self._attributes_set_getter (STATE,           self.get_state)
        self._attributes_set_getter (ID,              self.get_id)
        self._attributes_set_getter (NAME,            self.get_name)
        self._attributes_set_getter (EXIT_CODE,       self._get_exit_code)
        self._attributes_set_getter (CREATED,         self._get_created)
        self._attributes_set_getter (STARTED,         self._get_started)
        self._attributes_set_getter (FINISHED,        self._get_finished)
        self._attributes_set_getter (EXECUTION_HOSTS, self._get_execution_hosts)
        self._attributes_set_getter (SERVICE_URL    , self._get_service_url)

        self._valid = True
 

    # --------------------------------------------------------------------------
    #
    @rus.takes   ('Job')
    @rus.returns (basestring)
    def __str__  (self) :
        """
        __str__()

        String representation. Returns the job's ID.
        """

        if  not self._valid :
            return 'no job id'

        return str (self.id)


    # --------------------------------------------------------------------------
    #
    @rus.takes   ('Job',
                  rus.optional (rus.one_of (SYNC, ASYNC, TASK)))
    @rus.returns ((rus.nothing, basestring, st.Task))
    def get_id   (self, ttype=None) :
        """
        get_id()

        Return the job ID. 
        """
        id = self._adaptor.get_id (ttype=ttype)
        return id


    # --------------------------------------------------------------------------
    #
    @rus.takes   ('Job',
                  rus.optional (rus.one_of (SYNC, ASYNC, TASK)))
    @rus.returns ((rus.nothing, basestring, st.Task))
    def get_name (self, ttype=None) :
        """
        get_name()

        Return the job name. 
        """
        name = self._adaptor.get_name(ttype=ttype)
        return name

    # --------------------------------------------------------------------------
    #
    @rus.takes          ('Job',
                         rus.optional (rus.one_of (SYNC, ASYNC, TASK)))
    @rus.returns        ((basestring, st.Task))
    def get_description (self, ttype=None) :
        """
        get_description()
          
        Return the job description this job was created from.
        
        The returned description can be used to inspect job properties
        (executable name, arguments, etc.).  It can also be used to start
        identical job instances.

        The returned job description will in general reflect the actual state of
        the running job, and is not necessarily a simple copy of the job
        description which was used to create the job instance.  For example, the
        environment variables in the returned job description may reflect the
        actual environment of the running job instance.


        **Example**::


          service = saga.job.Service("fork://localhost")
          jd = saga.job.Description ()
          jd.executable = '/bin/date'

          j1 = service.create_job(jd)
          j1.run()

          j2 = service.create_job(j1.get_description())
          j2.run()

          service.close()
        """
        return self._adaptor.get_description (ttype=ttype)


    # --------------------------------------------------------------------------
    #
    @rus.takes    ('Job',
                   rus.optional (rus.one_of (SYNC, ASYNC, TASK)))
    @rus.returns  ((file, st.Task))
    def get_stdin (self, ttype=None) :
        """
        get_stdin()

        ttype:     saga.task.type enum
        ret:       string / saga.Task
    
        Return the job's STDIN as string. 
        """

        # FIXME: we have no means to set a stdin stream
        return self._adaptor.get_stdin (ttype=ttype)


    # --------------------------------------------------------------------------
    #
    @rus.takes     ('Job',
                    rus.optional (rus.one_of (SYNC, ASYNC, TASK)))
    @rus.returns   ((file, st.Task))
    def get_stdout (self, ttype=None) :
        """
        get_stdout()

        ttype:     saga.task.type enum
        ret:       string / saga.Task

        Return the job's STDOUT as string. 
        """
        return self._adaptor.get_stdout (ttype=ttype)


    # --------------------------------------------------------------------------
    #
    @rus.takes     ('Job',
                    rus.optional (rus.one_of (SYNC, ASYNC, TASK)))
    @rus.returns   ((str, st.Task))
    def get_stdout_string (self, ttype=None) :
        """
        get_stdout_string()

        Return the job's STDOUT as string.

        ttype:     saga.task.type enum
        ret:       string / saga.Task

        THIS METHOD IS DEPRECATED AND WILL BE REMOVED IN A FUTURE RELEASE.
        USE job.get_stdout() INSTEAD.
        """
        return self.get_stdout (ttype=ttype)


    # --------------------------------------------------------------------------
    #
    @rus.takes     ('Job',
                    rus.optional (rus.one_of (SYNC, ASYNC, TASK)))
    @rus.returns   ((file, st.Task))
    def get_stderr (self, ttype=None) :
        """
        get_stderr()

        Return the job's STDERR as string.

        ttype:     saga.task.type enum
        ret:       string / saga.Task
        """
        return self._adaptor.get_stderr (ttype=ttype)


    # --------------------------------------------------------------------------
    #
    @rus.takes     ('Job',
                    rus.optional (rus.one_of (SYNC, ASYNC, TASK)))
    @rus.returns   ((str, st.Task))
    def get_stderr_string (self, ttype=None) :
        """
        get_stderr_string()

        Return the job's STDERR.

        ttype:     saga.task.type enum
        ret:       string / saga.Task

        THIS METHOD IS DEPRECATED AND WILL BE REMOVED IN A FUTURE RELEASE.
        USE job.get_stderr() INSTEAD.
        """
        return self.get_stderr(ttype=ttype)


    # --------------------------------------------------------------------------
    #
    @rus.takes     ('Job',
                    rus.optional (rus.one_of (SYNC, ASYNC, TASK)))
    @rus.returns   ((str, st.Task))
    def get_log (self, ttype=None) :
        """
        get_log_string()

        Return the job's log information, ie. backend specific log messages
        which have been collected during the job execution.  Those messages also
        include stdout/stderr from the job's pre- and post-exec.  The returned
        string generally contains one log message per line, but the format of
        the string is ultimately undefined.

        ttype:     saga.task.type enum
        ret:       string / saga.Task
        """
        return self._adaptor.get_log (ttype=ttype)


    # --------------------------------------------------------------------------
    #
    @rus.takes     ('Job',
                    rus.optional (rus.one_of (SYNC, ASYNC, TASK)))
    @rus.returns   ((str, st.Task))
    def get_log_string (self, ttype=None) :
        """
        get_log_string()

        Return the job's log information, ie. backend specific log messages
        which have been collected during the job execution.  Those messages also
        include stdout/stderr from the job's pre- and post-exec.  The returned
        string generally contains one log message per line, but the format of
        the string is ultimately undefined.

        ttype:     saga.task.type enum
        ret:       string / saga.Task

        THIS METHOD IS DEPRECATED AND WILL BE REMOVED IN A FUTURE RELEASE.
        USE job.get_log() INSTEAD.
        """
        return self.get_log (ttype=ttype)


    # --------------------------------------------------------------------------
    #
    @rus.takes   ('Job',
                  rus.optional (rus.one_of (SYNC, ASYNC, TASK)))
    @rus.returns ((rus.nothing, st.Task))
    def suspend  (self, ttype=None) :
        """
        suspend()

        Suspend the job.
        """
        return self._adaptor.suspend (ttype=ttype)


    # --------------------------------------------------------------------------
    #
    @rus.takes   ('Job',
                  rus.optional (rus.one_of (SYNC, ASYNC, TASK)))
    @rus.returns ((rus.nothing, st.Task))
    def resume   (self, ttype=None) :
        """
        resume()

        Resume the job.
        """
        return self._adaptor.resume (ttype=ttype)


    # --------------------------------------------------------------------------
    #
    @rus.takes     ('Job',
                    rus.optional (rus.one_of (SYNC, ASYNC, TASK)))
    @rus.returns   ((rus.nothing, st.Task))
    def checkpoint (self, ttype=None) :
        """
        checkpoint()
    
        Checkpoint the job.
        """
        return self._adaptor.checkpoint (ttype=ttype)


    # --------------------------------------------------------------------------
    #
    @rus.takes   ('Job',
                  descr.Description,
                  rus.optional (rus.one_of (SYNC, ASYNC, TASK)))
    @rus.returns ((rus.nothing, st.Task))
    def migrate  (self, jd, ttype=None) :
        """
        jd:        saga.job.Description  
        ttype:     saga.task.type enum
        ret:       None / saga.Task
        """
        return self._adaptor.migrate (jd, ttype=ttype)


    # --------------------------------------------------------------------------
    #
    @rus.takes   ('Job',
                  int,
                  rus.optional (rus.one_of (SYNC, ASYNC, TASK)))
    @rus.returns ((rus.nothing, st.Task))
    def signal   (self, signum, ttype=None) :
        """
        signal(signum)

        Send a signal to the job.

        :param signum: signal to send
        :type  signum: int
        """
        return self._adaptor.signal (signum, ttype=ttype)


    #-----------------------------------------------------------------
    #
    id          = property (get_id)            # string
    description = property (get_description)   # Description
   #stdin       = property (get_stdin)         # File
    stdout      = property (get_stdout_string) # string
    stderr      = property (get_stderr_string) # string
    log         = property (get_log)           # string


    #-----------------------------------------------------------------
    #
    # task methods flattened into job :-/
    #
    @rus.takes   ('Job',
                  rus.optional (rus.one_of (SYNC, ASYNC, TASK)))
    @rus.returns ((rus.nothing, st.Task))
    def run      (self, ttype=None) :
        """
        run()

        Run (start) the job. 
        
        Request that the job is being executed by the backend.  If the backend
        is accepting this run request, the job will move to the 'Pending' or
        'Running' state -- otherwise this method will raise an error, and the
        job will be moved to 'Failed'.


        **Example**::

            js = saga.job.Service("fork://localhost")
            jd = saga.job.Description ()
            jd.executable = '/bin/date'
            j  = js.create_job(jd)

            if j.get_state() == saga.job.NEW : 
                print "new"
            else : 
                print "oops!"

            j.run()

            if   j.get_state() == saga.job.PENDING :
                print "pending"
            elif j.get_state() == saga.job.RUNNING :
                print "running"
            else :
                print "oops!"
          """

        return self._adaptor.run (ttype=ttype)


    # --------------------------------------------------------------------------
    #
    @rus.takes   ('Job',
                  float,
                  rus.optional (rus.one_of (SYNC, ASYNC, TASK)))
    @rus.returns ((rus.nothing, st.Task))
    def cancel   (self, timeout=None, ttype=None) :
        """
        cancel(timeout)

        Cancel the execution of the job.

        :param timeout: `cancel` will return after timeout
        :type  timeout: float

        **Example**::

          js = saga.job.Service("fork://localhost")
          jd = saga.job.Description ()
          jd.executable = '/bin/date'
          j  = js.create_job(jd)

          if   j.get_state() == saga.job.NEW :
              print "new"
          else :
              print "oops!"

          j.run()

          if   j.get_state() == saga.job.PENDING :
              print "pending"
          elif j.get_state() == saga.job.RUNNING :
              print "running"
          else :
              print "oops!"

          j.cancel()

          if   j.get_state() == saga.job.CANCELED :
              print "canceled"
          else :
              print "oops!"
        """
        return self._adaptor.cancel (timeout, ttype=ttype)


    # --------------------------------------------------------------------------
    #
    @rus.takes   ('Job',
                  rus.optional (float),
                  rus.optional (rus.one_of (SYNC, ASYNC, TASK)))
    @rus.returns ((bool, st.Task))
    def wait     (self, timeout=None, ttype=None) :
        """
        wait(timeout)

        :param timeout: `wait` will return after timeout
        :type  timeout: float
        
        Wait for a running job to finish execution.

        The optional timeout parameter specifies the time to wait, and accepts
        the following values::

          timeout <  0  : wait forever (block) -- same for 'None'
          timeout == 0  : wait not at all (non-blocking test)
          timeout >  0  : wait for 'timeout' seconds

        On a non-negative timeout, the call can thus return even if the job is
        not in final state, and the application should check the actual job
        state.  The default timeout value is 'None' (blocking).


        **Example**::

          js = saga.job.Service("fork://localhost")
          jd = saga.job.Description ()
          jd.executable = '/bin/date'
          j  = js.create_job(jd)

          if   j.get_state() == saga.job.NEW :
              print "new"
          else :
              print "oops!"

          j.run()

          if   j.get_state() == saga.job.PENDING :
              print "pending"
          elif j.get_state() == saga.job.RUNNING :
              print "running"
          else :
              print "oops!"

          j.wait(-1.0)

          if   j.get_state() == saga.job.DONE :
              print "done"
          elif j.get_state() == saga.job.FAILED :
              print "failed"
          else :
              print "oops!"
        """

        if  None == timeout :
            timeout = -1.0 # FIXME

        return self._adaptor.wait (timeout, ttype=ttype)


    # --------------------------------------------------------------------------
    #
    @rus.takes    ('Job',
                   rus.optional (rus.one_of (SYNC, ASYNC, TASK)))
    @rus.returns  ((rus.one_of (UNKNOWN, NEW, PENDING, RUNNING, SUSPENDED, DONE, FAILED, CANCELED), st.Task))
    def get_state (self, ttype=None) :
        """
        get_state()
        
        Return the current state of the job.
    
        **Example**::
    
          js = saga.job.Service("fork://localhost")
          jd = saga.job.Description ()
          jd.executable = '/bin/date'
          j  = js.create_job(jd)
    
          if   j.get_state() == saga.job.NEW : 
              print "new"
          else : 
              print "oops!"
    
          j.run()
    
          if   j.get_state() == saga.job.PENDING : 
              print "pending"
          elif j.get_state() == saga.job.RUNNING : 
              print "running"
          else :
              print "oops!"
        """
        return self._adaptor.get_state (ttype=ttype)


    # --------------------------------------------------------------------------
    #
    @rus.takes     ('Job',
                    rus.optional (rus.one_of (SYNC, ASYNC, TASK)))
    @rus.returns   ((rus.anything, st.Task))
    def get_result (self, ttype=None) :
        """
        get_result()
        """
        return self._adaptor.get_result (ttype=ttype)

    # --------------------------------------------------------------------------
    #
    @rus.takes     ('Job',
                    rus.optional (rus.one_of (SYNC, ASYNC, TASK)))
    @rus.returns   ((sb.Base, st.Task))
    def get_object (self, ttype=None) :
        """ :todo: describe me
            :note: this will return the job_service which created the job.
        """
        return self._adaptor.get_object (ttype=ttype)


    # --------------------------------------------------------------------------
    #
    @rus.takes     ('Job',
                    rus.optional (rus.one_of (SYNC, ASYNC, TASK)))
    @rus.returns   ((se.SagaException, st.Task))
    def get_exception (self, ttype=None) :
        """ :todo: describe me

            :note: if job failed, that will get an exception describing 
                   why, if that exists.  Otherwise, the call returns None.
        """
        
        if  self.state == FAILED :
            return se.NoSuccess ("job stderr: %s" % self.get_stderr_string ())
        else :
            return None


    # --------------------------------------------------------------------------
    #
    @rus.takes     ('Job')
    @rus.returns   (rus.nothing)
    def re_raise   (self) :
        """ :todo: describe me

            :note: if job failed, that will re-raise an exception describing 
                   why, if that exists.  Otherwise, the call does nothing.
        """

        if  self.state == FAILED :
            raise se.NoSuccess ("job stderr: %s" % self.get_stderr_string ())
        else :
            return 


    # ----------------------------------------------------------------
    # 
    # attribute getters
    #
    @rus.takes         ('Job',
                        rus.optional (rus.one_of (SYNC, ASYNC, TASK)))
    @rus.returns       ((rus.nothing, int, st.Task))
    def _get_exit_code (self, ttype=None) :
        ec = self._adaptor.get_exit_code(ttype=ttype)
        if ec in [None, ""]:
            return None
        else:
            # Exit code is always an int. If this 'cast' fails, 
            # the adaptor is doing something stupid.
            return ec

    # --------------------------------------------------------------------------
    #
    @rus.takes       ('Job',
                      rus.optional (rus.one_of (SYNC, ASYNC, TASK)))
    @rus.returns     ((rus.nothing, float, st.Task))
    def _get_created (self, ttype=None) :
        return self._adaptor.get_created (ttype=ttype)

    # --------------------------------------------------------------------------
    #
    @rus.takes       ('Job',
                      rus.optional (rus.one_of (SYNC, ASYNC, TASK)))
    @rus.returns     ((rus.nothing, float, st.Task))
    def _get_started (self, ttype=None) :
        return self._adaptor.get_started (ttype=ttype)

    # --------------------------------------------------------------------------
    #
    @rus.takes       ('Job',
                      rus.optional (rus.one_of (SYNC, ASYNC, TASK)))
    @rus.returns     ((rus.nothing, float, st.Task))
    def _get_finished (self, ttype=None) :
        return self._adaptor.get_finished (ttype=ttype)

    # --------------------------------------------------------------------------
    #
    @rus.takes       ('Job',
                      rus.optional (rus.one_of (SYNC, ASYNC, TASK)))
    @rus.returns     ((rus.nothing, rus.list_of (basestring), st.Task))
    def _get_execution_hosts (self, ttype=None) :
        return self._adaptor.get_execution_hosts (ttype=ttype)

    # --------------------------------------------------------------------------
    #
    @rus.takes       ('Job',
                      rus.optional (rus.one_of (SYNC, ASYNC, TASK)))
    @rus.returns     ((rus.nothing, surl.Url, st.Task))
    def _get_service_url (self, ttype=None) :
        return self._adaptor.get_service_url (ttype=ttype)

    name      = property (get_name)        # job name
    state     = property (get_state)       # state enum
    result    = property (get_result)      # result type    (None)
    object    = property (get_object)      # object type    (job_service)
    exception = property (re_raise)        # exception type


    # ----------------------------------------------------------------
    #
    ExitCode = property (doc = """
    The job's exitcode.

    this attribute is only meaningful if the job is in 'Done' or 'Final'
    state - for all other job states, this attribute value is undefined.

    **Example**::


      js = saga.job.Service("fork://localhost")
      jd = saga.job.Description ()
      jd.executable = '/bin/date'
      j  = js.create_job(jd)

      j.run()
      j.wait()

      if j.get_state() == saga.job.FAILED :
        if j.exitcode == "42" :
            print "Ah, galaxy bypass error!"
        else :
            print "oops!"

    """)


    # ----------------------------------------------------------------
    #
    JobID = property (doc = """
    The job's identifier.

    This attribute is equivalent to the value returned by job.get_job_id()
    """)


    # ----------------------------------------------------------------
    #
    ServiceURL = property (doc = """
    The URL of the :class:`saga.job.Service` instance managing this job.

    This attribute is represents the URL under where the job management
    service can be contacted which owns the job.  The value is equivalent to
    the service part of the job_id.

    **Example**::


      js = saga.job.Service("fork://localhost")
      jd = saga.job.Description ()
      jd.executable = '/bin/date'
      j  = js.create_job(jd)

      if j.serviceurl == "fork://localhost" :
          print "yes!"
      else :
          print "oops!"

    """)


# ------------------------------------------------------------------------------
#
# class Self (Job, monitoring.Steerable) :
#
class Self (Job) :

    # --------------------------------------------------------------------------
    #
    @rus.takes   ('Self',
                  rus.optional (basestring),
                  rus.optional (sab.Base),
                  rus.optional (dict),
                  rus.optional (rus.one_of (SYNC, ASYNC, TASK)))
    @rus.returns (rus.nothing)
    def __init__ (self, _method_type='run', _adaptor=None, _adaptor_state={}, _ttype=None) : 

    
        self._base = super  (Job, self)
        self._base.__init__ (_method_type, _adaptor, _adaptor_state, _ttype=_ttype)




