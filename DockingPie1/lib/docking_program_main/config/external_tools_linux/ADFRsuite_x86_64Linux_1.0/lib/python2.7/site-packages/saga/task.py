
__author__    = "Andre Merzky, Ole Weidner"
__copyright__ = "Copyright 2012-2013, The SAGA Project"
__license__   = "MIT"


""" Task interface
"""

import inspect
import Queue

import radical.utils.signatures as rus
import radical.utils            as ru

import base                     as sbase
import exceptions               as se
import attributes               as satt
import adaptors.cpi.base        as sacb

from   saga.constants       import *


# ------------------------------------------------------------------------------
#
class Task (sbase.SimpleBase, satt.Attributes) :

    # --------------------------------------------------------------------------
    #
    @rus.takes   ('Task', 
                  sacb.CPIBase, 
                  basestring,
                  dict, 
                  rus.one_of (SYNC, ASYNC, TASK))
    @rus.returns (rus.nothing)
    def __init__ (self, _adaptor, _method_type, _method_context, _ttype) :
        """ 
        This saga.Task constructor is private.

        ``_adaptor`` references the adaptor class instance from which this
        task was created via an asynchronous function.  Note that the API level
        object instance can be inferred via ``_adaptor.get_api ()``.  Further, the
        adaptor will reference an _adaptor._container class, which will be
        considered the target for bulk operations for this task.

        ``_method_type`` specifies the SAGA API method which task is
        representing.  For example, for the following code::

          d = saga.filesystem.Directory ("file:///")
          t = d.copy ('/etc/passwd', '/tmp/passwd.bak', saga.task.ASYNC)

        The resulting task ``t`` would represent the *'copy'* method.  This is
        required to forward :class:`saga.task.Container` calls to the correct
        bulk method, in this case ``container_copy()``.

        ``_method_context`` describes the context in which the task method is
        running.  It is up to the creator of the task to provide that context --
        in general, it will at least include method parameters.

        ``ttype`` determines in what state the constructor will leave the task:
        ``DONE`` for ``ttype=SYNC``, ``RUNNING`` for ``ttype=ASYNC`` and ``NEW``
        for ``ttype=TASK``.

        If the ``_method_context`` has *exactly* two elements, names ``_call``
        and ``args``, then the created task will wrap
        a :class:`ru.Thread` with that ``_call (_args)``.
        """
        
        self._base = super  (Task, self)
        self._base.__init__ ()

        self._thread         = None
        self._ttype          = _ttype
        self._adaptor        = _adaptor
        self._method_type    = _method_type
        self._method_context = _method_context

        # set attribute interface properties
        self._attributes_extensible    (False)
        self._attributes_allow_private (True)
        self._attributes_camelcasing   (True)

        # register properties with the attribute interface
        self._attributes_register   (RESULT,    None,    satt.ANY,  satt.SCALAR, satt.READONLY)
        self._attributes_set_getter (RESULT,    self.get_result)
        self._attributes_set_setter (RESULT,    self._set_result)

        self._attributes_register   (EXCEPTION, None,    satt.ANY,  satt.SCALAR, satt.READONLY)
        self._attributes_set_getter (EXCEPTION, self.get_exception)
        self._attributes_set_setter (EXCEPTION, self._set_exception)

        self._attributes_register   (STATE,     UNKNOWN, satt.ENUM, satt.SCALAR, satt.READONLY)
        self._attributes_set_enums  (STATE,    [UNKNOWN, NEW, RUNNING, DONE, FAILED, CANCELED])
        self._attributes_set_getter (STATE,     self.get_state)
        self._attributes_set_setter (STATE,     self._set_state)
              
        self._set_state (NEW)

        # check if this task is supposed to wrap a callable in a thread
        if  '_call'   in self._method_context :

            if not '_args'   in self._method_context :
                self._method_context['_args'] = ()

            if not '_kwargs' in self._method_context :
                self._method_context['_kwargs'] = {}

            if  3 !=  len (self._method_context) :
                raise se.BadParameter ("invalid call context for callable task")
            
            call   = self._method_context['_call']
            args   = self._method_context['_args']
            kwargs = self._method_context['_kwargs']

            # if the called function expects a task handle, provide it.
            if  '_from_task' in inspect.getargspec (call)[0] :
                if  not    '_from_task' in kwargs :
                    kwargs['_from_task'] = self

            self._thread = ru.Thread (call, *args, **kwargs)


        # ensure task goes into the correct state
        if self._ttype == SYNC :
            self.run  ()
            self.wait ()
        elif self._ttype == ASYNC :
            self.run  ()
        elif self._ttype == TASK :
            pass



    # --------------------------------------------------------------------------
    #
    @rus.takes   ('Task')
    @rus.returns (rus.nothing)
    def run      (self) :

        if  self._thread :
            self._thread.run ()

        else :
            # FIXME: make sure task_run exists.  Should be part of the CPI!
            self._adaptor.task_run (self)


    # --------------------------------------------------------------------------
    #
    @rus.takes   ('Task', 
                  rus.optional (float))
    @rus.returns (bool)
    def wait (self, timeout=None) :

        if  None == timeout :
            timeout = -1.0 # FIXME

        if self._thread :
            self._thread.wait ()  # FIXME: timeout?!
            self._set_state   (self._thread.state)

        else :
            # FIXME: make sure task_wait exists.  Should be part of the CPI!
            self._adaptor.task_wait (self, timeout)


    # ----------------------------------------------------------------
    #
    @rus.takes   ('Task', 
                  float)
    @rus.returns (rus.nothing)
    def cancel (self) :

        if self._thread :
            self._thread.cancel ()
            self._set_state (CANCELED)

        else :
            # FIXME: make sure task_cancel exists.  Should be part of the CPI!
            self._adaptor.task_cancel (self)


    # --------------------------------------------------------------------------
    #
    @rus.takes   ('Task', 
                  rus.one_of (UNKNOWN, NEW, RUNNING, DONE, FAILED, CANCELED))
    @rus.returns (rus.nothing)
    def _set_state (self, state) :

        if not state in [UNKNOWN, NEW, RUNNING, DONE, FAILED, CANCELED] :
            raise se.BadParameter ("attempt to set invalid task state '%s'" % state)

        self._attributes_i_set (self._attributes_t_underscore (STATE), state, force=True)


    # --------------------------------------------------------------------------
    #
    @rus.takes   ('Task')
    @rus.returns (rus.one_of (UNKNOWN, NEW, RUNNING, DONE, FAILED, CANCELED))
    def get_state (self) :

        if self._thread :
            self._set_state (self._thread.state)

        return self.state


    # --------------------------------------------------------------------------
    #
    @rus.takes   ('Task', 
                  rus.anything)
    @rus.returns (rus.nothing)
    def _set_result (self, result) :

        self._attributes_i_set (self._attributes_t_underscore (RESULT), result, force=True)
        self._attributes_i_set (self._attributes_t_underscore (STATE),  DONE,   force=True)


    # --------------------------------------------------------------------------
    #
    @rus.takes   ('Task')
    @rus.returns (rus.anything)
    def get_result (self) :
        
        if not self.state in [DONE, FAILED, CANCELED] :
            self.wait ()

        assert (self.state in [DONE, FAILED, CANCELED]) 
        
        if  self.state == FAILED :
            self.re_raise ()
            return

        if self.state == CANCELED :
            raise se.IncorrectState ("task.get_result() cannot be called on cancelled tasks")

        if  self.state == DONE :

            if  self._thread :
                self._set_result (self._thread.result)

            return self.result


    # --------------------------------------------------------------------------
    #
    @rus.takes   ('Task', 
                  basestring, 
                  rus.anything)
    @rus.returns (rus.nothing)
    def _set_metric (self, metric, value) :

        self._attributes_i_set (self._attributes_t_underscore (metric), value, force=True)


    # --------------------------------------------------------------------------
    #
    @rus.takes   ('Task', 
                  se.SagaException)
    @rus.returns (rus.nothing)
    def _set_exception (self, e) :
        self._attributes_i_set (self._attributes_t_underscore (EXCEPTION), e, force=True)

    # --------------------------------------------------------------------------
    #
    @rus.takes   ('Task')
    @rus.returns (se.SagaException)
    def get_exception (self) :

        if  self._thread :
            self._set_exception (self._thread.exception)

        return self.exception

    # --------------------------------------------------------------------------
    #
    @rus.takes   ('Task')
    @rus.returns (rus.nothing)
    def re_raise (self) :

        if self.exception :
            raise self.exception




# ------------------------------------------------------------------------------
#
class Container (sbase.SimpleBase, satt.Attributes) :


    # --------------------------------------------------------------------------
    #
    @rus.takes   ('Container')
    @rus.returns (rus.nothing)
    def __init__ (self) :


        self._base = super  (Container, self)
        self._base.__init__ ()

        # set attribute interface properties
        self._attributes_allow_private (True)
        self._attributes_extensible    (False)
        self._attributes_camelcasing   (True)

        # register properties with the attribute interface
        self._attributes_register   (SIZE,    0,  satt.INT,  satt.SCALAR, satt.READONLY)
        self._attributes_set_getter (SIZE,    self.get_size)

        self._attributes_register   (TASKS,   [], satt.ANY,  satt.VECTOR, satt.READONLY)
        self._attributes_set_getter (TASKS,   self.get_tasks)

        self._attributes_register   (STATES,  [], satt.ENUM, satt.VECTOR, satt.READONLY)
        self._attributes_set_getter (STATES,  self.get_states)

        self._attributes_set_enums  (STATES,  [UNKNOWN, NEW, RUNNING, DONE, FAILED, CANCELED])

        # cache for created container instances
        self._containers = {}


    # --------------------------------------------------------------------------
    #
    def __str__ (self) :

        ret  = "["
        for task in self.tasks :
            ret += "'%s', "  %  str(task)
        ret += "]"

        return ret


    # --------------------------------------------------------------------------
    #
    @rus.takes   ('Container', 
                  Task)
    @rus.returns (rus.nothing)
    def add      (self, task) :

        import saga.job as sjob

        if  not isinstance (task, Task) :
            
            raise se.BadParameter ("Container handles tasks, not %s" \
                                % (type(task)))

        if not task in self.tasks :
            self.tasks.append (task)



    # --------------------------------------------------------------------------
    #
    @rus.takes   ('Container', 
                  Task)
    @rus.returns (rus.nothing)
    def remove   (self, task) :

        if task in self.tasks :
            self.tasks.delete (task)


    # --------------------------------------------------------------------------
    #
    @rus.takes   ('Container')
    @rus.returns (rus.nothing)
    def run      (self) :

        if not len (self.tasks) :
            # nothing to do
            return None

        buckets = self._get_buckets ()
        threads = []  # threads running container ops

        # handle all container
        for c in buckets['bound'] :
        
            # handle all methods
            for m in buckets['bound'][c] :

                tasks    = buckets['bound'][c][m]
                m_name   = "container_%s" % m
                m_handle = None

                for (name, handle) in inspect.getmembers (c, predicate=inspect.ismethod) :
                    if name == m_name :
                        m_handle = handle
                        break

                if not handle :
                    # Hmm, the specified container can't handle the call after
                    # all -- fall back to the unbound handling
                    buckets['unbound'] += tasks

                else :
                    # hand off to the container function, in a separate task
                    threads.append (ru.Thread.Run (m_handle, tasks))


        # handle tasks not bound to a container
        for task in buckets['unbound'] :

            threads.append (ru.Thread.Run (task.run))
            

        # wait for all threads to finish
        for thread in threads :
            if  thread.isAlive () :
                thread.join ()

            if  thread.get_state () == FAILED :
                raise se.NoSuccess ("thread exception: %s" \
                                 % (thread.get_exception ()))


    # --------------------------------------------------------------------------
    #
    # --------------------------------------------------------------------------
    #
    @rus.takes   ('Container', 
                  rus.one_of   (ANY, ALL),
                  rus.optional (float))
    @rus.returns (rus.list_of (Task))
    def wait (self, mode=ALL, timeout=None) :

        if  None == timeout :
            timeout = -1.0 # FIXME

        if not mode in [ANY, ALL] :
            raise se.BadParameter ("wait mode must be saga.task.ANY or saga.task.ALL")

        if type (timeout) not in [int, long, float] : 
            raise se.BadParameter ("wait timeout must be a floating point number (or integer)")

        if not len (self.tasks) :
            # nothing to do
            return None

        if mode == ALL :
            return self._wait_all (timeout)
        else : 
            return self._wait_any (timeout)



    # --------------------------------------------------------------------------
    #
    @rus.takes   ('Container', 
                  float)
    @rus.returns (rus.list_of (Task))
    def _wait_any (self, timeout) :

        buckets = self._get_buckets ()
        threads = []  # threads running container ops

        # handle all tasks bound to containers
        for c in buckets['bound'] :

            # handle all methods -- all go to the same 'container_wait' though)
            tasks = []
            for m in buckets['bound'][c] :
                tasks += buckets['bound'][c][m]

            threads.append (ru.Thread.Run (c.container_wait, tasks, ANY, timeout))

        
        # handle all tasks not bound to containers
        for task in buckets['unbound'] :

            threads.append (ru.Thread.Run (task.wait, timeout))
            

        # mode == ANY: we need to watch our threads, and whenever one
        # returns, and declare success.  Note that we still need to get the
        # finished task from the 'winner'-thread -- we do that via a Queue
        # object.  Note also that looser threads are not canceled, but left
        # running (FIXME: consider sending a signal at least)

        timeout = 0.01 # seconds, heuristic :-/

        for thread in threads :
            thread.join (timeout)

            if thread.get_state () == FAILED :
                raise thread.get_exception ()

            if not thread.isAlive :
                # thread indeed finished -- dig return value from this
                # threads queue
                result = thread.get_result ()

                # ignore other threads, and simply declare success
                return result



    # --------------------------------------------------------------------------
    #
    @rus.takes   ('Container', 
                  float)
    @rus.returns (rus.list_of (Task))
    def _wait_all (self, timeout) :
        # this method should actually be symmetric to _wait_any, and could
        # almost be mapped to it, but the code below is a kind of optimization
        # (does not need threads, thus simpler code).

        buckets = self._get_buckets ()
        ret     = None

        # handle all tasks bound to containers
        for c in buckets['bound'] :

            # handle all methods -- all go to the same 'container_wait' though)
            tasks = []
            for m in buckets['bound'][c] :
                tasks += buckets['bound'][c][m]

            # TODO: this is semantically not correct: timeout is applied
            #       n times...
            c.container_wait (tasks, ALL, timeout)
            ret = tasks[0]
 
        # handle all tasks not bound to containers
        for task in buckets['unbound'] :
            task.wait ()
            ret = task

        # all done - return random task (first from last container, or last
        # unbound task)
        # FIXME: that task should be removed from the task container
        return ret


    # --------------------------------------------------------------------------
    #
    @rus.takes   ('Container', 
                  rus.optional (float))
    @rus.returns (rus.nothing)
    def cancel   (self, timeout=None) :

        if  None == timeout :
            timeout = -1.0 # FIXME

        buckets = self._get_buckets ()
        threads = []  # threads running container ops

        # handle all tasks bound to containers
        for c in buckets['bound'] :

            # handle all methods -- all go to the same 'container_cancel' though)
            tasks = []
            for m in buckets['bound'][c] :
                tasks += buckets['bound'][c][m]

            threads.append (ru.Thread.Run (c.container_cancel, tasks, timeout))

        
        # handle all tasks not bound to containers
        for task in buckets['unbound'] :

            threads.append (ru.Thread.Run (task.cancel, timeout))
            

        for thread in threads :
            thread.join ()


    # ----------------------------------------------------------------
    #
    @rus.takes   ('Container')
    @rus.returns (int)
    def get_size (self) :

        return len (self.tasks)


    # --------------------------------------------------------------------------
    #
    @rus.takes   ('Container', 
                  'basestring')
    @rus.returns (Task)
    def get_task (self, id) :

        # FIXME: this should not be a search, but a lookup
        if not id:
            raise se.NoSuccess ("Lookup requires non-empty id (not '%s')" % id)

        for t in self.tasks:
            if t.id == id:
                return t

        raise se.NoSuccess ("task '%s' not found in container" % id)


    # --------------------------------------------------------------------------
    #
    @rus.takes   ('Container')
    @rus.returns (rus.list_of (Task))
    def get_tasks (self) :

        return self.tasks


    # --------------------------------------------------------------------------
    #
    @rus.takes   ('Container')
    @rus.returns (rus.list_of (rus.one_of (UNKNOWN, NEW, RUNNING, DONE, FAILED, CANCELED)))
    def get_states (self) :

        buckets = self._get_buckets ()
        threads = []  # threads running container ops

        # handle all tasks bound to containers
        for c in buckets['bound'] :

            # handle all methods -- all go to the same 'container_get_states' though)
            tasks = []
            for m in buckets['bound'][c] :
                tasks += buckets['bound'][c][m]

            threads.append (ru.Thread.Run (c.container_get_states, tasks))

        
        # handle all tasks not bound to containers
        for task in buckets['unbound'] :

            threads.append (ru.Thread.Run (task.get_state))
            

        # We still need to get the states from all threads.
        # FIXME: order
        states  = []

        for thread in threads :
            thread.join ()

            if thread.get_state () == FAILED :
                raise thread.get_exception ()

            # FIXME: what about ordering tasks / states?
            res = thread.get_result ()

            if res != None :
                states += res

        return states


    # ----------------------------------------------------------------
    #
    @rus.takes   ('Container')
    @rus.returns (dict)
    def _get_buckets (self) :
        # collective container ops: walk through the task list, and sort into
        # buckets of tasks which have (a) the same task._container, or if that
        # is not set, the same class type (for which one container instance is
        # created).  All tasks were neither is available are handled one-by-one

        buckets = {}
        buckets['unbound'] = [] # no container adaptor for these [tasks]
        buckets['bound']   = {} # dict  of container adaptors [tasks]

        for task in self.tasks :

            if  task._adaptor and task._adaptor._container :

                # the task's adaptor has a valid associated container class 
                # which can handle the container ops - great!
                c = task._adaptor._container
                m = task._method_type

                if not c in buckets['bound'] :
                    buckets['bound'][c] = {}

                if not m in buckets['bound'][c] :
                    buckets['bound'][c][m] = []

                buckets['bound'][c][m].append (task)

            else :

                # we have no container to handle this task -- so
                # put it into the fallback list
                buckets['unbound'].append (task)

        return buckets


# FIXME: add get_apiobject



