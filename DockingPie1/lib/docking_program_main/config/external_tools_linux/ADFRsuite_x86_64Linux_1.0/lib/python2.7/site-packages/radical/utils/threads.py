
__author__    = "Radical.Utils Development Team (Andre Merzky)"
__copyright__ = "Copyright 2013, RADICAL@Rutgers"
__license__   = "MIT"


import os
import sys
import time
import signal
import thread
import threading
import traceback

import misc  as rumisc


_out_lock = threading.RLock()


# ------------------------------------------------------------------------------
#
NEW     = 'New'
RUNNING = 'Running'
DONE    = 'Done'
FAILED  = 'Failed'


# ------------------------------------------------------------------------------
#
def lout(txt, stream=sys.stdout):

    with _out_lock:
        stream.write(txt)
        stream.flush()


# ------------------------------------------------------------------------------
#
def Event(*args, **kwargs):
    return threading.Event(*args, **kwargs)


# ------------------------------------------------------------------------------
#
class RLock(object):
    """
    This threading.RLock wrapper is supportive of lock debugging.  The only
    semantic difference to threading.RLock is that a lock acquired via the
    'with' statement can be released within the 'with' scope, w/o penalty when
    leaving the locked scope.  This supports up-calling callback semantics, but
    should be used with utter care, and rarely (such as on close()).

    see http://stackoverflow.com/questions/6780613/
         is-it-possible-to-subclass-lock-objects-in-python-if-not-other-ways-to-debug
    """

    # --------------------------------------------------------------------------
    #
    def __init__(self, obj=None):

        self._lock = threading.RLock()

      # with self._lock:
      #     self._obj = obj
      #     self._cnt = 0


    # --------------------------------------------------------------------------
    #
    def acquire(self):

      # ind = (self._cnt)*' '+'>'+(30-self._cnt)*' '
      # lout("%s -- %-10s %50s acquire  - %s\n" % (ind, threading.current_thread().name, self, self._lock))

        self._lock.acquire()

      # self._cnt += 1
      # ind = (self._cnt)*' '+'|'+(30-self._cnt)*' '
      # lout("%s    %-10s %50s acquired - %s\n" % (ind, threading.current_thread().name, self, self._lock))


    # --------------------------------------------------------------------------
    #
    def release(self):

      # ind = (self._cnt)*' '+'-'+(30-self._cnt)*' '
      # lout("%s    %-10s %50s release  - %s\n" % (ind, threading.current_thread().name, self, self._lock))

        try:
            self._lock.release()
        except RuntimeError as e:
            # lock has been released meanwhile - we allow that
          # print 'ignore double lock release'
            pass

      # self._cnt -= 1
      # ind = (self._cnt)*' '+'<'+(30-self._cnt)*' '
      # lout("%s -- %-10s %50s released - %s\n" % (ind, threading.current_thread().name, self, self._lock))


    # --------------------------------------------------------------------------
    #
    def __enter__(self)                        : self.acquire() 
    def __exit__ (self, type, value, traceback): self.release()


# ------------------------------------------------------------------------------
#
class Thread(threading.Thread):
    """
    This `Thread` class is a thin wrapper around Python's native
    `threading.Thread` class, which adds some convenience methods.  
    """

    # --------------------------------------------------------------------------
    #
    def __init__(self, call, *args, **kwargs):

        if not callable(call):
            raise ValueError("Thread requires a callable to function, not %s" \
                            % (str(call)))

        threading.Thread.__init__(self)

        self._call      = call
        self._args      = args
        self._kwargs    = kwargs
        self._state     = NEW
        self._result    = None
        self._exception = None
        self._traceback = None
        self.daemon     = True


    # --------------------------------------------------------------------------
    #
    @classmethod
    def Run(self, call, *args, **kwargs):

        t = self(call, *args, **kwargs)
        t.start()
        return t


    # --------------------------------------------------------------------------
    #
    @property 
    def tid(self):
        return self.tid


    # --------------------------------------------------------------------------
    #
    def run(self):

        try:
            self._state     = RUNNING
            self._result    = self._call(*self._args, **self._kwargs)
            self._state     = DONE

        except Exception as e:
            tb = traceback.format_exc()
            self._traceback = tb
            self._exception = e
            self._state     = FAILED


    # --------------------------------------------------------------------------
    #
    def wait(self):

        if  self.isAlive():
            self.join()


    # --------------------------------------------------------------------------
    #
    def cancel(self):
        # FIXME: this is not really implementable generically, so we ignore 
        # cancel requests for now.
        pass


    # --------------------------------------------------------------------------
    #
    def get_state(self):
        return self._state 

    state = property(get_state)


    # --------------------------------------------------------------------------
    #
    def get_result(self):

        if  self._state == DONE:
            return self._result

        return None

    result = property(get_result)


    # --------------------------------------------------------------------------
    #
    def get_exception(self):

        return self._exception

    exception = property(get_exception)


    # --------------------------------------------------------------------------
    #
    def get_traceback(self):

        return self._traceback

    traceback = property(get_traceback)


# ------------------------------------------------------------------------------
#
def is_main_thread():

    return isinstance(threading.current_thread(), threading._MainThread)


# ------------------------------------------------------------------------------
#
_signal_lock = threading.Lock()
_signal_sent = dict()
def cancel_main_thread(signame=None, once=False):
    """
    This method will call thread.interrupt_main from any calling subthread.
    That will cause a 'KeyboardInterrupt' exception in the main thread.  This
    can be excepted via `except KeyboardInterrupt`
    
    The main thread MUST NOT have a SIGINT signal handler installed (other than
    the default handler or SIGIGN), otherwise this call will cause an exception
    in the core python signal handling (see http://bugs.python.org/issue23395).

    The thread will exit after this, via sys.exit(0), and can then be joined
    from the main thread.

    When being called *from* the main thread, no interrupt will be generated,
    but sys.exit() will still be called.  This can be excepted in the code via 
    `except SystemExit:`.

    Another way to avoid the SIGINT problem is to send a different signal to the
    main thread.  We do so if `signal` is specified.

    After sending the signal, any sub-thread will call sys.exit(), and thus
    finish.  We leave it to the main thread thogh to decide if it will exit at
    this point or not.  Either way, it will have to handle the signal first.

    If `once` is set to `True`, we will send the given signal at most once.
    This will mitigate races between multiple error causes, specifically during
    finalization.
    """

    global _signal_lock
    global _signal_sent


    if signame: signal = get_signal_by_name(signame)
    else      : signal = None


    with _signal_lock:

        if once:
            if signal in _signal_sent:
                # don't signal again
                return

        if signal:
            # send the given signal to the process to which this thread belongs
            os.kill(os.getpid(), signal)
        else:
            # this sends a SIGINT, resulting in a KeyboardInterrupt.
            # NOTE: see http://bugs.python.org/issue23395 for problems on using
            #       SIGINT in combination with signal handlers!
            thread.interrupt_main()

        # record the signal sending
        _signal_sent[signal] = True


    # the sub thread will at this point also exit.
    if not is_main_thread():
        sys.exit()


# ------------------------------------------------------------------------------
#
def get_signal_by_name(signame):
    """
    Translate a signal name into the respective signal number.  If `signame` is
    not found to be a valid signal name, this method will raise a `KeyError`
    exception.  Lookup is case insensitive.
    """

    table = {'abrt'    : signal.SIGABRT,
             'alrm'    : signal.SIGALRM,
             'bus'     : signal.SIGBUS,
             'chld'    : signal.SIGCHLD,
             'cld'     : signal.SIGCLD,
             'cont'    : signal.SIGCONT,
             'fpe'     : signal.SIGFPE,
             'hup'     : signal.SIGHUP,
             'ill'     : signal.SIGILL,
             'int'     : signal.SIGINT,
             'io'      : signal.SIGIO,
             'iot'     : signal.SIGIOT,
             'kill'    : signal.SIGKILL,
             'pipe'    : signal.SIGPIPE,
             'poll'    : signal.SIGPOLL,
             'prof'    : signal.SIGPROF,
             'pwr'     : signal.SIGPWR,
             'quit'    : signal.SIGQUIT,
             'rtmax'   : signal.SIGRTMAX,
             'rtmin'   : signal.SIGRTMIN,
             'segv'    : signal.SIGSEGV,
             'stop'    : signal.SIGSTOP,
             'sys'     : signal.SIGSYS,
             'term'    : signal.SIGTERM,
             'trap'    : signal.SIGTRAP,
             'tstp'    : signal.SIGTSTP,
             'ttin'    : signal.SIGTTIN,
             'ttou'    : signal.SIGTTOU,
             'urg'     : signal.SIGURG,
             'usr1'    : signal.SIGUSR1,
             'usr2'    : signal.SIGUSR2,
             'vtalrm'  : signal.SIGVTALRM,
             'winch'   : signal.SIGWINCH,
             'xcpu'    : signal.SIGXCPU,
             'xfsz'    : signal.SIGXFSZ,
             }
    
    return table[signame.lower()]


# ------------------------------------------------------------------------------
#
class ThreadExit(SystemExit):
    pass

class SignalRaised(SystemExit):

    def __init__(self, msg, signum=None):
        if signum:
            msg = '%s [signal: %s]' % (msg, signum)
        SystemExit.__init__(self, msg)


# ------------------------------------------------------------------------------
#
def get_thread_name():

    return threading.current_thread().name


# ------------------------------------------------------------------------------
#
def get_thread_id():

    return threading.current_thread().ident


# ------------------------------------------------------------------------------
#
def raise_in_thread(e=None, tname=None, tident=None):
    """
    This method uses an internal Python function to inject an exception 'e' 
    into any given thread.  That thread can be specified by its name ('tname')
    or thread id ('tid').  If not specified, the exception is sent to the
    MainThread.

    The target thread will receive the exception with some delay.  More
    specifically, it needs to call up to 100 op codes before the exception 
    is evaluated and raised.

    The default exception raised is 'radical.utils.ThreadExit' which inherits
    from 'SystemExit'.

    NOTE: this is not reliable: the exception is not raised immediately, but is
          *scheduled* for raising at some point, ie. in general after about 100
          opcodes (`sys.getcheckinterval()`).  Depending on when exactly the 
          exception is finally raised, the interpreter might silently swallow
          it, if that happens in a generic try/except clause.  Those exist in
          the Python core, even if discouraged by some PEP or the other.

          See https://bugs.python.org/issue1779233


    NOTE: we can only raise exception *types*, not exception *instances*

          See https://bugs.python.org/issue1538556
    """

    if not tident:
        if not tname:
            tname = 'MainThread'

        for th in threading.enumerate():
            if tname  == th.name:
                tident = th.ident
                break

    if not tident:
        raise ValueError('no target thread given/found')

    if not e:
        e = ThreadExit

    self_thread = threading.current_thread()
    if self_thread.ident == tident:
        # if we are in the target thread, we simply raise the exception.  
        # This specifically also applies to the main thread.
        # Alas, we don't have a decent message to use...
        raise e('raise_in_thread')

    else:
        import ctypes
        ctypes.pythonapi.PyThreadState_SetAsyncExc(ctypes.c_long(tident),
                                                   ctypes.py_object(e))


# ------------------------------------------------------------------------------
#
def fs_event_create(fname, msg=None):
    """
    Atomically create a file at the given `fname` (relative to pwd),
    with `msg` as content (or empty otherwise).

    NOTE: this expects a POSIX compliant file system to be accessible by the two
          entities which use this fs event mechanism.
    NOTE: this also assumes `os.rename()` to be an atomic operation.
    """

    pid = os.getpid()

    if not msg:
        msg = ''

    with open('%s.%d.in' % (fname, pid), 'w') as f:
        f.write(msg)

    os.rename('%s.%d.in' % (fname, pid), fname)


# ------------------------------------------------------------------------------
#
def fs_event_wait(fname, timeout=None):
    """
    Wait for a file ate the given `fname` to appear.  Return `None` at timeout,
    or the file content otherwise.
    """

    msg   = None
    pid   = os.getpid()
    start = time.time()

    while True:

        try:
            with open(fname, 'r') as f:
                msg = f.read()
        except Exception as e:
            print 'e: %s' % type(e)
            print 'e: %s' % e
            pass

        if msg != None:
            try:
                os.rename(fname, '%s.%d.out' % (fname, pid))
                os.unlink('%s.%d.out' % (fname, pid))
            except Exception as e:
                # there is not much we can do at this point...
                print 'unlink: %s' % e
                pass
            return msg

        if timeout and start + timeout <= time.time():
            return None

        time.sleep(0.1)


# ------------------------------------------------------------------------------

