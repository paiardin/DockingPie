
import os
import sys
import time
import pprint
import signal
import thread
import random
import threading
import traceback

from .misc import gettid

# --------------------------------------------------------------------
#
def get_trace():

    trace = sys.exc_info ()[2]

    if  trace :
        stack           = traceback.extract_tb  (trace)
        traceback_list  = traceback.format_list (stack)
        return "".join (traceback_list)

    else :
        stack           = traceback.extract_stack ()
        traceback_list  = traceback.format_list (stack)
        return "".join (traceback_list[:-1])


# ------------------------------------------------------------------------------
#
class DebugHelper (object) :
    """
    When instantiated, and when "RADICAL_DEBUG" is set in the environment, this
    class will install a signal handler for SIGUSR1.  When that signal is
    received, a stacktrace for all threads is printed to stdout.
    We also check if SIGINFO is available, which is generally bound to CTRL-T.

    Additionally, a call to 'dh.fs_block(info=None)' will create a file system
    based barrier: it will create a unique file in /tmp/ (based on 'name' if
    given), and dump the stack trace and any 'info' into it.  It then waits
    until that file has changed (touched or removed etc), and then returns.  The
    wait is a simple pull based 'os.stat()' (once per sec).
    """

    # --------------------------------------------------------------------------
    #
    def __init__ (self, name=None, info=None):
        """
        name: string to identify fs barriers
        info: static info to dump into fs barriers
        """

        self._name       = name
        self._info       = info

        if not self._name:
            self._name = str(id(self))

        if 'MainThread' not in threading.current_thread().name:
            # python only supports signals in main threads :/
            return

        if 'RADICAL_DEBUG' in os.environ :
            signal.signal(signal.SIGUSR1, print_stacktraces) # signum 30
            signal.signal(signal.SIGQUIT, print_stacktraces) # signum  3

            try:
                assert signal.SIGINFO
                signal.signal(signal.SIGINFO, print_stacktraces) # signum 29
            except AttributeError as e:
                pass


    # --------------------------------------------------------------------------
    #
    def fs_block(self, info=None):
        """
        Dump state, info in barrier file, and wait for it tou be touched or
        read or removed, then continue.  Leave no trace.
        """

        if not 'RADICAL_DEBUG' in os.environ:
            return

        try:
            pid = os.getpid()
            tid = threading.currentThread().ident

            fs_barrier = "/tmp/ru.dh.%s.%s.%s" % (self._name, pid, tid)
            fs_handle  = open(fs_barrier, 'w+')
            fs_handle.seek(0,0)
            fs_handle.write("\nSTACK TRACE:\n%s\n%s\n" % (time.time(), get_trace()))
            fs_handle.write("\nSTATIC INFO:\n%s\n\n" % pprint.pformat(self._info))
            fs_handle.write("\nINFO:\n%s\n\n" % pprint.pformat(info))
            fs_handle.flush()

            new = os.stat(fs_barrier)
            old = new

            while old == new:
                new = os.stat(fs_barrier)
                time.sleep(1)

            fs_handle.close()
            os.unlink(fs_barrier)

        except Exception as e:
            # we don't care (much)...
            print get_trace()
            print e
            pass


# ------------------------------------------------------------------------------
#
def print_stacktraces(signum=None, sigframe=None):
    """
    signum, sigframe exist to satisfy signal handler signature requirements
    """

    this_tid = threading.currentThread().ident

    # if multiple processes (ie. a process group) get the signal, then all
    # traces are mixed together.  Thus we waid 'pid%100' milliseconds, in
    # the hope that this will stagger the prints.
    pid = int(os.getpid())
    time.sleep((pid%100)/1000)

    out  = "===============================================================\n"
    out += "RADICAL Utils -- Debug Helper -- Stacktraces\n"
    out += os.popen('ps -efw --forest | grep " %s " | grep -v grep' % pid).read()
    try:
        info = get_stacktraces()
    except Exception as e:
        out += 'skipping frame (%s)' % e
        info = None

    if info:
        for tid,tname in info:

            if tid == this_tid : marker = '[active]'
            else               : marker = ''
            out += "---------------------------------------------------------------\n"
            out += "Thread: %s %s\n" % (tname, marker)
            out += "  PID : %s \n"   % os.getpid()
            out += "  TID : %s \n"   % tid
            for fname,lineno,method,code in info[tid,tname]:

                if code: code = code.strip()
                else   : code = '<no code>'

              # # [:-1]: .py vs. .pyc :/
              # if not (__file__[:-1] in fname and \
              #         method in ['get_stacktraces', 'print_stacktraces']):
                if method not in ['get_stacktraces', 'print_stacktraces']:
                    out += "  File: %s, line %d, in %s\n" % (fname, lineno, method)
                    out += "        %s\n" % code

    out += "==============================================================="

    sys.stdout.write("%s\n" % out)

    if 'RADICAL_DEBUG' in os.environ:
        with open('/tmp/ru.stacktrace.%s.log' % pid, 'w') as f:
            f.write ("%s\n" % out)


# --------------------------------------------------------------------------
#
def get_stacktraces():

    id2name = {}
    for th in threading.enumerate():
        id2name[th.ident] = th.name

    ret = dict()
    for tid, stack in sys._current_frames().items():
        if tid in id2name:
            ret[tid,id2name[tid]] = traceback.extract_stack(stack)
        else:
            ret[tid,'noname'] = traceback.extract_stack(stack)

    return ret


# --------------------------------------------------------------------------
#
def print_stacktrace(msg=None):

    if not msg:
        msg = ''

    pid  = int(os.getpid())
    out  = "--------------\n"
    out += "RADICAL Utils -- Stacktrace [%s] [%s]\n" % (pid, threading.currentThread().name)
    out += "%s\n" % msg
    out += os.popen('ps -efw --forest | grep " %s " | grep -v grep' % pid).read()
    for line in get_stacktrace():
        out += line.strip()
        out += "\n"
    out += "--------------\n"

    sys.stdout.write(out)


# --------------------------------------------------------------------------
#
def print_exception_trace(msg=None):

    if not msg:
        msg = ''

    pid  = int(os.getpid())
    out  = "--------------\n"
    out += "RADICAL Utils -- Stacktrace [%s] [%s]\n" % (pid, threading.currentThread().name)
    out += "%s\n" % msg
    out += os.popen('ps -efw --forest | grep " %s " | grep -v grep' % pid).read()
    for line in traceback.format_exc().split('\n'):
        out += line.strip()
        out += "\n"
    out += "--------------\n"

    sys.stdout.write(out)


# --------------------------------------------------------------------------
#
def get_stacktrace():

    return traceback.format_stack()[:-1]


# ------------------------------------------------------------------------------
#
def get_caller_name(skip=2):
    """
    Get a name of a caller in the format module.class.method

    `skip` specifies how many levels of stack to skip while getting caller
    name. skip=1 means "who calls me", skip=2 "who calls my caller" etc.

    An empty string is returned if skipped levels exceed stack height

    Kudos: http://stackoverflow.com/questions/2654113/ \
            python-how-to-get-the-callers-method-name-in-the-called-method
    """
    import inspect

    stack = inspect.stack()
    start = 0 + skip
    if len(stack) < start + 1:
        return ''

    parentframe = stack[start][0]

    name   = []
    module = inspect.getmodule(parentframe)

    # `modname` can be None when frame is executed directly in console
    # TODO(techtonik): consider using __main__
    if module:
        name.append(module.__name__)

    # detect classname
    if 'self' in parentframe.f_locals:
        name.append(parentframe.f_locals['self'].__class__.__name__)

    codename = parentframe.f_code.co_name

    if codename != '<module>':  # top level usually
        name.append(codename)   # function or a method

    del parentframe

    return ".".join(name)


# ------------------------------------------------------------------------------
#
_raise_on_state = dict()
_raise_on_lock  = threading.Lock()
def raise_on(tag, log=None):
    """
    The purpose of this method is to artificially trigger error conditions for
    testing purposes, for example when handling the n'th unit, getting the n'th
    heartbeat signal, etc.  

    The tag parameter is interpreted as follows: on the `n`'th invocation of
    this method with any given `tag`, an exception is raised, and the counter
    for that tag is reset.
    
    The limit `n` is set via an environment variable `RU_RAISE_ON_<tag>`, with
    `tag` in upper casing.  The environment will only be inspected during the
    first invocation of the method with any given tag.  The tag counter is
    process-local, but is shared amongst threads of that process.
    """

    global _raise_on_state
    global _raise_on_lock

    with _raise_on_lock:

        if tag not in _raise_on_state:
            env = os.environ.get('RU_RAISE_ON_%s' % tag.upper())
            if env and env.startswith('RANDOM_'):
                # env is rnd spec
                rate  = (float(env[7:]) / 100.0)
                limit = 1
            elif env:
                # env is int
                rate  = 1
                limit = int(env)
            else:
                # no env set
                rate  = 1
                limit = 0

            _raise_on_state[tag] = { 
                    'count' : 0,
                    'rate'  : rate,
                    'limit' : limit
                    }
            
        _raise_on_state[tag]['count'] += 1

        count = _raise_on_state[tag]['count']
        limit = _raise_on_state[tag]['limit']
        rate  = _raise_on_state[tag]['rate']

        if log: log.debug('raise_on checked   %s [%s / %s]' % (tag, count, limit))
      # else:   print    ('raise_on checked   %s [%s / %s]' % (tag, count, limit))

        if limit and count == limit:

            if rate == 1:
                val = limit
            else:
                val = random.random()
                if val > rate:
                    if log: log.warn('raise_on untriggered %s [%s / %s]' % (tag, val, rate))
                  # else:   print   ('raise_on untriggered %s [%s / %s]' % (tag, val, rate))
                    return

            if log: log.warn('raise_on triggered %s [%s]' % (tag, val))
          # else:   print    'raise_on triggered %s [%s]' % (tag, val)

            # reset counter and raise exception
            _raise_on_state[tag]['count'] = 0
            raise RuntimeError('raise_on for %s [%s]' % (tag, val))


# ------------------------------------------------------------------------------
#
def attach_pudb(logger=None):

    host = '127.0.0.1'
  # host = gethostip()
    tid  = gettid()
    port = tid + 10000

    if logger:
        logger.info('debugger open: telnet %s %d', host, port)
    else:
        print 'debugger open: telnet %s %d' % (host, port)

    import pudb
    import signal
    pudb.DEFAULT_SIGNAL = signal.SIGALRM

    from pudb.remote import set_trace
    set_trace(host=host, port=port, term_size=(200, 50))


# ------------------------------------------------------------------------------

