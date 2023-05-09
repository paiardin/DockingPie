
__copyright__ = "Copyright 2013-2014, http://radical.rutgers.edu"
__license__   = "MIT"

# ------------------------------------------------------------------------------
#
"""

Using the RU logging module can lead to deadlocks when used in multiprocess and
multithreaded environments.  This is due to the native python logging library
not being thread save.  See [1] for a Python bug report on the topic from 2009.
As of 2016, a patch is apparently submitted, but not yet accepted, for Python
3.5.  So we have to solve the problem on our own.

The fix is encoded in the 'after_fork()' method.  That method *must* be called
immediately after a process fork.  A fork will carry over all locks from the
parent process -- but it will not carry over any threads.  That means that, if
a fork happens while a thread owns the lock, the child process will start with
a locked lock -- but since there are no threads, there is actually nobody who
can unlock the lock, thus badaboom.

Since after the fork we *know* the logging locks should be unset (after all, we
are not in a logging call right now, are we?), we pre-emtively unlock them here.
But, to do so we need to know what locks exist in the first place.  For that
purpose, we create a process-singletone of all the loggers we hand out via
'get_logger()'.

Logging locks are 'threading.RLock' instances.  As such they can be locked
multiple times (from within the same thread), and we have to unlock them that
many times.  We use a shortcut, and create a new, unlocked lock.
"""

# ------------------------------------------------------------------------------
#
import os
import sys
import time
import threading

from .atfork import *

# monkeypatching can be disabled by setting RADICAL_UTILS_NOATFORK
if not 'RADICAL_UTILS_NOATFORK' in os.environ:
    stdlib_fixer.fix_logging_module()
    monkeypatch_os_fork_functions()

# ------------------------------------------------------------------------------
#
class _LoggerRegistry(object):

    from .singleton import Singleton
    __metaclass__ = Singleton

    def __init__(self):
        self._registry = list()

    def add(self, logger):
        self._registry.append(logger)

    def release_all(self):
        for logger in self._registry:
            while logger:
                for handler in logger.handlers:
                    handler.lock = threading.RLock()
                  # handler.reset()
                logger = logger.parent

# ------------------------------------------------------------------------------
#
_logger_registry = _LoggerRegistry()


# ------------------------------------------------------------------------------
#
def _after_fork():

    _logger_registry.release_all()
    logging._lock = threading.RLock()

# ------------------------------------------------------------------------------
#
def _atfork_prepare():
    pass

# ------------------------------------------------------------------------------
#
def _atfork_parent():
    pass

# ------------------------------------------------------------------------------
#
def _atfork_child():
    _after_fork()

if not 'RADICAL_UTILS_NOATFORK' in os.environ:
    atfork(_atfork_prepare, _atfork_parent, _atfork_child)

#
# ------------------------------------------------------------------------------


from   logging  import DEBUG, INFO, WARNING, WARN, ERROR, CRITICAL
import logging
import colorama

from .misc      import import_module
from .reporter  import Reporter


# The 'REPORT' level is used for demo output and the like.
# log.report.info(msg) style messages always go to stdout, and are enabled if:
#   - log level is set exactly to 'REPORT' (35)
#   - log level is < than 'REPORT' (ie. DEBUG, INFO, REPORT), AND 'stdout' is
#     not in the log target list
#
# To redirect debug logging to some file and have the reporter enabled, use:
#   - export RADICAL_PILOT_VERBOSE=DEBUG
#   - export RADICAL_PILOT_LOG_TGT=rp.log
#
# To interleave both log and reporter on screen, set the log target to 'stderr'.

_DEFAULT_LEVEL = 'ERROR'
REPORT = 35
OFF    = -1


# ------------------------------------------------------------------------------
#
class ColorStreamHandler(logging.StreamHandler):
    """
    A colorized output SteamHandler
    """

    colours = {'DEBUG'    : colorama.Fore.CYAN,
               'INFO'     : colorama.Fore.GREEN,
               'WARN'     : colorama.Fore.YELLOW,
               'WARNING'  : colorama.Fore.YELLOW,
               'ERROR'    : colorama.Fore.RED,
               'CRITICAL' : colorama.Back.RED + colorama.Fore.WHITE
    }


    # --------------------------------------------------------------------------
    #
    def __init__(self, target):

        logging.StreamHandler.__init__(self, target)
        self._tty  = self.stream.isatty()
        self._term = getattr(self, 'terminator', '\n')

    # --------------------------------------------------------------------------
    #
    def emit(self, record):

        # only write in color when using a tty
        if self._tty:
            self.stream.write(self.colours[record.levelname] \
                             + self.format(record)           \
                             + colorama.Style.RESET_ALL      \
                             + self._term)
        else:
            self.stream.write(self.format(record) + self._term)

      # self.flush()


# ------------------------------------------------------------------------------
#
class FSHandler(logging.FileHandler):

    def __init__(self, target):

        try:
            os.makedirs(os.path.abspath(os.path.dirname(target)))
        except:
            pass # exists
        logging.FileHandler.__init__(self, target)


# ------------------------------------------------------------------------------
#
def get_logger(name, target=None, level=None, path=None, header=True):
    """
    Get a logging handle.

    'name'   is used to identify log entries on this handle.
    'target' is a comma separated list (or Python list) of specifiers, where
             specifiers are:
             '-'      : stdout
             '1'      : stdout
             'stdout' : stdout
             '='      : stderr
             '2'      : stderr
             'stderr' : stderr
             '.'      : logfile named ./<name>.log
             <string> : logfile named <string>
    'level'  log level (DEBUG, INFO, WARNING, ERROR, CRITICAL)
    'header' print some version info on log open
    """

    if not name:
        name = 'radical'

    if not path:
        path = os.getcwd()

    if '/' in name:
        try:
            os.makedirs(os.path.normpath(os.path.dirname(name)))
        except:
            # dir exists
            pass

    logger = logging.getLogger(name)
    logger.propagate = False   # let messages not trickle upward
    logger.name = name

    if logger.handlers:
        # we already conifgured that logger in the past -- just reuse it
        return logger

    # --------------------------------------------------------------------------
    # unconfigured loggers get configured.  We try to get the log level and
    # target from the environment.  We try env vars like this:
    #
    #     name  : radical.saga.pty
    #
    #     level : RADICAL_SAGA_PTY_VERBOSE
    #             RADICAL_SAGA_VERBOSE
    #             RADICAL_VERBOSE
    #
    #     target: RADICAL_SAGA_PTY_LOG_TARGET
    #             RADICAL_SAGA_LOG_TARGET
    #             RADICAL_LOG_TARGET
    # whatever is found first is applied.  Well, we actually try the way around
    # and then apply the last one -- but whatever ;)
    # Default level  is 'CRITICAL',
    # default target is '-' (stdout).


    env_name = name.upper().replace('.', '_')
    if '_' in env_name:
        elems = env_name.split('_')
    else:
        elems = [env_name]

    if not level:
        level = _DEFAULT_LEVEL
        for i in range(0,len(elems)):
            env_test = '_'.join(elems[:i+1]) + '_VERBOSE'
            level    = os.environ.get(env_test, level)

    # backward compatible interpretation of SAGA_VERBOSE
    if env_name.startswith('RADICAL_SAGA'):
        level = os.environ.get('SAGA_VERBOSE', level)

    # translate numeric levels into symbolic ones
    level = {50 : 'CRITICAL',
             40 : 'ERROR',
             30 : 'WARNING',
             35 : 'REPORT',
             20 : 'INFO',
             10 : 'DEBUG',
              0 : _DEFAULT_LEVEL,
             -1 : 'OFF'}.get(level, str(level))

    # we want levels to be uppercase
    level = level.upper()

    level_warning = None
    if level not in ['OFF', 'DEBUG', 'INFO', 'WARN', 'WARNING', 'ERROR', 'CRITICAL', 'REPORT']:
        level_warning = "log level '%s' not supported -- reset to '%s'" % (level, _DEFAULT_LEVEL)
        level = _DEFAULT_LEVEL

    if level in ['REPORT']:
        level = REPORT

    if level in [OFF, 'OFF']:
        # we don't want logging actually
        target = 'null'

    if not target:
        target = 'stderr'
        for i in range(0,len(elems)):
            env_test = '_'.join(elems[:i+1]) + '_LOG_TARGET'
            target   = os.environ.get(env_test, target)
            env_test = '_'.join(elems[:i+1]) + '_LOG_TGT'
            target   = os.environ.get(env_test, target)

    if isinstance(target, list):
        targets = target
    else:
        targets = target.split(',')

    formatter = logging.Formatter('%(asctime)s: ' \
                                  '%(name)-20s: ' \
                                  '%(processName)-32s: ' \
                                  '%(threadName)-15s: ' \
                                  '%(levelname)-8s: ' \
                                  '%(message)s')

    # add a handler for each targets (using the same format)
    logger.targets = targets
    for t in logger.targets:
        if t in ['null']:
            continue
        if t in ['-', '1', 'stdout']:
            handle = ColorStreamHandler(sys.stdout)
        elif t in ['=', '2', 'stderr']:
            handle = ColorStreamHandler(sys.stderr)
        elif t in ['.']:
            handle = FSHandler("%s/%s.log" % (path, name))
        elif t.startswith('/'):
            handle = FSHandler(t)
        else:
            handle = FSHandler("%s/%s" % (path, t))
        handle.setFormatter(formatter)
        handle.name = '%s.%s' % (logger.name, str(t))
        logger.addHandler(handle)

    if level in [OFF, 'OFF']:
        logger.setLevel('CRITICAL')
    else:
        logger.setLevel(level)

    if level_warning:
        logger.warn(level_warning)

    # if the given name points to a version or version_detail, log those
    if header:
        try:
            logger.info("%-20s version: %s", 'python.interpreter', ' '.join(sys.version.split()))
            tmp = import_module(name)
            if hasattr(tmp, 'version_detail'):
                logger.info("%-20s version: %s", name, getattr(tmp, 'version_detail'))
            elif hasattr(tmp, 'version'):
                logger.info("%-20s version: %s", name, getattr(tmp, 'version'))
        except:
            pass

        try:
            logger.info("%-20s pid: %s", '', os.getpid())
            logger.info("%-20s tid: %s", '', threading.current_thread().name)
        except:
            pass

    # we also equip our logger with reporting capabilities, so that we can
    # report, for example, demo output whereever we have a logger.
    class _LogReporter(object):

        def __init__(self, logger):
            self._logger   = logger
            self._reporter = Reporter()

            # we always enable report if the log level is REPORT
            # otherwise, we enable report if log level  < REPORT and targets does
            # not contain stdout.
            if logger.getEffectiveLevel() == REPORT:
                self._enabled = True
            else:
                if  logger.getEffectiveLevel() < REPORT and \
                    '-'      not in logger.targets      and \
                    '1'      not in logger.targets      and \
                    'stdout' not in logger.targets      :
                    self._enabled = True
                else:
                    self._enabled = False

        def title(self, *args, **kwargs):
            if self._enabled:
                self._reporter.title(*args, **kwargs)

        def header(self, *args, **kwargs):
            if self._enabled:
                self._reporter.header(*args, **kwargs)

        def info(self, *args, **kwargs):
            if self._enabled:
                self._reporter.info(*args, **kwargs)

        def idle(self, *args, **kwargs):
            if self._enabled:
                self._reporter.idle(*args, **kwargs)

        def progress(self, *args, **kwargs):
            if self._enabled:
                self._reporter.progress(*args, **kwargs)

        def ok(self, *args, **kwargs):
            if self._enabled:
                self._reporter.ok(*args, **kwargs)

        def warn(self, *args, **kwargs):
            if self._enabled:
                self._reporter.warn(*args, **kwargs)

        def error(self, *args, **kwargs):
            if self._enabled:
                self._reporter.error(*args, **kwargs)

        def exit(self, *args, **kwargs):
            if self._enabled:
                self._reporter.exit(*args, **kwargs)

        def plain(self, *args, **kwargs):
            if self._enabled:
                self._reporter.plain(*args, **kwargs)

        def set_style(self, *args, **kwargs):
            self._reporter.set_style(args, kwargs)


    # we also equip our logger with reporting capabilities, so that we can
    # report, for example, demo output whereever we have a logger.
    logger.report = _LogReporter(logger)

    # keep the handle a round, for cleaning up on fork
    _logger_registry.add(logger)

    return logger


# -----------------------------------------------------------------------------
#
class LogReporter(object):
    """
    This class provides a wrapper around the Logger and (indirectly) Reporter
    classes, which adds uniform output filtering for the reporter while
    preserving the Reporter's API.
    """

    # --------------------------------------------------------------------------
    #
    def __init__(self, title=None, targets=['stdout'], name='radical',
                 level='REPORT'):

        self._logger = get_logger(name=name, target=targets, level=level)
        if title:
            self._logger.report.title(title)

        self.title     = self._logger.report.title
        self.header    = self._logger.report.header
        self.info      = self._logger.report.info
        self.progress  = self._logger.report.progress
        self.ok        = self._logger.report.ok
        self.warn      = self._logger.report.warn
        self.error     = self._logger.report.error
        self.exit      = self._logger.report.exit
        self.plain     = self._logger.report.plain
        self.set_style = self._logger.report.set_style



# ------------------------------------------------------------------------------

if __name__ == "__main__":

    import radical.utils as ru

    r = ru.Reporter(title='test')

    r.header  ('header  \n')
    r.info    ('info    \n')
    r.progress('progress\n')
    r.ok      ('ok      \n')
    r.warn    ('warn    \n')
    r.error   ('error   \n')
    r.plain   ('plain   \n')

    r.set_style('error', color='yellow', style='ELTTMLE', segment='X')
    r.error('error ')

    i = 0
    j = 0
    for cname,col in r.COLORS.items():
        if cname == 'reset':
            continue
        i+=1
        for mname,mod in r.COLOR_MODS.items():
            if mname == 'reset':
                continue
            j+=1
            sys.stdout.write("%s%s[%12s-%12s] " % (col, mod, cname, mname))
            sys.stdout.write("%s%s" % (r.COLORS['reset'], r.COLOR_MODS['reset']))
        sys.stdout.write("\n")
        j = 0

    r.exit    ('exit    \n', 2)

# ------------------------------------------------------------------------------

