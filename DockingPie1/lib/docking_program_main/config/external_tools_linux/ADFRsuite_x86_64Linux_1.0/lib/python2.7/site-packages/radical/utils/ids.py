
__author__    = "Radical.Utils Development Team (Andre Merzky)"
__copyright__ = "Copyright 2013, RADICAL@Rutgers"
__license__   = "MIT"


import os
import time
import uuid
import fcntl
import socket
import datetime
import threading
import singleton


# ------------------------------------------------------------------------------
#
class _IDRegistry(object):
    """
    This helper class (which is not exposed to any user of radical.utils)
    generates a sequence of continous numbers for each known ID prefix.  It is
    a singleton, and thread safe (assuming that the Singleton metaclass supports
    thread safe construction).
    """

    __metaclass__ = singleton.Singleton


    # --------------------------------------------------------------------------
    def __init__(self):
        """
        Initialized the registry dict and the threading lock
        """

        self._rlock    = threading.RLock()
        self._registry = dict()


    # --------------------------------------------------------------------------
    def get_counter(self, prefix):
        """
        Obtain the next number in the sequence for the given prefix.
        If the prefix is not known, a new registry counter is created.
        """

        with self._rlock:

            if not prefix in self._registry:
                self._registry[prefix] = 0

            ret = self._registry[prefix]

            self._registry[prefix] += 1

            return ret


    # --------------------------------------------------------------------------
    def reset_counter(self, prefix, reset_all_others=False):
        """
        Reset the given counter to zero.
        """

        with self._rlock:
            
            if reset_all_others:
                # reset all counters *but* the one given
                for p in self._registry:
                    if p != prefix:
                        self._registry[p] = 0
            else:
                self._registry[prefix] = 0


# ------------------------------------------------------------------------------
#
# we create on private singleton instance for the ID registry.
_id_registry = _IDRegistry()
_BASE        = "%s/.radical/utils" % os.environ.get("HOME", "/tmp")
os.system("mkdir -p %s" % _BASE)

# ------------------------------------------------------------------------------
#
ID_SIMPLE  = 'simple'
ID_UNIQUE  = 'unique'
ID_PRIVATE = 'private'
ID_CUSTOM  = 'custom'
ID_UUID    = 'uiud'

# ------------------------------------------------------------------------------
#
def generate_id(prefix, mode=ID_SIMPLE):
    """
    Generate a human readable, sequential ID for the given prefix.

    The ID is by default very simple and thus very readable, but cannot be
    assumed to be globally unique -- simple ID uniqueness is only guaranteed
    within the scope of one python instance.

    If `mode` is set to the non-default type `ID_UNIQUE`, an attempt is made to
    generate readable but globally unique IDs -- although the level of
    confidence for uniqueness is significantly smaller than for, say UUIDs.

    Examples::

        print radical.utils.generate_id('item.')
        print radical.utils.generate_id('item.')
        print radical.utils.generate_id('item.', mode=radical.utils.ID_SIMPLE)
        print radical.utils.generate_id('item.', mode=radical.utils.ID_SIMPLE)
        print radical.utils.generate_id('item.', mode=radical.utils.ID_UNIQUE)
        print radical.utils.generate_id('item.', mode=radical.utils.ID_UNIQUE)
        print radical.utils.generate_id('item.', mode=radical.utils.ID_PRIVATE)
        print radical.utils.generate_id('item.', mode=radical.utils.ID_PRIVATE)
        print radical.utils.generate_id('item.', mode=radical.utils.ID_UUID)

    The above will generate the IDs:

        item.0001
        item.0002
        item.0003
        item.0004
        item.2014.07.30.13.13.44.0001
        item.2014.07.30.13.13.44.0002
        item.cameo.merzky.021342.0001
        item.cameo.merzky.021342.0002
        item.23cacb7e-0b08-11e5-9f0f-08002716eaa9

    where 'cameo' is the (short) hostname, 'merzky' is the username, and '02134'
    is 'days since epoch'.  The last element, the counter is unique for each id
    type and item type, and restarts for each session (application process).  In
    the last case though (`ID_PRIVATE`), the counter is reset for every new day,
    and can thus span multiple applications.
    """

    if not prefix or \
        not isinstance(prefix, basestring):
        raise TypeError("ID generation expect prefix in basestring type")

    template = ""

    if   mode == ID_SIMPLE : template = "%(prefix)s.%(counter)04d"
    elif mode == ID_UNIQUE : template = "%(prefix)s.%(date)s.%(time)s.%(pid)06d.%(counter)04d"
    elif mode == ID_PRIVATE: template = "%(prefix)s.%(host)s.%(user)s.%(days)06d.%(day_counter)04d"
    elif mode == ID_CUSTOM : template = prefix
    elif mode == ID_UUID   : template = "%(prefix)s.%(uuid)s"
    else:
        raise ValueError("mode '%s' not supported for ID generation", mode)

    return _generate_id(template, prefix)

# ------------------------------------------------------------------------------
#
def _generate_id(template, prefix):

    # FIXME: several of the vars below are constants, and many of them are
    # rarely used in IDs.  They should be created only once per module instance,
    # and/or only if needed.

    import getpass

    # seconds since epoch(float), and timestamp
    seconds = time.time()
    now     = datetime.datetime.fromtimestamp(seconds)
    days    = int(seconds / (60*60*24))

    try:
        user = getpass.getuser()
    except Exception:
        user = 'nobody'

    info = dict()

    info['day_counter' ]  = 0
    info['item_counter']  = 0
    info['counter'     ]  = 0
    info['prefix'      ]  = prefix
    info['now'         ]  = now
    info['seconds'     ]  = int(seconds)              # full seconds since epoch
    info['days'        ]  = days                      # full days since epoch
    info['user'        ]  = user                      # local username
    info['date'        ]  = "%04d.%02d.%02d" % (now.year, now.month,  now.day)
    info['time'        ]  = "%02d.%02d.%02d" % (now.hour, now.minute, now.second)
    info['pid'         ]  = os.getpid()

    # the following ones are time consuming, and only done when needed
    if '%(host)' in template: info['host'] = socket.gethostname() # local hostname
    if '%(uuid)' in template: info['uuid'] = uuid.uuid1()         # pain old uuid

    if '%(day_counter)' in template:
        fd = os.open("%s/rp_%s_%s.cnt" % (_BASE, user, days), os.O_RDWR | os.O_CREAT)
        fcntl.flock(fd, fcntl.LOCK_EX)
        os.lseek(fd, 0, os.SEEK_SET )
        data = os.read(fd, 256)
        if not data: data = 0
        info['day_counter'] = int(data)
        os.lseek(fd, 0, os.SEEK_SET )
        os.write(fd, "%d\n" % (info['day_counter']+1))
        os.close(fd)

    if '%(item_counter)' in template:
        fd = os.open("%s/rp_%s_%s.cnt" % (_BASE, user, prefix), os.O_RDWR | os.O_CREAT)
        fcntl.flock(fd, fcntl.LOCK_EX)
        os.lseek(fd, 0, os.SEEK_SET)
        data = os.read(fd, 256)
        if not data: data = 0
        info['item_counter'] = int(data)
        os.lseek(fd, 0, os.SEEK_SET)
        os.write(fd, "%d\n" % (info['item_counter']+1))
        os.close(fd)

    if '%(counter)' in template:
        info['counter'] = _id_registry.get_counter(prefix.replace('%', ''))


    ret = template % info

    if '%(' in ret:
        import pprint
        pprint.pprint(info)
        print template
        print ret
        raise ValueError('unknown replacement pattern in template (%s)' % template)

    return ret


# ------------------------------------------------------------------------------
#
def reset_id_counters(prefix=None, reset_all_others=False):

    if not isinstance(prefix, list):
        prefix = [prefix]

    for p in prefix:
        _id_registry.reset_counter(p.replace('%', ''), reset_all_others)


# ------------------------------------------------------------------------------

