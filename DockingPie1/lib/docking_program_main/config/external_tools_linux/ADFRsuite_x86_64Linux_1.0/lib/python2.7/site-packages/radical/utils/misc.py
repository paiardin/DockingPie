import os
import sys
import time
import regex
import signal
import socket
import netifaces
import threading
import url as ruu

# ------------------------------------------------------------------------------
#
def split_dburl(dburl, default_dburl=None) :
    """
    we split the url into the base mongodb URL, and the path element, whose
    first element is the database name, and the remainder is interpreted as
    collection id.
    """

    # if the given URL does not contain schema nor host, the default URL is used
    # as base, and the given URL string is appended to the path element.
    
    url = ruu.Url (dburl)

    if not url.schema and not url.host :
        url      = ruu.Url (default_dburl)
        url.path = dburl

    # NOTE: add other data base schemes here...
    if 'mongodb' not in url.schema.split('+'):
        raise ValueError ("url must be a 'mongodb://' or 'mongodb+ssl://' url, not '%s'" % dburl)

    host = url.host
    port = url.port
    path = url.path
    user = url.username
    pwd  = url.password
    ssl  = False

    if 'ssl' in url.schema.split('+'):
        ssl = True
        url.schema = 'mongodb'

    if not host:
        host = 'localhost'

    if  path.startswith ('/') :
        path = path[1:]
    path_elems = path.split ('/')

    dbname = None
    cname  = None
    pname  = None

    if  len(path_elems)  >  0 :
        dbname = path_elems[0]

    if  len(path_elems)  >  1 :
        dbname = path_elems[0]
        cname  = path_elems[1]

    if  len(path_elems)  >  2 :
        dbname = path_elems[0]
        cname  = path_elems[1]
        pname  = '.'.join (path_elems[2:])

    if  dbname == '.' : 
        dbname = None

    return [host, port, dbname, cname, pname, user, pwd, ssl]


# ------------------------------------------------------------------------------
#
def mongodb_connect (dburl, default_dburl=None) :
    """
    connect to the given mongodb, perform auth for the database (if a database
    was given).
    """

    try :
        import pymongo
    except ImportError :
        msg  = " \n\npymongo is not available -- install radical.utils with: \n\n"
        msg += "  (1) pip install --upgrade -e '.[pymongo]'\n"
        msg += "  (2) pip install --upgrade    'radical.utils[pymongo]'\n\n"
        msg += "to resolve that dependency (or install pymongo manually).\n"
        msg += "The first version will work for local installation, \n"
        msg += "the second one for installation from pypi.\n\n"
        raise ImportError (msg)

    [host, port, dbname, cname, pname, user, pwd, ssl] = split_dburl(dburl, default_dburl)

    mongo = pymongo.MongoClient (host=host, port=port, ssl=ssl)
    db    = None

    if  dbname :
        db = mongo[dbname]

        if  user and pwd :
            db.authenticate (user, pwd)


    else :

        # if no DB is given, we try to auth against all databases.
        for dbname in mongo.database_names () :
            try :
                mongo[dbname].authenticate (user, pwd)
            except Exception as e :
                pass 


    return mongo, db, dbname, cname, pname


# ------------------------------------------------------------------------------
#
def parse_file_staging_directives (directives) :
    """
    staging directives

       [local_path] [operator] [remote_path]

    local path: 
        * interpreted as relative to the application's working directory
        * must point to local storage (localhost)
    
    remote path
        * interpreted as relative to the job's working directory

    operator :
        * >  : stage to remote target, overwrite if exists
        * >> : stage to remote target, append    if exists
        * <  : stage to local  target, overwrite if exists
        * << : stage to local  target, append    if exists

    This method returns a tuple [src, tgt, op] for each given directive.  This
    parsing is backward compatible with the simple staging directives used
    previously -- any strings which do not contain staging operators will be
    interpreted as simple paths (identical for src and tgt), operation is set to
    '=', which must be interpreted in the caller context.  
    """

    bulk = True
    if  not isinstance (directives, list) :
        bulk       = False
        directives = [directives]

    ret = list()

    for directive in directives :

        if  not isinstance (directive, basestring) :
            raise TypeError ("file staging directives muct by of type string, "
                             "not %s" % type(directive))

        rs = regex.ReString (directive)

        if  rs // '^(?P<one>.+?)\s*(?P<op><|<<|>|>>)\s*(?P<two>.+)$' :
            res = rs.get ()
            ret.append ([res['one'], res['two'], res['op']])

        else :
            ret.append ([directive, directive, '='])

    if  bulk : return ret
    else     : return ret[0]


# ------------------------------------------------------------------------------
#
def time_stamp (spec) :

    if  isinstance (spec, int)   or \
        isinstance (spec, long)  or \
        isinstance (spec, float) :

        import datetime
        return datetime.datetime.utcfromtimestamp (spec)

    return spec


# ------------------------------------------------------------------------------
#
def time_diff (dt_abs, dt_stamp) :
    """
    return the time difference bewteen  two datetime 
    objects in seconds (incl. fractions).  Exceptions (like on improper data
    types) fall through.
    """

    delta = dt_stamp - dt_abs

    # make it easy to use seconds since epoch instead of datetime objects
    if  isinstance (delta, int)   or \
        isinstance (delta, long)  or \
        isinstance (delta, float) :
        return delta

    import datetime
    if  not isinstance  (delta, datetime.timedelta) :
        raise TypeError ("difference between '%s' and '%s' is not a .timedelta" \
                      % (type(dt_abs), type(dt_stamp)))

    # get seconds as float 
    seconds = delta.seconds + delta.microseconds/1E6
    return seconds


# ------------------------------------------------------------------------------
#
def all_pairs (iterable, n) :
    """
    [ABCD] -> [AB], [AC], [AD], [BC], [BD], [CD]
    """

    import itertools 
    return list(itertools.combinations (iterable, n))


# ------------------------------------------------------------------------------
#
def cluster_list (iterable, n) :
    """
    s -> (s0,s1,s2,...sn-1), (sn,sn+1,sn+2,...s2n-1), (s2n,s2n+1,s2n+2,...s3n-1), ...
    """

    from itertools import izip
    return izip(*[iter(iterable)]*n)


# ------------------------------------------------------------------------------
# From https://docs.python.org/release/2.3.5/lib/itertools-example.html
#
def window (seq, n=2) :
    """
    Returns a sliding window (of width n) over data from the iterable"
    s -> (s0,s1,...s[n-1]), (s1,s2,...,sn), ... 
    """

    from itertools import islice

    it = iter(seq)
    result = tuple(islice(it, n))

    if len(result) == n :
        yield result

    for elem in it :
        result = result[1:] + (elem,)
        yield result

# ------------------------------------------------------------------------------
# 
def round_to_base (value, base=1):
    """
    This method expects an integer or float value, and will round it to any
    given integer base.  For example:

      1.5, 2 -> 2
      3.5, 2 -> 4
      4.5, 2 -> 4

      11.5, 20 -> 20
      23.5, 20 -> 20
      34.5, 20 -> 40

    The default base is '1'.
    """

    return int(base * round(float(value) / base))


# ------------------------------------------------------------------------------
# 
def round_upper_bound (value):
    """
    This method expects an integer or float value, and will return an integer upper
    bound suitable for example to define plot ranges.  The upper bound is the
    smallest value larger than the input value which is a multiple of 1, 2 or
    5 times the order of magnitude (10**x) of the value.
    """

    bound = 0
    order = 0
    check = [1, 2, 5]

    while True:

        for c in check :

            bound = c*(10**order)

            if value < bound: 
                return bound

        order += 1


# ------------------------------------------------------------------------------
#
def islist(thing):
    """
    return True if a thing is a list thing, False otherwise
    """

    return isinstance(thing, list)


# ------------------------------------------------------------------------------
#
def tolist(thing):
    """
    return a non-list thing into a list thing
    """

    if islist(thing):
        return thing
    return [thing]


# ------------------------------------------------------------------------------
#
# to keep RU 2.6 compatible, we provide import_module which works around some
# quirks of __import__ when being used with dotted names. This is what the
# python docs recommend to use.  This basically steps down the module path and
# loads the respective submodule until arriving at the target.
#
def import_module(name):

    mod = __import__(name)
    for s in name.split('.')[1:]:
        mod = getattr(mod, s)
    return mod


# ------------------------------------------------------------------------------
#
def gettid():
    """
    Python is not able to give us the native thread ID.  We thus use a syscall
    to do so.  Since this is not portable, we fall back to None in case of error
    (Hi MacOS).
    """
    try:
        import ctypes
        SYS_gettid = 186
        libc = ctypes.cdll.LoadLibrary('libc.so.6')
        return int(libc.syscall(SYS_gettid))
    except:
        return None


# ------------------------------------------------------------------------------
#
_hostname = None
def get_hostname():
    """
    Look up the hostname
    """

    global _hostname
    if not _hostname:

        if socket.gethostname().find('.')>=0:
            _hostname = socket.gethostname()

        else:
            _hostname = socket.gethostbyaddr(socket.gethostname())[0]

    return _hostname
    

# ------------------------------------------------------------------------------
#
_hostip = None
def get_hostip(req=None, logger=None):
    """
    Look up the ip number for a given requested interface name.
    If interface is not given, do some magic.
    """

    AF_INET = netifaces.AF_INET

    # We create a ordered preference list, consisting of:
    #   - given arglist
    #   - white list (hardcoded preferred interfaces)
    #   - black_list (hardcoded unfavorable interfaces)
    #   - all others (whatever is not in the above)
    # Then this list is traversed, we check if the interface exists and has an
    # IP address.  The first match is used.

    if req: 
        if not isinstance(req, list):
            req = [req]
    else:
        req = []

    white_list = [
            'ipogif0', # Cray's
            'br0',     # SuperMIC
            'eth0',    # desktops etc.
            'wlan0'    # laptops etc.
            ]

    black_list = [
            'lo',      # takes the 'inter' out of the 'net'
            'sit0'     # ?
            ]

    all  = netifaces.interfaces()
    rest = [iface for iface in all \
                   if iface not in req and \
                      iface not in white_list and \
                      iface not in black_list]

    preflist = req + white_list + black_list + rest

    for iface in preflist:

        if iface not in all:
            if logger:
                logger.debug('check iface %s: does not exist', iface)
            continue

        info = netifaces.ifaddresses(iface)
        if AF_INET not in info:
            if logger:
                logger.debug('check iface %s: no information', iface)
            continue

        if not len(info[AF_INET]):
            if logger:
                logger.debug('check iface %s: insufficient information', iface)
            continue

        if not info[AF_INET][0].get('addr'):
            if logger:
                logger.debug('check iface %s: disconnected', iface)
            continue

      
        ip = info[AF_INET][0].get('addr')
        if logger:
            logger.debug('check iface %s: ip is %s', iface, ip)
        break

    _hostip = ip
    return ip


# ------------------------------------------------------------------------------
#
def watch_condition(cond, target=None, timeout=None, interval=0.1):
    """
    Watch a given condition (a callable) until it returns the target value, and
    return that value.  Stop watching on timeout, in that case return None.  The
    condition is tested approximately every 'interval' seconds.
    """
    
    start = time.time()
    while True:
        ret = cond()
        if ret == target:
            return ret
        if timeout and time.time() > start+timeout:
            return None
        time.sleep(interval)


# ------------------------------------------------------------------------------
#
def name2env(name):
    """
    convert a name of the for 'radical.pilot' to an env vare base named
    'RADICAL_PILOT'.
    """
    
    return name.replace('.', '_').upper()


# ------------------------------------------------------------------------------

