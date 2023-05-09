
__author__    = "Andre Merzky"
__copyright__ = "Copyright 2012-2013, The SAGA Project"
__license__   = "MIT"


import re
import os
import sys
import time
import socket
import traceback

import saga


""" Provides an assortment of utilities """

_latencies = {}


# --------------------------------------------------------------------
#
def get_trace () :

    trace = sys.exc_info ()[2]

    if  trace :
        stack           = traceback.extract_tb  (trace)
        traceback_list  = traceback.format_list (stack)
        return "".join (traceback_list)

    else :
        stack           = traceback.extract_stack ()
        traceback_list  = traceback.format_list (stack)
        return "".join (traceback_list[:-1])

# --------------------------------------------------------------------
#
def host_is_local (host) :
    """ Returns True if the given host is the localhost
    """

    if  not host                   or \
        host == 'localhost'        or \
        host == socket.gethostname () :
        return True
    else :
        return False


# --------------------------------------------------------------------
#
def host_is_valid (host) :
    """ 
    Returns True if the given hostname can be resolved.
    We also test the reverse DNS lookup -- some seriously stupid and standard
    violating internet providers implement a DNS catchall -- the reverse lookup
    can catch that case in some cases (say that quickly 3 times!)
    """

    # FIXME: cache results so that further lookups are quick

    if  host_is_local (host) :
        return True
    
    try :
        ip   = socket.gethostbyname (host)
        return True
    except :
        return False


# --------------------------------------------------------------------
#
def get_host_latency (host_url) :
    """ 
    This call measures the base tcp latency for a connection to the target
    host.  Note that port 22 is used for connection tests, unless the URL
    explicitly specifies a different port.  If the used port is blocked, the
    returned latency can be wrong by orders of magnitude.
    """

    try :

        # FIXME see comments to #62bebc9 -- this breaks for some cases, or is at
        # least annoying.  Thus we disable latency checking for the time being,
        # and return a constant assumed latency of 250ms (which approximately 
        # represents a random WAN link).
        return 0.25


        global _latencies

        if  host_url in _latencies :
            return _latencies[host_url]

        u = saga.Url (host_url)

        if u.host : host = u.host
        else      : host = 'localhost'
        if u.port : port = u.port
        else      : port = 22  #  FIXME: we should guess by protocol 

        import socket
        import time

        # ensure host is valid
        ip   = socket.gethostbyname (host)

        start = time.time ()

        s = socket.socket (socket.AF_INET, socket.SOCK_STREAM)
        s.connect  ((host, port))
        s.shutdown (socket.SHUT_RDWR)

        stop = time.time ()

        latency = stop - start

        _latencies[host_url] = latency

        return latency

    except :

        raise



# --------------------------------------------------------------------
#
def url_is_local (arg) :
    """ 
    Returns True if the given url points to localhost.

    We consider all URLs which explicitly define a port as non-local, because it
    looks like someone wants to explicitly use some protocol --
    `ssh://localost:2222/` is likely to point at an ssh tunnel.  
    
    If, however, the port matches the default port for the given protocol, we
    consider it local again -- `ssh://localhost:22/` is most likely a real local
    activity.

    Note that the schema set operates on substring level, so that we will accept
    port 22 for `ssh` and also for `sge+ssh` -- this may break in corner cases
    (altough I can't think of any right now).
    """
    
    u = saga.Url (arg)

    if  not host_is_local (u.host) :
        return False

    # host is local, but what does the port indicate?
    if u.port and u.port > 0 :
        
        try :
            if  socket.getservbyport (u.port) in u.schema :
                # some non-default port is used -- consider remote
                return False
        except :
            # unknown service port --assume this is non-standard...
            return False

    # port is not set or points to default port for service
    return True



# --------------------------------------------------------------------
#
def url_is_relative (url_1) :
    """ an URL is considered relative if it only contains a path element, and
    that path element does not start with '/'.
    """

    u1 = saga.Url (url_1)

    if  str (u1) == str(u1.path) :
        if  not u1.path :
            return True 
        elif u1.path[0] != '/' :
            return True

    return False


# --------------------------------------------------------------------
#
def url_get_dirname (url_1) :
    """ 
    Extract the directory part of the given URL's path element.  We consider
    everything up to the last '/' as directory.  That also holds for relative
    paths.
    """

    u1 = saga.Url (url_1)
    p1 = u1.path

    return re.sub (r"[^/]*$", "", p1)


# --------------------------------------------------------------------
#
def url_get_filename (url_1) :
    """ 
    Extract the directory part of the given URL's path element.  We consider
    everything up to the last '/' as directory.  That also holds for relative
    paths.
    """

    u1 = saga.Url (url_1)
    p1 = u1.path

    if '/' in p1 :
        return re.sub (r"^.*/(.*)$", "\1", p1)
    else :
        return p1


# --------------------------------------------------------------------
#
def url_normalize (url_1) :
    """ 
    The path element of the URL is normalized
    """

    ret      = saga.Url (url_1)
    ret.path = os.path.normpath (ret.path)

    return ret


# --------------------------------------------------------------------
#
def url_make_absolute (url_1, url_2) :
    """ 
    URL1 is expected to only have a path
    Missing elements in url_1 are copied from url_2 -- but path stays the
    same.  Unless it is a relative path in the first place: then it is
    interpreted as relative to url_2.path, and is made absolute.
    protocol/port/user etc.
    """

    if  not isinstance(url_1, saga.Url):
        url_1 = saga.Url(url_1)

    if  not isinstance(url_2, saga.Url):
        url_2 = saga.Url(url_2)

    if  not url_is_compatible (url_1, url_2) :
        raise saga.BadParameter ("Cannot interpret url %s in the context of url %s" \
                              % (url_2, url_1))

    # re-interpret path of url_2, using url_1 as base directory
    ret = saga.Url (url_1)

    if  url_is_relative (url_2) :
        # note that we have no means if 'file://localhost/tmp/test.txt' refers
        # to a file or an directory -- so we have to trust that url_1 is
        # a dir...
        ret.path = url_1.path + '/' + url_2.path
    else :
        # absolute path, replace url path...
        ret.path = url_2.path


    return ret


# --------------------------------------------------------------------
#
def url_is_compatible (url_1, url_2) :
    """ 
    Returns True if the given urls point to the same host, using the same
    protocol/port/user etc.  If one of the URLs only contains a path, it is
    considered compatible with any other URL.
    """
    
    u1 = saga.Url (url_1)
    u2 = saga.Url (url_2)

    # if either one url only contains a path, it is compatible to anything.
    if os.path.normpath(u1.path) == os.path.normpath (str(u1)) : return True
    if os.path.normpath(u2.path) == os.path.normpath (str(u2)) : return True

    # more than path in both URLs -- check compatibility for all elements
    if u1.scheme   and     u2.scheme   and u1.scheme   != u2.scheme   : return False 
    if u1.host     and     u2.host     and u1.host     != u2.host     : return False
    if u1.port     and     u2.port     and u1.port     != u2.port     : return False
    if u1.username and     u2.username and u1.username != u2.username : return False
    if u1.password and     u2.password and u1.password != u2.password : return False

    if u1.scheme   and not u2.scheme                                  : return False 
    if u1.host     and not u2.host                                    : return False
    if u1.port     and not u2.port                                    : return False
    if u1.username and not u2.username                                : return False
    if u1.password and not u2.password                                : return False

    if u2.scheme   and not u1.scheme                                  : return False 
    if u2.host     and not u1.host                                    : return False
    if u2.port     and not u1.port                                    : return False
    if u2.username and not u1.username                                : return False
    if u2.password and not u1.password                                : return False

    # no differences detected (ignored fragments and query though)
    return True


# --------------------------------------------------------------------
#
def normalize_version (v) :
    """
    For a given version string (numeric only!), return an ordered tuple of
    integers, removing trailing zeros.  That tuple can then be used for
    comparison.
    """
    # parts = [int (x) for x in v.split (".")]
    # while parts[-1] == 0:
    #     parts.pop ()
    # return parts

    return tuple (v.split ("."))


# --------------------------------------------------------------------



