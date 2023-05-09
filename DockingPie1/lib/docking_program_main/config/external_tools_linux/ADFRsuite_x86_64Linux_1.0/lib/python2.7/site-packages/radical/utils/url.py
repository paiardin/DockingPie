
__author__    = "Radical.Utils Development Team (Andre Merzky, Ole Weidner)"
__copyright__ = "Copyright 2013, RADICAL@Rutgers"
__license__   = "MIT"


import os

import contrib.urlparse25 as urlparse
import signatures         as rus


# ------------------------------------------------------------------------------
#
class Url (object):
    """ The RADICAL Url class.

        URLs are used in several places in the RADICAL software projects: to
        specify service endpoints for job submission or resource management, for
        file or directory locations, etc.

        The URL class is designed to simplify URL management for these
        purposes -- it allows to manipulate individual URL elements, while
        ensuring that the resulting URL is well formatted. Example::

          # create a URL from a string
          location = radical.utils.Url ("file://localhost/tmp/file.dat")
          d = radical.utils.filesystem.Directory(location)

        A URL consists of the following components (where one ore more can
        be 'None')::

          <scheme>://<user>:<pass>@<host>:<port>/<path>?<query>#<fragment>

        Each of these components can be accessed via its property or
        alternatively, via getter / setter methods. Example::

          url = radical.utils.Url ("scheme://pass:user@host:123/path?query#fragment")

          # modify the scheme
          url.scheme = "anotherscheme"

          # above is equivalent with
          url.set_scheme("anotherscheme")
    """

    # --------------------------------------------------------------------------
    #
    @rus.takes   ('Url', 
                  rus.optional ((basestring, 'Url')))
    @rus.returns (rus.nothing)
    def __init__(self, url_in=''):
        """ 
        __init__(url_in='')

        Create a new Url object from a string or another Url object.
        """

        if  not url_in :
            url_in = ""

        self._urlobj = urlparse.urlparse (str (url_in), allow_fragments=True)
        self._renew_url ()

    # --------------------------------------------------------------------------
    #
    ##
    @rus.takes   ('Url')
    @rus.returns ((rus.nothing, basestring))
    def __str__  (self):
        """
        __str__()

        String representation.
        """
        return self._urlobj.geturl()

    # --------------------------------------------------------------------------
    #
    ##
    @rus.takes   ('Url')
    @rus.returns (basestring)
    def __unicode__(self):
        """ 
        __unicode__()

        Unicode representation.
        """
        return u'%s'  %  unicode(self._urlobj.geturl())


    # --------------------------------------------------------------------------
    #
    ##
    @rus.takes   ('Url', 
                  ('Url', dict))
    @rus.returns ('Url')
    def __deepcopy__(self, memo):
        """ 
        __deepcopy__(self, memo)

        Deep copy of a Url
        """
        return Url(self)


    # --------------------------------------------------------------------------
    #
    ##
    @rus.takes   ('Url')
    @rus.returns (bool)
    def __nonzero__(self) :

        if str(self):
            return 1
        else:
            return 0


    # --------------------------------------------------------------------------
    #
    ##
    @rus.takes   ('Url', 
                  rus.optional(basestring),
                  rus.optional(basestring),
                  rus.optional(basestring),
                  rus.optional((basestring, int)))
    @rus.returns (basestring)
    def _make_netloc (self, username, password, hostname, port):
        """ 
        _make_netloc(self, username, password, hostname, port)

        Private helper function to generate netloc string.
        """
        netloc = str()

        if  username :
            if  password : netloc += "%s:%s@" % (username, password)
            else :         netloc += "%s@"    % (username)
        if  hostname :     netloc += "%s"     % (hostname)
        if  port     :     netloc += ":%s"    % (port)

        return netloc


    def _renew_netloc (self, username='', password='', hostname='', port='') :

        newloc = self._make_netloc (username or self._urlobj.username,
                                    password or self._urlobj.password,
                                    hostname or self._urlobj.hostname,
                                    port     or self._urlobj.port)

        self._renew_url (netloc=newloc)


    def _renew_url (self, scheme='', netloc='', path='', 
                          params='', query='',  fragment='', force_path=False) :

        # always normalize the path.  
        path = self.normpath(path)

        if  force_path :
            forced_path = path or ''
        else :
            forced_path = path or self._urlobj.path


        newurl = urlparse.urlunparse ((scheme   or self._urlobj.scheme,
                                       netloc   or self._urlobj.netloc,
                                       forced_path,
                                       params   or self._urlobj.params,
                                       query    or self._urlobj.query,
                                       fragment or self._urlobj.fragment))

        self._urlobj = urlparse.urlparse (newurl, allow_fragments=True)


    def normpath(self, path):

        if  path :

            # Alas, os.path.normpath removes trailing slashes, 
            # so we re-add them.
            if len(path) > 1 and path.endswith('/'):
                trailing_slash=True
            else:
                trailing_slash=False

            path = os.path.normpath (path)

            if trailing_slash and not path.endswith('/'):
                path += '/'

        return path


    # --------------------------------------------------------------------------
    #
    # Scheme property
    #
    @rus.takes   ('Url', 
                  (rus.nothing, basestring))
    @rus.returns (rus.nothing)
    def set_scheme(self, scheme):
        """ 
        set_scheme(scheme)

        Set the URL 'scheme' component.

        :param scheme: The new scheme
        :type  scheme: str

        """
        self._renew_url (scheme=scheme)


    @rus.takes   ('Url')
    @rus.returns ((rus.nothing, basestring))
    def get_scheme(self):
        """
        get_scheme()

        Return the URL 'scheme' component.
        """
        return self._urlobj.scheme

    scheme = property(get_scheme, set_scheme)
    schema = scheme  # alias, as both terms are used...
    """ The scheme component.  """


    # --------------------------------------------------------------------------
    #
    # Host property
    #
    @rus.takes   ('Url', 
                  (rus.nothing, basestring))
    @rus.returns (rus.nothing)
    def set_host(self, hostname):
        """ 
        set_host(hostname)

        Set the 'hostname' component.

        :param hostname: The new hostname
        :type  hostname: str
        """
        netloc = self._renew_netloc (hostname=hostname)


    @rus.takes   ('Url')
    @rus.returns ((rus.nothing, basestring))
    def get_host(self):
        """ 
        get_host()

        Return the URL 'hostname' component.
        """
        return self._urlobj.hostname

    host = property(get_host, set_host)
    """ The hostname component.  """


    # --------------------------------------------------------------------------
    #
    # Port property
    #
    @rus.takes   ('Url', 
                  (rus.nothing, basestring, int))
    @rus.returns (rus.nothing)
    def set_port(self, port):
        """ 
        set_port(port)

        Set the URL 'port' component.

        :param port: The new port
        :type  port: int
        """

        self._renew_netloc (port=port)

    @rus.takes   ('Url')
    @rus.returns ((rus.nothing, int))
    def get_port(self):
        """ 
        get_port()

        Return the URL 'port' component.
        """
        if self._urlobj.port is not None:
            return int(self._urlobj.port)
        else:
            return None

    port = property(get_port, set_port)
    """ The port component.  """


    # --------------------------------------------------------------------------
    #
    # Username property
    #
    @rus.takes   ('Url', 
                  (rus.nothing, basestring))
    @rus.returns (rus.nothing)
    def set_username(self, username):
        """ 
        set_username(username)

        Set the URL 'username' component.

        :param username: The new username
        :type  username: str
        """
        self._renew_netloc (username=username)

    @rus.takes   ('Url')
    @rus.returns ((rus.nothing, basestring))
    def get_username(self):
        """ 
        get_username()

        Return the URL 'username' component.
        """
        return self._urlobj.username

    username = property(get_username, set_username)
    """ The username component.  """


    # --------------------------------------------------------------------------
    #
    # Password property
    #
    @rus.takes   ('Url', 
                  (rus.nothing, basestring))
    @rus.returns (rus.nothing)
    def set_password(self, password):
        """ 
        set_password(password)

        Set the URL 'password' component.

        :param password: The new password
        :type password:  str
        """
        self._renew_netloc (password=password)

    @rus.takes   ('Url')
    @rus.returns ((rus.nothing, basestring))
    def get_password(self):
        """ 
        get_password()

        Return the URL 'username' component.
        """
        return self._urlobj.password

    password = property(get_password, set_password)
    """ The password component.  """


    # --------------------------------------------------------------------------
    #
    # Fragment property
    #
    @rus.takes   ('Url', 
                  (rus.nothing, basestring))
    @rus.returns (rus.nothing)
    def set_fragment(self, fragment):
        """ 
        set_fragment(fragment)

        Set the URL 'fragment' component.

        :param fragment: The new fragment
        :type fragment:  str
        """
        self._renew_url (fragment=fragment)

    @rus.takes   ('Url')
    @rus.returns ((rus.nothing, basestring))
    def get_fragment(self):
        """ 
        get_fragment()

        Return the URL 'fragment' component.
        """
        return self._urlobj.fragment

    fragment = property(get_fragment, set_fragment)
    """ The fragment component.  """


    # --------------------------------------------------------------------------
    #
    # Path property
    #
    @rus.takes   ('Url', 
                  (rus.nothing, basestring))
    @rus.returns (rus.nothing)
    def set_path(self, path):
        """ 
        set_path(path)

        Set the URL 'path' component.

        :param path: The new path
        :type path:  str
        """
        self._renew_url (path=path, force_path=True)

    @rus.takes   ('Url')
    @rus.returns ((rus.nothing, basestring))
    def get_path(self):
        """ 
        get_path()

        Return the URL 'path' component.
        """

        path = self._urlobj.path

        if '?' in path:
            (path, query) = path.split('?')

        return self.normpath(path)

    path = property(get_path, set_path)
    """ The path component.  """


    # --------------------------------------------------------------------------
    #
    # Query property
    #
    @rus.takes   ('Url', 
                  (rus.nothing, basestring))
    @rus.returns (rus.nothing)
    def set_query(self, query):
        """ 
        set_query(query)

        Set the URL 'query' component.

        :param query: The new query
        :type query:  str
        """
        self._renew_url (query=query)

    @rus.takes   ('Url')
    @rus.returns ((rus.nothing, basestring))
    def get_query(self):
        """
        get_query()

        Return the URL 'query' component.
        """
        if not self._urlobj.query:
            if '?' in self._urlobj.path:
                (path, query) = self._urlobj.path.split('?')
                return query
        else:
            return self._urlobj.query

    query = property(get_query, set_query)
    """ The query component.  """
   

# --------------------------------------------------------------------
#


