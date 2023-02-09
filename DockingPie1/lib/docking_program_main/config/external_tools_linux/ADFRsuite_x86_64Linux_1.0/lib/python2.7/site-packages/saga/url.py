
__author__    = "Andre Merzky, Ole Weidner"
__copyright__ = "Copyright 2012-2013, The SAGA Project"
__license__   = "MIT"


import radical.utils            as ru
import radical.utils.signatures as rus


# ------------------------------------------------------------------------------
#
class Url (ru.Url) :
    """ 
    The SAGA Url class.

    URLs are used in several places in the SAGA API: to specify service
    endpoints for job submission or resource management, for file or
    directory locations, etc.

    The URL class is designed to simplify URL management for these
    purposes -- it allows to manipulate individual URL elements, while
    ensuring that the resulting URL is well formatted. Example::

      # create a URL from a string
      location = saga.Url ("file://localhost/tmp/file.dat")
      d = saga.filesystem.Directory(location)

    A URL consists of the following components (where one ore more can
    be 'None')::

      <scheme>://<user>:<pass>@<host>:<port>/<path>?<query>#<fragment>

    Each of these components can be accessed via its property or
    alternatively, via getter / setter methods. Example::

      url = saga.Url ("scheme://pass:user@host:123/path?query#fragment")

      # modify the scheme
      url.scheme = "anotherscheme"

      # above is equivalent with
      url.set_scheme("anotherscheme")
    """

    @rus.takes   ('Url', 
                  rus.optional (basestring, 'Url'))
    @rus.returns (rus.nothing)
    def __init__ (self, url_in=''):
        """ 
        __init__(url_string='')

        Create a new Url object from a string or another Url object.
        """

        self._super = super (Url, self)
        self._super.__init__ (url_in)


    # --------------------------------------------------------------------------
    #
    @rus.takes   ('Url', 
                 ('Url', dict))
    @rus.returns ('Url')
    def __deepcopy__ (self, memo) :
        """ 
        __deepcopy__(self, memo)

        Deep copy of a Url
        """

        return Url (self)


# ------------------------------------------------------------------------------
#


