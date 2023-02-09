
__author__    = "Andre Merzky"
__copyright__ = "Copyright 2012-2013, The SAGA Project"
__license__   = "MIT"


""" Provides the CPI base class and a number of call decorators. """

import weakref

import radical.utils         as ru
import radical.utils.config  as ruc

from   saga.exceptions import *


# ------------------------------------------------------------------------------
#
# CPI base class
#
class CPIBase (ruc.Configurable) :

    # --------------------------------------------------------------------------
    #
    def __init__ (self, api, adaptor) :

        self._session   = None
        self._adaptor   = adaptor
        self._cpi_cname = self.__class__.__name__
        self._logger    = ru.get_logger('radical.saga.cpi')

        # The API object must obviously keep an adaptor instance.  If we also
        # keep an API instance ref in this adaptor base, we create a ref cycle
        # which will annoy (i.e. disable) garbage collection.  We thus use weak
        # references to break that cycle.  The inheriting classes MUST use
        # get_api() to obtain the API reference.
      # if  api :
      #     self._api   = weakref.ref (api)
      # else :
      #     self._api   = None
        self._api   = weakref.ref (api)

        # by default, we assume that no bulk optimizations are supported by the
        # adaptor class.  Any adaptor class supporting bulks ops must overwrite
        # the ``_container`` attribute (via
        # ``self._set_container(container=None)``, and have it point to the
        # class which implements the respective ``container_*`` methods.
        self._container = None


    def _set_container (self, container=None) :
        self._container = container


    def get_cpi_cname (self) :
        return self._cpi_cname


    def get_api (self) :

        # get api from weakref.  We can be quite confident that the api
        # object has *not* been garbage collected, yet, as it obviously is
        # still binding this adaptor instance.
        return self._api ()

    #   if self._api :
    #       # get api from weakref.  We can be quite confident that the api
    #       # object has *not* been garbage collected, yet, as it obviously is
    #       # still binding this adaptor instance.
    #       return self._api ()
    #   else :
    #       # no need to de-weakref 'None'
    #       return self._api


    def get_adaptor_name (self) :
        return self._adaptor.get_name ()


    def _set_session (self, session) :
        self._session = session


    def get_session (self) :
        return self._session

    @property
    def session(self):
        return self._session

    @session.setter
    def session(self, session):
        self._session = session
