
__author__    = "Andre Merzky"
__copyright__ = "Copyright 2012-2013, The SAGA Project"
__license__   = "MIT"


import saga.adaptors.cpi.decorators as cpi_dec
import saga.adaptors.cpi.namespace  as cpi_ns
import saga.adaptors.cpi.attributes as cpi_att

SYNC  = cpi_dec.CPI_SYNC_CALL
ASYNC = cpi_dec.CPI_ASYNC_CALL


# keep order of inheritance!  super() below uses MRO
class Entry (cpi_ns.entry.Entry, 
             cpi_att.Attributes) :

    # ----------------------------------------------------------------
    #
    # initialization methods
    #
    def __init__ (self, api, adaptor) :

        self._cpi_nsentry = super  (Entry, self)
        self._cpi_nsentry.__init__ (api, adaptor)


    @SYNC
    def init_instance         (self, url, flags, session)      : pass
    @ASYNC
    def init_instance_async   (self, url, flags, session)      : pass


    # ----------------------------------------------------------------
    #
    # advert methods
    #
    @SYNC
    def set_ttl                 (self, ttl, ttype=None)              : pass
    @ASYNC
    def set_ttl_async           (self, ttl, ttype=None)              : pass

    @SYNC
    def get_ttl                 (self, ttype)                        : pass
    @ASYNC
    def get_ttl_async           (self, ttype)                        : pass

    @SYNC
    def store_object            (self, object, ttype)                : pass
    @ASYNC
    def store_object_async      (self, object, ttype)                : pass

    @SYNC
    def retrieve_object         (self, ttype)                        : pass
    @ASYNC
    def retrieve_object_async   (self, ttype)                        : pass

    @SYNC
    def delete_object           (self, ttype)                        : pass
    @ASYNC
    def delete_object_async     (self, ttype)                        : pass




