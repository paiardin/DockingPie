
__author__    = "Andre Merzky"
__copyright__ = "Copyright 2012-2013, The SAGA Project"
__license__   = "MIT"


import saga.adaptors.cpi.decorators as cpi_dec
import saga.adaptors.cpi.namespace  as cpi_ns
import saga.adaptors.cpi.attributes as cpi_att

SYNC  = cpi_dec.CPI_SYNC_CALL
ASYNC = cpi_dec.CPI_ASYNC_CALL


# keep order of inheritance!  super() below uses MRO
class Directory (cpi_ns.directory.Directory,
                 cpi_att.Attributes) :

    # ----------------------------------------------------------------
    #
    # initialization methods
    #
    def __init__ (self, api, adaptor) :

        self._cpi_nsdirec = super  (Directory, self)
        self._cpi_nsdirec.__init__ (api, adaptor)

    @SYNC
    def init_instance           (self, url, flags, session)          : pass
    @ASYNC
    def init_instance_async     (self, url, flags, session)          : pass


    # ----------------------------------------------------------------
    #
    # advert methods
    #
    @SYNC
    def set_ttl_self            (self, ttl, ttype=None)              : pass
    @ASYNC
    def set_ttl_self_async      (self, ttl, ttype=None)              : pass

    @SYNC
    def get_ttl_self            (self, ttype)                        : pass
    @ASYNC
    def get_ttl_self_async      (self, ttype)                        : pass

    @SYNC
    def set_ttl                 (self, tgt, ttl, ttype)              : pass
    @ASYNC
    def set_ttl_async           (self, tgt, ttl, ttype)              : pass

    @SYNC
    def get_ttl                 (self, tgt, ttype)                   : pass
    @ASYNC
    def get_ttl_async           (self, tgt, ttype)                   : pass

    @SYNC
    def find_adverts            (self, name_pattern, attr_pattern,
                                 obj_type, flags, ttype)             : pass
    @ASYNC
    def find_adverts_async      (self, name_pattern, attr_pattern,
                                 obj_type, flags, ttype)             : pass





