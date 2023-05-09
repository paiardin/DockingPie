
__author__    = "Andre Merzky"
__copyright__ = "Copyright 2012-2013, The SAGA Project"
__license__   = "MIT"


import saga.adaptors.cpi.decorators as cpi_dec
import saga.adaptors.cpi.namespace  as cpi_ns
import saga.adaptors.cpi.attributes as cpi_att

SYNC  = cpi_dec.CPI_SYNC_CALL
ASYNC = cpi_dec.CPI_ASYNC_CALL


# keep order of inheritance!  super() below uses MRO
class LogicalFile (cpi_ns.entry.Entry, 
                   cpi_att.Attributes) :

    # ----------------------------------------------------------------
    #
    # initialization methods
    #
    def __init__ (self, api, adaptor) :

        self._cpi_nsentry = super  (LogicalFile, self)
        self._cpi_nsentry.__init__ (api, adaptor)


    @SYNC
    def init_instance         (self, url, flags, session)      : pass
    @ASYNC
    def init_instance_async   (self, url, flags, session)      : pass


    # ----------------------------------------------------------------
    #
    # replica methods
    #
    @SYNC
    def is_file_self          (self, ttype)                    : pass
    @ASYNC
    def is_file_self_async    (self, ttype)                    : pass

    @SYNC
    def get_size_self         (self, ttype)                    : pass
    @ASYNC
    def get_size_self_async   (self, ttype)                    : pass

    @SYNC
    def add_location          (self, name, ttype)              : pass
    @ASYNC
    def add_location_async    (self, name, ttype)              : pass

    @SYNC
    def remove_location       (self, name, ttype)              : pass
    @ASYNC
    def remove_location_async (self, name, ttype)              : pass

    @SYNC
    def update_location       (self, old, new, ttype)          : pass
    @ASYNC
    def update_location_async (self, old, new, ttype)          : pass

    @SYNC
    def list_locations        (self, ttype)                    : pass
    @ASYNC
    def list_locations_async  (self, ttype)                    : pass

    @SYNC
    def replicate             (self, name, flags, ttype)       : pass
    @ASYNC
    def replicate_async       (self, name, flags, ttype)       : pass

    @SYNC
    def upload                (self, name, tgt, flags, ttype)  : pass
    @ASYNC
    def upload_async          (self, name, tgt, flags, ttype)  : pass

    @SYNC
    def download              (self, name, src, flags, ttype)  : pass
    @ASYNC
    def download_async        (self, name, src, flags, ttype)  : pass




