
__author__    = "Andre Merzky"
__copyright__ = "Copyright 2012-2013, The SAGA Project"
__license__   = "MIT"


import saga.adaptors.cpi.decorators as cpi_dec
import saga.adaptors.cpi.namespace  as cpi_ns
import saga.adaptors.cpi.attributes as cpi_att

SYNC  = cpi_dec.CPI_SYNC_CALL
ASYNC = cpi_dec.CPI_ASYNC_CALL


# keep order of inheritance!  super() below uses MRO
class LogicalDirectory (cpi_ns.directory.Directory, 
                        cpi_att.Attributes) :

    # ----------------------------------------------------------------
    #
    # initialization methods
    #
    def __init__ (self, api, adaptor) :

        self._cpi_nsdirec = super  (LogicalDirectory, self)
        self._cpi_nsdirec.__init__ (api, adaptor)

    @SYNC
    def init_instance           (self, url, flags, session)          : pass
    @ASYNC
    def init_instance_async     (self, url, flags, session)          : pass


    @SYNC
    def open                    (self, tgt, flags, ttype)            : pass
    @ASYNC
    def open_async              (self, tgt, flags, ttype)            : pass

    @SYNC
    def open_dir                (self, tgt, flags, ttype)            : pass
    @ASYNC
    def open_dir_async          (self, tgt, flags, ttype)            : pass

    # ----------------------------------------------------------------
    #
    # replica methods
    #
    @SYNC
    def get_size                (self, tgt, ttype)                   : pass
    @ASYNC
    def get_size_async          (self, tgt, ttype)                   : pass

    @SYNC
    def is_file                 (self, tgt, ttype)                   : pass
    @ASYNC
    def is_file_async           (self, tgt, ttype)                   : pass

    @SYNC
    def is_file_self            (self, ttype)                        : pass
    @ASYNC
    def is_file_self_async      (self, ttype)                        : pass

    @SYNC
    def find_replicas       (self, name_pattern, attr_pattern, flags, ttype)  : pass
    @ASYNC
    def find_replicas_async (self, name_pattern, attr_pattern, flags, ttype)  : pass




