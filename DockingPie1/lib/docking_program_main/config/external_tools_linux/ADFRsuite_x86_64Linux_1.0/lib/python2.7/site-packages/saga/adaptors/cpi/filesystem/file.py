
__author__    = "Andre Merzky"
__copyright__ = "Copyright 2012-2013, The SAGA Project"
__license__   = "MIT"


import saga.adaptors.cpi.decorators as cpi_dec
import saga.adaptors.cpi.namespace  as cpi_ns

SYNC  = cpi_dec.CPI_SYNC_CALL
ASYNC = cpi_dec.CPI_ASYNC_CALL


# keep order of inheritance!  super() below uses MRO
class File (cpi_ns.entry.Entry) :

    # ----------------------------------------------------------------
    #
    # initialization methods
    #
    def __init__ (self, api, adaptor) :

        self._cpi_nsentry = super  (File, self)
        self._cpi_nsentry.__init__ (api, adaptor)


    @SYNC
    def init_instance        (self, url, flags, session)      : pass
    @ASYNC
    def init_instance_async  (self, url, flags, session)      : pass


    # ----------------------------------------------------------------
    #
    # filesystem methods
    #
    @SYNC
    def is_file_self         (self, ttype)                    : pass
    @ASYNC
    def is_file_self_async   (self, ttype)                    : pass

    @SYNC
    def get_size_self        (self, ttype)                    : pass
    @ASYNC
    def get_size_self_async  (self, ttype)                    : pass

    @SYNC
    def read                 (self, size, ttype)              : pass
    @ASYNC
    def read_async           (self, size, ttype)              : pass

    @SYNC
    def write                (self, data, ttype)              : pass
    @ASYNC
    def write_async          (self, data, ttype)              : pass

    @SYNC
    def seek                 (self, off, whence, ttype)       : pass
    @ASYNC
    def seek_async           (self, off, whence, ttype)       : pass

    @SYNC
    def read_v               (self, iovecs, ttype)            : pass
    @ASYNC
    def read_v_async         (self, iovecs, ttype)            : pass

    @SYNC
    def write_v              (self, data, ttype)              : pass
    @ASYNC
    def write_v_async        (self, data, ttype)              : pass

    @SYNC
    def size_p               (self, pattern, ttype)           : pass
    @ASYNC
    def size_p_async         (self, pattern, ttype)           : pass

    @SYNC
    def read_p               (self, pattern, ttype)           : pass
    @ASYNC
    def read_p_async         (self, pattern, ttype)           : pass

    @SYNC
    def write_p              (self, pattern, data, ttype)     : pass
    @ASYNC
    def write_p_async        (self, pattern, data, ttype)     : pass

    @SYNC
    def modes_e              (self, ttype)                    : pass
    @ASYNC
    def modes_e_async        (self, ttype)                    : pass

    @SYNC
    def size_e               (self, emode, spec, ttype)       : pass
    @ASYNC
    def size_e_async         (self, emode, spec, ttype)       : pass

    @SYNC
    def read_e               (self, emode, spec, ttype)       : pass
    @ASYNC
    def read_e_async         (self, emode, spec, ttype)       : pass

    @SYNC
    def write_e              (self, emode, spec, data, ttype) : pass
    @ASYNC
    def write_e_async        (self, emode, spec, data, ttype) : pass




