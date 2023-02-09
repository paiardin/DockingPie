
__author__    = "Andre Merzky"
__copyright__ = "Copyright 2012-2013, The SAGA Project"
__license__   = "MIT"


import saga.adaptors.cpi.decorators as cpi_dec
import saga.adaptors.cpi.base       as cpi_base
import saga.adaptors.cpi.async      as cpi_async

SYNC  = cpi_dec.CPI_SYNC_CALL
ASYNC = cpi_dec.CPI_ASYNC_CALL


# keep order of inheritance!  super() below uses MRO
class Entry (cpi_base.CPIBase, cpi_async.Async) :

    # ----------------------------------------------------------------
    #
    # initialization methods
    #
    def __init__ (self, api, adaptor) :

        _cpi_base = super  (Entry, self)
        _cpi_base.__init__ (api, adaptor)

    @SYNC
    def init_instance        (self, url, flags, session)      : pass
    @ASYNC
    def init_instance_async  (self, url, flags, session)      : pass


    # ----------------------------------------------------------------
    #
    # namespace entry methods
    #
    @SYNC
    def close                (self, timeout, ttype)           : pass
    @ASYNC
    def close_async          (self, timeout, ttype)           : pass

    @SYNC
    def get_url              (self, ttype)                    : pass
    @ASYNC
    def get_url_async        (self, ttype)                    : pass

    @SYNC
    def get_cwd              (self, ttype)                    : pass
    @ASYNC
    def get_cwd_async        (self, ttype)                    : pass

    @SYNC
    def get_name             (self, ttype)                    : pass
    @ASYNC
    def get_name_async       (self, ttype)                    : pass

    @SYNC
    def is_dir_self          (self, ttype)                    : pass
    @ASYNC
    def is_dir_self_async    (self, ttype)                    : pass

    @SYNC
    def is_entry_self        (self, ttype)                    : pass
    @ASYNC
    def is_entry_self_async  (self, ttype)                    : pass

    @SYNC
    def is_link_self         (self, ttype)                    : pass
    @ASYNC
    def is_link_self_async   (self, ttype)                    : pass

    @SYNC
    def read_link_self       (self, ttype)                    : pass
    @ASYNC
    def read_link_self_async (self, ttype)                    : pass

    @SYNC
    def copy_self            (self, tgt, flags, ttype)        : pass
    @ASYNC
    def copy_self_async      (self, tgt, flags, ttype)        : pass

    @SYNC
    def link_self            (self, tgt, flags, ttype)        : pass
    @ASYNC
    def link_self_async      (self, tgt, flags, ttype)        : pass

    @SYNC
    def move_self            (self, tgt, flags, ttype)        : pass
    @ASYNC
    def move_self_async      (self, tgt, flags, ttype)        : pass

    @SYNC
    def remove_self          (self, flags, ttype)             : pass
    @ASYNC
    def remove_self_async    (self, flags, ttype)             : pass




