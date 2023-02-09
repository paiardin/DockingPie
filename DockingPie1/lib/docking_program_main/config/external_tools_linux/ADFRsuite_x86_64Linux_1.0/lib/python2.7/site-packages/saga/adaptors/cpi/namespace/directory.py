
__author__    = "Andre Merzky"
__copyright__ = "Copyright 2012-2013, The SAGA Project"
__license__   = "MIT"


import saga.adaptors.cpi.decorators as cpi_dec
import entry

SYNC  = cpi_dec.CPI_SYNC_CALL
ASYNC = cpi_dec.CPI_ASYNC_CALL


# keep order of inheritance!  super() below uses MRO
class Directory (entry.Entry) :

    # ----------------------------------------------------------------
    #
    # initialization methods
    #
    def __init__ (self, api, adaptor) :

        self._cpi_nsentry = super  (Directory, self)
        self._cpi_nsentry.__init__ (api, adaptor)

    @SYNC
    def init_instance           (self, url, flags, session)    : pass
    @ASYNC
    def init_instance_async     (self, url, flags, session)    : pass 


    # ----------------------------------------------------------------
    #
    # namespace directory methods
    #
    @SYNC
    def change_dir              (self, url, flags, ttype)      : pass
    @ASYNC
    def change_dir_async        (self, url, flags, ttype)      : pass

    @SYNC
    def list                    (self, npat, ttype)            : pass
    @ASYNC
    def list_async              (self, npat, ttype)            : pass

    @SYNC
    def find                    (self, npat, flags, ttype)     : pass
    @ASYNC
    def find_async              (self, npat, flags, ttype)     : pass

    @SYNC
    def exists                  (self, name, ttype)            : pass
    @ASYNC
    def exists_async            (self, name, ttype)            : pass


    @SYNC
    def is_dir                  (self, name, ttype)            : pass
    @ASYNC
    def is_dir_async            (self, name, ttype)            : pass

    @SYNC
    def is_entry                (self, name, ttype)            : pass
    @ASYNC
    def is_entry_async          (self, name, ttype)            : pass

    @SYNC
    def is_link                 (self, name, ttype)            : pass
    @ASYNC
    def is_link_async           (self, name, ttype)            : pass

    @SYNC
    def read_link               (self, name, ttype)            : pass
    @ASYNC
    def read_link_async         (self, name, ttype)            : pass

    @SYNC
    def get_num_entries         (self, ttype)                  : pass
    @ASYNC
    def get_num_entries_async   (self, ttype)                  : pass

    @SYNC
    def get_entry               (self, num, ttype)             : pass
    @ASYNC
    def get_entry_async         (self, num, ttype)             : pass

    @SYNC
    def copy                    (self, src, tgt, flags, ttype) : pass
    @ASYNC
    def copy_async              (self, src, tgt, flags, ttype) : pass

    @SYNC
    def link                    (self, src, tgt, flags, ttype) : pass
    @ASYNC
    def link_async              (self, src, tgt, flags, ttype) : pass

    @SYNC
    def move                    (self, src, tgt, flags, ttype) : pass
    @ASYNC
    def move_async              (self, src, tgt, flags, ttype) : pass

    @SYNC
    def remove                  (self, tgt, flags, ttype)      : pass
    @ASYNC
    def remove_async            (self, tgt, flags, ttype)      : pass

    @SYNC
    def make_dir                (self, tgt, flags, ttype)      : pass
    @ASYNC
    def make_dir_async          (self, tgt, flags, ttype)      : pass




