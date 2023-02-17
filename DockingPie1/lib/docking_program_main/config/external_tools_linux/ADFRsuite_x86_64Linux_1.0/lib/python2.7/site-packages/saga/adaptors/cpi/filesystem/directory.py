
__author__    = "Andre Merzky"
__copyright__ = "Copyright 2012-2013, The SAGA Project"
__license__   = "MIT"


import saga.adaptors.cpi.decorators as cpi_dec
import saga.adaptors.cpi.namespace  as cpi_ns

SYNC  = cpi_dec.CPI_SYNC_CALL
ASYNC = cpi_dec.CPI_ASYNC_CALL


# keep order of inheritance!  super() below uses MRO
class Directory (cpi_ns.directory.Directory) :

    # ----------------------------------------------------------------
    #
    # initialization methods
    #
    def __init__ (self, api, adaptor) :

        self._cpi_nsdirec = super  (Directory, self)
        self._cpi_nsdirec.__init__ (api, adaptor)

    @SYNC
    def init_instance         (self, url, flags, session) : pass
    @ASYNC
    def init_instance_async   (self, url, flags, session) : pass


    # ----------------------------------------------------------------
    #
    # add filesystem directory methods
    #
    @SYNC
    def get_size              (self, name, flags, ttype)  : pass
    @ASYNC
    def get_size_async        (self, name, flags, ttype)  : pass

    @SYNC
    def get_size_self         (self,       flags, ttype)  : pass
    @ASYNC                                
    def get_size__self_async  (self,       flags, ttype)  : pass

    @SYNC
    def is_file               (self, name,        ttype)  : pass
    @ASYNC
    def is_file_async         (self, name,        ttype)  : pass

    @SYNC
    def is_file_self          (self,              ttype)  : pass
    @ASYNC                                        
    def is_file_self_async    (self,              ttype)  : pass




