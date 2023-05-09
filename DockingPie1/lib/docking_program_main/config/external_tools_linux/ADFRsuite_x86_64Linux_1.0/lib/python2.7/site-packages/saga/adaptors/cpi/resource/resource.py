
__author__    = "Andre Merzky"
__copyright__ = "Copyright 2012-2013, The SAGA Project"
__license__   = "MIT"


""" Provides the SAGA Resource CPI """

import saga.adaptors.cpi.decorators as cpi_dec
import saga.adaptors.cpi.base       as cpi_base
import saga.adaptors.cpi.async      as cpi_async

SYNC  = cpi_dec.CPI_SYNC_CALL
ASYNC = cpi_dec.CPI_ASYNC_CALL


class Resource (cpi_base.CPIBase, cpi_async.Async) :

    # ----------------------------------------------------------------
    #
    # initialization methods
    #
    def __init__ (self, api, adaptor) :

        self._cpi_base = super  (Resource, self)
        self._cpi_base.__init__ (api, adaptor)

    @SYNC
    def init_instance             (self, info,    ttype) : pass
    @ASYNC                                        
    def init_instance_async       (self, info,    ttype) : pass


    # ----------------------------------------------------------------
    #
    # resource methods
    #
    @SYNC
    def reconfig           (self, descr,          ttype) : pass
    @ASYNC                                       
    def reconfig_async     (self, descr,          ttype) : pass
                           
    @SYNC                  
    def destroy            (self,                 ttype) : pass
    @ASYNC                                                     
    def destroy_async      (self,                 ttype) : pass
                           
    @SYNC                  
    def wait               (self, state, timeout, ttype) : pass
    @ASYNC                 
    def wait_async         (self, state, timeout, ttype) : pass


    #-----------------------------------------------------------------
    # attribute getters
    @SYNC
    def get_id                 (self, ttype)             : pass
    @ASYNC                                     
    def get_id_async           (self, ttype)             : pass
                                               
    @SYNC                                      
    def get_rtype              (self, ttype)             : pass
    @ASYNC                                     
    def get_rtype_async        (self, ttype)             : pass
                                               
    @SYNC                                      
    def get_state              (self, ttype)             : pass
    @ASYNC                                     
    def get_state_async        (self, ttype)             : pass
                                               
    @SYNC                                      
    def get_state_detail       (self, ttype)             : pass
    @ASYNC                                     
    def get_state_detail_async (self, ttype)             : pass
                                               
    @SYNC                                      
    def get_access             (self, ttype)             : pass
    @ASYNC                                     
    def get_access_async       (self, ttype)             : pass
                                               
    @SYNC                                      
    def get_manager            (self, ttype)             : pass
    @ASYNC                                     
    def get_manager_async      (self, ttype)             : pass
                                               
    @SYNC                                      
    def get_description        (self, ttype)             : pass
    @ASYNC                                     
    def get_description_async  (self, ttype)             : pass


class Compute (Resource) : pass
class Storage (Resource) : pass
class Network (Resource) : pass



