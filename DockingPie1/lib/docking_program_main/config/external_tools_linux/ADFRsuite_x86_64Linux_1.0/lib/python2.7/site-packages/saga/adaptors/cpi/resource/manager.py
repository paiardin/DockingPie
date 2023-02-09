
__author__    = "Andre Merzky"
__copyright__ = "Copyright 2012-2013, The SAGA Project"
__license__   = "MIT"


""" Provides the SAGA Resource CPI """

import saga.adaptors.cpi.decorators as cpi_dec
import saga.adaptors.cpi.base       as cpi_base
import saga.adaptors.cpi.async      as cpi_async

SYNC  = cpi_dec.CPI_SYNC_CALL
ASYNC = cpi_dec.CPI_ASYNC_CALL


class Manager (cpi_base.CPIBase, cpi_async.Async) :

    # ----------------------------------------------------------------
    #
    # initialization methods
    #
    def __init__ (self, api, adaptor) :

        self._cpi_base = super  (Manager, self)
        self._cpi_base.__init__ (api, adaptor)

    @SYNC
    def init_instance              (self, url, session)        : pass
    @ASYNC
    def init_instance_async        (self, url, session)        : pass



    # ----------------------------------------------------------------
    #
    # resource manager methods
    #
    @SYNC
    def list                       (self, type,      ttype)    : pass
    @ASYNC
    def list_async                 (self, type,      ttype)    : pass

    @SYNC
    def get_description            (self, id,        ttype)    : pass
    @ASYNC
    def get_description_async      (self, id,        ttype)    : pass

    @SYNC
    def list_templates             (self, type,      ttype)    : pass
    @ASYNC
    def list_templates_async       (self, type,      ttype)    : pass

    @SYNC
    def get_template               (self, name,      ttype)    : pass
    @ASYNC
    def get_template_async         (self, name,      ttype)    : pass

    @SYNC
    def list_images                (self, type,      ttype)    : pass
    @ASYNC
    def list_images_async          (self, type,      ttype)    : pass

    @SYNC
    def get_image                  (self, name,      ttype)    : pass
    @ASYNC
    def get_image_async            (self, name,      ttype)    : pass

    @SYNC
    def acquire                    (self, descr,     ttype)    : pass
    @ASYNC
    def acquire_async              (self, descr,     ttype)    : pass

    @SYNC
    def destroy                    (self, id,        ttype)    : pass
    @ASYNC                                           
    def destroy_async              (self, id,        ttype)    : pass





