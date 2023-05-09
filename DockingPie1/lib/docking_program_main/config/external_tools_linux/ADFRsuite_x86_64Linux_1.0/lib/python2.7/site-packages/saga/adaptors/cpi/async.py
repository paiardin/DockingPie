
__author__    = "Andre Merzky"
__copyright__ = "Copyright 2012-2013, The SAGA Project"
__license__   = "MIT"


import saga.adaptors.cpi.decorators as cpi_dec

SYNC  = cpi_dec.CPI_SYNC_CALL
ASYNC = cpi_dec.CPI_ASYNC_CALL


class Async (object) :
    
    @SYNC
    def task_run    (self, task)          : pass

    @SYNC
    def task_wait   (self, task, timeout) : pass

    @SYNC
    def task_cancel (self, task)          : pass




