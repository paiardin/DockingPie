
__author__    = "Andre Merzky"
__copyright__ = "Copyright 2012-2013, The SAGA Project"
__license__   = "MIT"


""" Provides the SAGA Job CPI """

import saga.adaptors.cpi.decorators as cpi_dec
import saga.adaptors.cpi.base       as cpi_base
import saga.adaptors.cpi.async      as cpi_async

SYNC  = cpi_dec.CPI_SYNC_CALL
ASYNC = cpi_dec.CPI_ASYNC_CALL


class Job (cpi_base.CPIBase, cpi_async.Async) :

    # ----------------------------------------------------------------
    #
    # initialization methods
    #
    def __init__ (self, api, adaptor) :

        _cpi_base = super  (Job, self)
        _cpi_base.__init__ (api, adaptor)

    @SYNC
    def init_instance             (self, info, ttype)    : pass
    @ASYNC
    def init_instance_async       (self, info, ttype)    : pass


    # ----------------------------------------------------------------
    #
    # job methods
    #
    @SYNC
    def get_id                    (self, ttype)          : pass
    @ASYNC
    def get_id_async              (self, ttype)          : pass

    @SYNC
    def get_name                  (self, ttype)          : pass
    @ASYNC
    def get_name_async            (self, ttype)          : pass

    @SYNC
    def get_description           (self, ttype)          : pass
    @ASYNC
    def get_description_async     (self, ttype)          : pass

    @SYNC
    def get_stdin                 (self, ttype)          : pass
    @ASYNC
    def get_stdin_async           (self, ttype)          : pass

    @SYNC
    def get_stdout                (self, ttype)          : pass
    @ASYNC
    def get_stdout_async          (self, ttype)          : pass

    @SYNC
    def get_stderr                (self, ttype)          : pass
    @ASYNC
    def get_stderr_async          (self, ttype)          : pass

    @SYNC
    def get_log                   (self, ttype)          : pass
    @ASYNC
    def get_log_async             (self, ttype)          : pass

    @SYNC
    def suspend                   (self, ttype)          : pass
    @ASYNC
    def suspend_async             (self, ttype)          : pass

    @SYNC
    def resume                    (self, ttype)          : pass
    @ASYNC
    def resume_async              (self, ttype)          : pass

    @SYNC
    def checkpoint                (self, ttype)          : pass
    @ASYNC
    def checkpoint_async          (self, ttype)          : pass

    @SYNC
    def migrate                   (self, jd, ttype)      : pass
    @ASYNC
    def migrate_async             (self, jd, ttype)      : pass

    @SYNC
    def signal                    (self, signum, ttype)  : pass
    @ASYNC
    def signal_async              (self, signum, ttype)  : pass


    #-----------------------------------------------------------------
    # task methods flattened into job :-/
    @SYNC
    def run                       (self, ttype)          : pass
    @ASYNC
    def run_async                 (self, ttype)          : pass

    @SYNC
    def cancel                    (self, timeout, ttype) : pass
    @ASYNC
    def cancel_async              (self, timeout, ttype) : pass

    @SYNC
    def wait                      (self, timeout, ttype) : pass
    @ASYNC
    def wait_async                (self, timeout, ttype) : pass

    @SYNC
    def get_state                 (self, ttype)          : pass
    @ASYNC
    def get_state_async           (self, ttype)          : pass

    @SYNC
    def get_result                (self, ttype)          : pass
    @ASYNC
    def get_result_async          (self, ttype)          : pass

    @SYNC
    def get_object                (self, ttype)          : pass
    @ASYNC
    def get_object_async          (self, ttype)          : pass

    @SYNC
    def re_raise                  (self, ttype)          : pass
    @ASYNC
    def re_raise_async            (self, ttype)          : pass


    #-----------------------------------------------------------------
    # attribute getters
    @SYNC
    def get_exit_code             (self, ttype)          : pass
    @ASYNC
    def get_exit_code_async       (self, ttype)          : pass

    @SYNC
    def get_created               (self, ttype)          : pass
    @ASYNC
    def get_created_async         (self, ttype)          : pass

    @SYNC
    def get_started               (self, ttype)          : pass
    @ASYNC
    def get_started_async         (self, ttype)          : pass

    @SYNC
    def get_finished              (self, ttype)          : pass
    @ASYNC
    def get_finished_async        (self, ttype)          : pass

    @SYNC
    def get_execution_hosts       (self, ttype)          : pass
    @ASYNC
    def get_execution_hosts_async (self, ttype)          : pass


# class Self (Job, monitoring.Steerable) :
class Self (Job) :

    @SYNC
    def init_instance             (self)                 : pass
    @ASYNC
    def init_instance_async       (self)                 : pass




