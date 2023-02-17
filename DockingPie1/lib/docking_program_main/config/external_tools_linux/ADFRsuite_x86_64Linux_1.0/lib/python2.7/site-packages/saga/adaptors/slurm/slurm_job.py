
__author__    = "Andre Merzky, Ashley Z, Ole Weidner"
__copyright__ = "Copyright 2012-2013, The SAGA Project"
__license__   = "MIT"


""" SLURM job adaptor implementation """

#TODO: Throw errors if a user does not specify the MINIMUM number of
#      attributes required for SLURM in a job description

import saga.utils.pty_shell

import saga.url as surl
import saga.adaptors.base
import saga.adaptors.cpi.job

import re
import os
import time
import textwrap
import string
import tempfile

SYNC_CALL  = saga.adaptors.cpi.decorators.SYNC_CALL
ASYNC_CALL = saga.adaptors.cpi.decorators.ASYNC_CALL

# --------------------------------------------------------------------
#
def log_error_and_raise(message, exception, logger):
    logger.error(message)
    raise exception(message)

# --------------------------------------------------------------------
# some private defs
#
_PTY_TIMEOUT = 2.0

# --------------------------------------------------------------------
# the adaptor name
#
_ADAPTOR_NAME          = "saga.adaptor.slurm_job"
_ADAPTOR_SCHEMAS       = ["slurm", "slurm+ssh", "slurm+gsissh"]
_ADAPTOR_OPTIONS       = []

# --------------------------------------------------------------------
# the adaptor capabilities & supported attributes
#
# TODO: FILL ALL IN FOR SLURM
_ADAPTOR_CAPABILITIES  = {
    "jdes_attributes"  : [saga.job.NAME, 
                          saga.job.EXECUTABLE,
                          saga.job.PRE_EXEC,
                          saga.job.POST_EXEC,
                          saga.job.ARGUMENTS,
                          saga.job.ENVIRONMENT,
                          saga.job.SPMD_VARIATION, #implement later, somehow
                          saga.job.TOTAL_CPU_COUNT, 
                          saga.job.NUMBER_OF_PROCESSES,
                          saga.job.PROCESSES_PER_HOST,
                          saga.job.THREADS_PER_PROCESS, 
                          saga.job.WORKING_DIRECTORY,
                          #saga.job.INTERACTIVE,
                          saga.job.INPUT,
                          saga.job.OUTPUT, 
                          saga.job.ERROR,
                          saga.job.FILE_TRANSFER,
                          saga.job.CLEANUP,
                          saga.job.JOB_START_TIME,
                          saga.job.WALL_TIME_LIMIT, 
                          saga.job.TOTAL_PHYSICAL_MEMORY, 
                          #saga.job.CPU_ARCHITECTURE, 
                          #saga.job.OPERATING_SYSTEM_TYPE, 
                          saga.job.CANDIDATE_HOSTS,
                          saga.job.QUEUE,
                          saga.job.PROJECT,
                          saga.job.JOB_CONTACT],

    "job_attributes"   : [saga.job.EXIT_CODE,
                          saga.job.EXECUTION_HOSTS,
                          saga.job.CREATED,
                          saga.job.STARTED,
                          saga.job.FINISHED],
    "metrics"          : [saga.job.STATE, 
                          saga.job.STATE_DETAIL],
    "contexts"         : {"ssh"      : "public/private keypair",
                          "x509"     : "X509 proxy for gsissh",
                          "userpass" : "username/password pair for simple ssh"}
}

# --------------------------------------------------------------------
# the adaptor documentation
#

        # General Notes
        # *************

        # On Stampede, returning a non-zero exit code results in the scheduler
        # putting the job into a FAILED state and assigning it an exit code of 127.

        # **Example:**

        # ::

        #   js = saga.job.Service("slurm+ssh://stampede")
        #   jd.executable  = '/bin/exit'
        #   jd.arguments   = ['3']
        #   job = js.create_job(jd)
        #   job.run()

        # Will return something similar to (personal account information
        # removed)::

        #   (saga-python-env)ashleyz@login1:~$ scontrol show job 309684
        #   JobId=309684 Name=SlurmJob
        #      UserId=_____ GroupId=__________
        #      Priority=3000 Account=_____ QOS=normal
        #      JobState=FAILED Reason=NonZeroExitCode Dependency=(null)
        #      Requeue=0 Restarts=0 BatchFlag=1 ExitCode=127:0
        #      RunTime=00:00:05 TimeLimit=00:01:00 TimeMin=N/A
        #      SubmitTime=2013-02-22T20:26:50 EligibleTime=2013-02-22T20:26:50
        #      StartTime=2013-02-22T20:26:50 EndTime=2013-02-22T20:26:55
        #      PreemptTime=None SuspendTime=None SecsPreSuspend=0
        #      Partition=development AllocNode:Sid=login1:12070
        #      ReqNodeList=(null) ExcNodeList=(null)
        #      NodeList=c557-401
        #      BatchHost=c557-401
        #      NumNodes=1 NumCPUs=16 CPUs/Task=1 ReqS:C:T=*:*:*
        #      MinCPUsNode=1 MinMemoryNode=0 MinTmpDiskNode=0
        #      Features=(null) Gres=(null) Reservation=(null)
        #      Shared=0 Contiguous=0 Licenses=(null) Network=(null)
        #      Command=/home1/01414/_______/.saga/adaptors/slurm_job/wrapper.sh
        #      WorkDir=/home1/01414/_______/

        # I'm not sure how to fix this for the time being.

        # Suspend/resume do not appear to be supported for regular
        # users on Stampede.

        # run_job is not supported, as there are many attributed (queues,
        # projects, etc) which need to be passed to the adaptor.  I could
        # add URL parsing so that you could pile on queue/project/jobname
        # information if this has any strong usecases, but I avoided doing
        # so for now to decrease complexity/possible confusion.

        # Cancelling a job with scontrol, puts it into a COMPLETING state, which
        # is parsed by the SLURM status parser as saga.job.RUNNING (see the
        # SLURM docs, COMPLETING is a state a job goes into when it is done
        # running but still flushing IO/etc).  Anyhow, I put some code in to
        # manually put the job into CANCELED state when the job is canceled,
        # but I'm not sure that this is reported correctly everywhere yet.

        # What exit code should be returned for a CANCELED job?
        

_ADAPTOR_DOC           = {
    "name"             : _ADAPTOR_NAME,
    "cfg_options"      : _ADAPTOR_OPTIONS, 
    "capabilities"     : _ADAPTOR_CAPABILITIES,
    "description"      : """ 
        The SLURM adaptor allows to run and manage jobs on a 
        `SLURM <https://computing.llnl.gov/linux/slurm/slurm.html>`_ HPC cluster.

        Implementation Notes
        ********************

         - If scontrol can't find an exit code, it returns None
           (see _job_get_exit_code)
         - If scancel can't cancel a job, we raise an exception
           (see _job_cancel)
         - If we can't suspend a job with scontrol suspend, we raise an exception
           (see _job_suspend).  scontrol suspend NOT supported on Stampede
         - I started to implement a dictionary to keep track of jobs locally.
           It works to the point where the unit tests are passed, but I have
           not gone over theis extensively...
         - Relating to the above, _job_get_info is written, but unused/untested
           (mostly from PBS adaptor)

        """,
    "example": "examples/jobs/slurmjob.py",
    "schemas": {"slurm":        "connect to a local cluster",
                "slurm+ssh":    "conenct to a remote cluster via SSH",
                "slurm+gsissh": "connect to a remote cluster via GSISSH"}
}

# --------------------------------------------------------------------
# the adaptor info is used to register the adaptor with SAGA

_ADAPTOR_INFO          = {
    "name"             : _ADAPTOR_NAME,
    "version"          : "v0.2",
    "schemas"          : _ADAPTOR_SCHEMAS,
    "capabilities"     : _ADAPTOR_CAPABILITIES,
    "cpis"             : [
        { 
        "type"         : "saga.job.Service",
        "class"        : "SLURMJobService"
        }, 
        { 
        "type"         : "saga.job.Job",
        "class"        : "SLURMJob"
        }
    ]
}

###############################################################################
# The adaptor class

class Adaptor (saga.adaptors.base.Base):
    """ 
    This is the actual adaptor class, which gets loaded by SAGA (i.e. by the
    SAGA engine), and which registers the CPI implementation classes which
    provide the adaptor's functionality.
    """


    # ----------------------------------------------------------------
    #
    def __init__ (self) :

        saga.adaptors.base.Base.__init__ (self, _ADAPTOR_INFO, _ADAPTOR_OPTIONS)

        self.id_re = re.compile ('^\[(.*)\]-\[(.*?)\]$')

    # ----------------------------------------------------------------
    #
    def sanity_check (self) :
        pass


    def parse_id (self, id) :
        # split the id '[rm]-[pid]' in its parts, and return them.

        match = self.id_re.match (id)

        if  not match or len (match.groups()) != 2 :
            raise saga.BadParameter ("Cannot parse job id '%s'" % id)

        return (match.group(1), match.group (2))


###############################################################################
#
class SLURMJobService (saga.adaptors.cpi.job.Service) :
    """ Implements saga.adaptors.cpi.job.Service """

    # ----------------------------------------------------------------
    #
    def __init__ (self, api, adaptor) :
        
        _cpi_base = super  (SLURMJobService, self)
        _cpi_base.__init__ (api, adaptor)

        # TODO make sure this formats properly and works right!
        self.exit_code_re            = re.compile(r"\bExitCode  \b=(\d*)", re.VERBOSE)
        self.scontrol_jobstate_re    = re.compile(r"\bJobState  \b=(\S*)", re.VERBOSE)
        self.scontrol_create_time_re = re.compile(r"\bSubmitTime\b=(\S*)", re.VERBOSE)
        self.scontrol_start_time_re  = re.compile(r"\bStartTime \b=(\S*)", re.VERBOSE)
        self.scontrol_end_time_re    = re.compile(r"\bEndTime   \b=(\S*)", re.VERBOSE)
        self.scontrol_comp_time_re   = re.compile(r"\bRunTime   \b=(\S*)", re.VERBOSE)
        self.scontrol_exec_hosts_re  = re.compile(r"\bNodeList  \b=(\S*)", re.VERBOSE)

        # these are the commands that we need in order to interact with SLURM
        # the adaptor will try to find them when it first opens the shell
        # connection, and bails out in case they are not available.
        self._commands = {'sbatch': None,
                          'squeue': None,
                          'scontrol': None,
                          'scancel': None}

    # ----------------------------------------------------------------
    #
    def __del__ (self) :
        try :
            if self.shell : del (self.shell)
        except :
            pass


    # ----------------------------------------------------------------
    #
    @SYNC_CALL
    def init_instance (self, adaptor_state, rm_url, session) :
        """ Service instance constructor """

        self.rm      = rm_url
        self.session = session

        self.jobs = {}
        self._open ()

        return self.get_api ()


    # ----------------------------------------------------------------
    #
    def close (self) :
        if  self.shell :
            self.shell.finalize (True)


    # # ----------------------------------------------------------------
    # #
    # def _alive (self) :

    #     if  not self.shell or not self.shell.alive () :
    #         self._logger.info ("shell is dead - long live the shell")
            
    #         try :
    #             self._close ()  # for cleanup...
    #             self._open  ()

    #         except Exception :
    #             # did not work for good - give up
    #             raise saga.IncorrectState ("job service is not connected, can't reconnect")


    # ----------------------------------------------------------------
    #
    def _open (self) :
        """
        Open our persistent shell for this job adaptor.  We use
        the pty_shell functionality for this.  
        """
        # check to see what kind of connection we will want to create
        if self.rm.schema   == "slurm":
            shell_schema = "fork://"
        elif self.rm.schema == "slurm+ssh":
            shell_schema = "ssh://"
        elif self.rm.schema == "slurm+gsissh":
            shell_schema = "gsissh://"
        else:
            raise saga.IncorrectURL("Schema %s not supported by SLURM adaptor."
                                    % self.rm.schema)

        #<scheme>://<user>:<pass>@<host>:<port>/<path>?<query>#<fragment>
        # build our shell URL
        shell_url = shell_schema 
        
        # did we provide a username and password?
        if self.rm.username and self.rm.password:
            shell_url += self.rm.username + ":" + self.rm.password + "@"

        # only provided a username
        if self.rm.username and not self.rm.password:
            shell_url += self.rm.username + "@"

        #add hostname
        shell_url += self.rm.host

        #add port
        if  self.rm.port:
            shell_url += ":" + str(self.rm.port)

        shell_url = saga.url.Url(shell_url)

        # establish shell connection
        self._logger.debug("Opening shell of type: %s" % shell_url)
        self.shell = saga.utils.pty_shell.PTYShell (shell_url, 
                                                    self.session, 
                                                    self._logger)

        # verify our SLURM environment contains the commands we need for this
        # adaptor to work properly
        self._logger.debug("Verifying existence of remote SLURM tools.")
        for cmd in self._commands.keys():
            ret, out, _ = self.shell.run_sync("which %s " % cmd)
            if ret != 0:
                message = "Error finding SLURM tool %s on remote server %s!\n" \
                          "Locations searched:\n%s\n" \
                          "Is SLURM installed on that machine? " \
                          "If so, is your remote SLURM environment "\
                          "configured properly? " % (cmd, self.rm, out)
                raise saga.NoSuccess._log (self._logger, message)
                
        self._logger.debug ("got cmd prompt (%s)(%s)" % (ret, out))
        
        self.rm.detected_username = self.rm.username
        # figure out username if it wasn't made explicit
        # important if .ssh/config info read+connected with 
        # a different username than what we expect
        if not self.rm.username:
            self._logger.debug ("No username provided in URL %s, so we are"
                                " going to find it with whoami" % self.rm)
            ret, out, _ = self.shell.run_sync("whoami")
            self.rm.detected_username = out.strip()
            self._logger.debug("Username detected as: %s",
                               self.rm.detected_username)

        return

    # ----------------------------------------------------------------
    #
    def _close (self) :
        """
        Close our shell connection
        """
        del (self.shell)
        self.shell = None


    # ----------------------------------------------------------------
    #
    #
    def _job_run (self, jd) :
        """ runs a job on the wrapper via pty, and returns the job id """
        
        #define a bunch of default args
        exe = jd.executable
        pre = jd.pre_exec
        post= jd.post_exec
        arg = ""
        env = ""
        cwd = ""
        job_name = jd.name
        if not job_name:
            job_name = "SAGAPythonSLURMJob"
        spmd_variation = None
        total_cpu_count = None
        number_of_processes = None
        processes_per_host = None
        output = "saga-python-slurm-default.out"
        error = None
        file_transfer = None
        job_start_time = None
        wall_time_limit = None
        queue = None
        project = None
        job_memory = None
        job_contact = None
        candidate_hosts = None

        # check to see what's available in our job description
        # to override defaults

        if jd.attribute_exists ("name"):
            job_name = jd.name

        if jd.attribute_exists ("arguments") :
            for a in jd.arguments :
                arg += " %s" % a

        if jd.attribute_exists ("spmd_variation"):
            spmd_variation = jd.spmd_variation

        if jd.attribute_exists ("total_cpu_count"):
            total_cpu_count = jd.total_cpu_count

        if jd.attribute_exists ("number_of_processes"):
            number_of_processes = jd.number_of_processes

        if jd.attribute_exists("processes_per_host"):
            processes_per_host = jd.processes_per_host

        if jd.attribute_exists ("working_directory"):
            cwd = jd.working_directory

        if jd.attribute_exists ("output"):
            output = jd.output

        if jd.attribute_exists("error"):
            error = jd.error

        if jd.attribute_exists("wall_time_limit"):
            wall_time_limit = jd.wall_time_limit

        if jd.attribute_exists("queue"):
            queue = jd.queue

        if jd.attribute_exists("project"):
            project = jd.project

        if jd.attribute_exists("total_physical_memory"):
            job_memory = jd.total_physical_memory

        if jd.attribute_exists("candidate_hosts"):
            if isinstance(jd.candidate_hosts, list):
                candidate_hosts = ','.join(jd.candidate_hosts)
            else:
                candidate_hosts = jd.candidate_hosts

        if jd.attribute_exists("job_contact"):
            job_contact = jd.job_contact[0]

        slurm_script = "#!/bin/sh\n\n"

        if  job_name:
            slurm_script += '#SBATCH -J "%s"\n' % job_name

        if spmd_variation:
            pass #TODO

        #### HANDLE NUMBER OF CORES
        # make sure we have something for total_cpu_count
        if not total_cpu_count:
            self._logger.warning("total_cpu_count not specified in submitted "
                                 "SLURM job description -- defaulting to 1!")
            total_cpu_count = 1

        # make sure we have something for number_of_processes
        if not number_of_processes:
            self._logger.debug("number_of_processes not specified in submitted "
                                 "SLURM job description -- defaulting to 1 per total_cpu_count! (%s)" % total_cpu_count)
            number_of_processes = total_cpu_count

        # make sure we aren't given more processes than CPUs
        if number_of_processes > total_cpu_count:
            log_error_and_raise("More processes (%s) requested than total number of CPUs! (%s)" % (number_of_processes, total_cpu_count), saga.NoSuccess, self._logger)

        #make sure we aren't doing funky math
        if  total_cpu_count % number_of_processes != 0:
            log_error_and_raise ("total_cpu_count (%s) must be evenly " \
                                 "divisible by number_of_processes (%s)" \
                              % (total_cpu_count, number_of_processes), 
                                 saga.NoSuccess, self._logger)

        slurm_script += "#SBATCH --ntasks=%s\n" % (number_of_processes)

        if total_cpu_count != number_of_processes:
            slurm_script += "#SBATCH --cpus-per-task=%s\n" % (total_cpu_count / number_of_processes)

        if processes_per_host:
            slurm_script += "#SBATCH --ntasks-per-node=%s\n" % processes_per_host

        if  cwd is not "":
            slurm_script += "#SBATCH --workdir %s\n" % cwd

        if  output:
            slurm_script += "#SBATCH --output %s\n" % output
        
        if  error:
            slurm_script += "#SBATCH --error %s\n" % error

        if  wall_time_limit:
            hours   = wall_time_limit / 60
            minutes = wall_time_limit % 60
            slurm_script += "#SBATCH --time %02d:%02d:00\n" % (hours, minutes)

        if  queue:
            slurm_script += "#SBATCH --partition %s\n" % queue

        if  project:
            if not ':' in project:
                account = project
            else:
                account, reservation = project.split(':')
                slurm_script += "#SBATCH --reservation %s\n" % reservation

            slurm_script += "#SBATCH --account %s\n" % account

        if  job_memory:
            slurm_script += "#SBATCH --mem=%s\n" % job_memory

        if  candidate_hosts:
            slurm_script += "#SBATCH --nodelist=%s\n" % candidate_hosts

        if  job_contact:
            slurm_script += "#SBATCH --mail-user=%s\n" % job_contact

        # make sure we are not missing anything important
        if  not queue:
            raise saga.BadParameter._log (self._logger, 
                                          "No queue has been specified, "
                                          "and the SLURM adaptor "
                                          "requires that a queue be "
                                          "specified.  Please specify "
                                          "a queue to submit the job to.")

        # add on our environment variables
        if jd.attribute_exists ("environment") :
            slurm_script += "\n## ENVIRONMENT\n"
            for key,val in jd.environment.iteritems() :
                slurm_script += "export %s=%s\n"  %  (key, val)

        if pre :
            slurm_script += "\n## PRE_EXEC\n"
            for cmd in pre :
                slurm_script += "%s\n" % cmd

        # create our commandline
        slurm_script += "\n## EXEC\n"
        slurm_script += exe
        if jd.attribute_exists ("arguments") :
            slurm_script += ' ' + ' '.join (jd.arguments)
        slurm_script += '\n'

        if post :
            slurm_script += "\n## POST_EXEC\n"
            for cmd in post :
                slurm_script += "%s\n" % cmd

        # try to create the working directory (if defined)
        # WRANING: this assumes a shared filesystem between login node and
        #           comnpute nodes.
        if jd.working_directory:
            self._logger.info("Creating working directory %s" % jd.working_directory)
            ret, out, _ = self.shell.run_sync("mkdir -p %s" % (jd.working_directory))
            if ret != 0:
                # something went wrong
                message = "Couldn't create working directory - %s" % (out)
                log_error_and_raise(message, saga.NoSuccess, self._logger)


        # write script into a tmp file for staging
        self._logger.info ("SLURM script generated:\n%s" % slurm_script)

        tgt = os.path.basename (tempfile.mktemp (suffix='.slurm', prefix='tmp_'))
        self.shell.write_to_remote (src=slurm_script, tgt=tgt)

        # submit the job
        ret, out, _ = self.shell.run_sync ("sbatch '%s'; rm -vf '%s'" % (tgt, tgt))

        self._logger.debug ("staged/submit SLURM script (%s) (%s)" % (tgt, ret))

        # find out what our job ID will be
        # TODO: Could make this more efficient
        found_id = False
        for line in out.split("\n"):
            if "Submitted batch job" in line:
                self.job_id = "[%s]-[%s]" % \
                    (self.rm, int(line.split()[-1:][0]))
                found_id = True

        # if we have no job ID, there's a failure...
        if not found_id:
            raise saga.NoSuccess._log(self._logger, 
                             "Couldn't get job id from submitted job!"
                              " sbatch output:\n%s" % out)

        self._logger.debug("started job %s" % self.job_id)
        self._logger.debug("Batch system output:\n%s" % out)

        # create local jobs dictionary entry
        self.jobs[self.job_id] = {
                'state': saga.job.PENDING,
                'create_time': None,
                'start_time': None,
                'end_time': None,
                'comp_time': None,
                'exec_hosts': None,
                'gone': False
            }

        return self.job_id

    # ----------------  
    # FROM STAMPEDE'S SQUEUE MAN PAGE
    # 
    # JOB STATE CODES
    #    Jobs typically pass through several states in the course of their execution.  The typical states are PENDING, RUNNING, SUSPENDED, COMPLETING, and COMPLETED.   An  explanation  of  each
    #    state follows.

    #    CA  CANCELED        Job was explicitly cancelled by the user or system administrator.  The job may or may not have been initiated.
    #    CD  COMPLETED       Job has terminated all processes on all nodes.
    #    CF  CONFIGURING     Job has been allocated resources, but are waiting for them to become ready for use (e.g. booting).
    #    CG  COMPLETING      Job is in the process of completing. Some processes on some nodes may still be active.
    #    F   FAILED          Job terminated with non-zero exit code or other failure condition.
    #    NF  NODE_FAIL       Job terminated due to failure of one or more allocated nodes.
    #    PD  PENDING         Job is awaiting resource allocation.
    #    PR  PREEMPTED       Job terminated due to preemption.
    #    R   RUNNING         Job currently has an allocation.
    #    S   SUSPENDED       Job has an allocation, but execution has been suspended.
    #    TO  TIMEOUT         Job terminated upon reaching its time limit.
 
    def _slurm_to_saga_jobstate(self, slurmjs):
        """ translates a slurm one-letter state to saga
        """
        if slurmjs == "CANCELLED" or slurmjs == 'CA':
            return saga.job.CANCELED
        elif slurmjs == "COMPLETED" or slurmjs == 'CD':
            return saga.job.DONE
        elif slurmjs == "CONFIGURING" or slurmjs == 'CF':
            return saga.job.PENDING
        elif slurmjs == "COMPLETING" or slurmjs == 'CG':
            return saga.job.RUNNING
        elif slurmjs == "FAILED" or slurmjs == 'F':
            return saga.job.FAILED
        elif slurmjs == "NODE_FAIL" or slurmjs == 'NF':
            return saga.job.FAILED
        elif slurmjs == "PENDING" or slurmjs == 'PD':
            return saga.job.PENDING
        elif slurmjs == "PREEMPTED" or slurmjs == 'PR':
            return saga.job.CANCELED
        elif slurmjs == "RUNNING" or slurmjs == 'R': 
            return saga.job.RUNNING
        elif slurmjs == "SUSPENDED" or slurmjs == 'S':
            return saga.job.SUSPENDED
        elif slurmjs == "TIMEOUT" or slurmjs == 'TO':
            return saga.job.CANCELED
        else:
            return saga.job.UNKNOWN

    def _job_get_exit_code (self, id) :
        """ get the job exit code from the wrapper shell """
        rm, pid     = self._adaptor.parse_id (id)
        ret, out, _ = self.shell.run_sync("scontrol show job %s" % pid)
        match       = self.exit_code_re.search (out)

        if match:
            self.exit_code = int(match.group(1))

        else:
            self.exit_code = None

        self._logger.debug("Returning exit code %s" % self.exit_code)
        return self.exit_code

        
        ### couldn't get the exitcode -- maybe should change this to just return
        ### None?  b/c we will lose the code if a program waits too
        ### long to look for the exitcode...
        ###raise saga.NoSuccess._log (self._logger, 
        ###                           "Could not find exit code for job %s" % id)
        # decided to have it just return none
        self._logger.warning("Couldn't find job %s in SLURM -- returning None"
                             " for job status." % id)
        return None

    def _job_cancel (self, id):
        """
        Given a job id, attempt to cancel it through use of commandline
        scancel.  Raises exception when unsuccessful.
        """
        rm, pid     = self._adaptor.parse_id (id)
        ret, out, _ = self.shell.run_sync("scancel %s" % pid)
        if ret == 0:
            return True
        else:
            raise saga.NoSuccess._log(self._logger,
                                      "Could not cancel job %s because: %s" % (pid,out))

    def _job_suspend (self, job_id):
        """
        Attempt to suspend a job with commandline scontrol.  Raise
        exception when unsuccessful.
        """
        rm, pid     = self._adaptor.parse_id (job_id)
        ret, out, _ = self.shell.run_sync("scontrol suspend %s" % pid)
        if ret == 0:
            return True

        # check to see if the error was a permission error
        elif "Access/permission denied" in out:
            raise saga.PermissionDenied._log(self._logger,
                                      "Could not suspend job %s because: %s" % (pid, out))
        
        # it's some other error
        else:
            raise saga.NoSuccess._log(self._logger,
                                      "Could not suspend job %s because: %s" % (pid, out))

    def _job_resume (self, job_id):
        """
        Attempt to resume a job with commandline scontrol.  Raise
        exception when unsuccessful.
        """
        rm, pid     = self._adaptor.parse_id (job_id)
        ret, out, _ = self.shell.run_sync("scontrol resume %s" % pid)
        if ret == 0:
            return True

        # check to see if the error was a permission error
        elif "Access/permission denied" in out:
            raise saga.PermissionDenied._log(self._logger,
                                      "Could not suspend job %s because: %s" % (pid, out))
        
        # it's some other error
        else:
            raise saga.NoSuccess._log(self._logger,
                                      "Could not resume job %s because: %s" % (pid, out))


    # ----------------------------------------------------------------
    #
    @SYNC_CALL
    def create_job (self, jd) :
        """ Implements saga.adaptors.cpi.job.Service.create_job()
        """
        # this dict is passed on to the job adaptor class -- use it to pass any
        # state information you need there.
        adaptor_state = { "job_service"     : self, 
                          "job_description" : jd,
                          "job_schema"      : self.rm.schema,
                          "reconnect" : False}

        return saga.job.Job (_adaptor=self._adaptor,
                             _adaptor_state=adaptor_state)

    # ----------------------------------------------------------------
    @SYNC_CALL
    def get_url (self) :
        """ Implements saga.adaptors.cpi.job.Service.get_url()
        """
        return self.rm

    @SYNC_CALL
    def run_job (self, cmd, host):
        raise saga.NotImplemented._log (self._logger, "run_job not implemented"
                                        " for SLURM jobs -- please construct"
                                        " a job description and create a job"
                                        " manually.")

    # ----------------------------------------------------------------
    #
    @SYNC_CALL
    def list(self):
        """ Implements saga.adaptors.cpi.job.Service.list()
        """

        # ashleyz@login1:~$ squeue -h -o "%i" -u ashleyz                                                                                                                                                                       
        # 255042
        # 255035
        # 255028
        # 255018

        # this line gives us a nothing but jobids for our user
        ret, out, _ = self.shell.run_sync('squeue -h -o "%%i" -u %s' 
                                          % self.rm.detected_username)

        # mangle our results into the proper id format
        output = ["[%s]-[%s]" % (self.rm, i) for i in out.strip().split("\n")]
        return output

  # # ----------------------------------------------------------------
  # #
    @SYNC_CALL
    def get_job (self, jobid):
        # try to get some information about this job and throw it into
        # our job dictionary.
        # self.jobs[jobid] = self._retrieve_job(jobid)

        # this dict is passed on to the job adaptor class -- use it to pass any
        # state information you need there.
        adaptor_state = {"job_service": self,
                         # TODO: fill job description
                         "job_description": saga.job.Description(),
                         "job_schema": self.rm.schema,
                         "reconnect": True,
                         "reconnect_jobid": jobid
                        }
        return_job =  saga.job.Job(_adaptor=self._adaptor,
                            _adaptor_state=adaptor_state)
        return return_job
        
  #     if jobid not in self._jobs.values():
  #         msg = "Service instance doesn't know a Job with ID '%s'" % (jobid)
  #         raise saga.BadParameter._log (self._logger, msg)
  #     else:
  #         for (job_obj, job_id) in self._jobs.iteritems():
  #             if job_id == jobid:
  #                 return job_obj.get_api ()
  #
  #
    # ----------------------------------------------------------------
    #
    def container_run (self, jobs) :
        self._logger.debug("container run: %s"  %  str(jobs))
        # TODO: this is not optimized yet
        for job in jobs:
            job.run()
    
    
    # ----------------------------------------------------------------
    #
    def container_wait (self, jobs, mode, timeout) :
        self._logger.debug("container wait: %s"  %  str(jobs))
        # TODO: this is not optimized yet
        for job in jobs:
            job.wait()
    
    
    # ----------------------------------------------------------------
    #
    def container_cancel (self, jobs) :
        self._logger.debug("container cancel: %s"  %  str(jobs))
        for job in jobs:
            job.cancel()
    

    # ----------------------------------------------------------------
    #
    @SYNC_CALL
    def container_get_states(self, jobs):

        states = list()
        for job in jobs:
            job.get_state()

        return states


###############################################################################
#
class SLURMJob (saga.adaptors.cpi.job.Job):
    """ Implements saga.adaptors.cpi.job.Job
    """
    # ----------------------------------------------------------------
    #
    def __init__ (self, api, adaptor) :

        _cpi_base = super  (SLURMJob, self)
        _cpi_base.__init__ (api, adaptor)

    # ----------------------------------------------------------------
    #
    @SYNC_CALL
    def init_instance (self, job_info):
        """ Implements saga.adaptors.cpi.job.Job.init_instance()
        """

        self.jd = job_info["job_description"]
        self.js = job_info["job_service"] 

        # the js is responsible for job bulk operations -- which
        # for jobs only work for run()
        self._container       = self.js
        self._method_type     = "run"

        # initialize job attribute values
        self._id              = None
        self._name            = self.jd.name
        self._state           = saga.job.NEW
        self._exit_code       = None
        self._exception       = None
        self._started         = None
        self._finished        = None

        # think "reconnect" in terms of "reloading" job id, _NOT_
        # physically creating a new network connection
        if job_info['reconnect'] is True:
            self._id = job_info['reconnect_jobid']
            self._started = True
        else:
            self._started = False

        return self.get_api ()

    def _job_get_info (self, job_id):
        """ 
        use scontrol to grab job info 
        NOT CURRENTLY USED/TESTED, here for later
        """
        # if we don't have the job in our dictionary, we don't want it
        # TODO: verify correctness, we should probably probe anyhow
        #       in case it was added by an external app
        if job_id not in self.js.jobs:
            message = "Unknown job ID: %s. Can't update state." % job_id
            raise saga.NoSuccess._log(self._logger, message)

        # prev. info contains the info collect when _job_get_info
        # was called the last time
        prev_info = self.js.jobs[job_id]

        # if the 'gone' flag is set, there's no need to query the job
        # state again. it's gone forever
        if prev_info['gone'] is True:
            self._logger.warning("Job information is not available anymore.")
            return prev_info

        # curr. info will contain the new job info collect. it starts off
        # as a copy of prev_info (don't use deepcopy because there is an API 
        # object in the dict -> recursion)
        curr_info = dict()
        curr_info['job_id'     ] = prev_info.get ('job_id'     )
        curr_info['state'      ] = prev_info.get ('state'      )
        curr_info['create_time'] = prev_info.get ('create_time')
        curr_info['start_time' ] = prev_info.get ('start_time' )
        curr_info['end_time'   ] = prev_info.get ('end_time'   )
        curr_info['comp_time'  ] = prev_info.get ('comp_time'  )
        curr_info['exec_hosts' ] = prev_info.get ('exec_hosts' )
        curr_info['gone'       ] = prev_info.get ('gone'       )

        rm, pid = self._adaptor.parse_id(job_id)

        # update current info with scontrol
        ret, out, _ = self.js.shell.run_sync('scontrol show job %s' % pid)
      # self._logger.debug("Updating job status using the following information:\n%s" % out) 

        # update the state
        curr_info['state'] = self._job_get_state(job_id)


        match = self.js.scontrol_create_time_re.search(out)
        if match:
            curr_info['create_time'] = match.group(1)
            self._logger.debug("create_time for job %s detected as %s" % \
                               (pid, curr_info['create_time']))

        match = self.js.scontrol_start_time_re.search(out)
        if match:
            curr_info['start_time'] = match.group(1)
            self._logger.debug("start_time for job %s detected as %s" % \
                               (pid, curr_info['start_time']))

        match = self.js.scontrol_end_time_re.search(out)
        if match:
            curr_info['end_time'] = match.group(1)
            self._logger.debug("end_time for job %s detected as %s" % \
                               (pid, curr_info['end_time']))

        match = self.js.scontrol_comp_time_re.search(out)
        if match:
            curr_info['comp_time'] = match.group(1)
            self._logger.debug("comp_time for job %s detected as %s" % \
                               (pid, curr_info['comp_time']))

        match = self.js.scontrol_exec_hosts_re.search(out)
        if match:
            curr_info['exec_hosts'] = match.group(1)
            self._logger.debug("exec_hosts for job %s detected as %s" % \
                               (pid, curr_info['exec_hosts']))

        return curr_info


    def _job_get_state (self, job_id) :
        """ get the job state from the wrapper shell """

        # if the state is NEW and we haven't sent out a run command, keep
        # it listed as NEW
        if self._state == saga.job.NEW and not self._started:
            return saga.job.NEW

        # if we don't even have an ID, state is unknown
        # TODO: VERIFY CORRECTNESS

        if job_id==None:
            return saga.job.UNKNOWN

        # if the state is DONE, CANCELED or FAILED, it is considered
        # final and we don't need to query the backend again
        if self._state == saga.job.CANCELED or self._state == saga.job.FAILED \
            or self._state == saga.job.DONE:
            return self._state

        rm, pid = self._adaptor.parse_id (job_id)

        try:
            ret, out, _ = self.js.shell.run_sync('scontrol show job %s' % pid)
            match       = self.js.scontrol_jobstate_re.search(out)

            if match:
                slurm_state = match.group(1)
            else:
                # no jobstate found from scontrol
                # the job may have finished a while back, use sacct to
                # look at the full slurm history
                slurm_state = self._sacct_jobstate_match(pid)
                if not slurm_state:
                    # no jobstate found in slurm
                    return saga.job.UNKNOWN

            return self.js._slurm_to_saga_jobstate(slurm_state)

        except Exception, ex:
            raise saga.NoSuccess("Error getting the job state for "
                                 "job %s:\n%s"%(pid,ex))
        
        raise saga.NoSuccess._log (self._logger,
                                   "Internal SLURM adaptor error"
                                   " in _job_get_state")

    def _sacct_jobstate_match (self, pid):
        """ get the job state from the slurm accounting data """
        ret, sacct_out, _ = self.js.shell.run_sync(
            "sacct --format=JobID,State --parsable2 --noheader --jobs=%s" % pid)
        # output will look like:
        # 500723|COMPLETED
        # 500723.batch|COMPLETED
        # or:
        # 500682|CANCELLED by 900369
        # 500682.batch|CANCELLED

        for line in sacct_out.strip().split('\n'):
            (slurm_id, slurm_state) = line.split('|', 1)
            if slurm_id == pid and slurm_state:
                return slurm_state.split()[0].strip()

        return None

    # ----------------------------------------------------------------
    #
    @SYNC_CALL
    def get_state(self):
        """ Implements saga.adaptors.cpi.job.Job.get_state()
        """
        self._state = self._job_get_state (self._id)
        return self._state

    # ----------------------------------------------------------------
    #
    @SYNC_CALL
    def get_description (self):
        return self.jd

    # ----------------------------------------------------------------
    #
    @SYNC_CALL
    def get_service_url(self):
        """ implements saga.adaptors.cpi.job.Job.get_service_url()
        """
        return self.js.rm
   

    # ----------------------------------------------------------------
    #
    @SYNC_CALL
    def wait(self, timeout):
        time_start = time.time()
        time_now   = time_start
        rm, pid    = self._adaptor.parse_id(self._id)

        while True:
            state = self._job_get_state(self._id)
            self._logger.debug("wait() state for job id %s:%s"%(self._id, state))

            if state == saga.job.UNKNOWN :
                log_error_and_raise("cannot get job state", saga.IncorrectState, self._logger)

            if state == saga.job.DONE or \
               state == saga.job.FAILED or \
               state == saga.job.CANCELED:
                    return True
            # avoid busy poll
            time.sleep(0.5)

            # check if we hit timeout
            if timeout >= 0:
                time_now = time.time()
                if time_now - time_start > timeout:
                    return False

        return True

    # ----------------------------------------------------------------
    #
    # Andre Merzky: In general, the job ID is something which is generated by the adaptor or by the backend, and the user should not interpret it.  So, you can do that.  Two caveats though:
    # (a) The ID MUST remain constant once it is assigned to a job (imagine an application hashes on job ids, for example.
    # (b) the ID SHOULD follow the scheme [service_url]-[backend-id] -- and in that case, you should make sure that the URL part of the ID can be used to create a new job service instance...

    @SYNC_CALL
    def get_id (self) :
        """ Implements saga.adaptors.cpi.job.Job.get_id() """        
        return self._id
   
    # ----------------------------------------------------------------
    #
    @SYNC_CALL
    def get_name (self):
        """ Implements saga.adaptors.cpi.job.Job.get_name() """        
        return self._name

    # ----------------------------------------------------------------
    #
    @SYNC_CALL
    def get_exit_code(self) :
        """ Implements saga.adaptors.cpi.job.Job.get_exit_code()
        """   
        return self.js._job_get_exit_code(self._id)

    #  ----------------------------------------------------------------
    #
    @SYNC_CALL
    def suspend(self) :
        """ Implements saga.adaptors.cpi.job.Job.get_exit_code()
        """ 
        return self.js._job_suspend(self._id)

    # ----------------------------------------------------------------
    @SYNC_CALL
    def resume(self) :
        """ Implements saga.adaptors.cpi.job.Job.get_exit_code()
        """ 
        return self.js._job_resume(self._id)

  
    # ----------------------------------------------------------------
    #
    @SYNC_CALL
    def get_created(self) :
        """ Implements saga.adaptors.cpi.job.Job.get_created()
        """     
        return self._job_get_info(self._id)['create_time']
        #raise saga.NotImplemented._log (self._logger, "get_created not"
        #                                " implemented for SLURM jobs.")
        #return
  
    # ----------------------------------------------------------------
    #
    @SYNC_CALL
    def get_started(self) :
        """ Implements saga.adaptors.cpi.job.Job.get_started()
        """        
        return self._job_get_info(self._id)['start_time']
        #raise saga.NotImplemented._log (self._logger, "get_started not"
        #                                " implemented for SLURM jobs.")
        #return self._started
  
    # ----------------------------------------------------------------
    #
    @SYNC_CALL
    def get_finished(self) :
        """ Implements saga.adaptors.cpi.job.Job.get_finished()
        """
        return self._job_get_info(self._id)['end_time']
        #raise saga.NotImplemented._log (self._logger, "get_finished not"
        #                                " implemented for SLURM jobs.")
        #return self._finished
   
    # ----------------------------------------------------------------
    #
    @SYNC_CALL
    def get_execution_hosts(self) :
        """ Implements saga.adaptors.cpi.job.Job.get_execution_hosts()
        """        

        return self._job_get_info(self._id)['exec_hosts']
        #raise saga.NotImplemented._log (self._logger, "get_execution_hosts not"
        #                                " implemented for SLURM jobs.")
        #return
  
    # ----------------------------------------------------------------
    #
    @SYNC_CALL
    def cancel(self, timeout):
        #scancel id
        self.js._job_cancel(self._id)
        self._state=saga.job.CANCELED

    # ----------------------------------------------------------------
    #
    @SYNC_CALL
    def run(self): 
        """ Implements saga.adaptors.cpi.job.Job.run()
        """
        self._id = self.js._job_run (self.jd)
        self._started = True



