
__author__    = "Hangi,Kim"
__copyright__ = "Copyright 2012-2013, The SAGA Project"
__license__   = "MIT"


""" IBM LoadLeveler job adaptor implementation
    reference for pbs job adaptor & sge job adaptor implementation
	Hangi, Kim hgkim@kisti.re.kr
"""

import saga.utils.pty_shell

import saga.url as surl
import saga.adaptors.cpi.base
import saga.adaptors.cpi.job

from saga.adaptors.sge.sgejob import SgeKeyValueParser

import os
import re
import time
from urlparse import parse_qs
from datetime import datetime


SYNC_CALL = saga.adaptors.cpi.decorators.SYNC_CALL


# --------------------------------------------------------------------
#
def log_error_and_raise(message, exception, logger):
    logger.error(message)
    raise exception(message)


# --------------------------------------------------------------------
#
def _ll_to_saga_jobstate(lljs):
    """ translates a loadleveler one-letter state to saga
        pbs_loadl_comparison.xlsx
    """
    if lljs == 'C':
        return saga.job.DONE
    elif lljs == 'S':
        return saga.job.PENDING
    elif lljs == 'ST':
        return saga.job.PENDING
    elif lljs == 'I':
        return saga.job.PENDING
    elif lljs == 'R':
        return saga.job.RUNNING
    else:
        return saga.job.UNKNOWN



def getId(out):

    t = out.split('\n')

    jobId = None

    for line in t:
        if line.startswith('Job'):
            tmpStr=line.split(' ')
            jobId=tmpStr[1]
            break

        elif re.search('The job ".+" has been submitted.', line):
            # Format: llsubmit: The job "srv03-ib.443336" has been submitted.
            jobId = re.findall(r'"(.*?)"', line)[0]
            break

    if not jobId:
        raise Exception("Failed to detect jobId.")

    return jobId


# --------------------------------------------------------------------
# some private defs
#
_PTY_TIMEOUT = 2.0

# --------------------------------------------------------------------
# the adaptor name
#
_ADAPTOR_NAME          = "saga.adaptor.loadljob"
_ADAPTOR_SCHEMAS       = ["loadl", "loadl+ssh", "loadl+gsissh"]
_ADAPTOR_OPTIONS       = [
    {
    'category'         : 'saga.adaptor.loadljob',
    'name'             : 'purge_on_start',
    'type'             : bool,
    'default'          : True,
    'valid_options'    : [True, False],
    'documentation'    : '''Purge temporary job information for all
                          jobs which are older than a number of days.
                          The number of days can be configured with <purge_older_than>.''',
    'env_variable'     : None
    },
    {
    'category'         : 'saga.adaptor.loadljob',
    'name'             : 'purge_older_than',
    'type'             : int,
    'default'          : 30,
    #'valid_options'    : [True, False],
    'documentation'    : '''When <purge_on_start> is enabled this specifies the number
                            of days to consider a temporary file older enough to be deleted.''',
    'env_variable'     : None
    },
]


# --------------------------------------------------------------------
# the adaptor capabilities & supported attributes
#
_ADAPTOR_CAPABILITIES = {
    "jdes_attributes":   [saga.job.NAME,
                          saga.job.EXECUTABLE,
                          saga.job.ARGUMENTS,
                          saga.job.ENVIRONMENT,
                          saga.job.INPUT,
                          saga.job.OUTPUT,
                          saga.job.ERROR,
                          saga.job.QUEUE,
                          saga.job.PROJECT,
                          saga.job.JOB_CONTACT,
                          saga.job.WALL_TIME_LIMIT,
                          saga.job.WORKING_DIRECTORY,
                          saga.job.TOTAL_PHYSICAL_MEMORY,
                          saga.job.PROCESSES_PER_HOST,
                          saga.job.TOTAL_CPU_COUNT],
    "job_attributes":    [saga.job.EXIT_CODE,
                          saga.job.EXECUTION_HOSTS,
                          saga.job.CREATED,
                          saga.job.STARTED,
                          saga.job.FINISHED],
    "metrics":           [saga.job.STATE],
    "contexts":          {"ssh": "SSH public/private keypair",
                          "x509": "GSISSH X509 proxy context",
                          "userpass": "username/password pair (ssh)"}
}

# --------------------------------------------------------------------
# the adaptor documentation
#
_ADAPTOR_DOC = {
    "name":          _ADAPTOR_NAME,
    "cfg_options":   _ADAPTOR_OPTIONS,
    "capabilities":  _ADAPTOR_CAPABILITIES,
    "description":  """
The LoadLeveler adaptor allows to run and manage jobs on ` IBM LoadLeveler<http://www-03.ibm.com/systems/software/loadleveler/>`_
controlled HPC clusters.
""",
    "example": "examples/jobs/loadljob.py",
    "schemas": {"loadl":        "connect to a local cluster",
                "loadl+ssh":    "conenct to a remote cluster via SSH",
                "loadl+gsissh": "connect to a remote cluster via GSISSH"}
}

# --------------------------------------------------------------------
# the adaptor info is used to register the adaptor with SAGA
#
_ADAPTOR_INFO = {
    "name":    _ADAPTOR_NAME,
    "version": "v0.1",
    "schemas": _ADAPTOR_SCHEMAS,
    "cpis": [
        {
        "type": "saga.job.Service",
        "class": "LOADLJobService"
        },
        {
        "type": "saga.job.Job",
        "class": "LOADLJob"
        }
    ]
}


###############################################################################
# The adaptor class
#class Adaptor (saga.adaptors.cpi.base.AdaptorBase):
class Adaptor (saga.adaptors.base.Base):
    """ this is the actual adaptor class, which gets loaded by SAGA (i.e. by 
        the SAGA engine), and which registers the CPI implementation classes 
        which provide the adaptor's functionality.
    """

    # ----------------------------------------------------------------
    #
    def __init__(self):

        saga.adaptors.base.Base.__init__(self, _ADAPTOR_INFO, _ADAPTOR_OPTIONS)

        self.id_re = re.compile('^\[(.*)\]-\[(.*?)\]$')
        self.opts = self.get_config(_ADAPTOR_NAME)

        self.purge_on_start = self.opts['purge_on_start'].get_value()
        self.purge_older_than = self.opts['purge_older_than'].get_value()

    # ----------------------------------------------------------------
    #
    def sanity_check(self):
        # FIXME: also check for gsissh
        pass

    # ----------------------------------------------------------------
    #
    def parse_id(self, id):
        # split the id '[rm]-[pid]' in its parts, and return them.

        match = self.id_re.match(id)

        if not match or len(match.groups()) != 2:
            raise saga.BadParameter("Cannot parse job id '%s'" % id)

        return (match.group(1), match.group(2))


###############################################################################
#
class LOADLJobService (saga.adaptors.cpi.job.Service):
    """ implements saga.adaptors.cpi.job.Service
    """

    # ----------------------------------------------------------------
    #
    def __init__(self, api, adaptor):

        self._cpi_base = super(LOADLJobService, self)
        self._cpi_base.__init__(api, adaptor)

        self._adaptor = adaptor

    # ----------------------------------------------------------------
    #
    def __del__(self):

        self.finalize(kill_shell=True)

    # ----------------------------------------------------------------
    #
    @SYNC_CALL
    def init_instance(self, adaptor_state, rm_url, session):
        """ service instance constructor
        """
        self.rm      = rm_url
        self.session = session
        self.ppn     = 0 # check for remove
        self.jobs    = dict()
        self.cluster_option = ''
        self.energy_policy_tag = None
        self.island_count = None
        self.node_usage = None
        self.network_mpi = None
        self.blocking = None
        self.job_type = 'MPICH' # TODO: Is this a sane default?
        self.enforce_resource_submission = False
        self.enforce_consumable_cpus = False
        self.enforce_consumable_memory = False
        self.enforce_consumable_virtual_memory = False
        self.enforce_consumable_large_page_memory = False
        self.temp_path = "$HOME/.saga/adaptors/loadl_job"

        # LoadLeveler has two ways of specifying the executable and arguments.
        # - Explicit: the executable and arguments are specified as parameters.
        # - Implicit: the (remainder of the) job script is the task.
        #
        # Currently we don't know how this policy can be detected at runtime.
        # We know that providing both will not work in all cases.
        #
        # As the IBM Red Book documents the explicit exec only,
        # we should use that as a default.
        # Currently we just use a hack to workaround Joule.
        #
        # Note: does this now simply become a Joule hack?
        #
        # TODO: Split script into submission file and script and use that for
        #       explicit exec?
        self.explicit_exec = False

        rm_scheme = rm_url.scheme
        pty_url   = surl.Url (rm_url)

        # this adaptor supports options that can be passed via the
        # 'query' component of the job service URL.
        if rm_url.query is not None:
            for key, val in parse_qs(rm_url.query).iteritems():
                if key == 'cluster':
                    self.cluster_option = " -X %s" % val[0]
                elif key == 'energy_policy_tag':
                    self.energy_policy_tag = val[0]
                elif key == 'island_count':
                    self.island_count = val[0]
                elif key == 'node_usage':
                    self.node_usage = val[0]
                elif key == 'network_mpi':
                    self.network_mpi = val[0]
                elif key == 'blocking':
                    self.blocking = val[0]
                elif key == 'job_type':
                    self.job_type = val[0]
                elif key == 'enforce_consumable_cpus':
                    self.enforce_consumable_cpus = True
                    self.enforce_resource_submission = True
                elif key == 'enforce_consumable_memory':
                    self.enforce_consumable_memory = True
                    self.enforce_resource_submission = True
                elif key == 'enforce_consumable_virtual_memory':
                    self.enforce_consumable_virtual_memory = True
                    self.enforce_resource_submission = True
                elif key == 'enforce_consumable_large_page_memory':
                    self.enforce_consumable_large_page_memory = True
                    self.enforce_resource_submission = True
                elif key == 'explicit_exec':
                    self.explicit_exec = True

        # we need to extract the scheme for PTYShell. That's basically the
        # job.Service Url without the loadl+ part. We use the PTYShell to execute
        # loadleveler commands either locally or via gsissh or ssh.
        if rm_scheme == "loadl":
            pty_url.scheme = "fork"
        elif rm_scheme == "loadl+ssh":
            pty_url.scheme = "ssh"
        elif rm_scheme == "loadl+gsissh":
            pty_url.scheme = "gsissh"

        # these are the commands that we need in order to interact with Load Leveler.
        # the adaptor will try to find them during initialize(self) and bail
        # out in case they are note avaialbe.
        self._commands = {'llq': None,
                          'llsubmit':     None,
                          'llcancel':     None}

        self.shell = saga.utils.pty_shell.PTYShell(pty_url, self.session)

        #self.shell.set_initialize_hook(self.initialize)
        #self.shell.set_finalize_hook(self.finalize)

        self.initialize()

        return self.get_api ()

    # ----------------------------------------------------------------
    #
    def close (self) :
        if  self.shell :
            self.shell.finalize (True)


    # ----------------------------------------------------------------
    #
    def initialize(self):
        # check if all required loadleveler tools are available
        for cmd in self._commands.keys():
            ret, out, _ = self.shell.run_sync("which %s " % cmd)
            self._logger.info(ret)
            self._logger.info(out)
            if ret != 0:
                message = "Error finding LoadLeveler tools: %s" % out
                log_error_and_raise(message, saga.NoSuccess, self._logger)
            else:
                path = out.strip()  # strip removes newline
                ret, out, _ = self.shell.run_sync("%s -v" % cmd)
                if ret != 0:
                    message = "Error finding LoadLeveler tools: %s" % out
                    log_error_and_raise(message, saga.NoSuccess,
                        self._logger)
                else:
                    # version is reported as: "version: x.y.z"
                    version = out.strip().split()[1]

                    # add path and version to the command dictionary
                    self._commands[cmd] = {"path":    path,
                                           "version": version}

        self._logger.info("Found LoadLeveler tools: %s" % self._commands)

        # see if we can get some information about the cluster, e.g.,
        # different queues, number of processes per node, etc.
        # TODO: this is quite a hack. however, it *seems* to work quite
        #       well in practice.
        # modi by hgkim

        # purge temporary files
        if self._adaptor.purge_on_start:
            cmd = "find $HOME/.saga/adaptors/loadl_job" \
                  " -type f -mtime +%d -print -delete | wc -l" % self._adaptor.purge_older_than
            ret, out, _ = self.shell.run_sync(cmd)
            if ret == 0 and out != "0":
                self._logger.info("Purged %s temporary files" % out)

    # ----------------------------------------------------------------
    #
    def finalize(self, kill_shell=False):
        if  kill_shell :
            if  self.shell :
                self.shell.finalize (True)

    def __remote_mkdir(self, path):
        """
        Creates a directory on the remote host.
        :param path: the remote directory to be created.
        """
        # check if the path exists
        ret, out, _ = self.shell.run_sync(
                        "(test -d %s && echo -n 0) || (mkdir -p %s && echo -n 1)" % (path, path))

        if ret == 0 and out == "1":
            self._logger.info("Remote directory created: %s" % path)
        elif ret != 0:
            # something went wrong
            message = "Couldn't create remote directory - %s\n%s" % (out, path)
            log_error_and_raise(message, saga.NoSuccess, self._logger)


    def __remote_job_info_path(self, loadl_job_id="$LOADL_JOB_NAME"):
        """
        Returns the path of the remote job info file.
        :param loadl_job_id: the LoadLeveler job id, if omitted an environment variable representing the job id will be used.
        :return: path to the remote job info file
        """

        return "%s/%s" % (self.temp_path, loadl_job_id)

    def __clean_remote_job_info(self, loadl_job_id):
        """
        Removes the temporary remote file containing job info.
        :param loadl_job_id: the LoadLeveler job id
        """

        path = self.__remote_job_info_path(loadl_job_id)
        ret, out, _ = self.shell.run_sync("rm %s" % path)
        if ret != 0:
            self._logger.debug("Remote job info couldn't be removed: %s" % path)

    def __get_remote_job_info(self, loadl_job_id):
        """
        Obtains the job info from a temporary remote file created by the llsubmit script.
        :param loadl_job_id: the LoadLeveler job id
        :return: a dictionary with the job info
        """
        ret, out, _ = self.shell.run_sync("cat %s" % self.__remote_job_info_path(loadl_job_id))
        if ret != 0:
            return None

        qres = SgeKeyValueParser(out, key_suffix=":").as_dict()

        if "signal" in qres:
            state = saga.job.CANCELED
        elif "exit_status" in qres:
            if int(qres.get("exit_status")) == 0:
                state = saga.job.DONE
            else:
                state = saga.job.FAILED
        else:
            state = saga.job.RUNNING

        job_info = dict(
                    state=state,
                    exec_hosts=qres.get("hostname"),
                    returncode=int(qres.get("exit_status", -1)),
                    create_time=qres.get("qsub_time"),
                    start_time=qres.get("start_time"),
                    end_time=qres.get("end_time"),
                    gone=False)

        return job_info

    def __generate_llsubmit_script(self, jd):
        """ 
        generates a IMB LoadLeveler script from a SAGA job description
        :param jd: job descriptor
        :return: the llsubmit script
        """
        loadl_params = ''
        exec_string = ''
        args_strings = ''

        if jd.executable is not None:
            exec_string = "%s" % (jd.executable)
        if jd.arguments is not None:
            for arg in jd.arguments:
                args_strings += "%s " % (arg)

        if jd.name is not None:
            loadl_params += "#@ job_name = %s \n" % jd.name

        if jd.environment is not None:
            variable_list = ''
            for key in jd.environment.keys():
                variable_list += "%s=%s;" % (key, jd.environment[key])
            loadl_params += "#@ environment = %s \n" % variable_list

        # Energy
        if self.energy_policy_tag:
            loadl_params += "#@ energy_policy_tag = %s\n" % self.energy_policy_tag
            loadl_params += "#@ minimize_time_to_solution = yes\n"

        if jd.working_directory is not None:
            loadl_params += "#@ initialdir = %s\n" % jd.working_directory
        if jd.output is not None:
            loadl_params += "#@ output = %s\n" % jd.output
        if jd.error is not None:
            loadl_params += "#@ error = %s\n" % jd.error
        if jd.wall_time_limit is not None:
            hours = jd.wall_time_limit / 60
            minutes = jd.wall_time_limit % 60
            loadl_params += "#@ wall_clock_limit = %s:%s:00\n" \
                % (str(hours), str(minutes))

        if jd.total_cpu_count is None:
            # try to come up with a sensible (?) default value
            jd.total_cpu_count = 1
        else:
            if jd.total_cpu_count > 1:
                if self.job_type not in ['bluegene']: # 'bluegene' and total_tasks dont live well together
                    loadl_params += "#@ total_tasks = %s\n" % jd.total_cpu_count

                loadl_params += "#@ job_type = %s\n" % self.job_type

        if self.job_type == 'bluegene':
            BGQ_CORES_PER_NODE = 16 # Only true for BG/Q
            if jd.total_cpu_count % BGQ_CORES_PER_NODE > 0:
                raise Exception("Number of cores requested is no multiple of 16.")
            loadl_params += "#@ bg_size = %d\n" % (jd.total_cpu_count / BGQ_CORES_PER_NODE)

        if self.blocking:
            loadl_params += "#@ blocking = %s\n" % self.blocking

        if self.enforce_resource_submission:

            loadl_params += "#@ resources ="

            if self.enforce_consumable_cpus:
                loadl_params += " ConsumableCpus(%d)" % jd.total_cpu_count

            if self.enforce_consumable_memory:
                if jd.total_physical_memory is None:
                    raise Exception("total_physical_memory is not set, but required by enforce_consumable_memory.")
                loadl_params += " ConsumableMemory(%dmb)" % jd.total_physical_memory

            if self.enforce_consumable_large_page_memory:
                # TODO: Not sure how to get a sensible value for this
                if jd.total_physical_memory is None:
                    raise Exception("total_physical_memory is not set, but required by enforce_consumable_large_page_memory.")
                loadl_params += " ConsumableLargePageMemory(%dmb)" % jd.total_physical_memory

            if self.enforce_consumable_virtual_memory:
                # TODO: Not sure how to get a sensible value for this
                if jd.total_physical_memory is None:
                    raise Exception("total_physical_memory is not set, but required by enforce_consumable_virtual_memory.")
                loadl_params += " ConsumableVirtualMemory(%dmb)" % jd.total_physical_memory

            loadl_params += "\n"

        # Number of islands to allocate resources on, can specify a number, or a min/max
        if self.island_count:
            loadl_params += "#@ island_count = %s\n" % self.island_count

        # Specify network configuration
        if self.network_mpi:
            loadl_params += "#@ network.MPI = %s\n" % self.network_mpi

        # Specify node usage policy
        if self.node_usage:
            loadl_params += "#@ node_usage = %s\n" % self.node_usage

        if jd.job_contact is not None:
            if len(jd.job_contact) > 1:
                raise Exception("Only one notify user supported.")
            loadl_params += "#@ notify_user = %s\n" % jd.job_contact[0]
            loadl_params += "#@ notification = always\n"

        # some default (?) parameter that seem to work fine everywhere... 
        if jd.queue is not None:
            loadl_params += "#@ class = %s\n" % jd.queue
        else:
            loadl_params += "#@ class = edison\n"

        # finally, we 'queue' the job
        loadl_params += "#@ queue\n"

        # Job info, executable and arguments
        job_info_path = self.__remote_job_info_path()

        script_body = [
            'function aborted() {',
            '  echo Aborted with signal $1.',
            '  echo "signal: $1" >>%s' % job_info_path,
            '  echo "end_time: $(LC_ALL=en_US.utf8 date \'+%%a %%b %%d %%H:%%M:%%S %%Y\')" >>%s' % job_info_path,
            '  exit -1',
            '}',
            'mkdir -p %s' % self.temp_path,
            'for sig in SIGHUP SIGINT SIGQUIT SIGTERM SIGUSR1 SIGUSR2; do trap "aborted $sig" $sig; done',
            'echo "hostname: $HOSTNAME" >%s' % job_info_path,
            'echo "qsub_time: %s" >>%s' % (datetime.now().strftime("%a %b %d %H:%M:%S %Y"), job_info_path),
            'echo "start_time: $(LC_ALL=en_US.utf8 date \'+%%a %%b %%d %%H:%%M:%%S %%Y\')" >>%s' % job_info_path
        ]

        script_body += ['%s %s' % (exec_string, args_strings)]

        script_body += [
            'echo "exit_status: $?" >>%s' % job_info_path,
            'echo "end_time: $(LC_ALL=en_US.utf8 date \'+%%a %%b %%d %%H:%%M:%%S %%Y\')" >>%s' % job_info_path
        ]

        # convert exec and args into an string and
        # escape all double quotes and dollar signs, otherwise 'echo |'
        # further down won't work.
        # only escape '$' in args and exe. not in the params
        script_body = "\n".join(script_body).replace('$', '\\$')

        # Dirty Trick for Joule: it expects an "executable" parameter,
        # but doesn't really need it, therefore we pass it after the queue
        # parameter, where it is not used anymore.
        if self.explicit_exec:
            loadl_params += "#@ executable = BOGUS\n"

        loadlscript = "\n%s%s" % (loadl_params, script_body)

        return loadlscript.replace('"', '\\"')

    # ----------------------------------------------------------------
    #
    def _job_run(self, jd):
        """ runs a job via llsubmit
        """

        try:
            # create a LoadLeveler job script from SAGA job description
            script = self.__generate_llsubmit_script(jd)

            self._logger.debug("Generated LoadLeveler script: %s" % script)
        except Exception, ex:
            log_error_and_raise(str(ex), saga.BadParameter, self._logger)

        # try to create the working/output/error directories (if defined)
        # WARNING: this assumes a shared filesystem between login node and
        #           compute nodes.
        if jd.working_directory is not None and len(jd.working_directory) > 0:
            self.__remote_mkdir(jd.working_directory)

        if jd.output is not None and len(jd.output) > 0:
            self.__remote_mkdir(os.path.dirname(jd.output))

        if jd.error is not None and len(jd.error) > 0:
            self.__remote_mkdir(os.path.dirname(jd.error))

        # submit the LoadLeveler script
        # Now we want to execute the script. This process consists of two steps:
        # (1) we create a temporary file with 'mktemp' and write the contents of
        #     the generated Load Leveler script into it
        # (2) we call 'llsubmit <tmpfile>' to submit the script to the queueing system
        cmdline = """SCRIPTFILE=`mktemp -t SAGA-Python-LOADLJobScript.XXXXXX` &&  echo "%s" > $SCRIPTFILE && %s%s $SCRIPTFILE && rm -f $SCRIPTFILE""" %  (script, self._commands['llsubmit']['path'], self.cluster_option)
        self._logger.info("cmdline: %r", cmdline)
        ret, out, _ = self.shell.run_sync(cmdline)

        if ret != 0:
            # something went wrong
            message = "Error running job via 'llsubmit': %s. Script was: %s" \
                % (out, script)
            log_error_and_raise(message, saga.NoSuccess, self._logger)
        else:
            # stdout contains the job id
            #job_id = "[%s]-[%s]" % (self.rm, out.strip().split('.')[0])
            job_id = "[%s]-[%s]" % (self.rm, getId(out))
            self._logger.info("Submitted LoadLeveler job with id: %s" % job_id)

            # add job to internal list of known jobs.
            self.jobs[job_id] = {
                'state':        saga.job.PENDING,
                'exec_hosts':   None,
                'returncode':   None,
                'create_time':  None,
                'start_time':   None,
                'end_time':     None,
                'gone':         False
            }

            return job_id

    # ----------------------------------------------------------------
    #
    def _retrieve_job(self, job_id, max_retries=10):
        """ see if we can get some info about a job that we don't
            know anything about
            refactoring by referencing sgejob.py
        """
        rm, pid = self._adaptor.parse_id(job_id)

        # run the LoadLeveler 'llq' command to get some info about our job
        ret, out, _ = self.shell.run_sync("%s -j %s -r %%st %%dd %%cc %%jt %%c %%Xs" % \
                                          (self._commands['llq']['path'], pid))
        # output is something like
        # R!03/25/2014 13:47!!Serial!normal!kisti.kim
        # OR
        # llq: There is currently no job status to report.
        if ret != 0:
            message = "Couldn't reconnect to job '%s': %s" % (job_id, out)
            log_error_and_raise(message, saga.NoSuccess, self._logger)

        else:
            # the job seems to exist on the backend. let's gather some data
            job_info = {
                'state':        saga.job.UNKNOWN,
                'exec_hosts':   None,
                'returncode':   None,
                'create_time':  None,
                'start_time':   None,
                'end_time':     None,
                'gone':         False
            }

            #lastStr=out.rstrip().split('\n')[-1]
            lastStr=out.rstrip()
            self._logger.debug(lastStr)
            if lastStr.startswith('llq:'): # llq: There is currently no job status to report
                job_info = None
                retries = 0
                delay = 1
                while job_info is None and retries < max_retries:
                    job_info = self.__get_remote_job_info(pid)
                    #print "llq:", job_info
                    if job_info == None and retries > 0:
                        message = "__get_remote_job_info get None, pid: %s and retries: %d" % (pid, retries)
                        self._logger.debug(message)
                        # Exponential back-off
                        time.sleep(2**retries)
                    retries += 1

                if job_info == None:
                    message = "__get_remote_job_info exceed %d times(s), pid: %s" % (max_retries, pid)
                    log_error_and_raise(message, saga.NoSuccess, self._logger)

                self._logger.info("_retrieve_job: %r", job_info)
            else: # job is still in the queue
                results = lastStr.split('!')
                self._logger.info("results: %r",results)

                job_info['state'] = _ll_to_saga_jobstate(results[0])
                job_info['returncode'] = -1 # still running
                job_info['start_time'] = results[1]
                #job_info['exec_hosts'] = results[5]

            return job_info

    # ----------------------------------------------------------------
    #
    def _job_get_info(self, job_id):
        """ get job attributes via llq
        """

        # if we don't have the job in our dictionary, we don't want it
        if job_id not in self.jobs:
            message = "Unknown job ID: %s. Can't update state." % job_id
            log_error_and_raise(message, saga.NoSuccess, self._logger)

        # prev. info contains the info collect when _job_get_info
        # was called the last time
        prev_info = self.jobs[job_id]

        # if the 'gone' flag is set, there's no need to query the job
        # state again. it's gone forever
        if prev_info['gone'] is True:
            self._logger.warning("Job information is not available anymore.")
            return prev_info

        # if the job is in a terminal state don't expect it to change anymore
        if prev_info["state"] in [saga.job.CANCELED, saga.job.FAILED, saga.job.DONE]:
            return prev_info

        # retrieve updated job information
        curr_info = self._retrieve_job(job_id)
        if curr_info is None:
            prev_info["gone"] = True
            return prev_info

        # update the job info cache and return it
        self.jobs[job_id] = curr_info
        return curr_info

    # ----------------------------------------------------------------
    #
    def _job_get_state(self, job_id):
        """ get the job's state
        """
        # check if we have already reach a terminal state
        if self.jobs[job_id]['state'] == saga.job.CANCELED \
        or self.jobs[job_id]['state'] == saga.job.FAILED \
        or self.jobs[job_id]['state'] == saga.job.DONE:
            return self.jobs[job_id]['state']

        # check if we can / should update
        if (self.jobs[job_id]['gone'] is not True):
            self.jobs[job_id] = self._job_get_info(job_id=job_id)

        return self.jobs[job_id]['state']

    # ----------------------------------------------------------------
    #
    def _job_get_exit_code(self, job_id):
        """ get the job's exit code
        """
        # check if we can / should update
        if (self.jobs[job_id]['gone'] is not True) \
        and (self.jobs[job_id]['returncode'] is None):
            self.jobs[job_id] = self._job_get_info(job_id=job_id)

        return self.jobs[job_id]['returncode']

    # ----------------------------------------------------------------
    #
    def _job_get_execution_hosts(self, job_id):
        """ get the job's exit code
        """
        # check if we can / should update
        if (self.jobs[job_id]['gone'] is not True) \
        and (self.jobs[job_id]['exec_hosts'] is None):
            self.jobs[job_id] = self._job_get_info(job_id=job_id)

        return self.jobs[job_id]['exec_hosts']

    # ----------------------------------------------------------------
    #
    def _job_get_create_time(self, job_id):
        """ get the job's creation time
        """
        # check if we can / should update
        if (self.jobs[job_id]['gone'] is not True) \
        and (self.jobs[job_id]['create_time'] is None):
            self.jobs[job_id] = self._job_get_info(job_id=job_id)

        return self.jobs[job_id]['create_time']

    # ----------------------------------------------------------------
    #
    def _job_get_start_time(self, job_id):
        """ get the job's start time
        """
        # check if we can / should update
        if (self.jobs[job_id]['gone'] is not True) \
        and (self.jobs[job_id]['start_time'] is None):
            self.jobs[job_id] = self._job_get_info(job_id=job_id)

        return self.jobs[job_id]['start_time']

    # ----------------------------------------------------------------
    #
    def _job_get_end_time(self, job_id):
        """ get the job's end time
        """
        # check if we can / should update
        if (self.jobs[job_id]['gone'] is not True) \
        and (self.jobs[job_id]['end_time'] is None):
            self.jobs[job_id] = self._job_get_info(job_id=job_id)

        return self.jobs[job_id]['end_time']

    # ----------------------------------------------------------------
    #
    def _job_cancel(self, job_id):
        """ cancel the job via 'llcancel'
        """
        rm, pid = self._adaptor.parse_id(job_id)

        ret, out, _ = self.shell.run_sync("%s%s %s\n" \
            % (self._commands['llcancel']['path'], self.cluster_option, pid))

        if ret != 0:
            message = "Error canceling job via 'llcancel': %s" % out
            log_error_and_raise(message, saga.NoSuccess, self._logger)

        #self.__clean_remote_job_info(pid)

        # assume the job was succesfully canceld
        self.jobs[job_id]['state'] = saga.job.CANCELED

    # ----------------------------------------------------------------
    #
    def _job_wait(self, job_id, timeout):
        """ wait for the job to finish or fail
        """

        time_start = time.time()
        time_now   = time_start
        rm, pid    = self._adaptor.parse_id(job_id)

        while True:
            state = self._job_get_state(job_id=job_id)

            if state == saga.job.UNKNOWN :
                log_error_and_raise("cannot get job state", saga.IncorrectState, self._logger)

            if state == saga.job.DONE or \
               state == saga.job.FAILED or \
               state == saga.job.CANCELED:
                    #self.__clean_remote_job_info(pid)
                    return True
            # avoid busy poll
            time.sleep(0.5)

            # check if we hit timeout
            if timeout >= 0:
                time_now = time.time()
                if time_now - time_start > timeout:
                    return False

    # ----------------------------------------------------------------
    #
    @SYNC_CALL
    def create_job(self, jd):
        """ implements saga.adaptors.cpi.job.Service.get_url()
        """
        # check that only supported attributes are provided
        for attribute in jd.list_attributes():
            if attribute not in _ADAPTOR_CAPABILITIES["jdes_attributes"]:
                message = "'jd.%s' is not supported by this adaptor" \
                    % attribute
                log_error_and_raise(message, saga.BadParameter, self._logger)

        # this dict is passed on to the job adaptor class -- use it to pass any
        # state information you need there.
        adaptor_state = {"job_service":     self,
                         "job_description": jd,
                         "job_schema":      self.rm.schema,
                         "reconnect":       False
                         }

        return saga.job.Job(_adaptor=self._adaptor,
                            _adaptor_state=adaptor_state)

    # ----------------------------------------------------------------
    #
    @SYNC_CALL
    def get_job(self, jobid):
        """ Implements saga.adaptors.cpi.job.Service.get_job()
        """

        self._logger.info("get_job: %r", jobid)
        # try to get some information about this job and throw it into
        # our job dictionary.
        self.jobs[jobid] = self._retrieve_job(jobid)

        # this dict is passed on to the job adaptor class -- use it to pass any
        # state information you need there.
        adaptor_state = {"job_service":     self,
                         # TODO: fill job description
                         "job_description": saga.job.Description(),
                         "job_schema":      self.rm.schema,
                         "reconnect":       True,
                         "reconnect_jobid": jobid
                         }

        return saga.job.Job(_adaptor=self._adaptor,
                            _adaptor_state=adaptor_state)

    # ----------------------------------------------------------------
    #
    @SYNC_CALL
    def get_url(self):
        """ implements saga.adaptors.cpi.job.Service.get_url()
        """
        return self.rm

    # ----------------------------------------------------------------
    #
    @SYNC_CALL
    def list(self):
        """ implements saga.adaptors.cpi.job.Service.list()
        """
        ids = []

        ret, out, _ = self.shell.run_sync("unset GREP_OPTIONS; %s | grep `whoami`" %
                                          self._commands['llq']['path'])

        if ret != 0 and len(out) > 0:
            message = "failed to list jobs via 'llq': %s" % out
            log_error_and_raise(message, saga.NoSuccess, self._logger)
        elif ret != 0 and len(out) == 0:
            # llq | grep `` exits with 1 if the list is empty
            pass
        else:
            for line in out.split("\n"):
                # output looks like this:
                # v4c064.8637.0            ydkim       3/27 13:33 R  50  normal       v4c064
                # v4c064.8638.0            ydkim       3/27 13:37 R  50  normal       v4c064
                # v4c064.8639.0            ydkim       3/27 13:37 R  50  normal       v4c065
                # v4c064.8640.0            ydkim       3/27 13:37 R  50  normal       v4c065
                # v4c064.8641.0            ydkim       3/27 13:37 I  50  normal
                lineArray=line.split()
                if len(lineArray) > 1:
                    # lineArray[0] : v4c064.8637.0
                    tmpStr=lineArray[0].split('.')
                    jobid = "[%s]-[%s]" % (self.rm, ".".join(tmpStr[:2]))
                    ids.append(str(jobid))

        return ids


  # # ----------------------------------------------------------------
  # #
  # def container_run (self, jobs) :
  #     self._logger.debug ("container run: %s"  %  str(jobs))
  #     # TODO: this is not optimized yet
  #     for job in jobs:
  #         job.run ()
  #
  #
  # # ----------------------------------------------------------------
  # #
  # def container_wait (self, jobs, mode, timeout) :
  #     self._logger.debug ("container wait: %s"  %  str(jobs))
  #     # TODO: this is not optimized yet
  #     for job in jobs:
  #         job.wait ()
  #
  #
  # # ----------------------------------------------------------------
  # #
  # def container_cancel (self, jobs) :
  #     self._logger.debug ("container cancel: %s"  %  str(jobs))
  #     raise saga.NoSuccess ("Not Implemented");


###############################################################################
#
class LOADLJob (saga.adaptors.cpi.job.Job):
    """ implements saga.adaptors.cpi.job.Job
    """

    def __init__(self, api, adaptor):

        # initialize parent class
        self._cpi_base = super(LOADLJob, self)
        self._cpi_base.__init__(api, adaptor)

    @SYNC_CALL
    def init_instance(self, job_info):
        """ implements saga.adaptors.cpi.job.Job.init_instance()
        """
        # init_instance is called for every new saga.job.Job object
        # that is created
        self.jd = job_info["job_description"]
        self.js = job_info["job_service"]

        if job_info['reconnect'] is True:
            self._id      = job_info['reconnect_jobid']
            self._name    = self.jd.get(saga.job.NAME)
            self._started = True
        else:
            self._id      = None
            self._name    = self.jd.get(saga.job.NAME)
            self._started = False

        return self.get_api()

    # ----------------------------------------------------------------
    #
    @SYNC_CALL
    def get_state(self):
        """ implements saga.adaptors.cpi.job.Job.get_state()
        """
        if self._started is False:
            # jobs that are not started are always in 'NEW' state
            return saga.job.NEW
        else:
            return self.js._job_get_state(self._id)

    # ----------------------------------------------------------------
    #
    @SYNC_CALL
    def wait(self, timeout):
        """ implements saga.adaptors.cpi.job.Job.wait()
        """
        if self._started is False:
            log_error_and_raise("Can't wait for job that hasn't been started",
                saga.IncorrectState, self._logger)
        else:
            self.js._job_wait(self._id, timeout)

    # ----------------------------------------------------------------
    #
    @SYNC_CALL
    def cancel(self, timeout):
        """ implements saga.adaptors.cpi.job.Job.cancel()
        """
        if self._started is False:
            log_error_and_raise("Can't wait for job that hasn't been started",
                saga.IncorrectState, self._logger)
        else:
            self.js._job_cancel(self._id)

    # ----------------------------------------------------------------
    #
    @SYNC_CALL
    def run(self):
        """ implements saga.adaptors.cpi.job.Job.run()
        """
        self._id = self.js._job_run(self.jd)
        self._started = True

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
    def get_id(self):
        """ implements saga.adaptors.cpi.job.Job.get_id()
        """
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
    def get_exit_code(self):
        """ implements saga.adaptors.cpi.job.Job.get_exit_code()
        """
        if self._started is False:
            return None
        else:
            return self.js._job_get_exit_code(self._id)

    # ----------------------------------------------------------------
    #
    @SYNC_CALL
    def get_created(self):
        """ implements saga.adaptors.cpi.job.Job.get_created()
        """
        if self._started is False:
            return None
        else:
            return self.js._job_get_create_time(self._id)

    # ----------------------------------------------------------------
    #
    @SYNC_CALL
    def get_started(self):
        """ implements saga.adaptors.cpi.job.Job.get_started()
        """
        if self._started is False:
            return None
        else:
            return self.js._job_get_start_time(self._id)

    # ----------------------------------------------------------------
    #
    @SYNC_CALL
    def get_finished(self):
        """ implements saga.adaptors.cpi.job.Job.get_finished()
        """
        if self._started is False:
            return None
        else:
            return self.js._job_get_end_time(self._id)

    # ----------------------------------------------------------------
    #
    @SYNC_CALL
    def get_execution_hosts(self):
        """ implements saga.adaptors.cpi.job.Job.get_execution_hosts()
        """
        if self._started is False:
            return None
        else:
            return self.js._job_get_execution_hosts(self._id)

# vim: tabstop=8 expandtab shiftwidth=4 softtabstop=4
