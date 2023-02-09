
__author__    = "Andre Merzky, Mark Santcroos, Ole Weidner"
__copyright__ = "Copyright 2012-2015, The SAGA Project"
__license__   = "MIT"


""" Condor job adaptor implementation
"""

import saga.utils.pty_shell

import saga.url as surl
import saga.adaptors.base
import saga.adaptors.cpi.job

from saga.job.constants import *
from transferdirectives import TransferDirectives

import re
import os
import time
from urlparse import parse_qs
from tempfile import NamedTemporaryFile

SYNC_CALL = saga.adaptors.cpi.decorators.SYNC_CALL
ASYNC_CALL = saga.adaptors.cpi.decorators.ASYNC_CALL


# --------------------------------------------------------------------
#
def log_error_and_raise(message, exception, logger):
    logger.error(message)
    raise exception(message)


# --------------------------------------------------------------------
#
def _condor_to_saga_jobstate(condorjs):
    """ translates a condor one-letter state to saga
    """
    # From: http://pages.cs.wisc.edu/~adesmet/status.html
    #
    # JobStatus in job ClassAds
    #
    # 0   Unexpanded  U
    # 1   Idle    I
    # 2   Running R
    # 3   Removed X
    # 4   Completed   C
    # 5   Held    H
    # 6   Submission_err  E

    if int(condorjs) == 0:
        return saga.job.PENDING
    elif int(condorjs) == 1:
        return saga.job.PENDING
    elif int(condorjs) == 2:
        return saga.job.RUNNING
    elif int(condorjs) == 3:
        return saga.job.CANCELED
    elif int(condorjs) == 4:
        return saga.job.DONE
    elif int(condorjs) == 5:
        return saga.job.PENDING
    elif int(condorjs) == 6:
        return saga.job.FAILED
    else:
        return saga.job.UNKNOWN


# --------------------------------------------------------------------
#
def _condorscript_generator(url, logger, jds, option_dict=None):
    """ generates a Condor script from a SAGA job description
    """
    condor_file = str()

    #
    # HTCondor quoting/escaping:
    # http://research.cs.wisc.edu/htcondor/manual/current/condor_submit.html#SECTION0012514000000000000000
    #

    # Assuming that most of the JDs are similar, we take the first one
    common_jd = jds[0]

    ##### OPTIONS PASSED VIA JOB SERVICE URL #####
    ##
    if option_dict is not None:
        condor_file += "\n##### OPTIONS PASSED VIA JOB SERVICE URL #####\n##"

        for (key, value) in option_dict.iteritems():
            condor_file += "\n%s = %s" % (key, value)

    ##### OPTIONS PASSED VIA JOB DESCRIPTION #####
    ##
    condor_file += "\n\n##### OPTIONS PASSED VIA JOB DESCRIPTION #####\n##"

    # special treatment for universe - defaults to 'vanilla'
    if common_jd.queue:
        condor_file += "\nuniverse = %s" % common_jd.queue
    else:
        condor_file += "\nuniverse = vanilla"

    # Condor doesn't expand environment variables in arguments.
    # To support this functionality, we wrap by calling /bin/env.
    condor_file += "\nexecutable = /bin/env"

    # environment -> environment
    # http://research.cs.wisc.edu/htcondor/manual/current/condor_submit.html#SECTION0012514000000000000000
    environment = "environment ="
    if common_jd.environment is not None:
        variable_list = str()
        for key in common_jd.environment.keys():
            variable_list += "%s=%s " % (key, common_jd.environment[key])
        environment += " \"%s\"" % variable_list.strip()
    condor_file += "\n%s" % environment

    # project -> +ProjectName
    if common_jd.project is not None:
        condor_file += "\n+ProjectName = \"%s\"" % str(common_jd.project)

    if common_jd.total_cpu_count:
        condor_file += "\nrequest_cpus = %d" % common_jd.total_cpu_count

    # arguments -> arguments
    common_arguments = "arguments = \"/bin/sh -c '"

    # The actual executable becomes the first argument.
    if common_jd.executable is not None:
        common_arguments += "%s" % common_jd.executable
    else:
        message = "Executable not set"
        log_error_and_raise(message, saga.NoSuccess, logger)

    if common_jd.candidate_hosts is not None:

        # Special Characters
        HOST_EXCLUSION_SYMBOL = '!'
        SPECIAL_REQ_SYMBOL = '~'

        requirements = ""

        # TODO: This is bound to GLIDEIN which it shouldnt be.
        resource_key = 'GLIDEIN_ResourceName'

        # Whitelist sites, filter out "special" entries from the candidate host lists
        incl_sites = [host for host in common_jd.candidate_hosts if not host.startswith((SPECIAL_REQ_SYMBOL, HOST_EXCLUSION_SYMBOL))]
        if incl_sites:
            requirements += '(' +  ' || '.join(['%s =?= "%s"' % (resource_key, site) for site in incl_sites]) + ')'

        # Blacklist sites
        excl_sites = [host for host in common_jd.candidate_hosts if host.startswith(HOST_EXCLUSION_SYMBOL)]
        if excl_sites:
            if incl_sites:
                # If there were sites, start with an AND operator again
                requirements += ' && '
            requirements += '(' +  ' && '.join(['%s =!= "%s"' % (resource_key, site.lstrip(' ' + HOST_EXCLUSION_SYMBOL)) for site in excl_sites]) + ')'

        # Generic requirements handling
        # Get the '~special_requirements' and strip leading ~ and possible spaces.
        special = ' && '.join([host.lstrip(' ' + SPECIAL_REQ_SYMBOL) for host in common_jd.candidate_hosts if host.startswith(SPECIAL_REQ_SYMBOL)])
        if special:
            if incl_sites or excl_sites:
                # If there were white and/or black sites, start with an AND operator again
                requirements += ' && '
            requirements += special

        if requirements:
            condor_file += "\nrequirements = %s\n" % requirements

    # Transfer Directives
    if common_jd.attribute_exists('transfer_directives'):
        # All checking is done in _handle_file_transfers() already.
        td = common_jd.transfer_directives

        if hasattr(td, 'transfer_input_files'):
            transfer_input_files = "transfer_input_files = "
            for source in td.transfer_input_files:
                transfer_input_files += "%s, " % source
            condor_file += "\n%s" % transfer_input_files

        if hasattr(td, 'transfer_output_files'):
            transfer_output_files = "transfer_output_files = "
            for target in td.transfer_output_files:
                transfer_output_files += "%s, " % target
            condor_file += "\n%s" % transfer_output_files

    #
    # Per Job in Bulk settings
    #

    for jd in jds:

        condor_file += "\n##### PER JOB OPTIONS #####\n##"

        arguments = common_arguments

        if jd.arguments is not None:

            for arg in jd.arguments:
                # Condor can't deal with multi-line arguments. yep, that's how
                # bad of a software it is.
                if len(arg.split("\n")) > 1:
                    message = "Condor doesn't support multi-line arguments: %s" % arg
                    log_error_and_raise(message, saga.NoSuccess, logger)

                # Condor HATES double quotes in the arguments. It'll return
                # some crap like: "Found illegal unescaped double-quote: ...
                # That's why we escape them by repeating.
                arg = arg.replace('"', '""')
                arg = arg.replace("'", "''")

                # Escape dollars (for environment variables)
                #arg = arg.replace('$', '\\$')

                arguments += " %s" % arg

        # close the quote opened earlier
        arguments += '\'\"'
        condor_file += "\n%s" % arguments

        condor_file += "\ninitialdir = %s " % jd.working_directory

        # TODO: what to do when working directory is not set?
        logname = "saga-condor-job-$(cluster)_$(process).log"
        condor_file += "\nlog = %s " % os.path.join(jd.working_directory, logname)

        # output -> output
        if jd.output is not None:
            condor_file += "\noutput = %s " % os.path.join(jd.working_directory, jd.output)

        # error -> error
        if jd.error is not None:
            condor_file += "\nerror = %s " % os.path.join(jd.working_directory, jd.error)


        condor_file += "\nqueue\n"

    condor_file += "\n##### END OF FILE #####\n"

    return condor_file


# --------------------------------------------------------------------
# some private defs
#
_PTY_TIMEOUT = 2.0

# --------------------------------------------------------------------
# the adaptor name
#
_ADAPTOR_NAME          = "saga.adaptor.condorjob"
_ADAPTOR_SCHEMAS       = ["condor", "condor+ssh", "condor+gsissh"]
_ADAPTOR_OPTIONS       = []

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
                          saga.job.QUEUE, # TODO: map jd.queue to universe or pool?!
                          saga.job.PROJECT,
                          saga.job.WALL_TIME_LIMIT,
                          saga.job.WORKING_DIRECTORY,
                          saga.job.CANDIDATE_HOSTS,
                          saga.job.TOTAL_CPU_COUNT,
                          saga.job.PROCESSES_PER_HOST,
                          saga.job.SPMD_VARIATION,
                          saga.job.FILE_TRANSFER],
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
The (HT)Condor(-G) adaptor allows to run and manage jobs on a 
`Condor <http://research.cs.wisc.edu/htcondor/>`_ gateway.
""",
    "example": "examples/jobs/condorjob.py",
    "schemas": {"condor"        : "connect to a local gateway",
                "condor+ssh"    : "connect to a remote gateway via SSH",
                "condor+gsissh ": "connect to a remote gateway via GSISSH"}
}

# --------------------------------------------------------------------
# the adaptor info is used to register the adaptor with SAGA
#
_ADAPTOR_INFO = {
    "name"             : _ADAPTOR_NAME,
    "version"          : "v0.2",
    "schemas"          : _ADAPTOR_SCHEMAS,
    "capabilities"     : _ADAPTOR_CAPABILITIES,
    "cpis": [
        {
        "type"         : "saga.job.Service",
        "class"        : "CondorJobService"
        },
        {
        "type"         : "saga.job.Job",
        "class"        : "CondorJob"
        }
    ]
}


###############################################################################
# The adaptor class
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
        self.opts  = self.get_config (_ADAPTOR_NAME)

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
class CondorJobService (saga.adaptors.cpi.job.Service):
    """ implements saga.adaptors.cpi.job.Service
    """

    # ----------------------------------------------------------------
    #
    def __init__(self, api, adaptor):

        _cpi_base = super(CondorJobService, self)
        _cpi_base.__init__(api, adaptor)

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
        self.rm            = rm_url
        self.session       = session
        self.ppn           = 0
        self.is_cray       = False
        self.jobs          = dict()
        self.query_options = dict()

        rm_scheme = rm_url.scheme
        pty_url   = surl.Url (rm_url)

        # this adaptor supports options that can be passed via the
        # 'query' component of the job service URL.
        if rm_url.query is not None:
            for key, val in parse_qs(rm_url.query).iteritems():
                self.query_options[key] = val[0]

        # we need to extract the scheme for PTYShell. That's basically the
        # job.Service Url without the condor+ part. We use the PTYShell to
        # execute condor commands either locally or via gsissh or ssh.
        if rm_scheme == "condor":
            pty_url.scheme = "fork"
        elif rm_scheme == "condor+ssh":
            pty_url.scheme = "ssh"
        elif rm_scheme == "condor+gsissh":
            pty_url.scheme = "gsissh"

        # these are the commands that we need in order to interact with Condor.
        # the adaptor will try to find them during initialize(self) and bail
        # out in case they are not available.
        self._commands = {'condor_version': None,
                          'condor_submit':  None,
                          'condor_q':       None,
                          'condor_history': None,
                          'condor_rm':      None}

        self.shell = saga.utils.pty_shell.PTYShell(pty_url, self.session)

      # self.shell.set_initialize_hook(self.initialize)
      # self.shell.set_finalize_hook(self.finalize)

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
        # check if all required condor tools are available
        for cmd in self._commands.keys():
            ret, out, _ = self.shell.run_sync("which %s " % cmd)
            if ret != 0:
                message = "Error finding Condor tools: %s" % out
                log_error_and_raise(message, saga.NoSuccess, self._logger)
            else:
                path = out.strip()  # strip removes newline
                if cmd == 'condor_version':
                    ret, out, _ = self.shell.run_sync("%s" % cmd)
                    if ret != 0:
                        message = "Error determining Condor version: %s" % out
                        log_error_and_raise(message, saga.NoSuccess,
                            self._logger)
                    else:
                        # version is reported as:
                        # $CondorVersion: 7.8.6 Oct 25 2012 $
                        # $CondorPlatform: X86_64-CentOS_5.7 $
                        lines = out.split('\n')
                        version = lines[0].replace("$CondorVersion: ", "")
                        version = version.strip(" $")

                        # add path and version to the command dictionary
                        # TODO: change indentation below?
                self._commands[cmd] = {"path":    path,
                                       "version": version}

        self._logger.info("Found Condor tools: %s" % self._commands)

    # ----------------------------------------------------------------
    #
    def finalize(self, kill_shell=False):
        pass

    # ----------------------------------------------------------------
    #
    def _job_run(self, jd):
        """ runs a job via condor_submit
        """

        # Because we do funky shit with env and sh, we need to explicitly add
        # the executable to the transfer list.
        # (Of course not if the executable is already on the target systems,
        # defined by the fact that it starts with ./
        if jd.executable.startswith('./') and os.path.exists(jd.executable):

            # TODO: Check if the executable is already in the file_transfer list,
            # because then we don't need to implicitly add it anymore.
            # (For example, if the executable is not in the local directory, it needs
            # to be explicitly added.)

            exe = jd.executable[len('./'):]
            exe_transfer = '%s > %s' % (exe, exe)

            if jd.file_transfer:
                jd.file_transfer.append(exe_transfer)
            else:
                jd.file_transfer = [exe_transfer]

        if jd.file_transfer is not None:
            jd.transfer_directives = TransferDirectives(jd.file_transfer)

            self._handle_file_transfers(jd)

            td = jd.transfer_directives
            if not (hasattr(td, 'transfer_output_files') or hasattr(td, 'transfer_input_files')):
                raise RuntimeError('One of them should be set!')

        # create a Condor job script from SAGA job description
        script = _condorscript_generator(url=self.rm, logger=self._logger,
                jds=[jd], option_dict=self.query_options)
        self._logger.info("Generated Condor script: %s" % script)

        submit_file = NamedTemporaryFile(mode='w', suffix='.condor',
                    prefix='tmp-saga-', delete=False)
        submit_file.write(script)
        submit_file.close()
        self._logger.info("Written Condor script locally: %s" % submit_file.name)

        if self.shell.url.scheme in ["ssh", "gsissh"]:
            self._logger.info("Transferring Condor script to: %s" % self.shell.url)
            submit_file_name = os.path.basename(submit_file.name)
            # TODO: the "-P" is a layer violation?
            self.shell.stage_to_remote(submit_file.name, submit_file_name, cp_flags="-P")

        else:
            raise NotImplementedError("%s support for Condor not implemented." % \
                    self.shell.url.scheme)


        ret, out, _ = self.shell.run_sync('%s -verbose %s' \
            % (self._commands['condor_submit']['path'], submit_file_name))


        if ret != 0:
            # something went wrong
            message = "Error running job via 'condor_submit': %s. Script was: %s" \
                % (out, script)
            log_error_and_raise(message, saga.NoSuccess, self._logger)
        else:
            # stdout contains the job id
            for line in out.split("\n"):
                if "** Proc" in line:
                    pid = line.split()[2][:-1]

            # we don't want the 'query' part of the URL to be part of the ID,
            # simply because it can get terribly long (and ugly). to get rid
            # of it, we clone the URL and set the query part to None.
            rm_clone = surl.Url (self.rm)
            rm_clone.query = ""
            rm_clone.path = ""

            job_id = "[%s]-[%s]" % (rm_clone, pid)
            self._logger.info("Submitted Condor job with id: %s" % job_id)

            # add job to internal list of known jobs.
            self.jobs[job_id] = {
                'state':        saga.job.PENDING,
                'exec_hosts':   None,
                'returncode':   None,
                'create_time':  None,
                'start_time':   None,
                'end_time':     None,
                'gone':         False,
                'transfers':    None,
                'stdout':       None,
                'stderr':       None,
                'name':         None
            }

            # remove submit file(s)
            # XXX: maybe leave them in case of debugging?
            self._logger.info("Submitted Condor job with scheme: '%s'" % self.shell.url.scheme)
            if self.shell.url.scheme in ['ssh', 'gsissh']:
                self._logger.info("1")
                ret, out, _ = self.shell.run_sync ('rm %s' % submit_file_name)
            else:
                self._logger.info("2")
                raise NotImplementedError("%s support for Condor not implemented." % \
                                 self.shell.url.scheme)
            os.remove(submit_file.name)

            return job_id

    # ----------------------------------------------------------------
    #
    def _retrieve_job(self, job_id):
        """ see if we can get some info about a job that we don't
            know anything about
        """
        rm, pid = self._adaptor.parse_id(job_id)

        # run the Condor 'condor_q' command to get some infos about our job
        ret, out, _ = self.shell.run_sync("unset GREP_OPTIONS; %s -long %s | \
            grep -E '(^JobStatus)|(ExitStatus)|(CompletionDate)'" \
            % (self._commands['condor_q']['path'], pid))
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
                'gone':         False,
                'transfers':    None,
                'stdout':       None,
                'stderr':       None,
                'name':         None
            }

            results = out.split('\n')
            for result in results:
                if len(result.split('=')) == 2:
                    key, val = result.split('=')
                    key = key.strip()  # strip() removes whitespaces at the
                    val = val.strip()  # beginning and the end of the string

                    if key == 'JobStatus':
                        job_info['state'] = _condor_to_saga_jobstate(val)
                    elif key == 'ExitStatus':
                        job_info['returncode'] = val
                    elif key == 'CompletionDate':
                        job_info['end_time'] = val

            return job_info

    def _handle_file_transfers(self, jd):

        td = jd.transfer_directives

        # Condor specific safety checks
        if td.in_append_dict:
            raise Exception('FileTransfer append syntax (>>) not supported by Condor: %s' % td.in_append_dict)
        if td.out_append_dict:
            raise Exception('FileTransfer append syntax (<<) not supported by Condor: %s' % td.out_append_dict)

        if td.in_overwrite_dict:
            td.transfer_input_files = []
            for (source, target) in td.in_overwrite_dict.iteritems():
                # make sure source is file an not dir
                (s_path, s_entry) = os.path.split(source)
                if len(s_entry) < 1:
                    raise Exception('Condor accepts only files (not directories) as FileTransfer sources: %s' % source)
                    # TODO: this doesn't seem true, but dont want to tackle it right now
                # make sure source and target file are the same
                (t_path, t_entry) = os.path.split(target)
                if s_entry != t_entry:
                    raise Exception('For Condor source file name and target file name have to be identical: %s != %s' % (s_entry, t_entry))

                # add for later use by job script generator
                td.transfer_input_files.append(target)

                if self.shell.url.scheme in ["ssh", "gsissh"]:
                    self._logger.info("Transferring file %s to %s" % (source, target))
                    self.shell.stage_to_remote(source, target, cp_flags=saga.filesystem.CREATE_PARENTS)

                else:
                    # TODO: this should just work
                    raise NotImplementedError("%s support for Condor not implemented." % \
                            self.shell.url.scheme)

        if td.out_overwrite_dict:
            td.transfer_output_files = []
            for (source, target) in td.out_overwrite_dict.iteritems():
                # make sure source is file and not dir
                (s_path, s_entry) = os.path.split(source)
                if len(s_entry) < 1:
                    raise Exception('Condor accepts only files (not directories) as FileTransfer sources: %s' % source)
                    # TODO: this doesn't seem true, but dont want to tackle it right now
                # make sure source and target file are the same
                (t_path, t_entry) = os.path.split(target)
                if s_entry != t_entry:
                    raise Exception('For Condor source file name and target file name have to be identical: %s != %s' % (s_entry, t_entry))

                # add for later use by job script generator
                td.transfer_output_files.append(target)


    # ----------------------------------------------------------------
    #
    def _job_get_info(self, job_id):
        """ get job attributes via condor_q
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

        # curr. info will contain the new job info collect. it starts off
        # as a copy of prev_info (don't use deepcopy because there is an API 
        # object in the dict -> recursion)
        curr_info = dict()
        #curr_info['job_id'     ] = prev_info.get ('job_id'     )
        curr_info['state'      ] = prev_info.get ('state'      )
        curr_info['exec_hosts' ] = prev_info.get ('exec_hosts' )
        curr_info['returncode' ] = prev_info.get ('returncode' )
        curr_info['create_time'] = prev_info.get ('create_time')
        curr_info['start_time' ] = prev_info.get ('start_time' )
        curr_info['end_time'   ] = prev_info.get ('end_time'   )
        curr_info['gone'       ] = prev_info.get ('gone'       )
        curr_info['transfers'  ] = prev_info.get ('transfers'  )
        curr_info['stdout'     ] = prev_info.get ('stdout'     )
        curr_info['stderr'     ] = prev_info.get ('stderr'     )
        curr_info['name'       ] = prev_info.get ('name'       )

        rm, pid = self._adaptor.parse_id(job_id)

        # run the Condor 'condor_q' command to get some infos about our job
        ret, out, _ = self.shell.run_sync("unset GREP_OPTIONS; %s -long %s | \
            grep -E '(^JobStatus)|(ExitStatus)|(CompletionDate)'" \
            % (self._commands['condor_q']['path'], pid))

        if ret != 0:
            if prev_info['state'] in [saga.job.RUNNING, saga.job.PENDING]:

                # run the Condor 'condor_history' command to get info about 
                # finished jobs
                ret, out, _ = self.shell.run_sync("unset GREP_OPTIONS; %s -long -match 1 %s | \
                    grep -E '(ExitCode)|(TransferOutput)|(CompletionDate)|(JobCurrentStartDate)|(QDate)|(Err)|(Out)'" \
                    % (self._commands['condor_history']['path'], pid))
                
                if ret != 0:
                    message = "Error getting job history via 'condor_history': %s" % out
                    log_error_and_raise(message, saga.NoSuccess, self._logger)

                # parse the egrep result. this should look something like this:
                # ExitCode = 0
                # TransferOutput = "radical.txt"
                results = out.split('\n')
                for result in results:
                    if len(result.split('=')) == 2:
                        key, val = result.split('=')
                        key = key.strip()  # strip() removes whitespaces at the
                        val = val.strip()  # beginning and the end of the string

                        if key == 'ExitCode':
                            curr_info['returncode'] = int(val)
                        elif key == 'TransferOutput':
                            curr_info['transfers'] = val
                        elif key == 'QDate':
                            curr_info['create_time'] = val
                        elif key == 'JobCurrentStartDate':
                            curr_info['start_time'] = val
                        elif key == 'CompletionDate':
                            curr_info['end_time'] = val
                        elif key == 'Out':
                            curr_info['stdout'] = val
                        elif key == 'Err':
                            curr_info['stderr'] = val

                if curr_info['returncode'] == 0:
                    curr_info['state'] = saga.job.DONE
                else:
                    curr_info['state'] = saga.job.FAILED

                curr_info['gone'] = True

            else:
                curr_info['gone'] = True
            #else:
            #    # something went wrong
            #    message = "Error retrieving job info via 'condor_q': %s" % out
            #    log_error_and_raise(message, saga.NoSuccess, self._logger)
        else:
            # parse the egrep result. this should look something like this:
            # JobStatus = 5
            # ExitStatus = 0
            # CompletionDate = 0
            results = out.split('\n')
            for result in results:
                if len(result.split('=')) == 2:
                    key, val = result.split('=')
                    key = key.strip()  # strip() removes whitespaces at the
                    val = val.strip()  # beginning and the end of the string

                    if key == 'JobStatus':
                        curr_info['state'] = _condor_to_saga_jobstate(val)
                    elif key == 'ExitStatus':
                        curr_info['returncode'] = val
                    elif key == 'CompletionDate':
                        curr_info['end_time'] = val

        if curr_info['gone'] is True:
            # If we are running over SSH, copy the output to our local system
            if self.shell.url.scheme == "ssh":
                files = []

                if curr_info['transfers']:
                    t = curr_info['transfers']
                    self._logger.debug("TransferOutput: %s" % t)

                    # Remove leading and ending double quotes
                    if t.startswith('"') and t.endswith('"'):
                        t = t[1:-1]

                    # Parse comma separated list
                    files += t.split(',')

                if curr_info['stdout']:
                    t = curr_info['stdout']
                    self._logger.debug("StdOut: %s" % t)

                    # Remove leading and ending double quotes
                    if t.startswith('"') and t.endswith('"'):
                        t = t[1:-1]

                    files.append(t)

                if curr_info['stderr']:
                    t = curr_info['stderr']
                    self._logger.debug("StdErr: %s" % t)

                    # Remove leading and ending double quotes
                    if t.startswith('"') and t.endswith('"'):
                        t = t[1:-1]

                    files.append(t)

                # Transfer list of files
                for f in files:
                    f = f.strip()
                    self._logger.info("Transferring file %s" % f)
                    self.shell.stage_from_remote(f, f)

        # return the new job info dict
        return curr_info


    # ----------------------------------------------------------------
    #
    def _job_get_info_bulk(self, cluster_id, cluster_size, rm):
        """ get job attributes via condor_q
        """

        prev_info = {}
        curr_info = {}

        proc_ids = map(str, range(cluster_size))

        for job_id in ['[%s]-[%s.%s]' % (rm, cluster_id, proc_id) for proc_id in proc_ids]:

            # if we don't have the job in our dictionary, we don't want it
            if job_id not in self.jobs:
                raise Exception("Unknown job ID: %s. Can't update state.", job_id)

            # prev. info contains the info collect when _job_get_info(_bulk)
            # was called the last time
            prev_info[job_id] = self.jobs[job_id]

            # if the 'gone' flag is set, there's no need to query the job
            # state again. it's gone forever
            # TODO: what and when to return?
            #if prev_info[job_id]['gone'] is True:
            #    self._logger.warning("Job information is not available anymore.")
            #    return prev_info[job_id]

            # curr. info will contain the new job info collect. it starts off
            # as a copy of prev_info (don't use deepcopy because there is an API
            # object in the dict -> recursion)
            curr_info[job_id] = {}
            #curr_info['job_id'] = prev_info.get ('job_id')
            curr_info[job_id]['state']       = prev_info[job_id].get('state')
            curr_info[job_id]['exec_hosts']  = prev_info[job_id].get('exec_hosts')
            curr_info[job_id]['returncode' ] = prev_info[job_id].get('returncode')
            curr_info[job_id]['create_time'] = prev_info[job_id].get('create_time')
            curr_info[job_id]['start_time']  = prev_info[job_id].get('start_time')
            curr_info[job_id]['end_time']    = prev_info[job_id].get('end_time')
            curr_info[job_id]['gone']        = prev_info[job_id].get('gone')
            curr_info[job_id]['transfers']   = prev_info[job_id].get('transfers')
            curr_info[job_id]['stdout']      = prev_info[job_id].get('stdout')
            curr_info[job_id]['stderr']      = prev_info[job_id].get('stderr')
            curr_info[job_id]['name']        = prev_info[job_id].get('name')


        # run the Condor 'condor_q' command to get some infos about our job
        ret, out, _ = self.shell.run_sync(
            "%s %s -autoformat:,v ProcId JobStatus ExitStatus CompletionDate" %
            (self._commands['condor_q']['path'], cluster_id))

        if ret != 0:
            raise Exception("condor_q failed")

        results = filter(None, out.split('\n'))
        procs_condor_q_found = set()

        for row in results:
            #
            # Some processes in cluster found with condor_q!
            #
            procid, jobstatus, exitstatus, completiondate = tuple([col.strip() for col in row.split(',')])

            job_id = "[%s]-[%s.%s]" % (rm, cluster_id, procid)
            procs_condor_q_found.update([procid])
            curr_info[job_id]['state'] = _condor_to_saga_jobstate(jobstatus)
            curr_info[job_id]['returncode'] = exitstatus
            curr_info[job_id]['end_time'] = completiondate

        procs_missing = set(proc_ids) - procs_condor_q_found

        if procs_missing:
            #
            # (Some) cluster processes not found with condor_q, trying condor_history now
            #
            #self._logger.debug('condor_q only returned (%s) of (%s) processes, missing (%s)', procs_condor_q_found, proc_ids, procs_missing)

            # TODO: Don't know why one would use this check ..., skipping for now
            # if prev_info['state'] in [saga.job.RUNNING, saga.job.PENDING]:

            # run the Condor 'condor_history' command to get info about finished jobs
            ret, out, _ = self.shell.run_sync(
                "%s %s -match %d -autoformat:, ProcId ExitCode TransferOutput CompletionDate JobCurrentStartDate QDate Err Out" %
                (self._commands['condor_history']['path'], cluster_id, len(procs_missing)))

            if ret != 0:
                raise Exception("Error getting job history via 'condor_history': %s" % out)

            procs_condor_history_found = set()

            results = filter(None, out.split('\n'))
            for row in results:

                procid, exitcode, transferoutput, completiondate, jobcurrentstartdate, qdate, stderr, stdout = \
                    tuple([col.strip() for col in row.split(',')])

                procs_condor_history_found.update([procid])

                job_id = "[%s]-[%s.%s]" % (rm, cluster_id, procid)

                curr_info[job_id]['returncode'] = int(exitcode)
                curr_info[job_id]['transfers'] = transferoutput
                curr_info[job_id]['create_time'] = qdate
                curr_info[job_id]['start_time'] = jobcurrentstartdate
                curr_info[job_id]['end_time'] = completiondate
                curr_info[job_id]['stdout'] = stdout
                curr_info[job_id]['stderr'] = stderr

                if curr_info[job_id]['returncode'] == 0:
                    curr_info[job_id]['state'] = saga.job.DONE
                else:
                    curr_info[job_id]['state'] = saga.job.FAILED

                curr_info[job_id]['gone'] = True

            procs_missing = set(proc_ids) - (procs_condor_q_found | procs_condor_history_found)
            if procs_missing:
                self._logger.debug('still missing procs: %s' % procs_missing)

        # Transfer stuff
        # if curr_info['gone'] is True:
        #     # If we are running over SSH, copy the output to our local system
        #     if self.shell.url.scheme == "ssh":
        #         files = []
        #
        #         if curr_info['transfers']:
        #             t = curr_info['transfers']
        #             self._logger.debug("TransferOutput: %s" % t)
        #
        #             # Remove leading and ending double quotes
        #             if t.startswith('"') and t.endswith('"'):
        #                 t = t[1:-1]
        #
        #             # Parse comma separated list
        #             files += t.split(',')
        #
        #         if curr_info['stdout']:
        #             t = curr_info['stdout']
        #             self._logger.debug("StdOut: %s" % t)
        #
        #             # Remove leading and ending double quotes
        #             if t.startswith('"') and t.endswith('"'):
        #                 t = t[1:-1]
        #
        #             files.append(t)
        #
        #         if curr_info['stderr']:
        #             t = curr_info['stderr']
        #             self._logger.debug("StdErr: %s" % t)
        #
        #             # Remove leading and ending double quotes
        #             if t.startswith('"') and t.endswith('"'):
        #                 t = t[1:-1]
        #
        #             files.append(t)
        #
        #         # Transfer list of files
        #         for f in files:
        #             f = f.strip()
        #             self._logger.info("Transferring file %s" % f)
        #             self.shell.stage_from_remote(f, f)

        # return the new job info dict
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

        ret = self.jobs[job_id]['returncode']

        # FIXME: 'None' should cause an exception
        if ret == None : return None
        else           : return int(ret)

    # ----------------------------------------------------------------
    #
    def _job_get_name(self, job_id):
        """ get the job's name
        """
        return self.jobs[job_id]['name']

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
        """ cancel the job via 'condor_rm'
        """
        rm, pid = self._adaptor.parse_id(job_id)

        ret, out, _ = self.shell.run_sync("%s %s\n" \
            % (self._commands['condor_rm']['path'], pid))

        if ret != 0:
            message = "Error canceling job via 'condor_rm': %s" % out
            log_error_and_raise(message, saga.NoSuccess, self._logger)

        # assume the job was successfully canceled
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

    # ----------------------------------------------------------------
    #
    @SYNC_CALL
    def create_job(self, jd):
        """ implements saga.adaptors.cpi.job.Service.get_url()
        """
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

        ret, out, _ = self.shell.run_sync("unset GREP_OPTIONS; %s | grep `whoami`"\
            % self._commands['condor_q']['path'])

        if ret != 0 and len(out) > 0:
            message = "failed to list jobs via 'condor_q': %s" % out
            log_error_and_raise(message, saga.NoSuccess, self._logger)
        elif ret != 0 and len(out) == 0:
            pass
        else:
            for line in out.split("\n"):
                # output looks like this:
                # 112059.svc.uc.futuregrid testjob oweidner 0 Q batch
                # 112061.svc.uc.futuregrid testjob oweidner 0 Q batch
                if len(line.split()) > 1:
                    rm_clone = surl.Url (self.rm)
                    rm_clone.query = ""
                    rm_clone.path = ""

                    jobid = "[%s]-[%s]" % (rm_clone, line.split()[0])
                    ids.append(str(jobid))

        return ids

    # ----------------------------------------------------------------
    #
    @SYNC_CALL
    def container_run(self, jobs):

        # TODO: Do we need to check that the bulks have at least something in common?
        # How about only allowing arguments to be different for the time being?
        # How about working_directory?

        executables = set([job.description.executable for job in jobs])
        if len(executables) != 1:
            raise Exception('Non-unique executables: %s' % executables)
        else:
            self._logger.info('Unique executable: %s' % executables)

        for job in jobs:

            jd = job.description

            # Because we do funky shit with env and sh, we need to explicitly add
            # the executable to the transfer list.
            # (Of course not if the executable is already on the target systems,
            # defined by the fact that it starts with ./
            if jd.executable.startswith('./') and \
                os.path.exists(jd.executable): # TODO: Why does it need to exist at this very path?

                # TODO: Check if the executable is already in the file_transfer list,
                # because then we don't need to implicitly add it anymore.
                # (For example, if the executable is not in the local directory, it needs
                # to be explicitly added.)

                exe = jd.executable[len('./'):]
                exe_transfer = '%s > %s' % (exe, exe)

                if jd.file_transfer:
                    jd.file_transfer.append(exe_transfer)
                else:
                    jd.file_transfer = [exe_transfer]

            if jd.file_transfer is not None:
                jd.transfer_directives = TransferDirectives(jd.file_transfer)

                # TODO: Given that input (and output) are likely similar for bulk tasks,
                # TODO: we probably don't want to transfer duplicates endlessly
                self._handle_file_transfers(jd)

                td = jd.transfer_directives
                if not (hasattr(td, 'transfer_output_files') or hasattr(td, 'transfer_input_files')):
                    raise RuntimeError('One of them should be set!')


        # create a Condor job script from SAGA job description
        script = _condorscript_generator(
            url=self.rm, logger=self._logger,
            jds=[job.description for job in jobs],
            option_dict=self.query_options)
        self._logger.info("Generated Condor script: %s" % script)

        submit_file = NamedTemporaryFile(mode='w', suffix='.condor',
                                         prefix='tmp-saga-', delete=False)
        submit_file.write(script)
        submit_file.close()
        self._logger.info("Written Condor script locally: %s" % submit_file.name)

        if self.shell.url.scheme in ["ssh", "gsissh"]:
            self._logger.info("Transferring Condor script to: %s" % self.shell.url)
            submit_file_name = os.path.basename(submit_file.name)
            # TODO: the "-P" is a layer violation?
            self.shell.stage_to_remote(submit_file.name, submit_file_name, cp_flags="-P")

        else:
            raise NotImplementedError("%s support for Condor not implemented." % \
                    self.shell.url.scheme)

        submit_file_name = submit_file.name

        ret, out, _ = self.shell.run_sync('%s -verbose %s' \
                                          % (self._commands['condor_submit']['path'], submit_file_name))

        if ret != 0:
            # condor_submit went wrong
            message = "Error running job via 'condor_submit': %s. Script was: %s" \
                      % (out, script)
            log_error_and_raise(message, saga.NoSuccess, self._logger)

        pids = []
        # stdout contains the job id
        for line in out.split("\n"):
            if line.startswith("** Proc "):
                pid = line.split()[2][:-1]
                if len(pid.split('.')) != 2:
                    raise Exception("Pid looks weird: %s" % pid)
                pids.append(pid)

        if len(pids) != len(jobs):
            raise Exception("Number of pids (%d) is different than number of jobs (%d)" % (len(pids), len(jobs)))

        # we don't want the 'query' part of the URL to be part of the ID,
        # simply because it can get terribly long (and ugly). to get rid
        # of it, we clone the URL and set the query part to None.
        rm_clone = surl.Url(self.rm)
        rm_clone.query = ""
        rm_clone.path = ""

        for job, pid in zip(jobs, pids):
            job_id = "[%s]-[%s]" % (rm_clone, pid)

            job._adaptor._id = job_id
            job._adaptor._started = True

            self._logger.info("Submitted Condor job with id: %s" % job_id)

            # add job to internal list of known jobs.
            self.jobs[job_id] = {
                'state':        saga.job.PENDING,
                'exec_hosts':   None,
                'returncode':   None,
                'create_time':  None,
                'start_time':   None,
                'end_time':     None,
                'gone':         False,
                'transfers':    None,
                'stdout':       None,
                'stderr':       None,
                'name':         None
            }

        # remove submit file(s)
        # XXX: maybe leave them in case of debugging?
        if self.shell.url.scheme in ['ssh', 'gsissh']:
            #ret, out, _ = self.shell.run_sync ('rm %s' % submit_file_name)
            pass
        else:
            raise NotImplementedError("%s support for Condor not implemented." % \
                    self.shell.url.scheme)

        os.remove(submit_file.name)

    # ----------------------------------------------------------------
    #
    @SYNC_CALL
    def container_wait(self, jobs, mode, timeout):

        if mode == saga.ANY:
            pass
        if mode == saga.ALL:
            pass

        while True:
            states = self.container_get_states(jobs)
            if saga.job.RUNNING not in states and saga.job.PENDING not in states:
                break
            time.sleep(5)

        return states[0]


    # ----------------------------------------------------------------
    #
    @SYNC_CALL
    def container_cancel(self, jobs):
        self._logger.debug("container cancel: %s" % str(jobs))

        # TODO: this is not optimized yet
        for job in jobs:
            job.cancel()

    # ----------------------------------------------------------------
    #
    @SYNC_CALL
    def container_get_states(self, jobs):

        # JobIds also include the rm, so we strip that out.
        job_ids = [job._adaptor._id for job in jobs]

        proc_ids = [self._adaptor.parse_id(job_id)[1] for job_id in job_ids]

        clusters = list(set([id.split('.', 1)[0] for id in proc_ids]))

        if len(clusters) > 1:
            raise Exception("More than one cluster in this container.")

        # Now that we assume we have just one ...
        cluster_id = clusters[0]
        cluster_size = len(jobs)
        self._logger.debug("wait for cluster: %s(%d)", cluster_id, cluster_size)

        rm, _ = self._adaptor.parse_id(job_ids[0])
        bulk_info = self._job_get_info_bulk(cluster_id, cluster_size, rm)

        states = [bulk_info[job_id]['state'] for job_id in job_ids]

        return states

        # TODO: check "cache" for final state jobs
        # check if we have already reach a terminal state
        # if self.jobs[job_id]['state'] == saga.job.CANCELED \
        #         or self.jobs[job_id]['state'] == saga.job.FAILED \
        #         or self.jobs[job_id]['state'] == saga.job.DONE:
        #     return self.jobs[job_id]['state']
        #
        # # check if we can / should update
        # if (self.jobs[job_id]['gone'] is not True):
        #     self.jobs[job_id] = self._job_get_info(job_id=job_id)
        #
        # return self.jobs[job_id]['state']

###############################################################################
#
class CondorJob (saga.adaptors.cpi.job.Job):
    """ implements saga.adaptors.cpi.job.Job
    """

    def __init__(self, api, adaptor):

        # initialize parent class
        _cpi_base = super(CondorJob, self)
        _cpi_base.__init__(api, adaptor)

    @SYNC_CALL
    def init_instance(self, job_info):
        """ implements saga.adaptors.cpi.job.Job.init_instance()
        """
        # init_instance is called for every new saga.job.Job object
        # that is created
        self.jd = job_info["job_description"]
        self.js = job_info["job_service"]

        # the js is responsible for job bulk operations -- which
        # for jobs only work for run()
        self._container = self.js

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
    def get_description(self):
        return self.jd

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
    def get_name(self):
        """ implements saga.adaptors.cpi.job.Job.get_name()
        """
        return self.js._job_get_name(self._id)

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

