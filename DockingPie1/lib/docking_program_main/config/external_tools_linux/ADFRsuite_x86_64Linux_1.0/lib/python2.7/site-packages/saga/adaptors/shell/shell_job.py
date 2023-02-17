
__author__    = "Andre Merzky, Ole Weidner"
__copyright__ = "Copyright 2012-2013, The SAGA Project"
__license__   = "MIT"


""" shell based job adaptor implementation """

import saga.utils.pty_shell

import saga.adaptors.base
import saga.adaptors.cpi.job

from   saga.job.constants import *
from   saga.utils.job     import TransferDirectives

import re
import time
import threading

import shell_wrapper

SYNC_CALL  = saga.adaptors.cpi.decorators.SYNC_CALL
ASYNC_CALL = saga.adaptors.cpi.decorators.ASYNC_CALL


# ------------------------------------------------------------------------------
#
class _job_state_monitor (threading.Thread) :
    """ thread that periodically monitors job states """

    # --------------------------------------------------------------------------
    #
    def __init__ (self, js, channel, rm, logger) :

        self.js      = js
        self.channel = channel
        self.rm      = rm 
        self.logger  = logger
        self.stop    = False
        self.events  = dict()

        super (_job_state_monitor, self).__init__ ()

        self.setDaemon (True)


    # --------------------------------------------------------------------------
    #
    def finalize (self) :

        self.stop = True


    # --------------------------------------------------------------------------
    #
    def run (self) :

        MONITOR_READ_TIMEOUT = 1.0   # check for stop signal now and then

        try:

            self.channel.run_async ("MONITOR")

            while self.channel.alive () :

                idx, out = self.channel.find (['\n'], timeout=MONITOR_READ_TIMEOUT)

                line = out.strip ()

                if  not line :

                    # just a read timeout, i.e. an opportiunity to check for
                    if  self.stop :
                        self.logger.debug ("stop monitoring")
                        return
                    pass


                elif line == 'EXIT' or line == "Killed" :
                    self.logger.error ("monitoring channel failed -- disable notifications")
                    return


                elif not ':' in line :
                    self.logger.warn ("monitoring channel noise: %s" % line)


                else :
                    job_pid, state, data = line.split (':', 2)
                    job_id = "[%s]-[%s]" % (self.rm, job_pid)

                    state = self.js._adaptor.string_to_state (state)

                    try :
                        job = self.js.get_job (job_id, no_reconnect=True)

                        if  not job :
                            # job not yet known -- keep event for later
                            if  not job_id in self.events :
                                self.events[job_id] = list()
                            self.events[job_id].append (state)

                        else :

                            # check for previous events :
                            if  job_id in self.events :
                                for event in self.events[job_id] :
                                    job._adaptor._set_state (event)
                                del (self.events[job_id])
                            job._adaptor._set_state (state)


                    except saga.DoesNotExist as e :
                        self.logger.error ("event for unknown job '%s'" % job_id)


        except Exception as e:

            self.logger.error ("Exception in job monitoring thread: %s" % e)
            self.logger.error ("Cancel job monitoring for %s" % self.rm)


# --------------------------------------------------------------------
#
# strip white space from a string, and hex-decode the remaining characters.
# This must be applied to stdout/stderr data returned from the shell wrapper.
#
def _decode (data) :

    elem = ""
    code = ""
    
    for c in data :
    
        if  c in ['0','1','2','3','4','5','6','7','8','9','a','b','c','d','e','f'] :
            code += c
        elif c not in [' ', '\n'] :
            raise saga.BadParameter("Cannot decode '%s' in data (%s)" % (c, data))

    return code.decode ('hex')


# --------------------------------------------------------------------
# the adaptor name
#
_ADAPTOR_NAME          = "saga.adaptor.shell_job"
_ADAPTOR_SCHEMAS       = ["fork", "local", "ssh", "gsissh"]
_ADAPTOR_OPTIONS       = [
    { 
    'category'         : 'saga.adaptor.shell_job',
    'name'             : 'enable_notifications', 
    'type'             : bool, 
    'default'          : False,
    'valid_options'    : [True, False],
    'documentation'    : '''Enable support for job state notifications.  Note that
                          enabling this option will create a local thread, a remote 
                          shell process, and an additional network connection.
                          In particular for ssh/gsissh where the number of
                          concurrent connections is limited to 10, this
                          effectively halfs the number of available job service
                          instances per remote host.''',
    'env_variable'     : None
    },
    {
    'category'         : 'saga.adaptor.shell_job',
    'name'             : 'purge_on_start',
    'type'             : bool,
    'default'          : True,
    'valid_options'    : [True, False],
    'documentation'    : '''Purge job information (state, stdio, ...) for all
                          jobs which are in final state when starting the job
                          service instance. Note that this will purge *all*
                          suitable jobs, including the ones managed by another,
                          live job service instance.''',
    'env_variable'     : None
    },
    {
    'category'         : 'saga.adaptor.shell_job',
    'name'             : 'base_workdir',
    'type'             : str,
    'default'          : "$HOME/.saga/adaptors/shell_job/",
    'documentation'    : '''The adaptor stores job state information on the
                          filesystem on the target resource.  This parameter
                          specified what location should be used.''',
    'env_variable'     : None
    }
]

# --------------------------------------------------------------------
# the adaptor capabilities & supported attributes
#
_ADAPTOR_CAPABILITIES  = {
    "jdes_attributes"  : [saga.job.NAME,
                          saga.job.EXECUTABLE,
                          saga.job.PRE_EXEC,
                          saga.job.POST_EXEC,
                          saga.job.ARGUMENTS,
                          saga.job.ENVIRONMENT,
                          saga.job.WORKING_DIRECTORY,
                          saga.job.FILE_TRANSFER,
                          saga.job.INPUT,
                          saga.job.OUTPUT,
                          saga.job.ERROR,
                          saga.job.NAME,
                          saga.job.WALL_TIME_LIMIT, # TODO: 'hot'-fix for BigJob - implement properly
                          saga.job.TOTAL_CPU_COUNT, # TODO: 'hot'-fix for BigJob - implement properly
                          saga.job.PROCESSES_PER_HOST,
                          saga.job.SPMD_VARIATION,  # TODO: 'hot'-fix for BigJob - implement properly
                         ],
    "job_attributes"   : [saga.job.EXIT_CODE,
                          saga.job.EXECUTION_HOSTS,
                          saga.job.CREATED,
                          saga.job.STARTED,
                          saga.job.FINISHED],
    "metrics"          : [saga.job.STATE, 
                          saga.job.STATE_DETAIL],
    "contexts"         : {"ssh"      : "public/private keypair",
                          "x509"     : "X509 proxy for gsissh",
                          "userpass" : "username/password pair for ssh"}
}

# --------------------------------------------------------------------
# the adaptor documentation
#
_ADAPTOR_DOC           = {
    "name"             : _ADAPTOR_NAME,
    "cfg_options"      : _ADAPTOR_OPTIONS, 
    "capabilities"     : _ADAPTOR_CAPABILITIES,
    "description"      : """ 
        The Shell job adaptor. This adaptor uses the sh command line tools (sh,
        ssh, gsissh) to run local and remote jobs.  The adaptor expects the
        login shell on the target host to be POSIX compliant.  However, one can
        also specify a custom POSIX shell via the resource manager URL, like::

          js = saga.job.Service ("ssh://remote.host.net/bin/sh")

        Note that custom shells in many cases will find a different environment
        than the users default login shell!


        Known Limitations:
        ******************

          * number of system pty's are limited:  each job.service object bound
            to this adaptor will use 2 pairs of pty pipes.  Systems usually
            limit the number of available pty's to 1024 .. 4096.  Given that
            other process also use pty's , that gives a hard limit to the number
            of object instances which can be created concurrently.  Hitting the
            pty limit will cause the following error message (or similar)::

              NoSuccess: pty_allocation or process creation failed (ENOENT: no more ptys)

            This limitation comes from saga.utils.pty_process.  On Linux
            systems, the utilization of pty's can be monitored::

               echo "allocated pty's: `cat /proc/sys/kernel/pty/nr`"
               echo "available pty's: `cat /proc/sys/kernel/pty/max`"


          * number of ssh connections are limited: sshd's default configuration,
            which is in place on many systems, limits the number of concurrent
            ssh connections to 10 per user -- beyond that, connections are
            refused with the following error::

              NoSuccess: ssh_exchange_identification: Connection closed by remote host

            As the communication with the ssh channel is unbuffered, the
            dropping of the connection will likely cause this error message to
            be lost.  Instead, the adaptor will just see that the ssh connection
            disappeared, and will issue an error message similar to this one::

              NoSuccess: read from pty process failed (Could not read line - pty process died)

 
          * number of processes are limited: the creation of an job.service
            object will create one additional process on the local system, and
            two processes on the remote system (ssh daemon clone and a shell
            instance).  Each remote job will create three additional processes:
            two for the job instance itself (double fork), and an additional
            process which monitors the job for state changes etc.  Additional
            temporary processes may be needed as well.  

            While marked as 'obsolete' by POSIX, the `ulimit` command is
            available on many systems, and reports the number of processes
            available per user (`ulimit -u`)

            On hitting process limits, the job creation will fail with an error
            similar to either of these::
            
              NoSuccess: failed to run job (/bin/sh: fork: retry: Resource temporarily unavailable)
              NoSuccess: failed to run job -- backend error

          * number of files are limited, as is disk space: the job.service will
            
            keep job state on the remote disk, in ``~/.saga/adaptors/shell_job/``.
            Quota limitations may limit the number of files created there,
            and/or the total size of that directory.  

            On quota or disk space limits, you may see error messages similar to
            the following ones::
            
              NoSuccess: read from pty process failed ([Errno 5] Quota exceeded)
              NoSuccess: read from pty process failed ([Errno 5] Input/output error)
              NoSuccess: find from pty process [Thread-5] failed (Could not read - pty process died)
              


          * Other system limits (memory, CPU, selinux, accounting etc.) apply as
            usual.

          
          * thread safety: it is safe to create multiple :class:`job.Service`
            instances to the same target host at a time -- they should not
            interfere with each other, but ``list()`` will list jobs created by
            either instance (if those use the same target host user account).

            It is **not** safe to use the *same* :class:`job.Service` instance
            from multiple threads concurrently -- the communication on the I/O
            channel will likely get screwed up.  This limitation may be removed
            in future versions of the adaptor.  Non-concurrent (i.e. serialized)
            use should work as expected though.

        """,
    "example": "examples/jobs/localjob.py",
    "schemas"          : {"fork"   :"use /bin/sh to run jobs", 
                          "local"  :"alias for fork://", 
                          "ssh"    :"use ssh to run remote jobs", 
                          "gsissh" :"use gsissh to run remote jobs"}
}

# --------------------------------------------------------------------
# the adaptor info is used to register the adaptor with SAGA

_ADAPTOR_INFO          = {
    "name"             : _ADAPTOR_NAME,
    "version"          : "v0.1",
    "schemas"          : _ADAPTOR_SCHEMAS,
    "capabilities"     : _ADAPTOR_CAPABILITIES,
    "cpis"             : [
        { 
        "type"         : "saga.job.Service",
        "class"        : "ShellJobService"
        }, 
        { 
        "type"         : "saga.job.Job",
        "class"        : "ShellJob"
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
        self.opts  = self.get_config (_ADAPTOR_NAME)

        self.notifications  = self.opts['enable_notifications'].get_value ()
        self.purge_on_start = self.opts['purge_on_start'      ].get_value ()
        self.base_workdir   = self.opts['base_workdir'        ].get_value ()


    # ----------------------------------------------------------------
    #
    def sanity_check (self) :

        # FIXME: also check for gsissh

        pass


    # ----------------------------------------------------------------
    #
    def parse_id (self, id) :
        """
        Split the id '[rm]-[pid]' in its parts, and return them.

        The callee makes sure that the ID is set and valid.
        """

        match = self.id_re.match (id)

        if  not match or len (match.groups()) != 2 :
            raise saga.BadParameter ("Cannot parse job id '%s'" % id)

        return (match.group(1), match.group (2))


    # ----------------------------------------------------------------
    #
    def string_to_state (self, state_str) :

        state_str = state_str.strip ().lower ()

        return {'new'       : saga.job.NEW,
                'running'   : saga.job.RUNNING,
                'suspended' : saga.job.SUSPENDED,
                'done'      : saga.job.DONE,
                'failed'    : saga.job.FAILED,
                'canceled'  : saga.job.CANCELED
               }.get (state_str, saga.job.UNKNOWN)


    # ----------------------------------------------------------------
    #
    def stage_input (self, shell, jd) :

        if not jd:
            return

        if  jd.file_transfer is not None:
            td = TransferDirectives (jd.file_transfer)

            if  len (td.in_append_dict)  > 0 or \
                len (td.out_append_dict) > 0 :
                raise saga.BadParameter('FileTransfer append (<</>>) not supported')

            if  td.in_overwrite_dict :
                for (source, target) in td.in_overwrite_dict.iteritems():
                    self._logger.info("Transferring file %s to %s" % (source, target))
                    shell.stage_to_remote(source, target)


    # ----------------------------------------------------------------
    #
    def stage_output (self, shell, jd) :

        if not jd:
            return

        if  jd.file_transfer is not None:
            td = TransferDirectives (jd.file_transfer)

            if  len (td.out_append_dict) > 0 :
                raise saga.BadParameter('FileTransfer append (<</>>) not supported')

            if  td.out_overwrite_dict :
                for (source, target) in td.out_overwrite_dict.iteritems():
                    self._logger.info("Transferring file %s to %s" % (source, target))
                    shell.stage_from_remote(source, target)


###############################################################################
#
class ShellJobService (saga.adaptors.cpi.job.Service) :
    """ Implements saga.adaptors.cpi.job.Service """

    # ----------------------------------------------------------------
    #
    def __init__ (self, api, adaptor) :

        _cpi_base = super  (ShellJobService, self)
        _cpi_base.__init__ (api, adaptor)

        self.opts = {}
        self.opts['shell'] = None  # default to login shell


    # ----------------------------------------------------------------
    #
    def __del__ (self) :

        try :
            self.close()

        except Exception as e :
          # print str(e)
            pass
    


    # ----------------------------------------------------------------
    #
    @SYNC_CALL
    def init_instance (self, adaptor_state, rm_url, session) :
        """ Service instance constructor """

        self.rm      = rm_url
        self.session = session
        self.jobs    = dict()
        self.njobs   = 0

        # Use `_set_session` method of the base class to set the session object.
        # `_set_session` and `get_session` methods are provided by `CPIBase`.
        self._set_session(session)

        # if the rm URL specifies a path, we interprete that as shell to run.
        # Otherwise, we default to running /bin/sh (for fork) or the user's
        # login shell (for ssh etc).
        if  self.rm.path and self.rm.path != '/' and self.rm.path != '.' :
            self.opts['shell'] = self.rm.path

        # create and initialize connection for starting jobs
        self.shell   = saga.utils.pty_shell.PTYShell (self.rm, self.get_session(), 
                                                      self._logger, opts=self.opts)
        self.channel = saga.utils.pty_shell.PTYShell (self.rm, self.get_session(), 
                                                      self._logger, opts=self.opts)
        self.initialize ()

        # the monitoring thread - one per service instance.  We wait for
        # initialize to finish to make sure that the shell_wrapper is set
        # up...
        self.monitor = _job_state_monitor (js      = self,
                                           channel = self.channel, 
                                           rm      = self.rm, 
                                           logger  = self._logger)
        self.monitor.start ()

        return self.get_api ()


    # ----------------------------------------------------------------
    #
    def close (self) :

        if  self.shell : 
            # FIXME: not sure if we should PURGE here -- that removes states which
            # might not be evaluated, yet.  Should we mark state evaluation
            # separately? 
            #   cmd_state () { touch $DIR/purgeable; ... }
            # When should that be done?

         #  self.shell.run_sync ("PURGE", iomode=None)
            self.shell.run_async("QUIT")
            self.shell.finalize(kill_pty=True)
            self.shell = None

        if  self.monitor : 
            self.monitor.finalize()
            # we don't care about join, really


    # ----------------------------------------------------------------
    #
    def initialize (self) :

        # start the shell, find its prompt.  If that is up and running, we can
        # bootstrap our wrapper script, and then run jobs etc.

        # -- now stage the shell wrapper script, and run it.  Once that is up
        # and running, we can requests job start / management operations via its
        # stdio.

        base = self._adaptor.base_workdir

        ret, out, _ = self.shell.run_sync (" mkdir -p %s" % base)
        if  ret != 0 :
            raise saga.NoSuccess ("host setup failed (%s): (%s)" % (ret, out))

        # TODO: replace some constants in the script with values from config
        # files, such as 'timeout' or 'purge_on_quit' ...
        tgt = "%s/wrapper.sh" % base

        # lets check if we actually need to stage the wrapper script.  We need
        # an adaptor lock on this one.
        with self._adaptor._lock :

            ret, out, _ = self.shell.run_sync (" test -f %s" % tgt)
            if  ret != 0 :
                # yep, need to stage...
              # src = shell_wrapper._WRAPPER_SCRIPT % ({ 'PURGE_ON_START' : str(self._adaptor.purge_on_start) })
                src = shell_wrapper._WRAPPER_SCRIPT
                src = src.replace('%(PURGE_ON_START)s', str(self._adaptor.purge_on_start))
                
                # If the target directory begins with $HOME or ${HOME} then we
                # need to remove this since scp won't expand the variable and
                # the copy will end up attempting to copy the file to 
                # /<path_to_home_dir>/$HOME/.....
                if tgt.startswith("$HOME") or tgt.startswith("${HOME}"):
                    tgt = tgt[tgt.find('/')+1:]
                self.shell.write_to_remote (src, tgt)

        # ----------------------------------------------------------------------
        # we run the script.  In principle, we should set a new / different
        # prompt -- but, due to some strange and very unlikely coincidence, the
        # script has the same prompt as the previous shell... - go figure ;-)
        #
        # Note that we use 'exec' - so the script replaces the shell process.
        # Thus, when the script times out, the shell dies and the connection
        # drops -- that will free all associated resources, and allows for
        # a clean reconnect.
        # ret, out, _ = self.shell.run_sync (" exec sh %s/wrapper.sh" % base)
      
        # Well, actually, we do not use exec, as that does not give us good
        # feedback on failures (the shell just quits) -- so we replace it with
        # this poor-man's version...
        ret, out, _ = self.shell.run_sync (" /bin/sh %s/wrapper.sh %s" % (base,
            base))

        # shell_wrapper.sh will report its own PID -- we use that to sync prompt
        # detection, too.
        if  ret != 0 :
            raise saga.NoSuccess ("failed to run bootstrap: (%s)(%s)" % (ret, out))

        id_pattern = re.compile ("\s*PID:\s+(\d+)\s*$")
        id_match   = id_pattern.search (out)

        if  not id_match :
            self.shell.run_async (" exit")
            self._logger.error   ("host bootstrap failed - no pid (%s)" % out)
            raise saga.NoSuccess ("host bootstrap failed - no pid (%s)" % out)

        # we actually don't care much about the PID :-P
        
        self._logger.debug ("got cmd prompt (%s)(%s)" % (ret, out.strip ()))


        # ----------------------------------------------------------------------
        # now do the same for the monitoring shell
        ret, out, _ = self.channel.run_sync (" /bin/sh %s/wrapper.sh %s" % (base, base))

        # shell_wrapper.sh will report its own PID -- we use that to sync prompt
        # detection, too.
        if  ret != 0 :
            raise saga.NoSuccess ("failed to run bootstrap: (%s)(%s)" % (ret, out))

        id_pattern = re.compile ("\s*PID:\s+(\d+)\s*$")
        id_match   = id_pattern.search (out)

        if  not id_match :
            self.channel.run_async (" exit")
            self._logger.error     ("host bootstrap failed - no pid (%s)" % out)
            raise saga.NoSuccess   ("host bootstrap failed - no pid (%s)" % out)

        # we actually don't care much about the PID :-P
        
        self._logger.debug ("got mon prompt (%s)(%s)" % (ret, out.strip ()))


    # ----------------------------------------------------------------
    #
    #
    def _jd2cmd (self, jd) :

        cmd = "true"

        if  jd.attribute_exists (ENVIRONMENT) :
            for e in jd.environment :
                cmd += " && export %s=%s"  %  (e, jd.environment[e])

        if  jd.attribute_exists (WORKING_DIRECTORY) :
            cmd += " && mkdir -p %s && cd %s" % (jd.working_directory, jd.working_directory)

        if  jd.attribute_exists (PRE_EXEC) :
            for p in jd.pre_exec :
                cmd += " && %s 2>&1 >> $SAGA_PWD/log"  %  p

        cmd += " && ("
        cmd += " %s" % jd.executable

        if  jd.attribute_exists (ARGUMENTS) :
            for a in jd.arguments :
                cmd += " %s" % a

        cmd += " )"

        if  jd.attribute_exists (INPUT) :
            cmd += " <%s" % jd.input

        if  jd.attribute_exists (OUTPUT) :
            cmd += " 1>%s" % jd.output

        if  jd.attribute_exists (ERROR) :
            cmd += " 2>%s" % jd.error

        if  jd.attribute_exists (POST_EXEC) :
            for p in jd.post_exec :
                cmd += " && %s 2>&1 >> $SAGA_PWD/log"  %  p

        return cmd

    # ----------------------------------------------------------------
    #
    #
    def _job_run (self, jd) :
        """ runs a job on the wrapper via pty, and returns the job id """

        # stage data, then run job
        self._adaptor.stage_input (self.shell, jd)

        # create command to run
        cmd = self._jd2cmd (jd)
        ret = 1
        out = ""

        run_cmd  = ""
        use_lrun = False

        # simple one-liners use RUN, otherwise LRUN
        if not "\n" in cmd :
            run_cmd = "RUN %s\n" % cmd
        else :
            use_lrun = True
            run_cmd  = "BULK\nLRUN\n%s\nLRUN_EOT\nBULK_RUN\n" % cmd

        run_cmd = run_cmd.replace ("\\", "\\\\\\\\") # hello MacOS

        ret, out, _ = self.shell.run_sync (run_cmd)
        if  ret != 0 :
            raise saga.NoSuccess ("failed to run Job '%s': (%s)(%s)" % (cmd, ret, out))

        lines = filter (None, out.split ("\n"))
        self._logger.debug (lines)

        if  len (lines) < 2 :
            raise saga.NoSuccess ("Failed to run job (%s)" % lines)
        
      # for i in range (0, len(lines)) :
      #     print "%d: %s" % (i, lines[i])

        if lines[-2] != "OK" :
            raise saga.NoSuccess ("Failed to run Job (%s)" % lines)

        # FIXME: verify format of returned pid (\d+)!
        pid    = lines[-1].strip ()
        job_id = "[%s]-[%s]" % (self.rm, pid)

        self._logger.debug ("started job %s" % job_id)

        self.njobs += 1

        # before we return, we need to clean the 'BULK COMPLETED message from lrun
        if use_lrun :
            ret, out = self.shell.find_prompt ()
            if  ret != 0 :
                raise saga.NoSuccess ("failed to run multiline job '%s': (%s)(%s)" % (run_cmd, ret, out))

        return job_id
        

    # ----------------------------------------------------------------
    #
    #
    def _job_get_stats (self, id) :
        """ get the job stats from the wrapper shell """

        rm, pid     = self._adaptor.parse_id (id)
        ret, out, _ = self.shell.run_sync ("STATS %s\n" % pid)

        if  ret != 0 :
            raise saga.NoSuccess ("failed to get job stats for '%s': (%s)(%s)" \
                               % (id, ret, out))

        # the filter removes also empty lines from stdout/stderr.  Oh well...
        lines = filter (None, out.split ("\n"))
        self._logger.debug (lines)

        if lines[0] != "OK" :
            raise saga.NoSuccess ("failed to get valid job state for '%s' (%s)" % (id, lines))

        ret = {}
        ret['STDOUT'] = ""
        ret['STDERR'] = ""
        in_stdout = False
        in_stderr = False

        for line in lines :

            if  in_stdout :
                if  'END_STDOUT' in line :
                    in_stdout = False
                else :
                    ret['STDOUT'] += line
                continue

            if  in_stderr :
                if  'END_STDERR' in line :
                    in_stderr = False
                else :
                    ret['STDERR'] += line
                continue

            if ':' in line :
                key, val = line.split (":", 1)
                ret[key.strip ().lower ()] = val.strip ()
                continue

            if 'START_STDOUT' in line :
                in_stdout = True
                in_stderr = False
                continue

            if 'START_STDERR' in line :
                in_stderr = True
                in_stdout = False
                continue

          # import pprint
          # pprint.pprint (ret)

        return ret

        

    # ----------------------------------------------------------------
    #
    # TODO: this should also fetch the (final) state, to safe a hop
    #
    def _job_get_exit_code (self, id) :
        """ get the job's exit code from the wrapper shell """

        rm, pid = self._adaptor.parse_id (id)

        ret, out, _ = self.shell.run_sync ("RESULT %s\n" % pid)
        if  ret != 0 :
            raise saga.NoSuccess ("failed to get exit code for '%s': (%s)(%s)" \
                               % (id, ret, out))
            self._logger.warning ("failed to get exit code for '%s': (%s)(%s)" \
                               % (id, ret, out))
            return None

        lines = filter (None, out.split ("\n"))
        self._logger.debug (lines)

        if  len (lines) == 3 :
            # shell did not manage to do 'stty -echo'?
            del (lines[0])

        if  len (lines) != 2 :
            raise saga.NoSuccess ("failed to get exit code for '%s': (%s)" % (id, lines))
            self._logger.warning ("failed to get exit code for '%s': (%s)" % (id, lines))
            return None

        if lines[0] != "OK" :
            raise saga.NoSuccess ("failed to get exit code for '%s' (%s)" % (id, lines))
            self._logger.warning ("failed to get exit code for '%s' (%s)" % (id, lines))
            return None

        exit_code = lines[1].strip ()

        if not exit_code.isdigit () :
            return None

        return int(exit_code)
        

    # ----------------------------------------------------------------
    #
    # TODO: this should also cache state
    #
    def _job_suspend (self, id) :

        rm, pid = self._adaptor.parse_id (id)

        ret, out, _ = self.shell.run_sync ("SUSPEND %s\n" % pid)
        if  ret != 0 :
            raise saga.NoSuccess ("failed to suspend job '%s': (%s)(%s)" \
                               % (id, ret, out))


    # ----------------------------------------------------------------
    #
    # TODO: this should also cache state
    #
    def _job_resume (self, id) :

        rm, pid = self._adaptor.parse_id (id)

        ret, out, _ = self.shell.run_sync ("RESUME %s\n" % pid)
        if  ret != 0 :
            raise saga.NoSuccess ("failed to resume job '%s': (%s)(%s)" \
                               % (id, ret, out))


    # ----------------------------------------------------------------
    #
    # TODO: this should also fetch the (final) state, to safe a hop
    #
    def _job_cancel (self, id) :

        rm, pid = self._adaptor.parse_id (id)

        ret, out, _ = self.shell.run_sync ("CANCEL %s\n" % pid)
        if  ret != 0 :
            raise saga.NoSuccess ("failed to cancel job '%s': (%s)(%s)" \
                               % (id, ret, out))

        lines = filter (None, out.split ("\n"))
        self._logger.debug (lines)

        if  len (lines) == 3 :
            # shell did not manage to do 'stty -echo'?
            del (lines[0])

        if  len (lines) != 2 :
            raise saga.NoSuccess ("failed to cancel job '%s': (%s)" % (id, lines))

        if lines[0] != "OK" :
            raise saga.NoSuccess ("failed to cancel job '%s' (%s)" % (id, lines))



    # ----------------------------------------------------------------
    #
    @SYNC_CALL
    def run_job (self, cmd, host) :
        """ Implements saga.adaptors.cpi.job.Service.run_job()
        """

        if not cmd :
            raise saga.BadParameter._log (self._logger, "run_job needs a command to run")

        if  host and host != self.rm.host :
            raise saga.BadParameter._log (self._logger, "Can only run jobs on %s, not on %s" \
                                       % (self.rm.host, host))

        cmd_quoted = cmd.replace ("'", "\\\\'")

        jd = saga.job.Description ()

        jd.executable = "/bin/sh"
        jd.arguments  = ["-c", "'%s'" % cmd_quoted]

        job = self.create_job (jd)
        job.run ()

        return job

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
                          "job_schema"      : self.rm.schema }

        return saga.job.Job (_adaptor=self._adaptor, _adaptor_state=adaptor_state)


    # ----------------------------------------------------------------
    @SYNC_CALL
    def get_url (self) :
        """ Implements saga.adaptors.cpi.job.Service.get_url()
        """
        return self.rm


    # ----------------------------------------------------------------
    #
    @SYNC_CALL
    def list (self):

        # FIXME: this should also fetch job state and metadata, and cache those

        ret, out, _ = self.shell.run_sync ("LIST\n")
        if  ret != 0 :
            raise saga.NoSuccess ("failed to list jobs: (%s)(%s)" \
                               % (ret, out))

        lines = filter (None, out.split ("\n"))
        self._logger.debug (lines)

        if lines[0] != "OK" :
            raise saga.NoSuccess ("failed to list jobs (%s)" % (lines))

        del lines[0]
        job_ids = list()

        for line in lines :

            try :
                pid    = line.strip ()
                job_id = "[%s]-[%s]" % (self.rm, pid)
                job_ids.append (job_id)

            except Exception as e:
                self._logger.debug ("Ignore ill-formatted job id (%s) (%s)" % (line, e))
                continue

        return job_ids
   
   
    # ----------------------------------------------------------------
    #
    @SYNC_CALL
    def get_job (self, job_id, no_reconnect=False):
        """ Implements saga.adaptors.cpi.job.Service.get_url()
        """

        if  job_id in self.jobs :
            # no need to reconnect
            return self.jobs[job_id]

        if  no_reconnect :
            return None

        known_jobs = self.list ()

        if  job_id not in known_jobs :
            # can't reconnect
            raise saga.BadParameter._log (self._logger, "job id '%s' unknown"
                                       % job_id)

        # this dict is passed on to the job adaptor class -- use it to pass any
        # state information you need there.
        adaptor_state = { "job_service"     : self, 
                          "job_id"          : job_id,
                          "job_schema"      : self.rm.schema }

        return saga.job.Job (_adaptor=self._adaptor, _adaptor_state=adaptor_state)


    # ----------------------------------------------------------------
    #
    @SYNC_CALL
    def container_run (self, jobs) :
        """
        From all the job descriptions in the container, build a bulk, and submit
        as async.  The read whaterver the wrapper returns, and sort through the
        messages, assigning job IDs etc.
        """

        # FIXME: this just assumes that all tasks are job creation tasks --
        #        which is not necessarily true...

        self._logger.debug ("container run: %s"  %  str(jobs))

        bulk = "BULK\n"

        for job in jobs :
            cmd   = self._jd2cmd (job.description)
            bulk += "RUN %s\n" % cmd

        # ------------------------------------------------------------
        # stage input data
        # FIXME: this is now blocking the run() method.  Ideally, this activity
        # should be passed to a data manager thread/process/service.
        for job in jobs :
            self._adaptor.stage_input (self.shell, job.description)
        # ------------------------------------------------------------

        bulk += "BULK_RUN\n"
        self.shell.run_async (bulk)

        for job in jobs :

            ret, out = self.shell.find_prompt ()

            if  ret != 0 :
                job._adaptor._set_state (saga.job.FAILED)
                job._adaptor._exception = saga.NoSuccess ("failed to run job: (%s)(%s)" % (ret, out))
                continue

            lines = filter (None, out.split ("\n"))

            if  len (lines) < 2 :
                job._adaptor._set_state (saga.job.FAILED)
                job._adaptor._exception = saga.NoSuccess ("failed to run job : (%s)(%s)" % (ret, out))
                continue

            if lines[-2] != "OK" :
                job._adaptor._set_state (saga.job.FAILED)
                job._adaptor._exception = saga.NoSuccess ("failed to run job : (%s)(%s)" % (ret, out))
                continue

            # FIXME: verify format of returned pid (\d+)!
            pid    = lines[-1].strip ()
            job_id = "[%s]-[%s]" % (self.rm, pid)

            self._logger.debug ("started job %s" % job_id)

            self.njobs += 1

            # FIXME: at this point we need to make sure that we actually created
            # the job.  Well, we should make sure of this *before* we run it.
            # But, actually, the container sorter should have done that already?
            # Check!
            job._adaptor._id = job_id

        # we also need to find the output of the bulk op itself
        ret, out = self.shell.find_prompt ()

        if  ret != 0 :
            self._logger.error ("failed to run (parts of the) bulk jobs: (%s)(%s)" % (ret, out))
            return

        lines = filter (None, out.split ("\n"))

        if  len (lines) < 2 :
            self._logger.error ("Cannot evaluate status of bulk job submission: (%s)(%s)" % (ret, out))
            return

        if lines[-2] != "OK" :
            self._logger.error ("failed to run (parts of the) bulk jobs: (%s)(%s)" % (ret, out))
            return

   
    # ----------------------------------------------------------------
    #
    @SYNC_CALL
    def container_wait (self, jobs, mode, timeout) :

        # FIXME: this just assumes that all tasks are job wait tasks --
        #        which is not necessarily true...
        # FIXME: we ignore the job wait mode (ALL/ANY), and always wait for all
        #        jobs...

        self._logger.debug ("container wait: %s"  %  str(jobs))

        bulk = "BULK\n"

        for job in jobs :
          # print type (job)
          # print type (job._adaptor)

            if  not isinstance (job._adaptor, ShellJob) :
                # this is not a job created by this adaptor.  Its probably
                # a task for a job operation where the job is owned by this
                # adaptor (FIXME: check).  Fall back to non-container wait.
                # FIXME: timeout handling is wrong
                job.wait (timeout)
            else :
                rm, pid = self._adaptor.parse_id (job.id)
                bulk   += "WAIT %s\n" % pid

        bulk += "BULK_RUN\n"
        self.shell.run_async (bulk)

        for job in jobs :

            ret, out = self.shell.find_prompt ()

            if  ret != 0 :
                job._adaptor._set_state (saga.job.FAILED)
                job._adaptor._exception = saga.NoSuccess ("failed to wait for job: (%s)(%s)" % (ret, out))
                continue

            lines = filter (None, out.split ("\n"))

            if  len (lines) < 2 :
                job._adaptor._set_state (saga.job.FAILED)
                job._adaptor._exception = saga.NoSuccess ("failed to wait for job : (%s)(%s)" % (ret, out))
                continue

            if lines[-2] != "OK" :
                job._adaptor._set_state (saga.job.FAILED)
                job._adaptor._exception = saga.NoSuccess ("failed to wait for job : (%s)(%s)" % (ret, out))
                continue

        # we also need to find the output of the bulk op itself
        ret, out = self.shell.find_prompt ()

        if  ret != 0 :
            self._logger.error ("failed to wait for (parts of the) bulk jobs: (%s)(%s)" % (ret, out))
            return

        lines = filter (None, out.split ("\n"))

        if  len (lines) < 2 :
            self._logger.error ("Cannot evaluate status of bulk job wait: (%s)(%s)" % (ret, out))
            return

        if lines[-2] != "OK" :
            self._logger.error ("failed to wait for (parts of the) bulk jobs: (%s)(%s)" % (ret, out))
            return

   
    # ----------------------------------------------------------------
    #
    @SYNC_CALL
    def container_cancel (self, jobs, timeout) :

        self._logger.debug ("container cancel: %s [%s]"  %  (str(jobs), timeout))

        bulk = "BULK\n"

        for job in jobs :
            rm, pid = self._adaptor.parse_id (job.id)
            bulk   += "CANCEL %s\n" % pid

        bulk += "BULK_RUN\n"
        self.shell.run_async (bulk)

        for job in jobs :

            ret, out = self.shell.find_prompt ()

            if  ret != 0 :
                job._adaptor._set_state (saga.job.FAILED)
                job._adaptor._exception = saga.NoSuccess ("failed to cancel job: (%s)(%s)" % (ret, out))
                continue

            lines = filter (None, out.split ("\n"))

            if  len (lines) < 2 :
                job._adaptor._set_state (saga.job.FAILED)
                job._adaptor._exception = saga.NoSuccess ("failed to cancel job : (%s)(%s)" % (ret, out))
                continue

            if lines[-2] != "OK" :
                job._adaptor._set_state (saga.job.FAILED)
                job._adaptor._exception = saga.NoSuccess ("failed to cancel job : (%s)(%s)" % (ret, out))
                continue

        # we also need to find the output of the bulk op itself
        ret, out = self.shell.find_prompt ()

        if  ret != 0 :
            self._logger.error ("failed to cancel (parts of the) bulk jobs: (%s)(%s)" % (ret, out))
            return

        lines = filter (None, out.split ("\n"))

        if  len (lines) < 2 :
            self._logger.error ("Cannot evaluate status of bulk job cancel: (%s)(%s)" % (ret, out))
            return

        if lines[-2] != "OK" :
            self._logger.error ("failed to cancel (parts of the) bulk jobs: (%s)(%s)" % (ret, out))
            return


    # ----------------------------------------------------------------
    #
    @SYNC_CALL
    def container_get_states (self, jobs) :

        self._logger.debug ("container get_state: %s"  %  str(jobs))

        bulk   = "BULK\n"
        states = []

        for job in jobs :

          # print job
          # job._attributes_dump ()

            rm, pid = self._adaptor.parse_id (job.id)
            bulk   += "STATE %s\n" % pid

        bulk += "BULK_RUN\n"
        self.shell.run_async (bulk)

        for job in jobs :

            ret, out = self.shell.find_prompt ()

            if  ret != 0 :
                job._adaptor._set_state (saga.job.FAILED)
                job._adaptor._exception = saga.NoSuccess ("failed to get job state: (%s)(%s)" % (ret, out))
                continue

            lines = filter (None, out.split ("\n"))

            if  len (lines) < 2 :
                job._adaptor._set_state (saga.job.FAILED)
                job._adaptor._exception = saga.NoSuccess ("failed to get job state : (%s)(%s)" % (ret, out))
                continue

            if lines[-2] != "OK" :
                job._adaptor._set_state (saga.job.FAILED)
                job._adaptor._exception = saga.NoSuccess ("failed to get job state : (%s)(%s)" % (ret, out))
                continue

            state = self._adaptor.string_to_state (lines[-1])

            job._adaptor._update_state (state)
            states.append (state)


        # we also need to find the output of the bulk op itself
        ret, out = self.shell.find_prompt ()

        if  ret != 0 :
            self._logger.error ("failed to get state for (parts of the) bulk jobs: (%s)(%s)" % (ret, out))
            return

        lines = filter (None, out.split ("\n"))

        if  len (lines) < 2 :
            self._logger.error ("Cannot evaluate status of bulk job get_state: (%s)(%s)" % (ret, out))
            return

        if lines[-2] != "OK" :
            self._logger.error ("failed to get state for (parts of the) bulk jobs: (%s)(%s)" % (ret, out))
            return

        return states


###############################################################################
#
class ShellJob (saga.adaptors.cpi.job.Job) :
    """ Implements saga.adaptors.cpi.job.Job
    """
    # ----------------------------------------------------------------
    #
    def __init__ (self, api, adaptor) :

        _cpi_base = super  (ShellJob, self)
        _cpi_base.__init__ (api, adaptor)


    # ----------------------------------------------------------------
    #
    @SYNC_CALL
    def init_instance (self, job_info):
        """ Implements saga.adaptors.cpi.job.Job.init_instance()
        """

        if  'job_description' in job_info :
            # comes from job.service.create_job()
            self.js = job_info["job_service"] 
            self.jd = job_info["job_description"]

            # the js is responsible for job bulk operations -- which
            # for jobs only work for run()
            self._container       = self.js

            # initialize job attribute values
            self._id              = None
            self._name            = self.jd.get(saga.job.NAME)
            self._log             = list()
            self._state           = None
            self._exit_code       = None
            self._exception       = None
            self._created         = time.time ()
            self._name            = self.jd.name
            self._started         = None
            self._finished        = None

            self._set_state (saga.job.NEW)

        elif 'job_id' in job_info :
            # initialize job attribute values
            self.js               = job_info["job_service"] 
            self.jd               = None
            self._id              = job_info['job_id']
            self._name            = job_info.get('job_name')
            self._log             = list()
            self._state           = None
            self._exit_code       = None
            self._exception       = None
            self._created         = None
            self._name            = None
            self._started         = None
            self._finished        = None

        else :
            # don't know what to do...
            raise saga.BadParameter ("Cannot create job, insufficient information")
        
        if self._created : self._created = float(self._created)

        return self.get_api ()


    # ----------------------------------------------------------------
    #
    @SYNC_CALL
    def get_description (self):
        return self.jd

    # --------------------------------------------------------------------------
    #
    def _update_state (self, state) :

        old_state = self._state

        if  state == saga.job.DONE and \
            old_state not in [saga.job.DONE, saga.job.FAILED, saga.job.CANCELED] :

            # stage output data
            # FIXME: _update_state blocks until data are staged.  That should not happen.
            self._adaptor.stage_output (self.js.shell, self.jd)
        
        # files are staged -- update state, and report to application
        self._state = state
        self._api ()._attributes_i_set ('state', self._state, self._api ()._UP)


    # ----------------------------------------------------------------
    #
    @SYNC_CALL
    def get_state (self):
        """ Implements saga.adaptors.cpi.job.Job.get_state() """

        # may not yet have backend representation, state is 'NEW'
        if  self._id == None :
            return self._state

        # no need to re-fetch final states
        if  self._state == saga.job.DONE      or \
            self._state == saga.job.FAILED    or \
            self._state == saga.job.CANCELED     :
                return self._state

        stats = self.js._job_get_stats (self._id)

        if 'start' in stats : self._started  = stats['start']
        if 'stop'  in stats : self._finished = stats['stop']

        if self._started  : self._started  = float(self._started)
        if self._finished : self._finished = float(self._finished)
        
        if  not 'state' in stats :
            raise saga.NoSuccess ("failed to get job state for '%s': (%s)" \
                               % (self._id, stats))

        self._update_state (self._adaptor.string_to_state (stats['state']))

        return self._state


    # ----------------------------------------------------------------
    #
    def _set_state (self, state) :

        old_state = self._state

        # on state changes, trigger notifications
        if  old_state != state :
            self._state  = state
            self._api ()._attributes_i_set ('state', state, self._api ()._UP)
        
        return self._state


    # ----------------------------------------------------------------
    #
    @SYNC_CALL
    def get_created (self) : 

        # no need to refresh stats -- this is set locally
        return self._created


    # ----------------------------------------------------------------
    #
    @SYNC_CALL
    def get_name (self) : 

        # no need to refresh stats -- this is set locally
        return self._name


    # ----------------------------------------------------------------
    #
    @SYNC_CALL
    def get_started (self) : 

        self.get_state () # refresh stats
        return self._started


    # ----------------------------------------------------------------
    #
    @SYNC_CALL
    def get_finished (self) : 

        self.get_state () # refresh stats
        return self._finished


    # ----------------------------------------------------------------
    #
    @SYNC_CALL
    def get_stdout (self) : 

        state = self.get_state () # refresh stats

        if  state == saga.job.NEW      or \
            state == saga.job.PENDING  :
            raise saga.IncorrectState ("Job output is only available after the job started")

        if  not self._id :
            raise saga.IncorrectState ("Job output is only available after the job started")

        rm, pid     = self._adaptor.parse_id (self._id)
        ret, out, _ = self.js.shell.run_sync ("STDOUT %s\n" % pid)

        if  ret != 0 :
            raise saga.NoSuccess ("failed to get job stdout for '%s': (%s)(%s)" \
                               % (self._id, ret, out))

        lines = filter (None, out.split ("\n"))

        if lines[0] != "OK" :
            raise saga.NoSuccess ("failed to get valid job stdout for '%s' (%s)" % (self._id, lines))

        return _decode('\n'.join (lines[1:]))


    # ----------------------------------------------------------------
    #
    @SYNC_CALL
    def get_stderr (self) : 

        state = self.get_state () # refresh stats

        if  state == saga.job.NEW      or \
            state == saga.job.PENDING  :
            raise saga.IncorrectState ("Job output is only available after the job started")

        if  not self._id :
            raise saga.IncorrectState ("Job output is only available after the job started")

        rm, pid     = self._adaptor.parse_id (self._id)
        ret, out, _ = self.js.shell.run_sync ("STDERR %s\n" % pid)

        if  ret != 0 :
            raise saga.NoSuccess ("failed to get job stderr for '%s': (%s)(%s)" \
                               % (self._id, ret, out))

        lines = filter (None, out.split ("\n"))

        if lines[0] != "OK" :
            raise saga.NoSuccess ("failed to get valid job stderr for '%s' (%s)" % (self._id, lines))

        return _decode('\n'.join (lines[1:]))


    # ----------------------------------------------------------------
    #
    @SYNC_CALL
    def get_log (self) : 

        state = self.get_state () # refresh stats

        if  state == saga.job.NEW      or \
            state == saga.job.PENDING  :
            raise saga.IncorrectState ("Job output is only available after the job started")

        if  not self._id :
            raise saga.IncorrectState ("Job output is only available after the job started")

        rm, pid     = self._adaptor.parse_id (self._id)
        ret, out, _ = self.js.shell.run_sync ("LOG %s\n" % pid)

        if  ret != 0 :
            raise saga.NoSuccess ("failed to get job log for '%s': (%s)(%s)" \
                               % (self._id, ret, out))

        lines = filter (None, out.split ("\n"))

        if lines[0] != "OK" :
            raise saga.NoSuccess ("failed to get valid job stderr for '%s' (%s)" % (self._id, lines))

        ret  = '\n'.join(self._log)  # pre-pend all local log messages
        ret += _decode('\n'.join(lines[1:]))

        return ret


    # ----------------------------------------------------------------
    #
    @SYNC_CALL
    def get_service_url (self):

        if not self.js :
            raise saga.IncorrectState ("Job Service URL unknown")
        else :
            return self.js.get_url ()


    # ----------------------------------------------------------------
    #
    # TODO: this should also fetch the (final) state, to safe a hop
    # TODO: implement via notifications
    #
    @SYNC_CALL
    def wait (self, timeout):
        """ 
        A call to the shell to do the WAIT would block the shell for any
        other interactions.  In particular, it would practically kill it if the
        Wait waits forever...

        So we implement the wait via a state pull.  The *real* solution is, of
        course, to implement state notifications, and wait for such
        a notification to arrive within timeout seconds...
        """

        time_start = time.time ()
        time_now   = time_start

        while True :

            state = self.get_state ()

            if  state == saga.job.DONE      or \
                state == saga.job.FAILED    or \
                state == saga.job.CANCELED     :
                    return True

            # avoid busy poll
            # FIXME: self-tune by checking call latency
            time.sleep (0.1)

            # check if we hit timeout
            if  timeout >= 0 :
                time_now = time.time ()
                if  time_now - time_start > timeout :
                    return False
   
    # ----------------------------------------------------------------
    #
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
    def get_exit_code (self) :
        """ Implements saga.adaptors.cpi.job.Job.get_exit_code() """

        if  self._exit_code != None :
            return self._exit_code

        if  self.get_state () not in [saga.job.DONE, 
                                      saga.job.FAILED, 
                                      saga.job.CANCELED] :
            raise saga.IncorrectState ("Cannot get exit code, job is not in final state")

        self._exit_code = self.js._job_get_exit_code (self._id)

        return self._exit_code
   
    # ----------------------------------------------------------------
    #
    # TODO: the values below should be fetched with every get_state...
    #
    @SYNC_CALL
    def get_execution_hosts (self) :
        """ Implements saga.adaptors.cpi.job.Job.get_execution_hosts()
        """        
        self._logger.debug ("this is the shell adaptor, reporting execution hosts")
        return [self.js.get_url ().host]
   
    # ----------------------------------------------------------------
    #
    @SYNC_CALL
    def run (self): 

        self._id = self.js._job_run (self.jd)
        self.js.jobs[self._id] = self._api ()

        self._set_state (saga.job.RUNNING)


    # ----------------------------------------------------------------
    #
    @SYNC_CALL
    def suspend (self):

        if  self.get_state () != saga.job.RUNNING :
            raise saga.IncorrectState ("Cannot suspend, job is not RUNNING")

        self.js._job_suspend (self._id)
   
   
    # ----------------------------------------------------------------
    #
    @SYNC_CALL
    def resume (self):

        if  self.get_state () != saga.job.SUSPENDED :
            raise saga.IncorrectState ("Cannot resume, job is not SUSPENDED")

        self.js._job_resume (self._id)
   
   
    # ----------------------------------------------------------------
    #
    @SYNC_CALL
    def cancel (self, timeout):

        if  self.get_state () not in [saga.job.RUNNING, 
                                      saga.job.SUSPENDED, 
                                      saga.job.CANCELED, 
                                      saga.job.DONE, 
                                      saga.job.FAILED] :
            raise saga.IncorrectState ("Cannot cancel, job is not running")

        if  self._state in [saga.job.CANCELED, 
                            saga.job.DONE, 
                            saga.job.FAILED] :
            self._set_state (saga.job.CANCELED)
            return

        self.js._job_cancel (self._id)
   
   
    # ----------------------------------------------------------------
    #
    @SYNC_CALL
    def re_raise (self):
        # nothing to do here actually, as run () is synchronous...
        return self._exception




