
__author__    = "Andre Merzky, Ole Weidner"
__copyright__ = "Copyright 2012-2013, The SAGA Project"
__license__   = "MIT"


import os
import sys
import pwd
import time
import string
import getpass

import radical.utils           as ru
import radical.utils.logger    as rul

import saga
import saga.exceptions         as se
import saga.utils.misc         as sumisc
import saga.utils.pty_process  as supp

import pty_exceptions               as ptye


# ------------------------------------------------------------------------------
#
# ssh options:
#   -e none         : no escape character
#   -M              : master mode for connection sharing
#   -S control_path : slave mode for connection sharing
#   -t              : force pty allocation
#   -x              : disable x11 forwarding
#
#   ServerAliveInterval
#   CheckHostIP no
#   ConnectTimeout
#   ControlMaster  yes | no | no ...
#   ControlPath    $BASE/ssh_control_%n_%p.$$.sock
#                  %r (remote id)? would need inspection
#   ControlPersist 100  : close master after 100 seconds idle
#   EscapeChar     none : transparent for binary data
#   TCPKeepAlive   yes  : detect connection failure
#
#   LoginGraceTime seconds : disconnect if no login after n seconds
#
# ------------------------------------------------------------------------------

# these arrays help to map requested client schemas to master schemas
_SCHEMAS_SH  = ['sh', 'fork', 'local', 'file']
_SCHEMAS_SSH = ['ssh', 'scp', 'sftp']
_SCHEMAS_GSI = ['gsissh', 'gsiscp', 'gsisftp']

_SCHEMAS = _SCHEMAS_SH + _SCHEMAS_SSH + _SCHEMAS_GSI

# FIXME: '-o ControlPersist' is only supported for newer ssh versions.  We
# should add detection, and enable that if available -- for now, just diable it.
#
# FIXME: we should use '%n' instead of '%h', but that is not supported by older
# ssh versions...

# ssh master/slave flag magic # FIXME: make timeouts configurable
_SSH_FLAGS_MASTER   = "-o ControlMaster=%(share_mode)s -o ControlPath=%(ctrl)s -o TCPKeepAlive=no  -o ServerAliveInterval=10 -o ServerAliveCountMax=20 %(connect_timeout)s"
_SSH_FLAGS_SLAVE    = "-o ControlMaster=%(share_mode)s -o ControlPath=%(ctrl)s -o TCPKeepAlive=no  -o ServerAliveInterval=10 -o ServerAliveCountMax=20 %(connect_timeout)s"
_SCP_FLAGS          = ""
_SFTP_FLAGS         = ""

# FIXME: right now, we create a shell connection as master --
# but a master does not actually need a shell, as it is never really
# used to run commands...
_SCRIPTS = {
    'ssh' : {
        'master'       : '%(ssh_env)s "%(ssh_exe)s" %(ssh_args)s %(m_flags)s %(host_str)s',
        'shell'        : '%(ssh_env)s "%(ssh_exe)s" %(ssh_args)s %(s_flags)s %(host_str)s',
        'copy_is_posix': True
    },
    'scp' : {
        'copy_to'      : '%(scp_env)s "%(scp_exe)s" %(scp_args)s %(s_flags)s %(cp_flags)s "%(src)s" "%(scp_root)s%(tgt)s"',
        'copy_from'    : '%(scp_env)s "%(scp_exe)s" %(scp_args)s %(s_flags)s %(cp_flags)s "%(scp_root)s%(src)s" "%(tgt)s"',
        'copy_to_in'   : '',
        'copy_from_in' : '',
        'copy_is_posix': False
    },
    'sftp' : {
        'copy_to'      : '%(sftp_env)s "%(sftp_exe)s" %(sftp_args)s %(s_flags)s %(host_str)s',
        'copy_from'    : '%(sftp_env)s "%(sftp_exe)s" %(sftp_args)s %(s_flags)s %(host_str)s',
        'copy_to_in'   : 'mput %(cp_flags)s "%(src)s" "%(tgt)s"',
        'copy_from_in' : 'mget %(cp_flags)s "%(src)s" "%(tgt)s"',
        'copy_is_posix': False
    },
    'sh' : {
        'master'       : '%(sh_env)s "%(sh_exe)s"  %(sh_args)s',
        'shell'        : '%(sh_env)s "%(sh_exe)s"  %(sh_args)s',
        'copy_to'      : '%(sh_env)s "%(sh_exe)s"  %(sh_args)s',
        'copy_from'    : '%(sh_env)s "%(sh_exe)s"  %(sh_args)s',
        'copy_to_in'   : 'cd ~ && "%(cp_exe)s" -v %(cp_flags)s "%(src)s" "%(tgt)s"',
        'copy_from_in' : 'cd ~ && "%(cp_exe)s" -v %(cp_flags)s "%(src)s" "%(tgt)s"',
        'copy_is_posix': True
    }
}

# ------------------------------------------------------------------------------
#
class PTYShellFactory (object) :
    """
    This is the place where all master and all client shell connections get
    created.  But also, this factory maintains a registry of master connections,
    to quickly spawn slave connections for any customer w/o repeated
    authorization overhead.  Masters are unique per
    a host/user/port/context/shell_type hash.

    Any ssh master connection in this registry can idle, and may thus shut down
    after ``ControlPersist`` seconds (see options).

    data model::


      self.registry
        |
        +-- "host[:port]"
        |   |
        |   +-- str(context)
        |   |   |
        |   |   +-- "shell_type" (ssh)
        |   |   |   |
        |   |   |   +-- pty_process (gc_timeout)
        |   |   |   +-- shell_initialize()
        |   |   |   +-- shell_finalize()
        |   |   |   +-- shell_alive()
        |   |   |
        |   |   +-- ...
        |   |
        |   +-- ...
        |
        +-- ...

    When Slave connections are requested, a suitable master connection is looked
    for and used.  'Suitable' means: ssh master for scp and sftp slaves; gsissh
    for gsiscp and gsisftp slaves; and sh master for file slaves

    """

    __metaclass__ = ru.Singleton


    # --------------------------------------------------------------------------
    #
    def __init__ (self) :

        self.logger     = ru.get_logger ('radical.saga.pty')
        self.registry   = {}
        self.rlock      = ru.RLock ('pty shell factory')


    # --------------------------------------------------------------------------
    #
    def initialize (self, url, session=None, prompt=None, logger=None,
            posix=True, interactive=True) :

        with self.rlock :

            # make sure we have a valid url type
            url = saga.Url (url)

            if  not prompt :
                prompt = "^(.*[\$#%>\]])\s*$"

            if  not logger :
                logger = self.logger

            # collect all information we have/need about the requested master
            # connection
            info = self._create_master_entry (url, session, prompt, logger,
                    posix, interactive)

            # we got master info - register the master, and create the instance!
            type_s = str(info['shell_type'])
            user_s = str(info['user'])
            host_s = str(info['host_str'])

            # Now, if we don't have that master, yet, we need to instantiate it
            if not host_s in self.registry                 : self.registry[host_s] = {}
            if not user_s in self.registry[host_s]         : self.registry[host_s][user_s] = {}
            if not type_s in self.registry[host_s][user_s] :

                # new master: create an instance, and register it
                m_cmd = info['scripts'][info['shell_type']]['master'] % info

                logger.debug ("open master pty for [%s] [%s] %s: %s'" \
                                % (type_s, host_s, user_s, m_cmd))

                info['pty'] = supp.PTYProcess (m_cmd, logger=logger)
                if not info['pty'].alive () :
                    raise se.NoSuccess._log (logger, \
                          "Shell not connected to %s" % info['host_str'])

                # authorization, prompt setup, etc.  Initialize as shell if not
                # explicitly marked as non-posix shell
                self._initialize_pty (info['pty'], info)

                # master was created - register it
                self.registry[host_s][user_s][type_s] = info


            else :
                # we already have a master: make sure it is alive, and restart as
                # needed
                info = self.registry[host_s][user_s][type_s]

                if  not info['pty'].alive (recover=True) :
                    raise se.IncorrectState._log (logger, \
                          "Lost shell connection to %s" % info['host_str'])

            return info


    # --------------------------------------------------------------------------
    #
    def _initialize_pty (self, pty_shell, info, posix=None) :

        # posix: only for posix shells we use prompt triggers.  sftp for example
        # does not deal well with triggers (no printf).

        with self.rlock :

            # import pprint
            # pprint.pprint (info)

            shell_pass = info['pass']
            key_pass   = info['key_pass']
            prompt     = info['prompt']
            logger     = info['logger']
            latency    = info['latency']
            timeout    = info['ssh_timeout']

            pty_shell.latency = latency

            if posix == None:
                posix = info['posix']

            # if we did not see a decent prompt within 'delay' time, something
            # went wrong.  Try to prompt a prompt (duh!)  Delay should be
            # minimum 0.1 second (to avoid flooding of local shells), and at
            # maximum 1 second (to keep startup time reasonable)
            # most one second.  We try to get within that range with 10*latency.
            delay = min (1.0, max (0.1, 10 * latency))

            try :
                prompt_patterns = ["[Pp]assword:\s*$",             # password   prompt
                                   "Enter passphrase for .*:\s*$", # passphrase prompt
                                   "Token_Response.*:\s*$",        # passtoken  prompt
                                   "Enter PASSCODE:$",             # RSA SecureID
                                   "want to continue connecting",  # hostkey confirmation
                                   ".*HELLO_\\d+_SAGA$",           # prompt detection helper
                                   prompt]                         # greedy native shell prompt

                # use a very aggressive, but portable prompt setting scheme.
                # Error messages may appear for tcsh and others.  Excuse
                # non-posix shells
                if posix:
                    pty_shell.write (" export PS1='$' ; set prompt='$'\n")

                # find a prompt
                n, match = pty_shell.find (prompt_patterns, delay)

                # this loop will run until we finally find the shell prompt, or
                # if we think we have tried enough and give up.  On success
                # we'll try to set a different prompt, and when we found that,
                # too, we exit the loop and are be ready to running shell

                # commands.
                retries       = 0
                retry_trigger = True
                used_trigger  = False
                found_trigger = ""
                time_start    = time.time()

                while True :

                    # --------------------------------------------------------------
                    if n == None :

                        # we found none of the prompts, yet, and need to try
                        # again.  But to avoid hanging on invalid prompts, we
                        # print 'HELLO_x_SAGA', and search for that one, too.
                        # We actually do 'printf HELLO_%d_SAGA x' so that the
                        # pattern only appears in the result, not in the
                        # command...

                        if time.time() - time_start > timeout:
                            raise se.NoSuccess ("Could not detect shell prompt (timeout)")

                        # make sure we retry a finite time...
                        retries += 1

                        if  not retry_trigger :
                            # just waiting for the *right* trigger or prompt,
                            # don't need new ones...
                            continue

                        if posix:
                            # use a very aggressive, but portable prompt setting scheme
                            pty_shell.write (" export PS1='$' > /dev/null 2>&1 || set prompt='$'\n")
                            pty_shell.write (" printf 'HELLO_%%d_SAGA\\n' %d\n" % retries)
                            used_trigger = True

                        # FIXME:  consider better timeout
                        n, match = pty_shell.find (prompt_patterns, delay)


                    # --------------------------------------------------------------
                    elif n == 0 :
                        logger.info ("got password prompt")
                        if  not shell_pass :
                            raise se.AuthenticationFailed ("prompted for unknown password (%s)" \
                                                          % match)

                        pty_shell.write ("%s\n" % shell_pass, nolog=True)
                        n, match = pty_shell.find (prompt_patterns, delay)


                    # --------------------------------------------------------------
                    elif n == 1 :
                        logger.info ("got passphrase prompt : %s" % match)

                        start = string.find (match, "'", 0)
                        end   = string.find (match, "'", start+1)

                        if start == -1 or end == -1 :
                            raise se.AuthenticationFailed ("could not extract key name (%s)" % match)

                        key = match[start+1:end]

                        if  not key in key_pass    :
                            raise se.AuthenticationFailed ("prompted for unknown key password (%s)" \
                                                          % key)

                        pty_shell.write ("%s\n" % key_pass[key], nolog=True)
                        n, match = pty_shell.find (prompt_patterns, delay)


                    # --------------------------------------------------------------
                    elif n == 2 or n == 3:
                        logger.info ("got token prompt")
                        import getpass
                        token = getpass.getpass ("enter token: ")
                        pty_shell.write ("%s\n" % token.strip(), nolog=True)
                        n, match = pty_shell.find (prompt_patterns, delay)


                    # --------------------------------------------------------------
                    elif n == 4:
                        logger.info ("got hostkey prompt")
                        pty_shell.write ("yes\n")
                        n, match = pty_shell.find (prompt_patterns, delay)


                    # --------------------------------------------------------------
                    elif n == 5:

                        # one of the trigger commands got through -- we can now
                        # hope to find the prompt (or the next trigger...)
                        logger.debug ("got shell prompt trigger (%s) (%s)" %  (n, match))

                        found_trigger = match
                        retry_trigger = False
                        n, match = pty_shell.find (prompt_patterns, delay)
                        continue


                    # --------------------------------------------------------------
                    elif n == 6 :

                        logger.debug ("got initial shell prompt (%s) (%s)" %  (n, match))

                        if  retries :
                            if  used_trigger :
                                # we already sent triggers -- so this match is only
                                # useful if saw the *correct* shell prompt trigger
                                # first
                                trigger = "HELLO_%d_SAGA" % retries

                                if  not trigger in found_trigger :

                                    logger.debug ("waiting for prompt trigger %s: (%s) (%s)" \
                                               % (trigger, n, match))
                                    # but more retries won't help...
                                    retry_trigger = False
                                    attempts      = 0
                                    n             = None

                                    while not n :

                                        attempts += 1
                                        n, match  = pty_shell.find (prompt_patterns, delay)

                                        if  not n :
                                            if  attempts == 1 :
                                                if posix:
                                                    pty_shell.write (" printf 'HELLO_%%d_SAGA\\n' %d\n" % retries)

                                            if  attempts > 100 :
                                                raise se.NoSuccess ("Could not detect shell prompt (timeout)")

                                    continue


                        logger.debug ("Got initial shell prompt (%s) (%s)" % (n, match))
                        # we are done waiting for a prompt
                        break

            except Exception as e :
                logger.exception(e)
                raise ptye.translate_exception (e)


    # --------------------------------------------------------------------------
    #
    def get_cp_slave (self, s_cmd, info, posix=None) :

        with self.rlock :

            if posix == None:
                posix = info.get('copy_is_posix')

          # print '> -- new cp  shell to %s' % s_cmd

            cp_slave = supp.PTYProcess (s_cmd, info['logger'])
            self._initialize_pty (cp_slave, info, posix)

            return cp_slave

    # --------------------------------------------------------------------------
    #
    def run_shell (self, info) :
        """
        This initiates a master connection.  If there is a suitable master
        connection in the registry, it is re-used, and no new master connection
        is created.  If needed, the existing master connection is revived.
        """

      # if True :
        with self.rlock :

            s_cmd = info['scripts'][info['shell_type']]['shell'] % info

            # at this point, we do have a valid, living master
            sh_slave = supp.PTYProcess (s_cmd, info['logger'])

            # authorization, prompt setup, etc
            self._initialize_pty (sh_slave, info)

            return sh_slave


    # --------------------------------------------------------------------------
    #
    def _which(self, cmd):

        ret = ru.which(cmd)
        if not ret:
            raise RuntimeError('cmd %s not found' % cmd)
        return ret


    # --------------------------------------------------------------------------
    #
    def _create_master_entry (self, url, session, prompt, logger, posix,
            interactive) :
        # FIXME: cache 'which' results, etc
        # FIXME: check 'which' results

        with self.rlock :


            info = {'posix' : posix}

            # get and evaluate session config
            if  not session :
                session = saga.Session (default=True)

            session_cfg = session.get_config ('saga.utils.pty')
            info['ssh_copy_mode']  = session_cfg['ssh_copy_mode'].get_value ()
            info['ssh_share_mode'] = session_cfg['ssh_share_mode'].get_value ()
            info['ssh_timeout']    = session_cfg['ssh_timeout'].get_value ()

            logger.info ("ssh copy  mode set to '%s'" % info['ssh_copy_mode' ])
            logger.info ("ssh share mode set to '%s'" % info['ssh_share_mode'])
            logger.info ("ssh timeout    set to '%s'" % info['ssh_timeout'])


            # fill the info dict with details for this master channel, and all
            # related future slave channels
            info['schema']    = url.schema.lower ()
            info['host_str']  = url.host
            info['prompt']    = prompt
            info['logger']    = logger
            info['url']       = url
            info['pass']      = ""
            info['key_pass']  = {}
            info['scripts']   = _SCRIPTS

            if  not info['schema'] :
                info['schema'] = 'local'


            # find out what type of shell we have to deal with
            if  info['schema'] in _SCHEMAS_SSH :
                info['shell_type'] = "ssh"
                info['copy_mode']  = info['ssh_copy_mode']
                info['share_mode'] = info['ssh_share_mode']
                info['ssh_exe']    = self._which ("ssh")
                info['scp_exe']    = self._which ("scp")
                info['sftp_exe']   = self._which ("sftp")

            elif info['schema'] in _SCHEMAS_GSI :
                info['shell_type'] = "ssh"
                info['copy_mode']  = info['ssh_copy_mode']
                info['share_mode'] = info['ssh_share_mode']
                info['ssh_exe']    = self._which ("gsissh")
                info['scp_exe']    = self._which ("gsiscp")
                info['sftp_exe']   = self._which ("gsisftp")

            elif info['schema'] in _SCHEMAS_SH :
                info['shell_type'] = "sh"
                info['copy_mode']  = "sh"
                info['share_mode'] = "auto"
                info['sh_env']     = "/usr/bin/env TERM=vt100 PS1='PROMPT-$?->'"
                info['cp_env']     = "/usr/bin/env TERM=vt100 PS1='PROMPT-$?->'"
                info['scp_root']   = "/"

                if interactive: info['sh_args'] = "-i"
                else          : info['sh_args'] = ""

                if  "SHELL" in os.environ :
                    info['sh_exe'] =  self._which (os.environ["SHELL"])
                    info['cp_exe'] =  self._which ("cp")
                else :
                    info['sh_exe'] =  self._which ("sh")
                    info['cp_exe'] =  self._which ("cp")

            else :
                raise se.BadParameter._log (self.logger, \
                          "cannot handle schema '%s://'" % url.schema)


            # If an SSH timeout has been specified set up the ConnectTimeout
            # string
            if info['ssh_timeout']:
                info['ssh_connect_timeout'] = ('-o ConnectTimeout=%s' 
                    % int(float(info['ssh_timeout'])))
            else:
                info['ssh_connect_timeout'] = ''

            # depending on type, create command line (args, env etc)
            #
            # We always set term=vt100 to avoid ansi-escape sequences in the prompt
            # and elsewhere.  Also, we have to make sure that the shell is an
            # interactive login shell, so that it interprets the users startup
            # files, and reacts on commands.

            try :
                info['latency'] = sumisc.get_host_latency (url)

                # FIXME: note that get_host_latency is considered broken (see
                # saga/utils/misc.py line 73), and will return a constant 250ms.

            except Exception  as e :
                info['latency'] = 1.0  # generic value assuming slow link
                info['logger'].warning ("Could not contact host '%s': %s" % (url, e))

            if  info['shell_type'] == "sh" :

                info['sh_env'] = "/usr/bin/env TERM=vt100 "  # avoid ansi escapes

                if not sumisc.host_is_local (url.host) :
                    raise se.BadParameter._log (self.logger, \
                            "expect local host for '%s://', not '%s'" % (url.schema, url.host))

                if  'user' in info and info['user'] :
                    pass
                else :
                    info['user'] = getpass.getuser ()

            else :
                info['ssh_env']   =  "/usr/bin/env TERM=vt100 "  # avoid ansi escapes
                info['scp_env']   =  "/usr/bin/env TERM=vt100 "  # avoid ansi escapes
                info['sftp_env']  =  "/usr/bin/env TERM=vt100 "  # avoid ansi escapes
                info['ssh_args']  =  "-t "                       # force pty
                info['scp_args']  =  _SCP_FLAGS
                info['sftp_args'] =  _SFTP_FLAGS

                if  session :

                    for context in session.contexts :

                        # ssh can also handle UserPass contexts, and ssh type contexts.
                        # gsissh can handle the same, but also X509 contexts.

                        if  context.type.lower () == "ssh" :
                            if  info['schema'] in _SCHEMAS_SSH + _SCHEMAS_GSI :

                                if  context.attribute_exists ("user_id") and context.user_id :
                                    info['user']  = context.user_id

                                if  context.attribute_exists ("user_key")  and  context.user_key  :
                                    info['ssh_args']  += "-o IdentityFile=%s " % context.user_key
                                    info['scp_args']  += "-o IdentityFile=%s " % context.user_key
                                    info['sftp_args'] += "-o IdentityFile=%s " % context.user_key

                                    if  context.attribute_exists ("user_pass") and context.user_pass :
                                        info['key_pass'][context.user_key] = context.user_pass

                        if  context.type.lower () == "userpass" :
                            if  info['schema'] in _SCHEMAS_SSH + _SCHEMAS_GSI :
                                if  context.attribute_exists ("user_id") and context.user_id :
                                    info['user']       = context.user_id
                                if  context.attribute_exists ("user_pass") and context.user_pass :
                                    info['pass']       = context.user_pass

                        if  context.type.lower () == "x509" :
                            if  info['schema'] in _SCHEMAS_GSI :

                                if  context.attribute_exists ("user_proxy")  and   context.user_proxy :
                                    info['ssh_env']   += "X509_USER_PROXY='%s' " % context.user_proxy
                                    info['scp_env']   += "X509_USER_PROXY='%s' " % context.user_proxy
                                    info['sftp_env']  += "X509_USER_PROXY='%s' " % context.user_proxy

                                if  context.attribute_exists ("user_cert")   and  context.user_cert :
                                    info['ssh_env']   += "X509_USER_CERT='%s' " % context.user_cert
                                    info['scp_env']   += "X509_USER_CERT='%s' " % context.user_cert
                                    info['sftp_env']  += "X509_USER_CERT='%s' " % context.user_cert

                                if  context.attribute_exists ("user_key")    and  context.user_key :
                                    info['ssh_env']   += "X509_USER_key='%s' "  % context.user_key
                                    info['scp_env']   += "X509_USER_key='%s' "  % context.user_key
                                    info['sftp_env']  += "X509_USER_key='%s' "  % context.user_key

                                if  context.attribute_exists ("cert_repository") and context.cert_repository :
                                    info['ssh_env']   += "X509_CERT_DIR='%s' "  % context.cert_repository
                                    info['scp_env']   += "X509_CERT_DIR='%s' "  % context.cert_repository
                                    info['sftp_env']  += "X509_CERT_DIR='%s' "  % context.cert_repository

                if url.port and url.port != -1 :
                    info['ssh_args']  += "-p %d " % int(url.port)
                    info['scp_args']  += "-p %d " % int(url.port)
                    info['sftp_args'] += "-P %d " % int(url.port)


                # all ssh based shells allow for user_id and user_pass from contexts
                # -- but the data given in the URL take precedence

                if url.username   :  info['user'] = url.username
                if url.password   :  info['pass'] = url.password

                ctrl_user = pwd.getpwuid (os.getuid ()).pw_name
                ctrl_base = "/tmp/saga_ssh_%s" % ctrl_user


                if  'user' in info and info['user'] :
                    info['host_str'] = "%s@%s"  % (info['user'], info['host_str'])
                    info['ctrl'] = "%s_%%h_%%p.%s.ctrl" % (ctrl_base, info['user'])
                else :
                    info['user'] = getpass.getuser ()
                    info['ctrl'] = "%s_%%h_%%p.ctrl" % (ctrl_base)

                info['m_flags']  = _SSH_FLAGS_MASTER % ({'share_mode' : info['share_mode'],
                                                         'ctrl'       : info['ctrl'],
                                                         'connect_timeout': info['ssh_connect_timeout']})
                info['s_flags']  = _SSH_FLAGS_SLAVE  % ({'share_mode' : info['share_mode'],
                                                         'ctrl'       : info['ctrl'],
                                                         'connect_timeout': info['ssh_connect_timeout']})

                logger.debug('SSH Connection M_FLAGS: %s' % info['m_flags'])
                logger.debug('SSH Connection S_FLAGS: %s' % info['s_flags'])
                
                # we want the userauth and hostname parts of the URL, to get the
                # scp-scope fs root.
                info['scp_root']  = ""
                has_auth          = False
                if  url.username :
                    info['scp_root'] += url.username
                    has_auth          = True
                if  url.password :
                    info['scp_root'] += ":"
                    info['scp_root'] += url.password
                    has_auth          = True
                if  has_auth :
                    info['scp_root'] += "@"
                info['scp_root']     += "%s:" % url.host

                # FIXME: port needs to be handled as parameter
              # if  url.port :
              #     info['scp_root'] += ":%d" % url.port


            # keep all collected info in the master dict, and return it for
            # registration
            return info


# ------------------------------------------------------------------------------

