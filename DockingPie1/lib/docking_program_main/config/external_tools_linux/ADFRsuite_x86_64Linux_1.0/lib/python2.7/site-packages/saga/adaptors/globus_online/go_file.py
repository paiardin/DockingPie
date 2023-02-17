
__author__    = "Andre Merzky, Mark Santcroos"
__copyright__ = "Copyright 2014-2015, The SAGA Project"
__license__   = "MIT"


""" (GSI)SSH based Globus Online Adaptor """

import os
import threading
import saga.utils.pty_shell as sups
import saga.utils.misc as sumisc

import saga.adaptors.base
import saga.adaptors.cpi.filesystem

from saga.adaptors.cpi.decorators import SYNC_CALL

# TODO: We could make this configurable,
# so people without gsissh can still perform some operations,
# as they could activate the endpoints through the globus web interface.
# go+ssh:// vs go+gsissh:// comes to mind
GO_DEFAULT_URL = "gsissh://cli.globusonline.org/"
#GO_DEFAULT_URL = "ssh://cli.globusonline.org/"


# --------------------------------------------------------------------
# the adaptor name
#
_ADAPTOR_NAME          = "saga.adaptor.globus_online_file"
_ADAPTOR_SCHEMAS       = ["go"] # TODO: also allow file:// ??
_ADAPTOR_OPTIONS       = [
    { 
    # fuck our config system!  I don't want these to be strings!  And its not
    # even using isinstance! :/
    'category'         : 'saga.adaptor.globus_online_file',
    'name'             : 'enable_notifications', 
    'type'             : str, 
    'default'          : 'None',
    'valid_options'    : ['True', 'False', 'None'],
    'documentation'    : '''Enable email notifications for all file transfers.
                            Note that 'True' and 'False' will result in
                            permanent changes to your GO notification settings.
                            'None' will leave your profile's settings upchanged.''',
    'env_variable'     : None
    },
    { 
    # fuck our config system!  I don't want these to be strings!  And its not
    # even using isinstance! :/
    'category'         : 'saga.adaptor.globus_online_file',
    'name'             : 'failure_mode', 
    'type'             : str, 
    'default'          : 'report',
    'valid_options'    : ['raise', 'report', 'ignore'],
    'documentation'    : '''Globus-Online seems to behave eratically.  This flag
                            defines how the adaptor should deal with intermittent
                            and fatal errors.  'raise' will cause exceptions on
                            all errors, 'report' will print error messages, but
                            otherwise continue, and 'ignore' will (duh!) ignore
                            errors.  'report' is the default, you should only use
                            'ignore' when you know what you are doing!''',
    'env_variable'     : None
    },
    # Set a GO specific prompt on the PTYShell
    {
        'category'      : 'saga.adaptor.globus_online_file',
        'name'          : 'prompt_pattern',
        'type'          : str,
        'default'       : '\$ $',
        'documentation' : 'GO always has a "$ " prompt',
        'env_variable'  : None
    },
    # Allow to specify which endpoint to use for "localhost"
    {
        'category'      : 'saga.adaptor.globus_online_file',
        'name'          : 'localhost_endpoint',
        'type'          : str,
        'default'       : 'None',
        'documentation' : '''Setting this value to the GO endpoint name of your
                             laptop/desktop/etc., allows you to refer to it by
                             "localhost" in your code. This enables consistency
                             with using your code with other adaptors.''',
        'env_variable'  : "SAGA_GLOBUS_ONLINE_LOCALHOST_ENDPOINT"
    }
]

# --------------------------------------------------------------------
# the adaptor capabilities & supported attributes
#
_ADAPTOR_CAPABILITIES  = {
    "metrics"          : [],
    "contexts"         : {"x509"     : "X509 proxy for Globus",
                          "userpass" : "username/password pair for GlobusOnline"}
}

# --------------------------------------------------------------------
# the adaptor documentation
#
_ADAPTOR_DOC           = {
    "name"             : _ADAPTOR_NAME,
    "cfg_options"      : _ADAPTOR_OPTIONS, 
    "capabilities"     : _ADAPTOR_CAPABILITIES,
    "description"      : """ 
        The globusonline file adaptor. This adaptor uses the GO file transfer
        service (https://www.globus.org/).
        """,
    "details"          : """ 
        """,
    "schemas"          : {"go": "use globus online for gridftp file transfer"}
}

# --------------------------------------------------------------------
# the adaptor info is used to register the adaptor with SAGA

_ADAPTOR_INFO          = {
    "name"             : _ADAPTOR_NAME,
    "version"          : "v0.2",
    "schemas"          : _ADAPTOR_SCHEMAS,
    "cpis"             : [
        { 
        "type"         : "saga.namespace.Directory",
        "class"        : "GODirectory"
        }, 
        { 
        "type"         : "saga.namespace.Entry",
        "class"        : "GOFile"
        },
        { 
        "type"         : "saga.filesystem.Directory",
        "class"        : "GODirectory"
        }, 
        { 
        "type"         : "saga.filesystem.File",
        "class"        : "GOFile"
        }
    ]
}

################################################################################
# The adaptor class

class Adaptor(saga.adaptors.base.Base):
    """ 
    This is the actual adaptor class, which gets loaded by SAGA (i.e. by the
    SAGA engine), and which registers the CPI implementation classes which
    provide the adaptor's functionality.
    """

    # --------------------------------------------------------------------------
    #
    def __init__(self):

        saga.adaptors.base.Base.__init__(self, _ADAPTOR_INFO, _ADAPTOR_OPTIONS)

        self.opts = self.get_config(_ADAPTOR_NAME)
        self.notify = self.opts['enable_notifications'].get_value()
        self.f_mode = self.opts['failure_mode'].get_value()
        self.prompt = self.opts['prompt_pattern'].get_value()
        self.localhost_ep = self.opts['localhost_endpoint'].get_value()
        self.shells = {}  # keep go shells for each session

        #
        # Lock to synchronize concurrent access to data structures
        #
        self.shell_lock = threading.RLock()

    # --------------------------------------------------------------------------
    #
    def sanity_check(self):
        pass

    # --------------------------------------------------------------------------
    #
    def get_go_shell(self, session, go_url=GO_DEFAULT_URL):

        # This returns a pty shell for: '[gsi]ssh username@cli.globusonline.org'
        #
        # X509 contexts are preferred, but ssh contexts, userpass and myproxy can
        # also be used.  If the given url has username / password encoded, we
        # create an userpass context out of it and add it to the (copy of) the
        # session.

        self._logger.debug("Acquiring lock")
        with self.shell_lock:
            self._logger.debug("Acquired lock")

            sid = session._id

            init = False
            create = False

            if sid in self.shells and self.shells[sid]['shell'].alive(recover=False):
                self._logger.debug("Shell in cache and alive, can reuse.")
            elif sid in self.shells:
                self._logger.debug("Shell in cache but not alive.")
                self.shells[sid]['shell'].finalize()
                self._logger.debug("Shell is finalized, need to recreate.")
                create = True
            else:
                self._logger.debug("Shell not in cache, create entry.")
                init = True
                create = True
                self.shells[sid] = {}

            # Acquire new shell
            if create:

                # deep copy URL (because of?)
                new_url = saga.Url(go_url)

                # GO specific prompt pattern
                opts = {'prompt_pattern': self.prompt}

                # create the shell.
                shell = sups.PTYShell(new_url, session=session, logger=self._logger, opts=opts, posix=False)
                self.shells[sid]['shell'] = shell

                # For this fresh shell, we get the list of public endpoints.
                # That list will contain the set of hosts we can potentially connect to.
                self.get_go_endpoint_list(session, shell, fetch=True)

            # Initialize other dict members and remote shell
            if init:
                shell = self.shells[sid]['shell']

                # Confirm the user ID for this shell
                self.shells[sid]['user'] = None
                _, out, _ = shell.run_sync('profile')
                for line in out.split('\n'):
                    if 'User Name:' in line:
                        self.shells[sid]['user'] = line.split(':', 2)[1].strip()
                        self._logger.debug("using account '%s'" % self.shells[sid]['user'])
                        break
                if not self.shells[sid]['user']:
                    raise saga.NoSuccess("Could not confirm user id")

                # Toggle notification
                if self.notify == 'True':
                    self._logger.debug("enable email notifications")
                    shell.run_sync('profile -n on')
                elif self.notify == 'False':
                    self._logger.debug("disable email notifications")
                    shell.run_sync('profile -n off')

            self._logger.debug("Release lock")

            # we have the shell for sure by now -- return it!
            return self.shells[sid]['shell']

    # ----------------------------------------------------------------
    #
    # Constructs the following from a SAGA/user URL:
    #
    #   - host: the original host element from the input URL
    #   - name: the EP name as used to interact with GO
    #   - url: the URL corrected with the EP name as the host element
    #
    def get_go_endpoint_ids(self, session, url):

        sid = session._id

        if not sid in self.shells:
            raise saga.IncorrectState("GO shell disconnected")

        ep_url = saga.Url()
        ep_url.schema = url.schema
        ep_url.port = url.port

        if '#' in url.host:
            # We assume an already existing EP.
            # We use the name for the host entry also.
            ep_name = url.host
            ep_url.host = ep_name
        elif url.host == 'localhost':
            # Special case to translate "localhost" to the users own machine.
            if self.localhost_ep != 'None':
                ep_name = self.localhost_ep
                ep_url.host = ep_name
            else:
                # TODO: could add more heuristics, like:
                # - looking for a "username#localhost"
                # - looking for special type of entries in the endpoint-list
                raise saga.BadParameter("localhost endpoint not configured")
        else:
            # Create an EP based on the username and hostname
            ep_name = "%s#%s" % (self.shells[sid]['user'], url.host)
            ep_url.host = url.host

        return ep_name, ep_url

    # ----------------------------------------------------------------
    #
    def get_path_spec(self, session, url, path=None, cwd_url=None, cwd_path=None):

        # we assume that, whenever we request a path spec, we also want to use
        # it, and we thus register and activate the endpoint, if needed.

        sid = session._id

        if not sid in self.shells:
            raise saga.IncorrectState("GO shell disconnected")

        shell = self.shells[sid]['shell']
        url   = saga.Url(url)

        if not path:
            path = url.path

        if not cwd_url:
            cwd_url = saga.Url(url)

            if not cwd_path:
                cwd_path = '.'
        else:
            if not cwd_path:
                cwd_path = cwd_url.path

        if not url.host:
            url.host = cwd_url.host
        if not url.schema:
            url.schema = cwd_url.schema

        if not url.host:
            raise saga.BadParameter('need host for GO ops')
        if not url.schema:
            raise saga.BadParameter('need schema for GO ops')

        ep_name, ep_url = self.get_go_endpoint_ids(session, url)

        # if both URLs point into the same namespace, and if the given path is
        # not absolute, then expand it relative to the cwd_path (if it exists).
        # Otherwise it is left to the unmodified path.
        ps_path = path
        # TODO: should the check be on cwd_url / url or on ep_url?
        if sumisc.url_is_compatible(cwd_url, url) and not path.startswith('/'):
            if cwd_path and path:
                ps_path = os.path.join(cwd_path, path)
            elif cwd_path:
                ps_path = cwd_path

        # the pathspec is the concatenation of ps_host and ps_path
        ps = "%s%s" % (ep_name, ps_path)

        # check if we know the endpoint in XXX, and create/activate as needed
        ep = self.get_go_endpoint(session, shell, ep_url)

        return ps

    # ----------------------------------------------------------------
    #
    def get_go_endpoint(self, session, shell, url):

        # for the given URL, derive the endpoint string.
        ep_name, ep_url = self.get_go_endpoint_ids(session, url)

        ep = self.get_go_endpoint_list(session, shell, ep_name, fetch=False)

        if not ep:

            # if not, check if it was created meanwhile (fetch again)
            ep = self.get_go_endpoint_list(session, shell, ep_name, fetch=True)

            if not ep:

                # Don't try to create endpoints that are supposed to be there
                # (plus, we can't!)
                if '#' in url.host:
                    raise saga.NoSuccess("# in hostname, not going to create!")

                # if not, create it, activate it, and refresh all entries
                shell.run_sync("endpoint-add %s -p %s" % (ep_name, ep_url))

                # refresh endpoint entries again
                ep = self.get_go_endpoint_list(session, shell, ep_name, fetch=True)

                if not ep:
                    # something above must have failed ...
                    raise saga.NoSuccess("endpoint initialization failed")

        # we have the endpoint now, for sure -- make sure its activated
        if not ep['Credential Status'] == 'ACTIVE':
            # TODO: I had an active endpoint, but still got an activation prompt,
            # probably because the remaining lifetime was not very long anymore.
            # or Credential Time Left    : 00:16:35 < ????
            # Answer: below 30 min there is a activation prompt,
            # but it does actually continue normally.
            # Need to capture that behavior.

            # Only Globus Connect Service Endpoints don't need -g?
            # Had contact on this with Globus Support, they couldn't suggest
            # anything better.
            if ep['MyProxy Server'] == 'myproxy.globusonline.org' and \
            '/C=US/O=Globus Consortium/OU=Globus Connect Service/CN=' in ep['Credential Subject']:
                shell.run_sync("endpoint-activate %s" % ep_name)
            else:
                shell.run_sync("endpoint-activate -g %s" % ep_name)

            # reload list to check status
            ep = self.get_go_endpoint_list(session, shell, ep_name, fetch=True)

            if not ep['Credential Status'] == 'ACTIVE':
                raise saga.AuthorizationFailed("endpoint activation failed")

        return ep
        
    # ----------------------------------------------------------------
    #
    def get_go_endpoint_list(self, session, shell, ep_name=None, fetch=False):

        # if 'fetch' is True, query the shell for an updated endpoint list.
        # then check if the given ep_name is a known endpoint name, and if so,
        # return that entry -- otherwise return None.  If no ep_name is given,
        # and fetch is True, we thus simply refresh the internal list.

        if fetch:

            with self.shell_lock:

                endpoints = {}
                name = None

                if ep_name:
                    endpoint_selection = ep_name
                else:
                    endpoint_selection = '-a'

                # Get the details of endpoints _OWNED_ by user
                _, out, _ = shell.run_sync("endpoint-details %s -f "
                                           "legacy_name,"        # Legacy Name
                                           "credential_status,"  # Credential Status
                                           "credential_subject," # Credential Subject
                                           "myproxy_server"     # MyProxy Server
                                           % endpoint_selection
                                           )

                for line in out.split('\n'):
                    elems = line.split(':', 1)

                    if len(elems) != 2:
                        continue

                    key = elems[0].strip()
                    val = elems[1].strip()

                    if not key or not val:
                        continue

                    if key == "Legacy Name":

                        # we now operate on a new entry -- initialize it
                        name = val

                        endpoints[name] = {}
                        endpoints[name]['Legacy Name']       = name

                    else:

                        # Continued passed of an entry, name should always exist
                        try:
                            endpoints[name][key] = val
                        except:
                            raise saga.NoSuccess("No entry to operate on for: %s[%s]" % (key,val))

                # replace the ep info dist with the new one, to clean out old entries.
                # TODO: merge and not replace(?)
                self.shells[session._id]['endpoints'] = endpoints

        if ep_name:
            # return the requested entry, or None
            return self.shells[session._id]['endpoints'].get(ep_name, None)

    # ----------------------------------------------------------------
    #
    def run_go_cmd(self, shell, cmd, mode=None):

        # available modes:
        #   raise : raise NoSuccess on error
        #   report: print error message, but continue
        #   ignore: do nothing

        if not mode:
            mode = self.f_mode

        _, out, err = shell.run_sync(cmd)

        # see if the second line starts with 'Error'.  Note that this assumes
        # that the command only is one line...
        lines = out.split('\n')

        if len(lines) > 1:
            if lines[1].startswith('Error:'):
                err = "%s\n%s" % (err, '\n'.join(lines))
                out = None

            else:
                # on success, we always remove the first line, as that is not
                # part of the output, really (this shell does not support 
                # 'stty -echo'...
                out = '\n'.join(lines[1:])

        if err:
            if mode == 'raise':
                if lines[3].startswith('Message:'):
                    cause = lines[3].split(':')[1].strip()
                    if cause == 'No such file or directory':
                        raise saga.DoesNotExist(cause)
                    elif cause == 'Could not connect to server':
                        raise saga.BadParameter(cause)
                    else:
                        raise saga.NoSuccess('Unknown error: %s' % cause)
                else:
                    raise saga.NoSuccess('Could not find find message in error: %s' % err)

                # TODO: Handle GO access to directories that are not allowed
                # TODO: by Globus Personal by default (e.g. /var/ , /tmp)
                # TODO: '''Message: Fatal FTP Response
                # TODO:    ---
                # TODO:    500 Command failed : Path not allowed.'''

            if mode == 'report':
                self._logger.error("Error in '%s': %s" % (cmd, err))

            if mode == 'silent':
                pass

        return out, err

    # ----------------------------------------------------------------
    #
    def mkparents(self, shell, tgt_ps):

        # GO does not support mkdir -p, so we need to split the dir into
        # elements and create one after the other, ignoring errors for already
        # existing elements.
        # TODO: Can't we check for existence? The errors are confusing.

        host_ps, path_ps = tgt_ps.split('/', 1)
        path_ps = '/' + path_ps

        self._logger.info('mkparents %s' % path_ps)

        if path_ps.startswith('/'):
            cur_path = ''
        else:
            cur_path = '.'

        error = {}
        path_elems = filter(None, path_ps.split('/'))

        for path_elem in path_elems :

            cur_path = "%s/%s" % (cur_path, path_elem)
            out, err = self.run_go_cmd(shell, "mkdir %s%s" % (host_ps, cur_path))

            if err:
                error[cur_path] = err

        if len(error):

            # some mkdir gave an error.  Check if the error occurred on the last
            # dir (the tgt), and if that is not a ignorable report that it
            # already exists -- anything else will raise an exception though...
            if cur_path in error:

                if not 'Path already exists' in error[cur_path]:

                    if self.f_mode == 'raise':
                        # TODO: a 'translate_exception' call would be useful here...
                        # TODO: We can use the exceptions from run_go_cmd?
                        raise saga.NoSuccess("Could not make dir hierarchy: %s" % str(error))

                    if self.f_mode == 'report':
                        self._logger.error("Could not make dir hierarchy: %s" % str(error))

                    if self.f_mode == 'silent':
                        pass

    def go_transfer(self, shell, flags, source, target):

        # TODO: I dont think we handle relative targets yet

        self._logger.debug('Adaptor:go_transfer(%s, %s)' % (source, target))

        # 0: Copy files that do not exist at the destination
        sync_level = 0

        # Create parents
        if flags & saga.filesystem.CREATE_PARENTS:
            self.mkparents(shell, os.path.dirname(target))

        #if flags & saga.filesystem.OVERWRITE:
            # 1: Copy files if the size of the destination does not match the
            #    size of the source
            # 2: Copy files if the timestamp of the destination is older than
            #    the timestamp of the source
            # 3: Copy files if checksums of the source and destination do not match
        #    sync_level = 3

        # Set recursive flag
        cmd_flags = ""
        if flags & saga.filesystem.RECURSIVE:
            cmd_flags += "-r"

        # Initiate background copy
        # TODO: Should we use a deadline?
        cmd = "transfer %s -s %d -- '%s' '%s'" % (cmd_flags, sync_level, source, target)
        out, _ = self.run_go_cmd(shell, cmd)
        # 'Task ID: 8c6f989d-b6aa-11e4-adc6-22000a97197b'
        key, value = out.split(':')
        if key != 'Task ID':
            raise Exception("Expected Task ID: <id>, got %s" % out)
        task_id = value.strip()

        # Wait until background copy has finished
        cmd = "wait -q %s" % task_id
        self.run_go_cmd(shell, cmd)

        # Retrieve task status
        cmd = "status -f status %s" % task_id
        out, _ = self.run_go_cmd(shell, cmd)
        # Status: SUCCEEDED
        key, value = out.split(':')
        if key != 'Status':
            raise Exception("Expected Status: <status>, got %s" % out)
        status = value.strip()

        # Validate task status
        if status == 'SUCCEEDED':
            # The task completed successfully.
            return
        elif status == 'ACTIVE':
            # The task is in progress.
            raise Exception('Task still active, this should not happen after wait')
        elif status == 'INACTIVE':
            # The task has been suspended and will not continue without intervention.
            # Currently, only credential expiration will cause this state.
            raise Exception('Task is inactive, probably credentials have expired')
        elif status == 'FAILED':
            # The task or one of its subtasks failed, expired, or was canceled.
            raise Exception('Task failed')
        else:
            raise Exception('Unknown status: %s' % status)

    ################################################################################
    #
    # Helper function to test for existence and type.
    # Will raise DoesNotExist for non-existing entries.
    def stat(self, shell, ps):

        # TODO: stat() would probably benefit from some caching

        out, err = self.run_go_cmd(shell, "ls -la '%s'" % ps, mode='raise')

        mode = out.split('\n')[0].split()[0][0]
        if mode == '-':
            mode = 'file'
        elif mode == 'd':
            mode = 'dir'
        elif mode == 'l':
            mode = 'link'
        else:
            raise saga.NoSuccess("stat() unknown mode: '%s' (%s)" % (mode, out))

        size = int(out.split('\n')[0].split()[3])

        return {
            'mode': mode,
            'size': size
        }


################################################################################
#
class GODirectory(saga.adaptors.cpi.filesystem.Directory):
    """ Implements saga.adaptors.cpi.filesystem.Directory """

    # ----------------------------------------------------------------
    #
    def __init__(self, api, adaptor):

        _cpi_base = super(GODirectory, self)
        _cpi_base.__init__(api, adaptor)

    # ----------------------------------------------------------------
    #
    def _is_valid(self):

        if not self.valid:
            raise saga.IncorrectState("this instance was closed or removed")

    # ----------------------------------------------------------------
    #
    @SYNC_CALL
    def init_instance(self, adaptor_state, url, flags, session):
        """ Directory instance constructor """

        # TODO: eval flags!
        if flags == None:
            flags = 0

        self.orig     = saga.Url(url) # deep copy
        self.url      = saga.Url(url) # deep copy
        self.path     = url.path       # keep path separate
        self.url.path = None

        self.flags    = flags
        self.session  = session
        self.valid    = False # will be set by initialize

        self.initialize()

        return self.get_api()

    # ----------------------------------------------------------------
    #
    def initialize(self):

        # GO shell got started, found its prompt.  Now, change
        # to the initial (or later current) working directory.

        self.shell = self._adaptor.get_go_shell(self.session)
        self.ep    = self._adaptor.get_go_endpoint(self.session, self.shell, self.url)

        ps = self.get_path_spec()

        if not self.ep:
            raise saga.BadParameter("invalid dir '%s'" % ps)

        if self.flags & saga.filesystem.CREATE_PARENTS:
            self._adaptor.mkparents(self.shell, ps)

        elif self.flags & saga.filesystem.CREATE:
            # TODO: check for errors?
            self._adaptor.run_go_cmd(self.shell, "mkdir '%s'" % ps)

        else:
            stat = self._adaptor.stat(self.shell, ps)
            if stat['mode'] not in ['dir', 'link']:
                # TODO: if link, check the target
                raise saga.IncorrectState('Is not a directory')

        self._logger.debug("Initialized directory %s/%s" % (self.url, self.path))

        self.valid = True

    # ----------------------------------------------------------------
    #
    def finalize(self, kill=False):

        if kill and self.shell:
            self.shell.finalize(True)
            self.shell = None

        self.valid = False

    # ----------------------------------------------------------------
    #
    def get_path_spec(self, url=None, path=None):

        return self._adaptor.get_path_spec(session  = self.session,
                                           url      = url,
                                           path     = path,
                                           cwd_url  = self.url,
                                           cwd_path = self.path)

    # ----------------------------------------------------------------
    #
    @SYNC_CALL
    def open(self, url, flags):

        self._is_valid()

        adaptor_state = { "from_open" : True,
                          "url"       : saga.Url(self.url),   # deep copy
                          "path"      : self.path}

        if sumisc.url_is_relative(url):
            url = sumisc.url_make_absolute(self.get_url(), url)

        return saga.filesystem.File(url=url, flags=flags, session=self.session,
                                     _adaptor=self._adaptor, _adaptor_state=adaptor_state)

    # ----------------------------------------------------------------
    #
    @SYNC_CALL
    def open_dir(self, url, flags):

        self._is_valid()

        adaptor_state = {"from_open": True,
                         "url"      : saga.Url(self.url),   # deep copy
                         "path"     : self.path}

        if sumisc.url_is_relative(url):
            url = sumisc.url_make_absolute(self.get_url(), url)

        return saga.filesystem.Directory(url=url, flags=flags, session=self.session,
                _adaptor=self._adaptor, _adaptor_state=adaptor_state)

    # ----------------------------------------------------------------
    #
    @SYNC_CALL
    def change_dir(self, tgt, flags):

        tgt_url = saga.Url(tgt)

        # TODO: attempt to get new EP
        if not sumisc.url_is_compatible(self.url, tgt_url):
            raise saga.BadParameter("Target dir outside of namespace '%s': %s" \
                                  % (self.url, tgt_url))

        if sumisc.url_is_relative(tgt_url):

            self.path      = tgt_url.path
            self.orig.path = self.path

        else:
            self.orig      = saga.Url(tgt_url)
            self.url       = tgt_url
            self.path      = self.url.path
            self.url.path  = None

        self.initialize()

        self._logger.debug("changed directory (%s)" % (tgt))

    # ----------------------------------------------------------------
    #
    @SYNC_CALL
    def close(self, timeout=None):

        if timeout:
            raise saga.BadParameter("timeout for close not supported")

        self.finalize(kill=True)

    # ----------------------------------------------------------------
    @SYNC_CALL
    def get_url(self):

        self._is_valid()

        return saga.Url(self.orig) # deep copy

    # ----------------------------------------------------------------
    #
    @SYNC_CALL
    def list(self, npat, flags):

        self._is_valid()

        npat_ps  = self.get_path_spec(url=npat)
        # TODO: catch errors?
        out, err = self._adaptor.run_go_cmd(self.shell, "ls '%s'" % (npat_ps))
        lines = filter(None, out.split("\n"))
        self._logger.debug(lines)

        self.entries = []
        for line in lines:
            self.entries.append(saga.Url(line.strip()))

        return self.entries
   
    # ----------------------------------------------------------------
    #
    @SYNC_CALL
    def copy_self(self, tgt, flags):

        self._is_valid()

        self._logger.debug('Directory:copy_self(%s)' % tgt)

        return self.copy(src_in=None, tgt_in=tgt, flags=flags)

    # ----------------------------------------------------------------
    #
    @SYNC_CALL
    def copy(self, src_in, tgt_in, flags, _from_task=None):

        self._logger.debug('Directory:copy(%s, %s)' % (src_in, tgt_in))
        self._is_valid()

        src_ps = self.get_path_spec(url=src_in)
        tgt_ps = self.get_path_spec(url=tgt_in)

        # Check for existence of source
        self._adaptor.stat(self.shell, src_ps)

        self._adaptor.go_transfer(self.shell, flags, src_ps, tgt_ps)

    # ----------------------------------------------------------------
    #
    @SYNC_CALL
    def move_self(self, tgt, flags):

        return self.move(src_in=None, tgt_in=tgt, flags=flags)
   
    # ----------------------------------------------------------------
    #
    @SYNC_CALL
    def move(self, src_in, tgt_in, flags):

        # if src and target are on the same endpoint, we might get away with an
        # actual 'rename' command -- in all other cases (or if rename failed),
        # we fake move as non-atomic copy/remove...

        src_ps = self.get_path_spec(url=src_in)
        tgt_ps = self.get_path_spec(url=tgt_in)

        src_ep_str = src_ps.split('/', 1)[0]
        tgt_ep_str = tgt_ps.split('/', 1)[0]

        if src_ep_str == tgt_ep_str:

            try:
                self._adaptor.run_go_cmd(self.shell, "rename '%s' '%s'" % (src_ps, tgt_ps))
                return
            except:
                self._logger.warn("rename op failed -- retry as copy/remove")

        # either the op spans endpoints, or the 'rename' op failed
        self.copy(src_in, tgt_in, flags)
        self.remove(src_in, flags)
   
    # ----------------------------------------------------------------
    #
    @SYNC_CALL
    def remove_self(self, flags):

        self._is_valid()

        self.remove(self.url, flags)
        self.invalid = True
   
    # ----------------------------------------------------------------
    #
    @SYNC_CALL
    def remove(self, tgt_in, flags):

        self._is_valid()

        tgt_ps = self.get_path_spec(url=tgt_in)

        cmd_flags = ""
        if flags & saga.filesystem.RECURSIVE:
            cmd_flags += "-r"

        # TODO: check for errors

        # TODO: a dir only gets removed (on some endpoints)
        # if the trailing '/' is specified -- otherwise the op *silently fails*!
        # Oh well, since we don't really (want to) know if the target is a dir
        # or not, we remove both versions... :/
        cmd      = "rm %s -f '%s/'" % (cmd_flags, tgt_ps)
        out, err = self._adaptor.run_go_cmd(self.shell, cmd)

        cmd      = "rm %s -f '%s'"  % (cmd_flags, tgt_ps)
        out, err = self._adaptor.run_go_cmd(self.shell, cmd, mode='ignore')

    # ----------------------------------------------------------------
    #
    @SYNC_CALL
    def make_dir(self, tgt_in, flags):

        self._is_valid()

        tgt_ps = self.get_path_spec(url=tgt_in)

        # TODO: check for errors

        if flags & saga.filesystem.CREATE_PARENTS:
            self._adaptor.mkparents(self.shell, tgt_ps)

        else:
            cmd = "mkdir '%s'" % tgt_ps
            self._adaptor.run_go_cmd(self.shell, cmd)

    # ----------------------------------------------------------------
    #
    @SYNC_CALL
    def get_size_self(self):

        self._is_valid()

        return self.get_size(self.url)

    # ----------------------------------------------------------------
    #
    @SYNC_CALL
    def get_size(self, tgt_in):

        self._is_valid()

        tgt_ps = self.get_path_spec(url=tgt_in)

        stat = self._adaptor.stat(self.shell, tgt_ps)

        return stat['size']

    # ----------------------------------------------------------------
    #
    @SYNC_CALL
    def is_dir_self(self):

        self._is_valid()

        return self.is_dir(self.url)

    # ----------------------------------------------------------------
    #
    @SYNC_CALL
    def is_dir(self, tgt_in):

        self._is_valid()

        tgt_ps = self.get_path_spec(url=tgt_in)

        stat = self._adaptor.stat(self.shell, tgt_ps)

        if stat['mode'] == 'dir':
            return True
        else:
            return False

    # ----------------------------------------------------------------
    #
    @SYNC_CALL
    def is_link_self(self):

        self._is_valid()

        return self.is_link(self.url)

    # ----------------------------------------------------------------
    #
    @SYNC_CALL
    def is_link(self, tgt_in):

        self._is_valid()

        tgt_ps = self.get_path_spec(url=tgt_in)

        stat = self._adaptor.stat(self.shell, tgt_ps)

        if stat['mode'] == 'link':
            return True
        else:
            return False

    # ----------------------------------------------------------------
    #
    @SYNC_CALL
    def is_file_self(self):

        self._is_valid()

        return self.is_link(self.url)

    # ----------------------------------------------------------------
    #
    @SYNC_CALL
    def is_file(self, tgt_in):

        self._is_valid()

        tgt_ps = self.get_path_spec(url=tgt_in)

        stat = self._adaptor.stat(self.shell, tgt_ps)

        if stat['mode'] == 'file':
            return True
        else:
            return False

    # ----------------------------------------------------------------
    #
    @SYNC_CALL
    def exists(self, tgt_in):

        self._is_valid()

        tgt_ps = self.get_path_spec(url=tgt_in)

        try:
            self._adaptor.stat(self.shell, tgt_ps)
        except saga.DoesNotExist:
            return False

        return True


###############################################################################
#
class GOFile(saga.adaptors.cpi.filesystem.File):
    """ Implements saga.adaptors.cpi.filesystem.File
    """
    # ----------------------------------------------------------------
    #
    def __init__(self, api, adaptor):

        _cpi_base = super(GOFile, self)
        _cpi_base.__init__(api, adaptor)

    # ----------------------------------------------------------------
    #
    def _is_valid(self):

        if not self.valid:
            raise saga.IncorrectState("this instance was closed or removed")

    # ----------------------------------------------------------------
    #
    @SYNC_CALL
    def init_instance(self, adaptor_state, url, flags, session):
        """ File instance constructor """

        # TODO: eval flags!
        if flags == None:
            flags = 0

        self.orig     = saga.Url(url) # deep copy
        self.url      = saga.Url(url) # deep copy
        self.path     = url.path       # keep path separate
        self.cwd      = sumisc.url_get_dirname(self.url)
        self.url.path = None

        self.flags    = flags
        self.session  = session
        self.valid    = False # will be set by initialize

        self.initialize()

        return self.get_api()

    # ----------------------------------------------------------------
    #
    def initialize(self):

        # GO shell got started, found its prompt.  Now, change
        # to the initial (or later current) working directory.

        self.shell = self._adaptor.get_go_shell(self.session)
        self.ep    = self._adaptor.get_go_endpoint(self.session, self.shell, self.url)
        ps         = self.get_path_spec()
        cwd_ps     = self.get_path_spec(path=self.cwd)

        if not self.ep:
            raise saga.BadParameter("invalid file '%s'" % ps)

        if self.flags & saga.filesystem.CREATE_PARENTS:
            self._adaptor.mkparents(self.shell, cwd_ps)

        elif self.flags & saga.filesystem.CREATE:
            self._logger.error("CREATE not supported for files via globus online")

        else:
            stat = self._adaptor.stat(self.shell, ps)
            if stat['mode'] not in ['file', 'link']:
                # TODO: if link, check the target
                raise saga.IncorrectState('Is not a (regular) file')

        self._logger.debug("Initialized file %s/%s" % (self.url, self.path))

        self.valid = True

    # ----------------------------------------------------------------
    #
    def get_path_spec(self, url=None, path=None):

        return self._adaptor.get_path_spec(
            session  = self.session,
            url      = url,
            path     = path,
            cwd_url  = self.url,
            cwd_path = self.path
        )

    # ----------------------------------------------------------------
    #
    def finalize(self, kill=False):

        if kill and self.shell:
            self.shell.finalize(True)
            self.shell = None

        self.valid = False

    # ----------------------------------------------------------------
    #
    @SYNC_CALL
    def close(self, timeout=None):

        if timeout:
            raise saga.BadParameter("timeout for close not supported")

        self.finalize(kill=True)

    # ----------------------------------------------------------------
    #
    @SYNC_CALL
    def get_url(self):

        self._is_valid()

        return saga.Url(self.orig) # deep copy

    # ----------------------------------------------------------------
    #
    @SYNC_CALL
    def copy_self(self, tgt_in, flags):

        self._is_valid()

        src_ps = self.get_path_spec()
        tgt_ps = self.get_path_spec(url=tgt_in)

        self._adaptor.go_transfer(self.shell, flags, src_ps, tgt_ps)

    # ----------------------------------------------------------------
    #
    @SYNC_CALL
    def move_self(self, tgt_in, flags):

        # if src and target are on the same endpoint, we might get away with an
        # actual 'rename' command -- in all other cases (or if rename failed),
        # we fake move as non-atomic copy/remove...

        src_ps = self.get_path_spec()
        tgt_ps = self.get_path_spec(url=tgt_in)

        src_ep_str = src_ps.split('/', 1)[0]
        tgt_ep_str = tgt_ps.split('/', 1)[0]

        # TODO: check for errors

        # TODO: should me extracted and merged with operations on Directory?

        if src_ep_str == tgt_ep_str:

            try:
                self._adaptor.run_go_cmd(self.shell, "rename '%s' '%s'" % (src_ps, tgt_ps))
                return
            except:
                self._logger.warn("Rename op failed -- retry as copy/remove")

        # either the op spans endpoints, or the 'rename' op failed
        self.copy_self(tgt_in, flags)
        self.remove_self(flags)

    # ----------------------------------------------------------------
    #
    @SYNC_CALL
    def remove_self(self, flags):

        self._is_valid()

        cmd_flags = ""
        if flags & saga.filesystem.RECURSIVE:
            cmd_flags += "-r"

        # TODO: check for errors

        # TODO: should me extracted and merged with operations on Directory?

        tgt_ps   = self.get_path_spec()
        cmd      = "rm %s -f '%s'"  % (cmd_flags, tgt_ps)
        out, err = self._adaptor.run_go_cmd(self.shell, cmd, mode='ignore')

    # ----------------------------------------------------------------
    #
    @SYNC_CALL
    def get_size_self(self):

        self._is_valid()

        tgt_ps = self.get_path_spec(url=self.url)

        stat = self._adaptor.stat(self.shell, tgt_ps)

        return stat['size']

    # ----------------------------------------------------------------
    #
    @SYNC_CALL
    def is_link_self(self):
        self._is_valid()

        tgt_ps = self.get_path_spec(url=self.url)

        stat = self._adaptor.stat(self.shell, tgt_ps)
        if stat['mode'] == 'link':
            return True
        else:
            return False

    # ----------------------------------------------------------------
    #
    @SYNC_CALL
    def is_file_self(self):
        self._is_valid()

        tgt_ps = self.get_path_spec(url=self.url)

        stat = self._adaptor.stat(self.shell, tgt_ps)
        if stat['mode'] == 'file':
            return True
        else:
            return False

    # ----------------------------------------------------------------
    #
    @SYNC_CALL
    def is_dir_self(self):

        # TODO: can probably be further extracted and merged with Directory calls

        self._is_valid()

        tgt_ps = self.get_path_spec(url=self.url)

        stat = self._adaptor.stat(self.shell, tgt_ps)
        if stat['mode'] == 'dir':
            return True
        else:
            return False
