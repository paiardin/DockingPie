
__author__    = "Andre Merzky, Ole Weidner, Alexander Grill"
__copyright__ = "Copyright 2012-2013, The SAGA Project"
__license__   = "MIT"


""" shell based file adaptor implementation """

import saga.utils.pty_shell as sups
import saga.utils.misc      as sumisc

import saga.adaptors.base
import saga.adaptors.cpi.filesystem

from   saga.filesystem.constants import *

import radical.utils as ru

import re
import os
import errno

import shell_wrapper

SYNC_CALL  = saga.adaptors.cpi.decorators.SYNC_CALL
ASYNC_CALL = saga.adaptors.cpi.decorators.ASYNC_CALL


# ------------------------------------------------------------------------------
# the adaptor name
#
_ADAPTOR_NAME          = "saga.adaptor.shell_file"
_ADAPTOR_SCHEMAS       = ["file", "local", "sftp", "gsisftp", "ssh", "gsissh"]
_ADAPTOR_OPTIONS       = [
  # { 
  # 'category'         : 'saga.adaptor.shell_file',
  # 'name'             : 'enable_notifications', 
  # 'type'             : bool, 
  # 'default'          : False,
  # 'valid_options'    : [True, False],
  # 'documentation'    : '''Enable support for filesystem notifications.  Note that
  #                       enabling this option will create a local thread and a remote 
  #                       shell process.''',
  # 'env_variable'     : None
  # }
]

# ------------------------------------------------------------------------------
# the adaptor capabilities & supported attributes
#
_ADAPTOR_CAPABILITIES  = {
    "metrics"          : [],
    "contexts"         : {"ssh"      : "public/private keypair",
                          "x509"     : "X509 proxy for gsissh",
                          "userpass" : "username/password pair for ssh"}
}

# ------------------------------------------------------------------------------
# the adaptor documentation
#
_ADAPTOR_DOC           = {
    "name"             : _ADAPTOR_NAME,
    "cfg_options"      : _ADAPTOR_OPTIONS, 
    "capabilities"     : _ADAPTOR_CAPABILITIES,
    "description"      : """ 
        The shell file adaptor. This adaptor uses the sh command line tools (sh,
        ssh, gsissh) to access remote filesystems.
        """,
    "details"          : """ 
        A more elaborate description....

        Known Limitations:
        ------------------

          * number of system pty's are limited:  each object bound
            to this adaptor will use 2 pairs of pty pipes.  Systems usually
            limit the number of available pty's to 1024 .. 4096.  Given that
            other processes also use pty's , that gives a hard limit to the number
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

 
          * Other system limits (memory, CPU, selinux, accounting etc.) apply as
            usual.


          * thread safety: it is safe to create multiple ``filesystem.*``
            instances to the same target host at a time -- they should not
            interfere with each other.

        """,
    "schemas"          : {"file"    :"use /bin/sh to access local filesystems",
                          "local"   :"alias for file://",
                          "ssh"     :"use sftp to access remote filesystems",
                          "sftp"    :"use sftp to access remote filesystems",
                          "gsissh"  :"use gsisftp to access remote filesystems",
                          "gsisftp" :"use gsisftp to access remote filesystems"}
}

# ------------------------------------------------------------------------------
# the adaptor info is used to register the adaptor with SAGA

_ADAPTOR_INFO          = {
    "name"             : _ADAPTOR_NAME,
    "version"          : "v0.1",
    "schemas"          : _ADAPTOR_SCHEMAS,
    "cpis"             : [
        { 
        "type"         : "saga.namespace.Directory",
        "class"        : "ShellDirectory"
        }, 
        { 
        "type"         : "saga.namespace.Entry",
        "class"        : "ShellFile"
        },
        { 
        "type"         : "saga.filesystem.Directory",
        "class"        : "ShellDirectory"
        }, 
        { 
        "type"         : "saga.filesystem.File",
        "class"        : "ShellFile"
        }
    ]
}

# ------------------------------------------------------------------------------
#
# use the native python facilities to create local directories
#
def _mkdir_p(path):

    try:
        os.makedirs(path)

    except OSError as exc:  # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass

    except Exception as e:
        raise saga.NoSuccess("mkdir failed for '%s': %s" \
                              % (path, e))


###############################################################################
# The adaptor class

class Adaptor (saga.adaptors.base.Base):
    """ 
    This is the actual adaptor class, which gets loaded by SAGA (i.e. by the
    SAGA engine), and which registers the CPI implementation classes which
    provide the adaptor's functionality.
    """


    # --------------------------------------------------------------------------
    #
    def __init__ (self) :

        saga.adaptors.base.Base.__init__ (self, _ADAPTOR_INFO, _ADAPTOR_OPTIONS)

        self.opts  = self.get_config (_ADAPTOR_NAME)


    # --------------------------------------------------------------------------
    #
    def sanity_check (self) :

        # FIXME: also check for gsissh

        pass


    # --------------------------------------------------------------------------
    #
    def get_lease_target (self, tgt) :

        """
        return a URL with empty path which can be used to identify leased copy
        shells (we don't care about the path while leasing copy shells)
        """

        lease_tgt = saga.Url (tgt)
        lease_tgt.path = '/shell_file_adaptor_command_shell/'

        return lease_tgt




###############################################################################
#
class ShellDirectory (saga.adaptors.cpi.filesystem.Directory) :
    """ Implements saga.adaptors.cpi.filesystem.Directory """

    # --------------------------------------------------------------------------
    #
    def __init__ (self, api, adaptor) :

        _cpi_base = super  (ShellDirectory, self)
        _cpi_base.__init__ (api, adaptor)


    # --------------------------------------------------------------------------
    #
    def __del__ (self) :

        self.finalize (kill=True)


    # --------------------------------------------------------------------------
    #
    def _is_valid (self) :

        if  not self.valid :
            raise saga.IncorrectState ("this instance was closed or removed")

    # --------------------------------------------------------------------------
    #
    def _create_parent (self, cwdurl, tgt) :

        dirname = sumisc.url_get_dirname (tgt)
        ret     = None
        out     = None

        if  sumisc.url_is_compatible (cwdurl, tgt) :

            ret, out, _ = self._command (" mkdir -p '%s'\n" % (dirname), make_location=True)
            if  ret != 0 :
                raise saga.NoSuccess ("failed at mkdir '%s': (%s) (%s)" \
                                   % (dirname, ret, out))

        elif sumisc.url_is_local (tgt) :

            if  tgt.scheme and not tgt.scheme.lower () in _ADAPTOR_SCHEMAS :
                raise saga.BadParameter ("schema of mkdir target is not supported (%s)" \
                                      % (tgt))
            _mkdir_p(dirname)

        else :
            lease_tgt = self._adaptor.get_lease_target (tgt)
            with self.lm.lease (lease_tgt, self.shell_creator, tgt) as tmp_shell :
                tmp_shell.run_sync ('mkdir -p %s' % dirname)


    # --------------------------------------------------------------------------
    #
    @SYNC_CALL
    def init_instance (self, adaptor_state, url, flags, session) :
        """ Directory instance constructor """

        if  flags == None :
            flags = 0

        self.url         = saga.Url (url) # deep copy
        self.flags       = flags
        self.session     = session
        self.valid       = False # will be set by initialize
        self.lm          = session._lease_manager

        # Use `_set_session` method of the base class to set the session object.
        # `_set_session` and `get_session` methods are provided by `CPIBase`.
        self._set_session(session)

        def _shell_creator (url) :
            return sups.PTYShell (url, self.get_session(), self._logger)
        self.shell_creator = _shell_creator

        # The dir command shell is leased, as the dir seems to be used
        # extensively in some cases.  Note that before each command, we need to
        # perform a 'cd' to the target location, to make sure we operate in the
        # right location (see self._command())

        self.initialize ()

        return self.get_api ()

    # --------------------------------------------------------------------------
    #
    def _command (self, command, location=None, make_location=False) :

        if  not location :
            location = self.url
        else :
            location = saga.Url (location)

        lease_tgt = self._adaptor.get_lease_target (location)
        with self.lm.lease (lease_tgt, self.shell_creator, location) \
             as cmd_shell :

            if  make_location and location.path :
                pre_cmd = "mkdir -p %s &&" % location.path
            else :
                pre_cmd = ""
             
            return cmd_shell.run_sync ("%s cd %s && %s" % (pre_cmd, location.path, command))


    # --------------------------------------------------------------------------
    #
    def initialize (self) :

        # shell got started, found its prompt.  Now, change
        # to the initial (or later current) working directory.

        cmd = ""
        mkl = False    # make location on command execution

        if  self.flags & saga.filesystem.CREATE_PARENTS :
            cmd = " mkdir -p '%s' ;  cd '%s'" % (self.url.path, self.url.path)
            mkl = True
        elif self.flags & saga.filesystem.CREATE :
            cmd = " mkdir    '%s' ;  cd '%s'" % (self.url.path, self.url.path)
            mkl = False
        else :
            cmd = " test -d  '%s' && cd '%s'" % (self.url.path, self.url.path)
            mkl = False

        # we can't yet 'cd' to self.url if we are to create it, so we start at
        # root on the target FS.
        if mkl:
            root = saga.Url(self.url)
            root.path = '/'
            ret, out, _ = self._command(cmd, location=root, make_location=mkl)
        else:
            ret, out, _ = self._command(cmd, make_location=mkl)

        if  ret != 0 :
            raise saga.BadParameter ("invalid dir '%s': %s" % (self.url.path, out))

        self._logger.debug ("initialized directory (%s)(%s)" % (ret, out))

        self.valid = True


    # --------------------------------------------------------------------------
    #
    def finalize (self, kill = False) :

        self.valid = False


    # --------------------------------------------------------------------------
    #
    @SYNC_CALL
    def open (self, url, flags) :

        self._is_valid ()

        adaptor_state = { "from_open" : True,
                          "cwd"       : saga.Url(self.url) }  # deep copy

        if  sumisc.url_is_relative (url) :
            url = sumisc.url_make_absolute (self.get_url (), url)

        return saga.filesystem.File (url=url, flags=flags, session=self.session, 
                                     _adaptor=self._adaptor, _adaptor_state=adaptor_state)


    # --------------------------------------------------------------------------
    #
    @SYNC_CALL
    def open_dir (self, url, flags) :

        self._is_valid ()

        adaptor_state = { "from_open" : True,
                          "cwd"       : saga.Url(self.url) }  # deep copy

        if  sumisc.url_is_relative (url) :
            url = sumisc.url_make_absolute (self.get_url (), url)

        return saga.filesystem.Directory (url=url, flags=flags, session=self.session, 
                                          _adaptor=self._adaptor, _adaptor_state=adaptor_state)

    # --------------------------------------------------------------------------
    #
    @SYNC_CALL
    def change_dir (self, tgt, flags) :

        cwdurl = saga.Url (self.url)
        tgturl = saga.Url (tgt)

        if  not sumisc.url_is_compatible (cwdurl, tgturl) :
            raise saga.BadParameter ("target dir outside of namespace '%s': %s" \
                                  % (cwdurl, tgturl))

        cmd = None

        if  tgturl.path == '.' or \
            tgturl.path == cwdurl.path :
            self._logger.debug ("change directory optimized away (%s) == (%s)" % (cwdurl, tgturl))


        if  flags & saga.filesystem.CREATE_PARENTS :
            cmd = " mkdir -p '%s' ;  cd '%s'" % (tgturl.path, tgturl.path)
        elif flags & saga.filesystem.CREATE :
            cmd = " mkdir    '%s' ;  cd '%s'" % (tgturl.path, tgturl.path)
        else :
            cmd = " test -d  '%s' && cd '%s'" % (tgturl.path, tgturl.path)

        ret, out, _ = self._command (cmd)

        if  ret != 0 :
            raise saga.BadParameter ("invalid dir '%s': %s" % (cwdurl, tgturl))

        self._logger.debug ("changed directory (%s)(%s)" % (ret, out))

        self.valid = True


    # --------------------------------------------------------------------------
    #
    @SYNC_CALL
    def close (self, timeout=None):

        if  timeout :
            raise saga.BadParameter ("timeout for close not supported")

        self.finalize (kill=True)


    # --------------------------------------------------------------------------
    @SYNC_CALL
    def get_url (self) :

        self._is_valid ()

        return saga.Url (self.url) # deep copy


    # --------------------------------------------------------------------------
    #
    @SYNC_CALL
    def list (self, npat, flags):

        self._is_valid ()

        # FIXME: eval flags

        if  None == npat :
            npat = "."
        else :
            npat = '-d %s' % npat

        ret, out, _ = self._command (" /bin/ls -C1 %s\n" % npat)

        if  ret != 0 :
            raise saga.NoSuccess ("failed to list(): (%s)(%s)" \
                               % (ret, out))

        lines = filter (None, out.split ("\n"))
        self._logger.debug (lines)

        self.entries = []
        for line in lines :
            # FIXME: convert to absolute URLs?
            self.entries.append (saga.Url (line.strip ()))

        return self.entries
   
   
    # --------------------------------------------------------------------------
    #
    @SYNC_CALL
    def copy_self (self, tgt, flags):

        self._is_valid ()

        # FIXME: eval flags

        return self.copy (self.url, tgt, flags)
   
   
    # --------------------------------------------------------------------------
    #
    @SYNC_CALL
    def copy (self, src_in, tgt_in, flags, _from_task=None):

        self._is_valid ()

        # FIXME: eval flags

        cwdurl = saga.Url (self.url) # deep copy
        src    = saga.Url (src_in)   # deep copy
        tgt    = saga.Url (tgt_in)   # deep copy

        if  sumisc.url_is_relative (src) : src = sumisc.url_make_absolute (cwdurl, src)
        if  sumisc.url_is_relative (tgt) : tgt = sumisc.url_make_absolute (cwdurl, tgt)
    
        rec_flag = ""
        if  flags & saga.filesystem.RECURSIVE : 
            rec_flag  += "-r "

        files_copied = list()

        # if cwd, src and tgt point to the same host, we just run a shell cp
        # command on that host
        if  sumisc.url_is_compatible (cwdurl, src) and \
            sumisc.url_is_compatible (cwdurl, tgt) :

            if  flags & saga.filesystem.CREATE_PARENTS : 
                self._create_parent (cwdurl, tgt)

            # print "shell cp"
            ret, out, err = self._command (" cp %s '%s' '%s'\n" % (rec_flag, src.path, tgt.path))
            if  ret != 0 :
                raise saga.NoSuccess ("copy (%s -> %s) failed (%s): (out: %s) (err: %s)" \
                                   % (src, tgt, ret, out, err))


        # src and tgt are on different hosts, we need to find out which of them
        # is local (stage_from vs. stage_to).
        else :
            # print "! shell cp"

            # for a remote target, we need to manually create the target dir (if
            # needed)
            if  flags & saga.filesystem.CREATE_PARENTS : 
                self._create_parent (cwdurl, tgt)


            # if cwd is remote, we use stage from/to
            if  not sumisc.url_is_local (cwdurl) :

                # print "cwd remote"
                if  sumisc.url_is_local (src)          and \
                    sumisc.url_is_compatible (cwdurl, tgt) :

                  # print "from local to remote: %s -> %s" % (src.path, tgt.path)
                    lease_tgt = self._adaptor.get_lease_target (self.url)
                    with self.lm.lease (lease_tgt, self.shell_creator, self.url) \
                        as copy_shell :
                        files_copied = copy_shell.stage_to_remote (src.path, tgt.path, rec_flag)

                elif sumisc.url_is_local (tgt)          and \
                     sumisc.url_is_compatible (cwdurl, src) :

                  # print "from remote to local: %s -> %s" % (src.path, tgt.path)
                    lease_tgt = self._adaptor.get_lease_target (self.url)
                    with self.lm.lease (lease_tgt, self.shell_creator, self.url) \
                        as copy_shell :
                        files_copied = copy_shell.stage_from_remote (src.path, tgt.path, rec_flag)

                else :
                    # print "from remote to other remote -- fail"
                    # we cannot support the combination of URLs
                    raise saga.BadParameter ("copy from %s to %s is not supported" \
                                          % (src, tgt))
   

            # if cwd is local, and src or tgt are remote, we need to actually
            # create a new pipe to the target host.  note that we may not have
            # a useful session for that!
            else : # sumisc.url_is_local (cwdurl) :

                # print "cwd local"

                if  sumisc.url_is_local (src) :

                    # need a compatible target scheme
                    if  tgt.scheme and not tgt.scheme.lower () in _ADAPTOR_SCHEMAS :
                        raise saga.BadParameter ("schema of copy target is not supported (%s)" \
                                              % (tgt))

                    # print "from local to remote"
                    lease_tgt = self._adaptor.get_lease_target (tgt)
                    with self.lm.lease (lease_tgt, self.shell_creator, tgt) \
                        as copy_shell :
                        files_copied = copy_shell.stage_to_remote (src.path, tgt.path, rec_flag)

                elif sumisc.url_is_local (tgt) :

                    # need a compatible source scheme
                    if  src.scheme and not src.scheme.lower () in _ADAPTOR_SCHEMAS :
                        raise saga.BadParameter ("schema of copy source is not supported (%s)" \
                                              % (src))

                    # print "from remote to local"
                    lease_tgt = self._adaptor.get_lease_target (tgt)
                    with self.lm.lease (lease_tgt, self.shell_creator, tgt) \
                        as copy_shell :
                        files_copied = copy_shell.stage_from_remote (src.path, tgt.path, rec_flag)

                else :

                    # print "from remote to other remote -- fail"
                    # we cannot support the combination of URLs
                    raise saga.BadParameter ("copy from %s to %s is not supported" \
                                          % (src, tgt))

   
        if  _from_task :
            _from_task._set_metric ('files_copied', files_copied)


    # --------------------------------------------------------------------------
    #
    @SYNC_CALL
    def link_self (self, tgt, flags):

        self._is_valid ()

        # FIXME: eval flags

        return self.link (self.url, tgt, flags)
   
   
    # --------------------------------------------------------------------------
    #
    @SYNC_CALL
    def link (self, src_in, tgt_in, flags, _from_task=None):

        # link will *only* work if src and tgt are on the same resource (and
        # even then may fail)
        self._is_valid ()

        cwdurl = saga.Url (self.url) # deep copy
        src    = saga.Url (src_in)   # deep copy
        tgt    = saga.Url (tgt_in)   # deep copy

        rec_flag = ""
        if  flags & saga.filesystem.RECURSIVE : 
            raise saga.BadParameter ("'RECURSIVE' flag not  supported for link()")

        if  flags & saga.filesystem.CREATE_PARENTS : 
            self._create_parent (cwdurl, tgt)

        # if src and tgt point to the same host, we just run a shell link
        # on that host
        if  sumisc.url_is_compatible (cwdurl, src) and \
            sumisc.url_is_compatible (cwdurl, tgt) :

            # print "shell ln"
            ret, out, err = self._command (" ln -s '%s' '%s'\n" % (src.path, tgt.path))
            if  ret != 0 :
                raise saga.NoSuccess ("link (%s -> %s) failed (%s): (out: %s) (err: %s)" \
                                   % (src, tgt, ret, out, err))


        # src and tgt are on different hosts, this is not supported
        else :
            raise saga.BadParameter ("link is only supported on same file system as cwd")



    # --------------------------------------------------------------------------
    #
    @SYNC_CALL
    def move_self (self, tgt, flags):

        self._is_valid ()

        # FIXME: eval flags

        self.move (self.url, tgt, flags)

        # need to re-initialize for new location
        self.url   = tgt
        self.flags = flags
        self.initialize ()
   
   
    # --------------------------------------------------------------------------
    #
    @SYNC_CALL
    def move (self, src_in, tgt_in, flags):

        # we handle move non-atomically, i.e. as copy/remove
        self.copy   (src_in, tgt_in, flags)
        self.remove (src_in, flags)
   
   
    # --------------------------------------------------------------------------
    #
    @SYNC_CALL
    def remove_self (self, flags):

        self._is_valid ()

        # FIXME: eval flags

        self.remove (self.url, flags)
        self.invalid = True
   
   
    # --------------------------------------------------------------------------
    #
    @SYNC_CALL
    def remove (self, tgt_in, flags):

        self._is_valid ()

        # FIXME: eval flags
        # FIXME: check if tgt remove happens to affect cwd... :-P

        cwdurl = saga.Url (self.url) # deep copy
        tgt    = saga.Url (tgt_in)   # deep copy

        rec_flag = ""
        if  flags & saga.filesystem.RECURSIVE : 
            rec_flag  += "-r "

        if  sumisc.url_is_compatible (cwdurl, tgt) :

            ret, out, err = self._command (" rm -f %s '%s'\n" % (rec_flag, tgt.path))
            if  ret != 0 :
                raise saga.NoSuccess ("remove (%s) failed (%s): (out: %s) (err: %s)" \
                                   % (tgt, ret, out, err))


        # we cannot support the URL
        else :
            raise saga.BadParameter ("remove of %s is not supported" % tgt)
   
   
    # --------------------------------------------------------------------------
    #
    @SYNC_CALL
    def make_dir (self, tgt_in, flags):

        self._is_valid ()

        cwdurl = saga.Url (self.url) # deep copy
        tgt    = saga.Url (tgt_in)   # deep copy

        if  flags & saga.filesystem.EXCLUSIVE : 
            # FIXME: this creates a race condition between testing for exclusive
            # mkdir and creating the dir.
            ret, out, _ = self._command (" test -d '%s' " % tgt.path)

            if  ret != 0 :
                raise saga.AlreadyExists ("make_dir target (%s) exists (%s)" \
                    % (tgt_in, out))


        options = ""

        if  flags & saga.filesystem.CREATE_PARENTS : 
            ret, out, _ = self._command (" mkdir -p '%s'" % tgt.path, make_location=True)
        else :
            ret, out, _ = self._command (" mkdir '%s'" % tgt.path)

        if  ret != 0 :
            raise saga.NoSuccess ("make_dir (%s) failed: %s" % (tgt_in, out))

   
    # --------------------------------------------------------------------------
    #
    @SYNC_CALL
    def get_size_self (self) :

        self._is_valid ()

        return self.get_size (self.url)

    # --------------------------------------------------------------------------
    #
    @SYNC_CALL
    def get_size (self, tgt_in) :

        self._is_valid ()

        cwdurl = saga.Url (self.url) # deep copy
        tgt    = saga.Url (tgt_in)   # deep copy

        ret, out, _ = self._command (" du -ks '%s'  | xargs | cut -f 1 -d ' '\n" % tgt.path)
        if  ret != 0 :
            raise saga.NoSuccess ("get size for (%s) failed (%s): (%s)" \
                               % (tgt, ret, out))

        size = None
        try :
            size = int (out) * 1024 # see '-k' option to 'du'
        except Exception as e :
            raise saga.NoSuccess ("could not get file size: %s" % out)

        return size
   

    # --------------------------------------------------------------------------
    #
    @SYNC_CALL
    def exists_self (self):

        self._is_valid ()

        return self.exists (self.url)
   
   
    # ----------------------------------------------------------------
    #
    @SYNC_CALL
    def exists (self, tgt_in):

        self._is_valid ()

        cwdurl = saga.Url (self.url) # deep copy
        tgt    = saga.Url (tgt_in)   # deep copy

        ret, out, _ = self._command (" test -e '%s'" % tgt.path)

        return True if ret == 0 else False
   
   
    # ----------------------------------------------------------------
    #
    @SYNC_CALL
    def is_dir_self (self):

        self._is_valid ()

        return self.is_dir (self.url)
   
   
    # --------------------------------------------------------------------------
    #
    @SYNC_CALL
    def is_dir (self, tgt_in):

        self._is_valid ()

        cwdurl = saga.Url (self.url) # deep copy
        tgt    = saga.Url (tgt_in)   # deep copy

        ret, out, _ = self._command (" test -d '%s' && test ! -h '%s'" % (tgt.path, tgt.path))

        return True if ret == 0 else False
   
   
    # --------------------------------------------------------------------------
    #
    @SYNC_CALL
    def is_entry_self (self):

        self._is_valid ()

        return self.is_entry (self.url)
   
   
    # --------------------------------------------------------------------------
    #
    @SYNC_CALL
    def is_entry (self, tgt_in):

        self._is_valid ()

        cwdurl = saga.Url (self.url) # deep copy
        tgt    = saga.Url (tgt_in)   # deep copy

        ret, out, _ = self._command (" test -f '%s' && test ! -h '%s'" % (tgt.path, tgt.path))

        return True if ret == 0 else False
   
   
    # --------------------------------------------------------------------------
    #
    @SYNC_CALL
    def is_link_self (self):

        self._is_valid ()

        return self.is_link (self.url)
   
   
    # --------------------------------------------------------------------------
    #
    @SYNC_CALL
    def is_link (self, tgt_in):

        self._is_valid ()

        cwdurl = saga.Url (self.url) # deep copy
        tgt    = saga.Url (tgt_in)   # deep copy

        ret, out, _ = self._command (" test -h '%s'" % tgt.path)

        return True if ret == 0 else False
   
   
    # --------------------------------------------------------------------------
    #
    @SYNC_CALL
    def is_file_self (self):

        return self.is_entry_self ()
   
   
    # --------------------------------------------------------------------------
    #
    @SYNC_CALL
    def is_file (self, tgt_in):

        return self.is_entry (tgt_in)
   
   
###############################################################################
#
class ShellFile (saga.adaptors.cpi.filesystem.File) :
    """ Implements saga.adaptors.cpi.filesystem.File
    """
    # --------------------------------------------------------------------------
    #
    def __init__ (self, api, adaptor) :

        _cpi_base = super  (ShellFile, self)
        _cpi_base.__init__ (api, adaptor)


    # --------------------------------------------------------------------------
    #
    def __del__ (self) :

        self.finalize (kill=True)


    # --------------------------------------------------------------------------
    #
    def _is_valid (self) :

        if  not self.valid :
            raise saga.IncorrectState ("this instance was closed or removed")

    # --------------------------------------------------------------------------
    #
    def _create_parent (self, cwdurl, tgt) :

        dirname = sumisc.url_get_dirname (tgt)

        if  sumisc.url_is_compatible (cwdurl, tgt) :

            ret, out, _ = self._run_sync (" mkdir -p '%s'\n" % (dirname))
            if  ret != 0 :
                raise saga.NoSuccess ("failed at mkdir '%s': (%s) (%s)" \
                                   % (dirname, ret, out))

        elif sumisc.url_is_local (tgt) :

            if  tgt.scheme and not tgt.scheme.lower () in _ADAPTOR_SCHEMAS :
                raise saga.BadParameter ("schema of mkdir target is not supported (%s)" \
                                      % (tgt))
            _mkdir_p(dirname)

        else :

            lease_tgt = self._adaptor.get_lease_target (tgt)
            with self.lm.lease (lease_tgt, self.shell_creator, tgt) as tmp_shell :
                tmp_shell.run_sync ('mkdir -p %s' % dirname)


    # --------------------------------------------------------------------------
    #
    def _run_sync(self, cmd):
        """
        lease the shell, set $PWD, run command, and release shell again
        """

        lease_tgt = self._adaptor.get_lease_target(self.cwdurl)
        cwd_path  = self.cwdurl.path

        with self.lm.lease(lease_tgt, self.shell_creator, self.cwdurl) as shell:

            return shell.run_sync("cd %s && %s\n" % (cwd_path, cmd))


    # --------------------------------------------------------------------------
    #
    @SYNC_CALL
    def init_instance (self, adaptor_state, url, flags, session):

        # FIXME: eval flags!
        if  flags == None :
            flags = 0

        self._logger.info ("init_instance %s" % url)

        if  'from_open' in adaptor_state and adaptor_state['from_open'] :

            # comes from job.service.create_job()
            self.url         = saga.Url(url) # deep copy
            self.flags       = flags
            self.session     = session
            self.valid       = False  # will be set by initialize
            self.lm          = session._lease_manager

            self.cwdurl      = saga.Url (adaptor_state["cwd"])
            self.cwd         = self.cwdurl.path

            if  sumisc.url_is_relative (self.url) :
                self.url = sumisc.url_make_absolute (self.cwd, self.url)

        else :

            if  sumisc.url_is_relative (url) :
                raise saga.BadParameter ("cannot interprete relative URL in this context ('%s')" % url)

            self.url         = url
            self.flags       = flags
            self.session     = session
            self.valid       = False  # will be set by initialize
            self.lm          = session._lease_manager

            self.cwd         = sumisc.url_get_dirname (url)
            self.cwdurl      = saga.Url (url) # deep copy
            self.cwdurl.path = self.cwd

            if  not self.flags :
                self.flags = 0

        # Use `_set_session` method of the base class to set the session object
        # `_set_session` and `get_session` methods are provided by `CPIBase`.
        self._set_session(session)

        def _shell_creator (url) :
            return sups.PTYShell (url, self.get_session(), self._logger)
        self.shell_creator = _shell_creator

     #  # self.shell is also a leased shell -- for File, it does not have any
     #  # state, really.
     #  # FIXME: get ssh Master connection from _adaptor dict
     #  lease_tgt  = self._adaptor.get_lease_target (self.url)
     #  self.shell = self.lm.lease (lease_tgt, self.shell_creator, self.url) 
     #  # TODO : release shell
     #
     ## self.shell.set_initialize_hook (self.initialize)
     ## self.shell.set_finalize_hook   (self.finalize)

        self.initialize ()

        return self.get_api ()


    # --------------------------------------------------------------------------
    #
    def initialize (self) :

        # shell got started, found its prompt.  Now, change
        # to the initial (or later current) working directory.

        cmd = ""
        dirname = sumisc.url_get_dirname  (self.url)

        if  self.flags & saga.filesystem.CREATE_PARENTS :
            cmd = " mkdir -p '%s'; touch '%s'" % (dirname, self.url.path)
            self._logger.info ("mkdir '%s'; touch '%s'" % (dirname, self.url.path))

        elif self.flags & saga.filesystem.CREATE :
            cmd = " touch '%s'" % (self.url.path)
            self._logger.info ("touch %s" % self.url.path)

        else :
            cmd = " true"


        if  self.flags & saga.filesystem.READ :
            cmd += "; test -r '%s'" % (self.url.path)

        if  self.flags & saga.filesystem.WRITE :
            cmd += "; test -w '%s'" % (self.url.path)

        ret, out, _ = self._run_sync(cmd)

        if  ret != 0 :
            if  self.flags & saga.filesystem.CREATE_PARENTS :
                raise saga.BadParameter ("cannot open/create: '%s' - %s" % (self.url.path, out))
            elif self.flags & saga.filesystem.CREATE :
                raise saga.BadParameter ("cannot open/create: '%s' - %s" % (self.url.path, out))
            else :
                raise saga.DoesNotExist("File does not exist: '%s' - %s" % (self.url.path, out))

        self._logger.info ("file initialized (%s)(%s)" % (ret, out))

        self.valid = True


    # --------------------------------------------------------------------------
    #
    def finalize (self, kill=False) :

        self.valid = False


    # --------------------------------------------------------------------------
    #
    @SYNC_CALL
    def close (self, timeout=None):

        if  timeout :
            raise saga.BadParameter ("timeout for close not supported")

        self.finalize (kill=True)


    # --------------------------------------------------------------------------
    #
    @SYNC_CALL
    def get_url (self):

        self._is_valid ()

        return saga.Url (self.url) # deep copy


    # --------------------------------------------------------------------------
    #
    @SYNC_CALL
    def copy_self (self, tgt_in, flags):

        self._is_valid ()

        # FIXME: eval flags

        # print "copy self %s -> %s" % (self.url, tgt_in)

        cwdurl = saga.Url (self.cwdurl) # deep copy
        src    = saga.Url (self.url)    # deep copy
        tgt    = saga.Url (tgt_in)      # deep copy

        if sumisc.url_is_relative (src) : src = sumisc.url_make_absolute (cwdurl, src)
        if sumisc.url_is_relative (tgt) : tgt = sumisc.url_make_absolute (cwdurl, tgt)

        rec_flag = ""
        if  flags & saga.filesystem.RECURSIVE :
            rec_flag  += "-r "

        if  flags & saga.filesystem.CREATE_PARENTS :
            self._create_parent (cwdurl, tgt)

        # if cwd, src and tgt point to the same host, we just run a shell cp
        # command on that host
        if  sumisc.url_is_compatible (cwdurl, src) and \
            sumisc.url_is_compatible (cwdurl, tgt) :

            # print "shell cp"
            ret, out, _ = self._run_sync (" cp %s '%s' '%s'\n" % (rec_flag, src.path, tgt.path))
            if  ret != 0 :
                raise saga.NoSuccess ("copy (%s -> %s) failed (%s): (%s)" \
                                   % (src, tgt, ret, out))


        # src and tgt are on different hosts, we need to find out which of them
        # is local (stage_from vs. stage_to).
        else :
            # print "! shell cp"

            # if cwd is remote, we use stage from/to on the existing pipe
            if  not sumisc.url_is_local (cwdurl) :

                # print "cwd remote"

                if  sumisc.url_is_local (src)          and \
                    sumisc.url_is_compatible (cwdurl, tgt) :

                    # print "from local to remote"
                    lease_tgt = self._adaptor.get_lease_target (self.url)
                    with self.lm.lease (lease_tgt, self.shell_creator, self.url) \
                        as copy_shell :
                        files_copied = copy_shell.stage_to_remote (src.path, tgt.path, rec_flag)

                elif sumisc.url_is_local (tgt)          and \
                     sumisc.url_is_compatible (cwdurl, src) :

                    # print "from remote to local"
                    lease_tgt = self._adaptor.get_lease_target (self.url)
                    with self.lm.lease (lease_tgt, self.shell_creator, self.url) \
                        as copy_shell :
                        files_copied = copy_shell.stage_from_remote (src.path, tgt.path, rec_flag)

                else :
                    # print "from remote to other remote -- fail"
                    # we cannot support the combination of URLs
                    raise saga.BadParameter ("copy from %s to %s is not supported" \
                                          % (src, tgt))


            # if cwd is local, and src or tgt are remote, we need to actually
            # create a new pipe to the target host.  note that we may not have
            # a useful session for that!
            else : # sumisc.url_is_local (cwdurl) :

                # print "cwd local"

                if  sumisc.url_is_local (src) :

                    # need a compatible target scheme
                    if  tgt.scheme and not tgt.scheme.lower () in _ADAPTOR_SCHEMAS :
                        raise saga.BadParameter ("schema of copy target is not supported (%s)" \
                                              % (tgt))

                    # print "from local to remote"
                    lease_tgt = self._adaptor.get_lease_target (tgt)
                    with self.lm.lease (lease_tgt, self.shell_creator, tgt) \
                        as copy_shell :
                        files_copied = copy_shell.stage_to_remote (src.path, tgt.path, rec_flag)

                elif sumisc.url_is_local (tgt) :

                    # need a compatible source scheme
                    if  src.scheme and not src.scheme.lower () in _ADAPTOR_SCHEMAS :
                        raise saga.BadParameter ("schema of copy source is not supported (%s)" \
                                              % (src))

                    # print "from remote to local"
                    lease_tgt = self._adaptor.get_lease_target (tgt)
                    with self.lm.lease (lease_tgt, self.shell_creator, tgt) \
                        as copy_shell :
                        files_copied = copy_shell.stage_from_remote (src.path, tgt.path, rec_flag)

                else :

                    # we cannot support two remote URLs
                    raise saga.BadParameter ("copy from %s to %s is not supported" \
                                          % (src, tgt))


    # --------------------------------------------------------------------------
    #
    @SYNC_CALL
    def link_self (self, tgt_in, flags, _from_task=None):

        # link will *only* work if src and tgt are on the same resource (and
        # even then may fail)
        self._is_valid ()

        cwdurl = saga.Url (self.url) # deep copy
        src    = saga.Url (self.url)    # deep copy
        tgt    = saga.Url (tgt_in)   # deep copy

        rec_flag = ""
        if  flags & saga.filesystem.RECURSIVE :
            raise saga.BadParameter ("'RECURSIVE' flag not  supported for link()")

        if  flags & saga.filesystem.CREATE_PARENTS :
            self._create_parent (cwdurl, tgt)

        # if src and tgt point to the same host, we just run a shell link
        # on that host
        if  sumisc.url_is_compatible (cwdurl, src) and \
            sumisc.url_is_compatible (cwdurl, tgt) :

            # print "shell ln"
            ret, out, err = self._run_sync (" ln -s '%s' '%s'\n" % (src.path, tgt.path))
            if  ret != 0 :
                raise saga.NoSuccess ("link (%s -> %s) failed (%s): (out: %s) (err: %s)" \
                                   % (src, tgt, ret, out, err))


        # src and tgt are on different hosts, this is not supported
        else :
            raise saga.BadParameter ("link is only supported on same file system as cwd")



    # --------------------------------------------------------------------------
    #
    @SYNC_CALL
    def move_self (self, tgt_in, flags):

        # we handle move non-atomically, i.e. as copy/remove
        self.copy_self   (tgt_in, flags)
        self.remove_self (flags)

        # however, we are not closed at this point, but need to re-initialize
        self.url   = tgt_in
        self.flags = flags
        self.initialize ()

    # --------------------------------------------------------------------------
    #
    @SYNC_CALL
    def write (self, string, flags=None):
        """
        This call is intended to write a string to a local or remote file.
        Since write() uses file staging calls, it cannot be used to randomly
        write certain parts of a file (i.e. seek()). Together with read(),
        it was designed to manipulate template files and write them back to
        the remote directory. Be aware, that writing large files will
        be very slow compared to native read(2) and write(2) calls.
        """
        self._is_valid ()
        if  flags==None:
            flags = self.flags
        else:
            self.flags=flags

        tgt = saga.Url (self.url)  # deep copy, is absolute

        if  flags==saga.filesystem.APPEND:
            string = self.read()+string
        # FIXME: eval flags

        lease_tgt = self._adaptor.get_lease_target(self.cwdurl)

        with self.lm.lease(lease_tgt, self.shell_creator, self.cwdurl) as shell:
            shell.write_to_remote(string,tgt.path)

    # --------------------------------------------------------------------------
    #
    @SYNC_CALL
    def read (self,size=None):
        """
        This call is intended to read a string wit length size from a local
        or remote file.	Since read() uses file staging calls, it cannot be
        used to randomly read certain parts of a file (i.e. seek()).
        Together with write(), it was designed to manipulate template files
        and write them back to the remote directory. Be aware, that reading
        large files will be very slow compared to native read(2) and write(2)
        calls.
        """

        self._is_valid ()

        lease_tgt = self._adaptor.get_lease_target(self.cwdurl)
        tgt       = saga.Url(self.url)  # deep copy, is absolute

        with self.lm.lease(lease_tgt, self.shell_creator, self.cwdurl) as shell:
            out = shell.read_from_remote(tgt.path)

        if size != None:
            return out[0:size-1]
        else:
            return out



    # --------------------------------------------------------------------------
    #
    @SYNC_CALL
    def remove_self (self, flags):

        self._is_valid ()

        # FIXME: eval flags
        # FIXME: check if tgt remove happens to affect cwd... :-P

        tgt = saga.Url (self.url)  # deep copy, is absolute

        rec_flag = ""
        if  flags & saga.filesystem.RECURSIVE :
            rec_flag  += "-r "

        ret, out, _ = self._run_sync (" rm -f %s '%s'\n" % (rec_flag, tgt.path))
        if  ret != 0 :
            raise saga.NoSuccess ("remove (%s) failed (%s): (%s)" \
                               % (tgt, ret, out))


    # --------------------------------------------------------------------------
    #
    @SYNC_CALL
    def get_size_self (self) :

        self._is_valid ()
        size      = None
        size_mult = 1
        ret       = None
        out       = None

        if  self.is_dir_self () :
            size_mult   = 1024   # see '-k' option to 'du'
            ret, out, _ = self._run_sync (" du -ks '%s'  | xargs | cut -f 1 -d ' '\n" \
                                            % self.url.path)
        else :
            ret, out, _ = self._run_sync (" wc -c '%s' | xargs | cut -f 1 -d ' '\n" \
                                            % self.url.path)

        if  ret != 0 :
            raise saga.NoSuccess ("get size for (%s) failed (%s): (%s)" \
                               % (self.url, ret, out))

        try :
            size = int (out) * size_mult
        except Exception as e :
            raise saga.NoSuccess ("could not get file size: %s" % out)


        return size


    # --------------------------------------------------------------------------
    #
    @SYNC_CALL
    def is_dir_self (self):

        self._is_valid ()

        cwdurl = saga.Url (self.url) # deep copy

        ret, out, _ = self._run_sync (" test -d '%s' && test ! -h '%s'" % (cwdurl.path, cwdurl.path))

        return True if ret == 0 else False


    # --------------------------------------------------------------------------
    #
    @SYNC_CALL
    def is_entry_self (self):

        self._is_valid ()

        cwdurl = saga.Url (self.url) # deep copy

        ret, out, _ = self._run_sync (" test -f '%s' && test ! -h '%s'" % (cwdurl.path, cwdurl.path))

        return True if ret == 0 else False


    # --------------------------------------------------------------------------
    #
    @SYNC_CALL
    def is_link_self (self):

        self._is_valid ()

        cwdurl = saga.Url (self.url) # deep copy

        ret, out, _ = self._run_sync (" test -h '%s'" % cwdurl.path)

        return True if ret == 0 else False


    # --------------------------------------------------------------------------
    #
    @SYNC_CALL
    def is_file_self (self):

        self._is_valid ()

        cwdurl = saga.Url (self.url) # deep copy

        ret, out, _ = self._run_sync (" test -f '%s' && test ! -h '%s'" % (cwdurl.path, cwdurl.path))

        return True if ret == 0 else False


# ------------------------------------------------------------------------------

