
__author__    = "Andre Merzky, Ole Weidner"
__copyright__ = "Copyright 2012-2013, The SAGA Project"
__license__   = "MIT"


import os

import saga.context
import saga.exceptions as se
import saga.adaptors.base
import saga.adaptors.cpi.context

SYNC_CALL  = saga.adaptors.cpi.decorators.SYNC_CALL
ASYNC_CALL = saga.adaptors.cpi.decorators.ASYNC_CALL

######################################################################
#
# adaptor meta data
#
_ADAPTOR_NAME          = 'saga.adaptor.ssh'
_ADAPTOR_SCHEMAS       = ['ssh']
_ADAPTOR_OPTIONS       = []

# FIXME: complete attribute list
_ADAPTOR_CAPABILITIES  = {
     'ctx_attributes'   : {saga.context.TYPE      : "This MUST be set to ssh",
                           saga.context.USER_ID   : "user name on target machine",
                           saga.context.USER_KEY  : "maps to public ssh key",
                           saga.context.USER_CERT : "maps to private ssh key",
                           saga.context.USER_PASS : "passphrase for encrypted keys"}
}

_ADAPTOR_DOC           = {
    'name'             : _ADAPTOR_NAME,
    'cfg_options'      : _ADAPTOR_OPTIONS, 
    'capabilities'     : _ADAPTOR_CAPABILITIES,
    'description'      : """ 
    
This SSH :class:`saga.Context` adaptor points to an ssh keypair and a user_id
to be used for ssh based backend connections.  For example, an ssh context can
be used to start jobs (:class:`saga.job.Job`) via ssh, to copy files
(:class:`saga.filesystem.File`) via sftp, etc.

Not all supported attributes have to be defined when using an ssh context
adaptor -- unspecified attributes will have sensible default values.  For
example, the ``c.user_id`` will default to the local user id, and the default
passphrase in ``c.user_pass`` will be empty.

The `UserKey` and `UserCert` attributes can point to either the public or
private key of the ssh keypair -- the SAGA-Python implementation will internally
complete the respective other key (public key file names are expected to be
derived from the private key, by appending the suffix `.pub` -- `.pem` files are
expected to contain both public and private key.).
    """,

    'schemas'          : {'ssh' : 'ssh key and userid information.'},
    'example'          : "examples/context/context_ssh.py"
}

_ADAPTOR_INFO          = {
    'name'             : _ADAPTOR_NAME,
    'version'          : 'v0.1',
    'schemas'          : _ADAPTOR_SCHEMAS,
    'cpis'             : [{ 
        'type'         : 'saga.Context',
        'class'        : 'ContextSSH'
        }
    ]
}


# ------------------------------------------------------------------------------
# 
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

        # there are no default myproxy contexts
        self._default_contexts = []
        self._have_defaults    = False


    # --------------------------------------------------------------------------
    #
    def sanity_check (self) :
        pass


    # --------------------------------------------------------------------------
    #
    def _get_default_contexts (self) :

        if not self._have_defaults :
 
            import glob
            candidate_certs = glob.glob ("%s/.ssh/*" % os.environ['HOME'])

            for key in candidate_certs :

                if  os.path.isdir (key) :
                    # don't want directories -- only keys
                    continue

                elif  key.endswith ('.pub') :
                    # don't want public keys in this loop
                    continue

                elif  key.endswith ('.pem') :
                    # we can't handle passwd protected PEMs well, since we can't
                    # detect when they are passwd protected in the first place,
                    # so we do not pick them up by default
                    continue
                else :
                    pub  = "%s.pub" % key

                # the private and public keys must exist
                if  not os.path.exists (key ) or \
                    not os.path.isfile (key )    :
                    self._logger.info ("ignore ssh key at %s (no private key: %s)" %  (key, key))
                    continue

                if  not os.path.exists (pub) or \
                    not os.path.isfile (pub)    :
                    self._logger.info ("ignore ssh key at %s (no public key: %s)" %  (key, pub))
                    continue


                try :
                    fh_key  = open (key )
                except Exception as e:
                    self._logger.info ("ignore ssh key at %s (key not readable: %s)" %  (key, e))
                    continue
                else :
                    fh_key .close ()


                try :
                    fh_pub  = open (pub )
                except Exception as e:
                    self._logger.info ("ignore ssh key at %s (public key %s not readable: %s)" %  (key, pub, e))
                    continue
                else :
                    fh_pub .close ()


                import subprocess
                if  not subprocess.call (["sh", "-c", "grep ENCRYPTED %s > /dev/null" % key]) :
                    # needs passphrase.  Great, actually, but won't work for
                    # default contexts as long as we can't prompt for pass
                    # phrases...
                    self._logger.warn ("ignore  ssh key at %s (requires passphrase)" %  key)
                    continue

                c = saga.Context ('ssh')
                c.user_key  = key
                c.user_cert = pub

                self._default_contexts.append (c)

                self._logger.info ("default ssh key at %s" %  key)

            self._have_defaults = True


        # have defaults, and can return them...
        return self._default_contexts



# ------------------------------------------------------------------------------
#
class ContextSSH (saga.adaptors.cpi.context.Context) :

    # --------------------------------------------------------------------------
    #
    def __init__ (self, api, adaptor) :

        _cpi_base = super  (ContextSSH, self)
        _cpi_base.__init__ (api, adaptor)


    # --------------------------------------------------------------------------
    #
    @SYNC_CALL
    def init_instance (self, adaptor_state, type) :

        if not type.lower () in (schema.lower() for schema in _ADAPTOR_SCHEMAS) :
            raise se.BadParameter \
                    ("the ssh context adaptor only handles ssh contexts - duh!")

        self.get_api ().type = type

        return self.get_api ()


    # --------------------------------------------------------------------------
    #
    @SYNC_CALL
    def _initialize (self, session) :

        # make sure we have can access the key
        api = self.get_api ()

        unexpanded_key = None
        unexpanded_pub = None
        pwd = None

        
        if api.attribute_exists (saga.context.USER_KEY ) :
            unexpanded_key  = api.get_attribute    (saga.context.USER_KEY )
        if api.attribute_exists (saga.context.USER_CERT) :
            unexpanded_pub  = api.get_attribute    (saga.context.USER_CERT)
        if api.attribute_exists (saga.context.USER_PASS) :
            pwd  = api.get_attribute    (saga.context.USER_PASS)

        # Expand any environment variables in the key/pub paths
        if unexpanded_key:
            key = os.path.expandvars(unexpanded_key)
        else:
            key = None
        if unexpanded_pub:
            pub = os.path.expandvars(unexpanded_pub)
        else:
            pub = None

        # either user_key or user_cert should be specified (or both), 
        # then we complement the other, and convert to/from private 
        # from/to public keys
        if  pub  and not key :
            key  = pub

        if  not key :
            # nothing to do, really.  This likely means that ssh setup is
            # done out-of-band.
            return

        # convert public key into private key
        if  key.endswith ('.pub') :
            if  not pub :
                pub = key
            key = key[:-4]
        elif key.endswith ('.pem') :
            if  not pub :
                pub = key
        else :
            if  not pub :
                pub = key+'.pub'

        # update the context with these setting
        api.set_attribute (saga.context.USER_KEY , key)
        api.set_attribute (saga.context.USER_CERT, pub)


        # the private and public keys must exist
        if  not os.path.exists (key) or \
            not os.path.isfile (key)    :
            raise se.BadParameter ("ssh key inaccessible: %s" % (key))

        if  not os.path.exists (pub) or \
            not os.path.isfile (pub)    :
            raise se.BadParameter ("ssh public key inaccessible: %s" % (pub))


        try :
            fh_key = open (key)
        except Exception as e:
            raise se.PermissionDenied ("ssh key '%s' not readable: %s" % (key, e))
        else :
            fh_key.close ()


        try :
            fh_pub = open (pub)
        except Exception as e:
            raise se.PermissionDenied ("ssh public key '%s' not readable: %s" % (pub, e))
        else :
            fh_pub.close ()


        import subprocess
        if  not subprocess.call (["sh", "-c", "grep ENCRYPTED %s > /dev/null" % key]) :
            if  pwd  :
                if  subprocess.call (["sh", "-c", "ssh-keygen -y -f %s -P '%s' > /dev/null" % (key, pwd)]) :
                    raise se.PermissionDenied ("ssh key '%s' is encrypted, incorrect password" % (key))
            else :
                self._logger.error ("ssh key '%s' is encrypted, unknown password" % (key))


        self._logger.info ("init SSH context for key  at '%s' done" % key)




