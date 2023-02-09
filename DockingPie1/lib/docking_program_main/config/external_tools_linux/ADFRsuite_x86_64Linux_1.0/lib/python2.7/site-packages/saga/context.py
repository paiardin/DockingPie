
__author__    = "Andre Merzky, Ole Weidner"
__copyright__ = "Copyright 2012-2013, The SAGA Project"
__license__   = "MIT"

import copy

import radical.utils.signatures as rus

import saga.adaptors.base    as sab
import saga.attributes       as sa
import saga.base             as sb

from   saga.constants import TYPE,       SERVER,    USER_CERT,   CERT_REPOSITORY
from   saga.constants import USER_PROXY, USER_KEY,  USER_ID,     USER_PASS,   USER_VO
from   saga.constants import LIFE_TIME,  REMOTE_ID, REMOTE_HOST, REMOTE_PORT, TOKEN

# ------------------------------------------------------------------------------
#
class Context (sb.Base, sa.Attributes) :
    '''A SAGA Context object as defined in GFD.90.

    A security context is a description of a security token.  It is important to
    understand that, in general, a context really just *describes* a token, but
    that a context *is not* a token (*). For example, a context may point to
    a X509 certificate -- but it will in general not hold the certificate
    contents.

    Context classes are used to inform the backends used by SAGA on what
    security tokens are expected to be used.  By default, SAGA will be able to
    pick up such tokens from their default location, but in some cases it might
    be necessary to explicitly point to them - then use Session with
    context instances to do so.

    The usage example for contexts is below::

        # define an ssh context
        ctx = saga.Context("SSH")
        ctx.user_cert = '$HOME/.ssh/special_id_rsa'
        ctx.user_key  = '$HOME/.ssh/special_id_rsa.pub'

        # add the context to a session
        session = saga.Session()
        session.add_context(ctx)

        # create a job service in this session -- that job service can now
        # *only* use that ssh context. 
        j = saga.job.Service('ssh://remote.host.net/', session=session)


    The Session argument to the job.Service constructor is fully optional
    -- if left out, SAGA will use default session, which picks up some default
    contexts as described above -- that will suffice for the majority of use
    cases.

    ----

    (*) The only exception to this rule is the 'UserPass' key, which is used to
    hold plain-text passwords.  Use this key with care -- it is not good
    practice to hard-code passwords in the code base, or in config files.
    Also, be aware that the password may show up in log files, when debugging or
    analyzing your application.

    '''

    # --------------------------------------------------------------------------
    #
    @rus.takes   ('Context', 
                  basestring, 
                  rus.optional (sab.Base),
                  rus.optional (dict))
    @rus.returns (rus.nothing)
    def __init__ (self, ctype, _adaptor=None, _adaptor_state={}) : 
        '''
        ctype: string
        ret:   None
        '''

        sb.Base.__init__ (self, ctype.lower(), _adaptor, _adaptor_state, ctype, ttype=None)


        import saga.attributes as sa

        # set attribute interface propertiesP
        self._attributes_extensible  (False)
        self._attributes_camelcasing (True)

        # register properties with the attribute interface
        self._attributes_register  (TYPE,            None, sa.STRING, sa.SCALAR, sa.WRITEABLE)
        self._attributes_register  (SERVER,          None, sa.STRING, sa.SCALAR, sa.WRITEABLE)
        self._attributes_register  (TOKEN,           None, sa.STRING, sa.SCALAR, sa.WRITEABLE)
        self._attributes_register  (CERT_REPOSITORY, None, sa.STRING, sa.SCALAR, sa.WRITEABLE)
        self._attributes_register  (USER_PROXY,      None, sa.STRING, sa.SCALAR, sa.WRITEABLE)
        self._attributes_register  (USER_CERT,       None, sa.STRING, sa.SCALAR, sa.WRITEABLE)
        self._attributes_register  (USER_KEY,        None, sa.STRING, sa.SCALAR, sa.WRITEABLE)
        self._attributes_register  (USER_ID,         None, sa.STRING, sa.SCALAR, sa.WRITEABLE)
        self._attributes_register  (USER_PASS,       None, sa.STRING, sa.SCALAR, sa.WRITEABLE)
        self._attributes_register  (USER_VO,         None, sa.STRING, sa.SCALAR, sa.WRITEABLE)
        self._attributes_register  (LIFE_TIME,       -1,   sa.INT,    sa.SCALAR, sa.WRITEABLE)
        self._attributes_register  (REMOTE_ID,       None, sa.STRING, sa.SCALAR, sa.WRITEABLE)
        self._attributes_register  (REMOTE_HOST,     None, sa.STRING, sa.SCALAR, sa.WRITEABLE)
        self._attributes_register  (REMOTE_PORT,     None, sa.STRING, sa.VECTOR, sa.WRITEABLE)

        self.type = ctype


    # --------------------------------------------------------------------------
    #
    @rus.takes   ('Context')
    @rus.returns (basestring)
    def __str__  (self) :

        d = self.as_dict ()
        s = "{"

        for key in sorted (d.keys ()) :
            if  key == 'UserPass' and d[key] :
                s += "'UserPass' : '%s'" % ('x'*len(d[key]))
            else :
                s += "'%s' : '%s'" % (key, d[key])
            s += ', '

        return "%s}" % s[0:-2]


    # --------------------------------------------------------------------------
    #
    @rus.takes   ('Context')
    @rus.returns (basestring)
    def __repr__ (self) :

        return str(self)


    # --------------------------------------------------------------------------
    #
    @rus.takes      ('Context', 
                     ('Session', '_DefaultSession'))
    @rus.returns    (rus.nothing)
    def _initialize (self, session) :
        '''
        ret:  None
        '''
        self._adaptor._initialize (session)


    # --------------------------------------------------------------------------
    #
    def __deepcopy__(self, memo):

        print self
        
        ret = Context(self.type)

        for a in self.list_attributes():
            ret.set_attribute(a, self.get_attribute(a))

        return ret


# ------------------------------------------------------------------------------

