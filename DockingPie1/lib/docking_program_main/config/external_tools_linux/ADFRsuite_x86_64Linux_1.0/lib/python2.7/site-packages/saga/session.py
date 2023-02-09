
__author__    = "Andre Merzky, Ole Weidner"
__copyright__ = "Copyright 2012-2013, The SAGA Project"
__license__   = "MIT"

import copy

import radical.utils            as ru
import radical.utils.signatures as rus

import saga.exceptions          as se

import saga.engine.engine
import saga.context
import saga.base



# ------------------------------------------------------------------------------
#
class _ContextList (list) :
    """
    The `saga.Session` class has a 'contexts' member which is a mutable list of
    `saga.Context` instances.  Adding contexts to that list is semantically
    equivalent to calling `saga.Session.add_context (ctx)`, which (a) creates
    a deep copy of the context, and (b) initializes the context (i.e. calls
    `ctx._initialize (self)`, with `self` being the session instance).  We thus
    create our own provate `_ContextList` class which inherits from the native
    Python `list` class, and overload the `append()` call with said semantics.
    """

    # --------------------------------------------------------------------------
    #
    def __init__ (self, session=None, *args, **kwargs) :

        self._session = session

        if  session : 
            self._logger  = session._logger
        else :
            self._logger  = ru.get_logger ('radical.saga')

        base_list = super  (_ContextList, self)
        base_list.__init__ (*args, **kwargs)


    # --------------------------------------------------------------------------
    #
    def append (self, ctx, session=None) :

        ctx_clone = self._initialise_context(ctx, session)

        # context initialized ok, add it to the list of known contexts
        super (_ContextList, self).append (ctx_clone)
    

    # --------------------------------------------------------------------------
    #
    def insert(self, index, ctx, session=None) :

        ctx_clone = self._initialise_context(ctx, session)

        # context initialized ok, add it to the list of known contexts
        super (_ContextList, self).insert (0, ctx_clone)
    

    # --------------------------------------------------------------------------
    # Initialise a context to be added to the list of known contexts
    # Returns a cloned, initialised context that can be added to the context
    # list. 
    def _initialise_context(self, ctx, session=None):

        if  not isinstance (ctx, saga.Context) :
            raise TypeError, "item to add is not a saga.Context instance"

        # create a deep copy of the context (this keeps _adaptor etc)
        ctx_clone = saga.Context  (ctx.type)
        ctx._attributes_deep_copy (ctx_clone)

        if  not session :
            session = self._session
            logger  = self._logger
        else :
            logger  = session._logger


        # try to initialize that context, i.e. evaluate its attributes and
        # infer additional runtime information as needed
      # logger.debug ("adding  context : %s" % (ctx_clone))

        if  not session :
            logger.warning ("cannot initialize context - no session: %s" \
                       % (ctx_clone))
        else :
            try :
                ctx_clone._initialize (session)
            except se.SagaException as e :
                msg = "Cannot add context, initialization failed (%s)"  %  str(e)
                raise se.BadParameter (msg)

        return ctx_clone


    # --------------------------------------------------------------------------
    #
    def __deepcopy__(self, memo):

        ret = _ContextList()

        for c in self:
            ret.append(c)

        return ret


# ------------------------------------------------------------------------------
#
class Session (saga.base.SimpleBase) :
    """A SAGA Session object as defined in GFD.90.

    A SAGA session has the purpose of scoping the use of security credentials
    for remote operations.  In other words, a session instance acts as
    a container for security Context instances -- SAGA objects (such as
    job.Service or filesystem.File) created in that session will then use
    exactly the security contexts from that session (and no others).

    That way, the session serves two purposes:  (1) it helps SAGA to decide
    which security mechanism should be used for what interaction, and (2) it
    helps SAGA to find security credentials which would be difficult to pick up
    automatically.
    
    The use of a session is as follows:


    Example::


        # define an ssh context
        c = saga.Context('ssh')
        c.user_cert = '$HOME/.ssh/special_id_rsa.pub'
        c.user_key  = '$HOME/.ssh/special_id_rsa'

        # add it to a session
        s = saga.Session
        s.add_context(c)

        # create a job service in this session -- that job service can now
        # *only* use that ssh context. 
        j = saga.job.Service('ssh://remote.host.net/', s)


    The session argument to the L{job.Service} constructor is fully optional --
    if left out, SAGA will use default session, which picks up some default
    contexts as described above -- that will suffice for the majority of use
    cases.

    A session instance exposes a `context` property, which is a list of
    authentication contexts managed by this session.  As the contexts and the
    session are stateless, it is safe to modify this list as needed.  
    """

    # FIXME: session deep copy not implemented

    # --------------------------------------------------------------------------
    #
    @rus.takes   ('Session', 
                  rus.optional(bool))
    @rus.returns (rus.nothing)
    def __init__ (self, default=True):
        """
        default: bool
        ret:     None
        """

        simple_base = super  (Session, self)
        simple_base.__init__ ()

        self._logger = ru.get_logger ('radical.saga')

        # if the default session is expected, we point our context list to the
        # shared list of the default session singleton.  Otherwise, we create
        # a private list which is not populated.

        # a session also has a lease manager, for adaptors in this session to use.

        if  default :
            default_session     = DefaultSession ()
            self.contexts       = copy.deepcopy(default_session.contexts)
            self._lease_manager = default_session._lease_manager
        else :
            self.contexts       = _ContextList (session=self)

            # FIXME: at the moment, the lease manager is owned by the session.  
            # Howevwer, the pty layer is the main user of the lease manager,
            # and we thus keep the lease manager options in the pty subsection.  
            # So here we are, in the session, evaluating the pty config options...
            config = self.get_config ('saga.utils.pty')
            self._lease_manager = ru.LeaseManager (
                    max_pool_size = config['connection_pool_size'].get_value (),
                    max_pool_wait = config['connection_pool_wait'].get_value (),
                    max_obj_age   = config['connection_pool_ttl'].get_value ()
                    )


    # ----------------------------------------------------------------
    #
    @rus.takes   ('Session')
    @rus.returns (basestring)
    def __str__  (self):
        """String represenation."""

        return "Registered contexts: %s" % (str(self.contexts))


    # ----------------------------------------------------------------
    #
    @rus.takes      ('Session', 
                     saga.context.Context)
    @rus.returns    (rus.nothing)
    def add_context (self, ctx) :
        """
        ctx:     saga.Context
        ret:     None

        Add a security L{Context} to the session.
        It is encouraged to use the L{contexts} property instead. 
        """

        return self.contexts.insert (0, ctx=ctx, session=self)


    # ----------------------------------------------------------------
    #
    @rus.takes   ('Session', 
                  saga.context.Context)
    @rus.returns (rus.nothing)
    def remove_context (self, ctx) :
        """
        ctx:     saga.Context
        ret:     None

        Remove a security L{Context} from the session.
        It is encouraged to use the L{contexts} property instead.
        """

        if ctx in self.contexts :
            self.contexts.remove (ctx)


    # ----------------------------------------------------------------
    #
    @rus.takes   ('Session')
    @rus.returns (rus.list_of (saga.context.Context))
    def list_contexts  (self) :
        """
        ret:     list[saga.Context]
        
        Retrieve all L{Context} objects attached to the session.
        It is encouraged to use the L{contexts} property instead.
        """

        return self.contexts


    # ----------------------------------------------------------------
    #
    @rus.takes   ('Session')
    @rus.returns (dict)
    def get_config (self, section="") :
        """
        ret:     radical.utils.Configuration
        
        Return the session configuration (optional a specific section).
        """

        return saga.engine.engine.Engine ().get_config (section)



# ------------------------------------------------------------------------------
#
class DefaultSession(Session):

    __metaclass__ = ru.Singleton


    # --------------------------------------------------------------------------
    #
    @rus.takes   ('DefaultSession')
    @rus.returns (rus.nothing)
    def __init__ (self):

        # the default session picks up default contexts, from all context
        # adaptors.  To implemented, we have to do some legwork: get the engine,
        # dig through the registered context adaptors, and ask each of them for
        # default contexts.
        
        super(DefaultSession, self).__init__(default=False)

        _engine = saga.engine.engine.Engine()

        if not 'saga.Context' in _engine._adaptor_registry :
            self._logger.warn ("no context adaptors found")
            return

        for schema   in _engine._adaptor_registry['saga.Context'] :
            for info in _engine._adaptor_registry['saga.Context'][schema] :

                default_ctxs = []

                try : 
                    default_ctxs = info['adaptor_instance']._get_default_contexts ()

                except se.SagaException as e :
                    self._logger.debug   ("adaptor %s failed to provide default" \
                                          "contexts: %s" % (info['adaptor_name'], e))
                    continue
                    

                for default_ctx in default_ctxs :

                    try :
                        self.contexts.append (ctx=default_ctx, session=self)
                        self._logger.debug   ("default context [%-20s] : %s" \
                                         %   (info['adaptor_name'], default_ctx))

                    except se.SagaException as e :
                        self._logger.debug   ("skip default context [%-20s] : %s : %s" \
                                         %   (info['adaptor_name'], default_ctx, e))
                        continue


# ------------------------------------------------------------------------------

