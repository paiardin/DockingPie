
__author__    = "Andre Merzky, Ole Weidner"
__copyright__ = "Copyright 2012-2013, The SAGA Project"
__license__   = "MIT"


import radical.utils.signatures as rus

import saga.adaptors.base       as sab
import saga.async               as async
import saga.task                as st
import saga.base                as sb
import saga.session             as ss
import saga.exceptions          as se
import saga.attributes          as sa
import saga.constants           as sc
import saga.url                 as surl
import constants                as const
import description              as descr
import resource                 as resrc

from   saga.resource.constants  import *
from   saga.constants           import SYNC, ASYNC, TASK


# ------------------------------------------------------------------------------
# 
class Manager (sb.Base, async.Async) :
    """
    In the context of SAGA-Python, a *ResourceManager* is a service which asserts
    control over a set of resources.  That manager can, on request, render
    control over subsets of those resources (resource slices) to an application.

    This :class:`Manager` class represents the contact point to such
    ResourceManager instances -- the application can thus acquire compute, data
    or network resources, according to some resource specification, for a bound
    or unbound amount of time. 
    """
    
    # --------------------------------------------------------------------------
    # 
    @rus.takes   ('Manager', 
                  rus.optional (basestring, surl.Url), 
                  rus.optional (ss.Session),
                  rus.optional (sab.Base), 
                  rus.optional (dict), 
                  rus.optional (rus.one_of (SYNC, ASYNC, TASK)))
    @rus.returns (rus.nothing)
    def __init__ (self, url=None, session=None,
                  _adaptor=None, _adaptor_state={}, _ttype=None) : 
        """
        __init__(url)

        Create a new Manager instance. Connect to a remote resource management endpoint.

        :type  url: :class:`saga.Url`
        :param url: resource management endpoint
        """

        # param checks
        _url   = surl.Url(url)
        scheme = _url.scheme.lower()

        if not session :
            session = ss.Session (default=True)

        self._base = super  (Manager, self)
        self._base.__init__ (scheme, _adaptor, _adaptor_state, 
                             _url, session, ttype=_ttype)


    # --------------------------------------------------------------------------
    #
    @classmethod
    @rus.takes   ('Manager', 
                  rus.optional ((surl.Url, basestring)), 
                  rus.optional (ss.Session),
                  rus.optional (rus.one_of (SYNC, ASYNC, TASK)))
    @rus.returns (st.Task)
    def create   (cls, url_in=None, session=None, ttype=sc.SYNC) :
        """ 
        This is the asynchronous class constructor, returning
        a :class:`saga:Task` instance.  For details on the accepted parameters,
        please see the description of :func:`__init__`.
        """

        return cls (url_in, session, _ttype=ttype)._init_task


    # --------------------------------------------------------------------------
    # 
    @rus.takes   ('Manager', 
                  rus.optional (rus.one_of (COMPUTE, STORAGE, NETWORK)),
                  rus.optional (rus.one_of (SYNC, ASYNC, TASK)))
    @rus.returns ((rus.list_of (basestring), st.Task))
    def list     (self, rtype=None, ttype=None) :
        """ 
        list(rtype=None)

        List known resource instances (which can be acquired). 
        Returns a list of IDs.  

        :type  rtype: None or enum (COMPUTE | STORAGE | NETWORK)
        :param rtype: filter for one or more resource types

        """
        return self._adaptor.list (rtype, ttype=ttype)


    # --------------------------------------------------------------------------
    # 
    @rus.takes   ('Manager', 
                  rus.optional (basestring),
                  rus.optional (rus.one_of (SYNC, ASYNC, TASK)))
    @rus.returns ((descr.Description, st.Task))
    def get_description (self, rid, ttype=None) :
        """ 
        get_description(rid)

        Get the resource :class:`Description` for the specified resource.

        :type  rid: str
        :param rid: identifies the resource to be described.
        """

        # TODO / NOTE: if rid is None, should we return a description of 
        # the managed resources?  
        return self._adaptor.get_description (id, ttype=ttype)

    # --------------------------------------------------------------------------
    # 
    @rus.takes   ('Manager', 
                  rus.optional (rus.one_of (COMPUTE, STORAGE, NETWORK)),
                  rus.optional (rus.one_of (SYNC, ASYNC, TASK)))
    @rus.returns ((rus.list_of (basestring), st.Task))
    def list_templates (self, rtype=None, ttype=None) :
        """
        list_templates(rtype=None)

        List template names available for the specified resource type(s).
        Returns a list of strings.

        :type  rtype: None or enum (COMPUTE | STORAGE | NETWORK)
        :param rtype: filter for one or more resource types

        """
        return self._adaptor.list_templates (rtype, ttype=ttype)


    # --------------------------------------------------------------------------
    # 
    @rus.takes   ('Manager', 
                  rus.optional (basestring),
                  rus.optional (rus.one_of (SYNC, ASYNC, TASK)))
    @rus.returns ((descr.Description, st.Task))
    def get_template (self, name, ttype=None) :
        """
        get_template(name)

        Get a :class:`Description` for the specified template.

        :type  name: str
        :param name: specifies the name of the template

        The returned resource description instance may not have all attributes
        filled, and may in fact not sufficiently complete to allow for
        successful resource acquisition.  The only guaranteed attribute in the
        returned description is `TEMPLATE`, containing the very template id
        specified in the call parameters.
        """

        return self._adaptor.get_template (name, ttype=ttype)


    # --------------------------------------------------------------------------
    # 
    @rus.takes   ('Manager', 
                  rus.optional (rus.one_of (COMPUTE, STORAGE, NETWORK)),
                  rus.optional (rus.one_of (SYNC, ASYNC, TASK)))
    @rus.returns ((rus.list_of (basestring), st.Task))
    def list_images (self, rtype=None, ttype=None) :
        """
        list_images(rtype=None)

        List image names available for the specified resource type(s).
        Returns a list of strings.

        :type  rtype: None or enum (COMPUTE | STORAGE | NETWORK)
        :param rtype: filter for one or more resource types
        """
        return self._adaptor.list_images (rtype, ttype=ttype)


    # --------------------------------------------------------------------------
    # 
    @rus.takes   ('Manager', 
                  basestring,
                  rus.optional (rus.one_of (SYNC, ASYNC, TASK)))
    @rus.returns ((dict, st.Task))
    def get_image (self, name, ttype=None) :
        """
        get_image(name)

        Get a description string for the specified image.

        :type  name: str
        :param name: specifies the image name
        """

        return self._adaptor.get_image (name, ttype=ttype)


    # --------------------------------------------------------------------------
    # 
    @rus.takes   ('Manager', 
                  (basestring, surl.Url, descr.Description),
                  rus.optional (rus.one_of (SYNC, ASYNC, TASK)))
    @rus.returns ((resrc.Resource, st.Task))
    def acquire  (self, spec, ttype=None) :
        """
        acquire(desc)

        Create a new :class:`saga.resource.Resource` handle for a 
        resource specified by the description.

        :type  spec: :class:`Description` or Url
        :param spec: specifies the resource

        Depending on the `RTYPE` attribute in the description, the
        returned resource may be a :class:`saga.resource.Compute`,
        :class:`saga.resource.Storage` or :class:`saga.resource.Network`
        instance. 

        If the `spec` parameter is 
        """

        if  isinstance (spec, surl.Url) or \
            isinstance (spec, basestring)  :

            return self._adaptor.acquire_by_id(spec, ttype=ttype)

        else :

            # make sure at least 'executable' is defined
            if spec.rtype is None:
                raise se.BadParameter ("No resource type defined in resource description")
    
            # spec_copy = descr.Description ()
            # spec._attributes_deep_copy (spec_copy)

            return self._adaptor.acquire (spec, ttype=ttype)


    # --------------------------------------------------------------------------
    # 
    @rus.takes   ('Manager', 
                  basestring,
                  rus.optional (rus.one_of (SYNC, ASYNC, TASK)))
    @rus.returns ((rus.nothing, st.Task))
    def destroy  (self, rid, ttype=None) :
        """
        destroy(rid)

        Destroy / release a resource. 

        :type  rid   : string
        :param rid   : identifies the resource to be released
        """

        return self._adaptor.destroy (rid, ttype=ttype)

  # FIXME: add
  # templates = property (list_templates, get_template)    # dict {string : Description}
  # images    = property (list_images,    get_image)       # dict {string : dict}
  # resources = property (list,           get_description) # dict {string : Description}




