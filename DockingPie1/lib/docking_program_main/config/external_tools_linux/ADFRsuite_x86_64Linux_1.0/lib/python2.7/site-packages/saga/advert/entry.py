
__author__    = "Andre Merzky"
__copyright__ = "Copyright 2012-2013, The SAGA Project"
__license__   = "MIT"


import radical.utils.signatures as rus

import saga.adaptors.base       as sab
import saga.attributes          as sa
import saga.session             as ss
import saga.task                as st
import saga.url                 as surl
import saga.namespace.entry     as nsentry

from   saga.advert.constants    import *
from   saga.constants           import SYNC, ASYNC, TASK


# ------------------------------------------------------------------------------
#
class Entry (nsentry.Entry, sa.Attributes) :


    # --------------------------------------------------------------------------
    #
    @rus.takes   ('Entry', 
                  rus.optional ((surl.Url, basestring)), 
                  rus.optional (int, rus.nothing),
                  rus.optional (ss.Session), 
                  rus.optional (sab.Base),
                  rus.optional (dict),
                  rus.optional (rus.one_of (SYNC, ASYNC, TASK)))
    @rus.returns (rus.nothing)
    def __init__ (self, url=None, flags=READ, session=None, 
                  _adaptor=None, _adaptor_state={}, _ttype=None) : 
        '''
        url:       saga.Url
        flags:     flags enum
        session:   saga.Session
        ret:       obj
        '''

        # param checks
        url = surl.Url (url)

        self._nsentry = super  (Entry, self)
        self._nsentry.__init__ (url, flags, session, 
                                _adaptor, _adaptor_state, _ttype=_ttype)


        # set attribute interface properties
        self._attributes_allow_private (True)
        self._attributes_camelcasing   (True)
        self._attributes_extensible    (True, getter=self._attribute_getter, 
                                              setter=self._attribute_setter,
                                              lister=self._attribute_lister,
                                              caller=self._attribute_caller)

        # register properties with the attribute interface 
        self._attributes_register   (ATTRIBUTE, None, sa.STRING, sa.SCALAR, sa.READONLY)
        self._attributes_register   (OBJECT,    None, sa.ANY,    sa.SCALAR, sa.READONLY)
        self._attributes_register   (TTL,       None, sa.INT,    sa.SCALAR, sa.WRITEABLE)

        self._attributes_set_setter (TTL,    self.set_ttl)
        self._attributes_set_getter (TTL,    self.get_ttl)

        self._attributes_set_setter (OBJECT, self.store_object)
        self._attributes_set_getter (OBJECT, self.retrieve_object)



    # --------------------------------------------------------------------------
    #
    @classmethod
    @rus.takes   ('Entry', 
                  rus.optional ((surl.Url, basestring)), 
                  rus.optional (int, rus.nothing), 
                  rus.optional (ss.Session), 
                  rus.optional (rus.one_of (SYNC, ASYNC, TASK)))
    @rus.returns (st.Task)
    def create (cls, url=None, flags=READ, session=None, ttype=None) :
        '''
        url:       saga.Url
        flags:     saga.advert.flags enum
        session:   saga.Session
        ttype:     saga.task.type enum
        ret:       saga.Task
        '''

        if not flags : flags = 0
        _nsentry = super (Entry, cls)
        return _nsentry.create (url, flags, session, ttype=ttype)



    # --------------------------------------------------------------------------
    #
    # attribute methods
    #
    # NOTE: we do not yet pass ttype, as async calls are not yet supported by
    # the attribute interface
    #
    @rus.takes   ('Entry', 
                  basestring,
                  rus.optional (rus.one_of (SYNC, ASYNC, TASK)))
    @rus.returns ((rus.anything, st.Task))
    def _attribute_getter (self, key, ttype=None) :

        return self._adaptor.attribute_getter (key)


    # --------------------------------------------------------------------------
    #
    @rus.takes   ('Entry', 
                  basestring,
                  rus.anything,
                  rus.optional (rus.one_of (SYNC, ASYNC, TASK)))
    @rus.returns ((rus.nothing, st.Task))
    def _attribute_setter (self, key, val, ttype=None) :

        return self._adaptor.attribute_setter (key, val)


    # --------------------------------------------------------------------------
    #
    @rus.takes   ('Entry', 
                  rus.optional (rus.one_of (SYNC, ASYNC, TASK)))
    @rus.returns ((rus.list_of (rus.anything), st.Task))
    def _attribute_lister (self, ttype=None) :

        return self._adaptor.attribute_lister ()


    # --------------------------------------------------------------------------
    #
    @rus.takes   ('Entry', 
                  basestring, 
                  int, 
                  callable, 
                  rus.optional (rus.one_of (SYNC, ASYNC, TASK)))
    @rus.returns ((rus.anything, st.Task))
    def _attribute_caller (self, key, id, cb, ttype=None) :

        return self._adaptor.attribute_caller (key, id, cb)



    # --------------------------------------------------------------------------
    #
    # advert methods
    #
    @rus.takes   ('Entry', 
                  float, 
                  rus.optional (rus.one_of (SYNC, ASYNC, TASK)))
    @rus.returns ((rus.nothing, st.Task))
    def set_ttl  (self, ttl=-1.0, ttype=None) : 
        """
        ttl :           int
        ttype:          saga.task.type enum
        ret:            None / saga.Task
        """

        return self._adaptor.set_ttl (ttl, ttype=ttype)

     
    # --------------------------------------------------------------------------
    #
    @rus.takes   ('Entry', 
                  rus.optional (rus.one_of (SYNC, ASYNC, TASK)))
    @rus.returns ((float, st.Task))
    def get_ttl  (self, ttype=None) : 
        """
        ttype:          saga.task.type enum
        ret:            int / saga.Task
        """

        return self._adaptor.get_ttl (ttype=ttype)


     
    # --------------------------------------------------------------------------
    #
    @rus.takes   ('Entry', 
                  object,
                  rus.optional (rus.one_of (SYNC, ASYNC, TASK)))
    @rus.returns ((rus.nothing, st.Task))
    def store_object (self, object, ttype=None) : 
        """
        object :        <object type>
        ttype:          saga.task.type enum
        ret:            None / saga.Task
        """
        return self._adaptor.store_object (object, ttype=ttype)

  
    # --------------------------------------------------------------------------
    #
    @rus.takes   ('Entry', 
                  rus.optional (rus.one_of (SYNC, ASYNC, TASK)))
    @rus.returns ((object, st.Task))
    def retrieve_object (self, ttype=None) : 
        """
        ttype:          saga.task.type enum
        ret:            any / saga.Task
        """
        return self._adaptor.retrieve_object (ttype=ttype)

     
    # --------------------------------------------------------------------------
    #
    @rus.takes   ('Entry', 
                  rus.optional (rus.one_of (SYNC, ASYNC, TASK)))
    @rus.returns ((rus.nothing, st.Task))
    def delete_object (self, ttype=None) : 
        """
        ttype:          saga.task.type enum
        ret:            None / saga.Task
        """
        return self._adaptor.delete_object (ttype=ttype)


