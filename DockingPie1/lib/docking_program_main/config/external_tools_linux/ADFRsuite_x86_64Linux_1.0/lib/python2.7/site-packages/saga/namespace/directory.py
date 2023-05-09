
__author__    = "Andre Merzky, Ole Weidner, Alexander Grill"
__copyright__ = "Copyright 2012-2013, The SAGA Project"
__license__   = "MIT"


import radical.utils.signatures  as rus

import saga.adaptors.base        as sab
import saga.session              as ss
import saga.task                 as st
import saga.url                  as surl

import entry

from   saga.namespace.constants  import *
from   saga.constants            import SYNC, ASYNC, TASK


# ------------------------------------------------------------------------------
#
class Directory (entry.Entry) :
    '''
    Represents a SAGA directory as defined in GFD.90
    
    The saga.namespace.Directory class represents, as the name indicates,
    a directory on some (local or remote) namespace.  That class offers
    a number of operations on that directory, such as listing its contents,
    copying entries, or creating subdirectories::
    
        # get a directory handle
        dir = saga.namespace.Directory("sftp://localhost/tmp/")
    
        # create a subdir
        dir.make_dir ("data/")
    
        # list contents of the directory
        entries = dir.list ()
    
        # copy *.dat entries into the subdir
        for f in entries :
            if f ^ '^.*\.dat$' :
                dir.copy (f, "sftp://localhost/tmp/data/")


    Implementation note:
    ^^^^^^^^^^^^^^^^^^^^

    The SAGA API Specification (GFD.90) prescribes method overloading on method
    signatures, but that is not supported by Python (Python only does method
    overwriting).  So we implement one generic method version here, and do the
    case switching based on the provided parameter set.
    '''

    # --------------------------------------------------------------------------
    #
    @rus.takes   ('Directory', 
                  rus.optional ((surl.Url, basestring)), 
                  rus.optional (int, rus.nothing), 
                  rus.optional (ss.Session),
                  rus.optional (sab.Base), 
                  rus.optional (dict), 
                  rus.optional (rus.one_of (SYNC, ASYNC, TASK)))
    @rus.returns (rus.nothing)
    def __init__ (self, url=None, flags=None, session=None, 
                  _adaptor=None, _adaptor_state={}, _ttype=None) : 
        '''
        :param url: Url of the (remote) entry system directory.
        :type  url: :class:`saga.Url` 

        flags:     flags enum
        session:   saga.Session
        ret:       obj
        
        Construct a new directory object

        The specified directory is expected to exist -- otherwise
        a DoesNotExist exception is raised.  Also, the URL must point to
        a directory (not to an entry), otherwise a BadParameter exception is
        raised.

        Example::

            # open some directory
            dir = saga.namespace.Directory("sftp://localhost/tmp/")

            # and list its contents
            entries = dir.list ()

        '''

        if  not flags : flags = 0
        self._nsentry = super  (Directory, self)
        self._nsentry.__init__ (url, flags, session, 
                                _adaptor, _adaptor_state, _ttype=_ttype)


    # --------------------------------------------------------------------------
    #
    @classmethod
    @rus.takes   ('Directory', 
                  rus.optional ((surl.Url, basestring)), 
                  rus.optional (int, rus.nothing), 
                  rus.optional (ss.Session),
                  rus.optional (rus.one_of (SYNC, ASYNC, TASK)))
    @rus.returns (st.Task)
    def create (cls, url=None, flags=None, session=None, ttype=None) :
        '''
        url:       saga.Url
        flags:     saga.namespace.flags enum
        session:   saga.Session
        ttype:     saga.task.type enum
        ret:       saga.Task
        '''

        if  not flags : flags = 0
        _nsentry = super (Directory, cls)
        return _nsentry.create (url, flags, session, ttype=ttype)


    # --------------------------------------------------------------------------
    #
    @rus.takes   ('Directory', 
                  (surl.Url, basestring), 
                  rus.optional (int, rus.nothing),
                  rus.optional (rus.one_of (SYNC, ASYNC, TASK)))
    @rus.returns ((entry.Entry, st.Task))
    def open (self, name, flags=None, ttype=None) :
        '''
        name:     saga.Url
        flags:    saga.namespace.flags enum
        ttype:    saga.task.type enum
        ret:      saga.namespace.Entry / saga.Task
        '''
        if  not flags : flags = 0
        url = surl.Url(name)
        return self._adaptor.open (url, flags, ttype=ttype)


    # --------------------------------------------------------------------------
    #
    @rus.takes   ('Directory', 
                  (surl.Url, basestring), 
                  rus.optional (int, rus.nothing),
                  rus.optional (rus.one_of (SYNC, ASYNC, TASK)))
    @rus.returns (('Directory', st.Task))
    def open_dir (self, path, flags=None, ttype=None) :
        '''
        :param path: name/path of the directory to open
        :param flags: directory creation flags

        ttype:    saga.task.type enum
        ret:      saga.namespace.Directory / saga.Task
        
        Open and return a new directoy

           The call opens and returns a directory at the given location.

           Example::

               # create a subdir 'data' in /tmp
               dir = saga.namespace.Directory("sftp://localhost/tmp/")
               data = dir.open_dir ('data/', saga.namespace.Create)
        '''
        if  not flags : flags = 0
        return self._adaptor.open_dir (surl.Url(path), flags, ttype=ttype)

    
    # ----------------------------------------------------------------
    #
    # namespace directory methods
    #
    @rus.takes   ('Directory', 
                  (surl.Url, basestring), 
                  rus.optional (int, rus.nothing),
                  rus.optional (rus.one_of (SYNC, ASYNC, TASK)))
    @rus.returns ((rus.nothing, st.Task))
    def make_dir (self, tgt, flags=0, ttype=None) :
        '''
        :param tgt:   name/path of the new directory
        :param flags: directory creation flags

        ttype:         saga.task.type enum
        ret:           None / saga.Task
        
        Create a new directoy

        The call creates a directory at the given location.

        Example::

            # create a subdir 'data' in /tmp
            dir = saga.namespace.Directory("sftp://localhost/tmp/")
            dir.make_dir ('data/')
        '''
        if  not flags : flags = 0
        return self._adaptor.make_dir (surl.Url (tgt), flags, ttype=ttype)
  
    
    # --------------------------------------------------------------------------
    #
    @rus.takes   ('Directory', 
                  (surl.Url, basestring), 
                  rus.optional (int, rus.nothing),
                  rus.optional (rus.one_of (SYNC, ASYNC, TASK)))
    @rus.returns ((rus.nothing, st.Task))
    def change_dir (self, url, flags=0, ttype=None) :
        '''
        url:           saga.Url
        flags:         flags enum
        ttype:         saga.task.type enum
        ret:           None / saga.Task
        '''
        if  not flags : flags = 0
        return self._adaptor.change_dir (url, flags=flags, ttype=ttype)
  
    
    # --------------------------------------------------------------------------
    #
    @rus.takes   ('Directory', 
                  rus.optional (basestring),
                  rus.optional (int, rus.nothing),
                  rus.optional (rus.one_of (SYNC, ASYNC, TASK)))
    @rus.returns ((rus.list_of (surl.Url), st.Task))
    def list (self, pattern=None, flags=0, ttype=None) :
        '''
        :param pattern: Entry name pattern (like POSIX 'ls', e.g. '\*.txt')

        flags:         flags enum
        ttype:         saga.task.type enum
        ret:           list [saga.Url] / saga.Task
        
        List the directory's content

        The call will return a list of entries and subdirectories within the
        directory::

            # list contents of the directory
            for f in dir.list() :
                print f
        '''
        if  not flags : flags = 0
        return self._adaptor.list (pattern, flags, ttype=ttype)


    # --------------------------------------------------------------------------
    #
    @rus.takes   ('Directory', 
                  (surl.Url, basestring),
                  rus.optional (rus.one_of (SYNC, ASYNC, TASK)))
    @rus.returns ((bool, st.Task))
    def exists (self, path, ttype=None) :

        '''
        :param path: path of the entry to check

        ttype:         saga.task.type enum
        ret:           bool / saga.Task
        
        Returns True if path exists, False otherwise. 


        Example::

            # inspect an entry
            dir  = saga.namespace.Directory("sftp://localhost/tmp/")
            if dir.exists ('data'):
                # do something
        '''
        return self._adaptor.exists (path, ttype=ttype)
  
  
    # --------------------------------------------------------------------------
    #
    @rus.takes   ('Directory', 
                  rus.optional (basestring),
                  rus.optional (int, rus.nothing),
                  rus.optional (rus.one_of (SYNC, ASYNC, TASK)))
    @rus.returns ((rus.list_of (surl.Url), st.Task))
    def find (self, pattern, flags=RECURSIVE, ttype=None) :
        '''
        pattern:       string
        flags:         flags enum
        ttype:         saga.task.type enum
        ret:           list [saga.Url] / saga.Task
        '''
        if  not flags : flags = 0
        return self._adaptor.find (pattern, flags, ttype=ttype)
  
    
    # --------------------------------------------------------------------------
    #
    @rus.takes   ('Directory', 
                  rus.optional (rus.one_of (SYNC, ASYNC, TASK)))
    @rus.returns ((int, st.Task))
    def get_num_entries (self, ttype=None) :
        '''
        ttype:         saga.task.type enum
        ret:           int / saga.Task
        '''
        return self._adaptor.get_num_entries (ttype=ttype)
  
    
    # --------------------------------------------------------------------------
    #
    @rus.takes   ('Directory', 
                  int, 
                  rus.optional (rus.one_of (SYNC, ASYNC, TASK)))
    @rus.returns ((surl.Url, st.Task))
    def get_entry (self, num, ttype=None) :
        '''
        num:           int 
        ttype:         saga.task.type enum
        ret:           saga.Url / saga.Task
        '''
        return self._adaptor.get_entry (num, ttype=ttype)


    # ----------------------------------------------------------------
    #
    # methods overloaded from namespace.Entry
    #
    @rus.takes   ('Directory', 
                  (surl.Url, basestring), 
                  rus.optional ((surl.Url, basestring)), 
                  rus.optional (int, rus.nothing),
                  rus.optional (rus.one_of (SYNC, ASYNC, TASK)))
    @rus.returns ((rus.nothing, st.Task))
    def copy (self, url_1, url_2=None, flags=0, ttype=None) :
        '''
        :param src: path of the entry to copy
        :param tgt: absolute URL of target name or directory
        
        url_1:         saga.Url
        url_2:         saga.Url / None
        flags:         flags enum / None
        ttype:         saga.task.type enum / None
        ret:           None / saga.Task
        
        Copy an entry from source to target

        The source is copied to the given target directory.  The path of the
        source can be relative::

            # copy an entry
            dir = saga.namespace.Directory("sftp://localhost/tmp/")
            dir.copy ("./data.bin", "sftp://localhost/tmp/data/")
        '''

        # FIXME: re-implement the url switching (commented out below)

        if  not flags : flags = 0
        if url_2  :  return self._adaptor.copy (url_1, url_2, flags, ttype=ttype) 
        else      :  return self._nsentry.copy (url_1,        flags, ttype=ttype)

    
    # --------------------------------------------------------------------------
    #
    @rus.takes   ('Directory', 
                  (surl.Url, basestring), 
                  rus.optional ((surl.Url, basestring)), 
                  rus.optional (int, rus.nothing),
                  rus.optional (rus.one_of (SYNC, ASYNC, TASK)))
    @rus.returns ((rus.nothing, st.Task))
    def link (self, url_1, url_2=None, flags=0, ttype=None) :
        '''
        src:           saga.Url
        tgt:           saga.Url
        flags:         flags enum
        ttype:         saga.task.type enum
        ret:           None / saga.Task
        '''
        if  not flags : flags = 0
        if url_2  :  return self._adaptor.link (url_1, url_2, flags, ttype=ttype)
        else      :  return self._nsentry.link (url_1,        flags, ttype=ttype)
  
    
    # --------------------------------------------------------------------------
    #
    @rus.takes   ('Directory', 
                  (surl.Url, basestring), 
                  rus.optional ((surl.Url, basestring)), 
                  rus.optional (int, rus.nothing),
                  rus.optional (rus.one_of (SYNC, ASYNC, TASK)))
    @rus.returns ((rus.nothing, st.Task))
    def move (self, url_1, url_2=None, flags=0, ttype=None) :
        '''
        :param src: path of the entry to copy
        :param tgt: absolute URL of target directory

        flags:         flags enum
        ttype:         saga.task.type enum
        ret:           None / saga.Task
        
        Move an entry from source to target

        The source is moved to the given target directory.  The path of the
        source can be relative::

            # copy an entry
            dir = saga.namespace.Directory("sftp://localhost/tmp/")
            dir.move ("./data.bin", "sftp://localhost/tmp/data/")

        '''
        if  not flags : flags = 0
        if url_2  :  return self._adaptor.move (url_1, url_2, flags, ttype=ttype)
        else      :  return self._nsentry.move (url_1,        flags, ttype=ttype)
  
    
    # --------------------------------------------------------------------------
    #
    @rus.takes   ('Directory', 
                  (surl.Url, basestring), 
                  rus.optional (int, rus.nothing),
                  rus.optional (rus.one_of (SYNC, ASYNC, TASK)))
    @rus.returns ((rus.nothing, st.Task))
    def remove (self, tgt=None, flags=0, ttype=None) :
        '''
        tgt:           saga.Url
        flags:         flags enum
        ttype:         saga.task.type enum
        ret:           None / saga.Task
        '''
        if  not flags : flags = 0
        if tgt    :  return self._adaptor.remove (tgt, flags, ttype=ttype)
        else      :  return self._nsentry.remove (     flags, ttype=ttype)
  

    # --------------------------------------------------------------------------
    #
    @rus.takes   ('Directory', 
                  rus.optional ((surl.Url, basestring)),
                  rus.optional (rus.one_of (SYNC, ASYNC, TASK)))
    @rus.returns ((bool, st.Task))
    def is_dir (self, tgt=None, ttype=None) :
        '''
        tgt:           saga.Url / None
        ttype:         saga.task.type enum
        ret:           bool / saga.Task
        
        Returns True if path is a directory, False otherwise. 

        Example::

            # inspect an entry
            dir  = saga.namespace.Directory("sftp://localhost/tmp/")
            if dir.is_dir ('data'):
                # do something
        '''
        if tgt    :  return self._adaptor.is_dir (tgt, ttype=ttype)
        else      :  return self._nsentry.is_dir (     ttype=ttype)
  
    
    # --------------------------------------------------------------------------
    #
    @rus.takes   ('Directory', 
                  rus.optional ((surl.Url, basestring)),
                  rus.optional (rus.one_of (SYNC, ASYNC, TASK)))
    @rus.returns ((bool, st.Task))
    def is_entry (self, tgt=None, ttype=None) :
        '''
        tgt:           saga.Url / None
        ttype:         saga.task.type enum
        ret:           bool / saga.Task
        '''
        if tgt    :  return self._adaptor.is_entry (tgt, ttype=ttype)
        else      :  return self._nsentry.is_entry (     ttype=ttype)
  
    
    # --------------------------------------------------------------------------
    #
    @rus.takes   ('Directory', 
                  rus.optional ((surl.Url, basestring)),
                  rus.optional (rus.one_of (SYNC, ASYNC, TASK)))
    @rus.returns ((bool, st.Task))
    def is_link (self, tgt=None, ttype=None) :
        '''
        tgt:           saga.Url / None
        ttype:         saga.task.type enum
        ret:           bool / saga.Task
        '''
        if tgt    :  return self._adaptor.is_link (tgt, ttype=ttype)
        else      :  return self._nsentry.is_link (     ttype=ttype)


    # --------------------------------------------------------------------------
    #
    @rus.takes   ('Directory', 
                  rus.optional ((surl.Url, basestring)),
                  rus.optional (rus.one_of (SYNC, ASYNC, TASK)))
    @rus.returns ((surl.Url, st.Task))
    def read_link (self, tgt=None, ttype=None) :
        '''
        tgt:           saga.Url / None
        ttype:         saga.task.type enum
        ret:           saga.Url / saga.Task
        '''

        if tgt    :  return self._adaptor.read_link (tgt, ttype=ttype)
        else      :  return self._nsentry.read_link (     ttype=ttype)
  



