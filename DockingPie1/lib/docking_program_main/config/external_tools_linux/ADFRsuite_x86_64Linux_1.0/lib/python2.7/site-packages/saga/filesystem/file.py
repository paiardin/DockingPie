
__author__    = "Andre Merzky, Ole Weidner, Alexander Grill"
__copyright__ = "Copyright 2012-2013, The SAGA Project"
__license__   = "MIT"


import radical.utils.signatures  as rus

import saga.adaptors.base        as sab
import saga.session              as ss
import saga.task                 as st
import saga.url                  as surl
import saga.namespace.entry      as nsentry

from   saga.filesystem.constants import *
from   saga.constants            import SYNC, ASYNC, TASK

# ------------------------------------------------------------------------------
#
class File (nsentry.Entry) :
    """
    Represents a local or remote file.

    The saga.filesystem.File class represents, as the name indicates,
    a file on some (local or remote) filesystem.  That class offers
    a number of operations on that file, such as copy, move and remove::
    
        # get a file handle
        file = saga.filesystem.File("sftp://localhost/tmp/data/data.bin")
    
        # copy the file
        file.copy ("sftp://localhost/tmp/data/data.bak")

        # move the file
        file.move ("sftp://localhost/tmp/data/data.new")
    """

    # --------------------------------------------------------------------------
    #
    @rus.takes   ('File', 
                  rus.optional ((surl.Url, basestring)), 
                  rus.optional (int, rus.nothing), 
                  rus.optional (ss.Session),
                  rus.optional (sab.Base), 
                  rus.optional (dict), 
                  rus.optional (rus.one_of (SYNC, ASYNC, TASK)))
    @rus.returns (rus.nothing)
    def __init__ (self, url=None, flags=READ, session=None, 
                  _adaptor=None, _adaptor_state={}, _ttype=None) : 
        """
        __init__(url, flags=READ, session)

        Construct a new file object

        :param url:     Url of the (remote) file
        :type  url:     :class:`saga.Url` 

        :fgs:   :ref:`filesystemflags`
        :param session: :class:`saga.Session`

        The specified file is expected to exist -- otherwise a DoesNotExist
        exception is raised.  Also, the URL must point to a file (not to
        a directory), otherwise a BadParameter exception is raised.

        Example::

            # get a file handle
            file = saga.filesystem.File("sftp://localhost/tmp/data/data.bin")
    
            # print the file's size
            print file.get_size ()
        """

        # param checks
        if  not flags : flags = 0
        url = surl.Url (url)

        if  not url.schema :
            url.schema = 'file'

        if  not url.host :
            url.host = 'localhost'

        self._nsentry = super  (File, self)
        self._nsentry.__init__ (url, flags, session, 
                                _adaptor, _adaptor_state, _ttype=_ttype)

    # --------------------------------------------------------------------------
    #
    @classmethod
    @rus.takes   ('File', 
                  rus.optional ((surl.Url, basestring)), 
                  rus.optional (int, rus.nothing), 
                  rus.optional (ss.Session),
                  rus.optional (rus.one_of (SYNC, ASYNC, TASK)))
    @rus.returns (st.Task)
    def create (cls, url=None, flags=READ, session=None, ttype=None) :
        """
        create(url, flags, session)

        url:       saga.Url
        flags:     saga.replica.flags enum
        session:   saga.Session
        ttype:     saga.task.type enum
        ret:       saga.Task
        """
        if  not flags : flags = 0
        _nsentry = super (File, cls)
        return _nsentry.create (url, flags, session, ttype=ttype)


    # ----------------------------------------------------------------
    #
    # filesystem methods
    #
    # --------------------------------------------------------------------------
    #
    @rus.takes   ('File', 
                  rus.optional (rus.one_of (SYNC, ASYNC, TASK)))
    @rus.returns ((bool, st.Task))
    def is_file (self, ttype=None) :
        """
        is_file()

        Returns `True` if instance points to a file, `False` otherwise. 
        """
        return self._adaptor.is_file_self (ttype=ttype)

  
    # --------------------------------------------------------------------------
    #
    @rus.takes   ('File', 
                  rus.optional (rus.one_of (SYNC, ASYNC, TASK)))
    @rus.returns ((int, st.Task))
    def get_size (self, ttype=None) :
        '''
        get_size()
        
        Returns the size (in bytes) of a file.

           Example::

               # get a file handle
               file = saga.filesystem.File("sftp://localhost/tmp/data/data.bin")
    
               # print the file's size
               print file.get_size ()

        '''
        return self._adaptor.get_size_self (ttype=ttype)

  
    # --------------------------------------------------------------------------
    #
    @rus.takes   ('File', 
                  rus.optional (int),
                  rus.optional (rus.one_of (SYNC, ASYNC, TASK)))
    @rus.returns ((basestring, st.Task))
    def read     (self, size=None, ttype=None) :
        '''
        size :    int
        ttype:    saga.task.type enum
        ret:      string / bytearray / saga.Task
        '''
        return self._adaptor.read (size, ttype=ttype)

    # --------------------------------------------------------------------------
    #
    @rus.takes   ('File', 
                  rus.optional (bool))
    @rus.returns (st.Task)
    def close     (self, kill=True, ttype=None) :
        '''
        kill :    bool
        ttype:    saga.task.type enum
        ret:      string / bytearray / saga.Task
        '''
        return self._adaptor.close ()
  
    # --------------------------------------------------------------------------
    #
    @rus.takes   ('File', 
                  basestring,
                  rus.optional (rus.one_of (SYNC, ASYNC, TASK)))
    @rus.returns ((int, st.Task))
    def write    (self, data, ttype=None) :
        '''
        data :    string / bytearray
        ttype:    saga.task.type enum
        ret:      int / saga.Task
        '''
        return self._adaptor.write (data, ttype=ttype)

  
    # --------------------------------------------------------------------------
    #
    @rus.takes   ('File', 
                  int,
                  rus.optional (rus.one_of (START, CURRENT, END )),
                  rus.optional (rus.one_of (SYNC,  ASYNC,   TASK)))
    @rus.returns ((int, st.Task))
    def seek     (self, offset, whence=START, ttype=None) :
        '''
        offset:   int
        whence:   seek_mode enum
        ttype:    saga.task.type enum
        ret:      int / saga.Task
        '''
        return self._adaptor.seek (offset, whence, ttype=ttype)

  
    # --------------------------------------------------------------------------
    #
    @rus.takes   ('File', 
                  rus.list_of  (rus.tuple_of (int)),
                  rus.optional (rus.one_of (SYNC, ASYNC, TASK)))
    @rus.returns ((basestring, st.Task))
    def read_v   (self, iovecs, ttype=None) :
        '''
        iovecs:   list [tuple (int, int)]
        ttype:    saga.task.type enum
        ret:      list [bytearray] / saga.Task
        '''
        return self._adaptor.read_v (iovecs, ttype=ttype)

  
    # --------------------------------------------------------------------------
    #
    @rus.takes   ('File', 
                  rus.list_of  (rus.tuple_of ((int, basestring))),
                  rus.optional (rus.one_of (SYNC, ASYNC, TASK)))
    @rus.returns ((rus.list_of (int), st.Task))
    def write_v (self, data, ttype=None) :
        '''
        data:     list [tuple (int, string / bytearray)]
        ttype:    saga.task.type enum
        ret:      list [int] / saga.Task
        '''
        return self._adaptor.write_v (data, ttype=ttype)

  
    # --------------------------------------------------------------------------
    #
    @rus.takes   ('File', 
                  basestring,
                  rus.optional (rus.one_of (SYNC, ASYNC, TASK)))
    @rus.returns ((int, st.Task))
    def size_p (self, pattern, ttype=None) :
        '''
        pattern:  string 
        ttype:    saga.task.type enum
        ret:      int / saga.Task
        '''
        return self._adaptor.size_p (pattern, ttype=ttype)
  

    # --------------------------------------------------------------------------
    #
    @rus.takes   ('File', 
                  basestring,
                  rus.optional (rus.one_of (SYNC, ASYNC, TASK)))
    @rus.returns ((basestring, st.Task))
    def read_p (self, pattern, ttype=None) :
        '''
        pattern:  string
        ttype:    saga.task.type enum
        ret:      string / bytearray / saga.Task
        '''
        return self._adaptor.read_p (pattern, ttype=ttype)

  
    # --------------------------------------------------------------------------
    #
    @rus.takes   ('File', 
                  basestring,
                  basestring,
                  rus.optional (rus.one_of (SYNC, ASYNC, TASK)))
    @rus.returns ((int, st.Task))
    def write_p (self, pattern, data, ttype=None) :
        '''
        pattern:  string
        data:     string / bytearray
        ttype:    saga.task.type enum
        ret:      int / saga.Task
        '''
        return self._adaptor.write_p (pattern, data, ttype=ttype)

  
    # --------------------------------------------------------------------------
    #
    @rus.takes   ('File', 
                  rus.optional (rus.one_of (SYNC, ASYNC, TASK)))
    @rus.returns ((rus.list_of (basestring), st.Task))
    def modes_e (self, ttype=None) :
        '''
        ttype:    saga.task.type enum
        ret:      list [string] / saga.Task
        '''
        return self._adaptor.modes_e (ttype=ttype)

  
    # --------------------------------------------------------------------------
    #
    @rus.takes   ('File', 
                  basestring,
                  basestring,
                  rus.optional (rus.one_of (SYNC, ASYNC, TASK)))
    @rus.returns ((int, st.Task))
    def size_e (self, emode, spec, ttype=None) :
        '''
        emode:    string
        spec:     string
        ttype:    saga.task.type enum
        ret:      int / saga.Task
        '''
        return self._adaptor.size_e (emode, spec, ttype=ttype)

  
    # --------------------------------------------------------------------------
    #
    @rus.takes   ('File', 
                  basestring,
                  basestring,
                  rus.optional (rus.one_of (SYNC, ASYNC, TASK)))
    @rus.returns ((basestring, st.Task))
    def read_e (self, emode, spec, ttype=None) :
        '''
        emode:    string
        spec:     string
        ttype:    saga.task.type enum
        ret:      bytearray / saga.Task
        '''
        return self._adaptor.read_e (emode, spec, ttype=ttype)

  
    # --------------------------------------------------------------------------
    #
    @rus.takes   ('File', 
                  basestring,
                  basestring,
                  basestring,
                  rus.optional (rus.one_of (SYNC, ASYNC, TASK)))
    @rus.returns ((int, st.Task))
    def write_e (self, emode, spec, data, ttype=None) :
        '''
        emode:    string
        spec:     string
        data:     string / bytearray
        ttype:    saga.task.type enum
        ret:      int / saga.Task
        '''
        return self._adaptor.read_e (emode, spec, data, ttype=ttype)

  
    size    = property (get_size)  # int
    modes_e = property (modes_e)   # list [string]
  
  


