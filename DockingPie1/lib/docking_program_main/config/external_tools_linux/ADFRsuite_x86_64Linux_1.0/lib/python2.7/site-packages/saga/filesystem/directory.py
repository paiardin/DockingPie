
__author__    = "Andre Merzky, Ole Weidner, Alexander Grill"
__copyright__ = "Copyright 2012-2013, The SAGA Project"
__license__   = "MIT"


import radical.utils.signatures  as rus

import saga.adaptors.base        as sab
import saga.session              as ss
import saga.task                 as st
import saga.url                  as surl
import saga.namespace.directory  as nsdir

from   saga.filesystem.constants import *
from   saga.constants            import SYNC, ASYNC, TASK


# ------------------------------------------------------------------------------
#
class Directory (nsdir.Directory) :
    """
    Represents a (remote) directory.
    
    The saga.filesystem.Directory class represents, as the name indicates,
    a directory on some (local or remote) filesystem.  That class offers
    a number of operations on that directory, such as listing its contents,
    copying files, or creating subdirectories::
    
        # get a directory handle
        dir = saga.filesystem.Directory("sftp://localhost/tmp/")
    
        # create a subdir
        dir.make_dir ("data/")
    
        # list contents of the directory
        files = dir.list ()
    
        # copy *.dat files into the subdir
        for f in files :
            if f ^ '^.*\.dat$' :
                dir.copy (f, "sftp://localhost/tmp/data/")
    """

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
    def __init__ (self, url=None, flags=READ, session=None, 
                  _adaptor=None, _adaptor_state={}, _ttype=None) : 
        """
        __init__(url, flags=READ, session)

        Construct a new directory object

        :param url:     Url of the (remote) directory
        :type  url:     :class:`saga.Url` 

        :param flags:   :ref:`filesystemflags`
        :param session: :class:`saga.Session`
        
        The specified directory is expected to exist -- otherwise
        a DoesNotExist exception is raised.  Also, the URL must point to
        a directory (not to a file), otherwise a BadParameter exception is
        raised.

        Example::

            # open some directory
            dir = saga.filesystem.Directory("sftp://localhost/tmp/")

            # and list its contents
            files = dir.list ()

        """

        # param checks
        if  not flags : flags = 0
        url = surl.Url (url)

        if  not url.schema :
            url.schema = 'file'

        if  not url.host :
            url.host = 'localhost'

        self._nsdirec = super  (Directory, self)
        self._nsdirec.__init__ (url, flags, session, 
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
    def create (cls, url=None, flags=READ, session=None, ttype=None) :
        """
        url:       saga.Url
        flags:     saga.replica.flags enum
        session:   saga.Session
        ttype:     saga.task.type enum
        ret:       saga.Task
        """

        if  not flags : flags = 0
        _nsdir = super (Directory, cls)
        return _nsdir.create (url, flags, session, ttype=ttype)

    # ----------------------------------------------------------------
    #
    # filesystem directory methods
    #
    @rus.takes   ('Directory', 
                  (surl.Url, basestring),
                  rus.optional (int, rus.nothing),
                  rus.optional (rus.one_of (SYNC, ASYNC, TASK)))
    @rus.returns (('File', st.Task))
    def open (self, path, flags=READ, ttype=None) :
        """
        open(path, flags=READ)

        Open a file in the directory instance namespace. Returns
        a new file object.

        :param path:     The name/path of the file to open
        :type path:      str()
        :param flags:    :ref:`filesystemflags`
        """
        if  not flags : flags = 0
        
        return self._adaptor.open (path, flags, ttype=ttype)


    # --------------------------------------------------------------------------
    #
    @rus.takes   ('Directory', 
                  (surl.Url, basestring),
                  rus.optional (int, rus.nothing),
                  rus.optional (rus.one_of (SYNC, ASYNC, TASK)))
    @rus.returns (('Directory', st.Task))
    def open_dir (self, path, flags=READ, ttype=None) :
        """
        open_dir(path, flags=READ)

        Open a directory in the directory instance namespace. Returns 
        a new directory object.

        :param path:     The name/path of the directory to open
        :type path:      str()
        :param flags:    :ref:`filesystemflags`        

        Example::

            # create a subdir 'data' in /tmp
            dir = saga.namespace.Directory("sftp://localhost/tmp/")
            data = dir.open_dir ('data/', saga.namespace.Create)
        """
        if  not flags : flags = 0

        return self._adaptor.open_dir (path, flags, ttype=ttype)


    # --------------------------------------------------------------------------
    #
    @rus.takes   ('Directory', 
                  rus.optional ((surl.Url, basestring)),
                  rus.optional (rus.one_of (SYNC, ASYNC, TASK)))
    @rus.returns ((int, st.Task))
    def get_size (self, path=None, ttype=None) :
        """
        get_size(path=None)

        Return the size of the directory itself or the entry pointed to by `path`. 
        
        :param path:     (Optional) name/path of an entry
        :type path:      str()

        Returns the size of a file or directory (in bytes)

        Example::

            # inspect a file for its size
            dir  = saga.filesystem.Directory("sftp://localhost/tmp/")
            size = dir.get_size ('data/data.bin')
            print size
        """
        if path   :  return self._adaptor.get_size      (path, ttype=ttype)
        else      :  return self._adaptor.get_size_self (      ttype=ttype)


    # --------------------------------------------------------------------------
    #
    @rus.takes   ('Directory', 
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
    @rus.takes   ('Directory', 
                  rus.optional ((surl.Url, basestring)),
                  rus.optional (rus.one_of (SYNC, ASYNC, TASK)))
    @rus.returns ((bool, st.Task))
    def is_file (self, path=None, ttype=None) :
        """
        is_file(path=None)

        Returns `True` if entry points to a file, `False` otherwise. If `path`
        is not none, the entry pointed to by `path` is inspected instead of the
        directory object itself. 

        :param path:     (Optional) name/path of an entry
        :type path:      str()
        """
        if path   :  return self._adaptor.is_file      (path, ttype=ttype)
        else      :  return self._adaptor.is_file_self (      ttype=ttype)


    size  = property (get_size)  # int

    
# ------------------------------------------------------------------------------

