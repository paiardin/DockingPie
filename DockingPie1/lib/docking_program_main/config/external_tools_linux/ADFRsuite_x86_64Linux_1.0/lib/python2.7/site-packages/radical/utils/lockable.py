
__author__    = "Radical.Utils Development Team (Andre Merzky)"
__copyright__ = "Copyright 2013, RADICAL@Rutgers"
__license__   = "MIT"


import threading 


def Lockable (cls) :
    """ 
    This class decorator will add lock/unlock methods to the thusly decorated
    classes, which will be enacted via an also added `threading.RLock` member
    (`self._rlock`)::

        @Lockable
        class A (object) :
        
            def call (self) :
                print 'locked: %s' % self._locked
        
    The class instance can then be used like this::

        a = A    ()
        a.call   ()
        a.lock   ()
        a.call   ()
        a.lock   ()
        a.call   ()
        a.unlock ()
        a.call   ()
        with a :
          a.call ()
        a.call   ()
        a.unlock ()
        a.call   ()
    
    which will result in::

        locked: 0
        locked: 1
        locked: 2
        locked: 1
        locked: 2
        locked: 1
        locked: 0

    The class A can also internally use the lock, and can, for example, use:

        @Lockable
        class A (object) :
            ...
            def work (self) :
                with self :
                    # locked code section
                    ...
    """

    if  hasattr (cls, '__enter__') :
        raise RuntimeError ("Cannot make '%s' lockable -- has __enter__" % cls)

    if  hasattr (cls, '__exit__') :
        raise RuntimeError ("Cannot make '%s' lockable -- has __exit__" % cls)

    if  hasattr (cls, '_rlock') :
        raise RuntimeError ("Cannot make '%s' lockable -- has _rlock" % cls)

    if  hasattr(cls, '_locked') :
        raise RuntimeError ("Cannot make '%s' lockable -- has _locked" % cls)

    if  hasattr(cls, 'locked') :
        raise RuntimeError ("Cannot make '%s' lockable -- has locked" % cls)

    if  hasattr (cls, 'lock') :
        raise RuntimeError ("Cannot make '%s' lockable -- has lock()" % cls)

    if  hasattr (cls, 'unlock') :
        raise RuntimeError ("Cannot make '%s' lockable -- has unlock()" % cls)


    def locked      (self)        : return self._locked
    def locker      (self)        : self._rlock.acquire (); self._locked += 1
    def unlocker    (self, *args) : self._rlock.release (); self._locked -= 1

    cls._rlock    = threading.RLock ()
    cls._locked   = 0
    cls.locked    = locked
    cls.is_locked = locked
    cls.lock      = locker
    cls.unlock    = unlocker
    cls.__enter__ = locker
    cls.__exit__  = unlocker

    return cls


# ------------------------------------------------------------------------------

# @Lockable
# class A (object) :
# 
#     def call (self) :
#         print 'locked 1: %s' % self.locked ()
# 
#         with self :
#             print 'locked 2: %s' % self.locked ()
# 
#         print 'locked 3: %s\n' % self.locked ()
# 
# print
# a = A()
# a.call ()
# a.lock ()
# a.call ()
# a.lock ()
# a.call ()
# a.unlock ()
# a.call ()
# with a :
#   a.call ()
# a.call ()
# a.unlock ()
# a.call ()
# try :
#     a.unlock ()
#     print 'oops\n'
# except :
#     print 'ok\n'
# a.lock ()
# a.call ()
# with a :
#   a.call ()
# a.call ()
# a.unlock ()
# a.call ()
# try :
#     a.unlock ()
#     print 'oops\n'
# except :
#     print 'ok\n'

# ------------------------------------------------------------------------------


