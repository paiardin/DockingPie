
import sys
import os
import time
import atexit
import signal
import multiprocessing

# from 
# http://stackoverflow.com/questions/1417631/python-code-to-daemonize-a-process


class Daemon (object) :
    """
    A generic daemon class.

    Usage: subclass the Daemon class and override the run() method
    """

    # --------------------------------------------------------------------------
    #
    def __init__ (self, pidfile=None, 
                  stdin='/dev/null', 
                  stdout='/dev/null', 
                  stderr='/dev/null'):

        self.stdin   = stdin
        self.stdout  = stdout
        self.stderr  = stderr
        self.pidfile = pidfile
        self.pid     = None
        self.queue   = multiprocessing.Queue ()


    # --------------------------------------------------------------------------
    #
    def daemonize (self, debug=False) :
        """
        do the UNIX double-fork magic, see Stevens' "Advanced 
        Programming in the UNIX Environment" for details (ISBN 0201563177)
        http://www.erlenstar.demon.co.uk/unix/faq_2.html#SEC16
        """

        try :

            # first fork
            try : 
                f1_pid = os.fork () 
                if  f1_pid > 0:
                    # wait for daemon pid from second parent
                    self.pid = self.queue.get ()

                    # we are done...
                    return False

            except OSError as e : 
                raise RuntimeError ("fork #1 failed: %d (%s)\n" % (e.errno, e.strerror))


            # decouple from parent environment
            os.chdir  ("/") 
          # os.setsid () 
            os.umask  (0) 

            # second fork
            try : 
                f2_pid = os.fork () 
                
                if  f2_pid > 0 :
                    # communicate pid to first parent
                    self.queue.put (f2_pid)

                    # exit from second parent
                    sys.exit(0) 

            except OSError as e : 
                # no use rasing exceptions at this point
                sys.stderr.write ("fork #2 failed: %d (%s)\n" % (e.errno, e.strerror))
                sys.exit (1) 

            if  not debug :
                # redirect standard file descriptors
                sys.stdout.flush ()
                sys.stderr.flush ()
                
                si = file (self.stdin,  'r' )
                so = file (self.stdout, 'a+')
                se = file (self.stderr, 'a+', 0)
                
                os.dup2 (si.fileno (), sys.stdin.fileno  ())
                os.dup2 (so.fileno (), sys.stdout.fileno ())
                os.dup2 (se.fileno (), sys.stderr.fileno ())

            # write pidfile
            if  self.pidfile :
                atexit.register (self.delpid)
                pid = str (os.getpid ())
                file (self.pidfile,'w+').write ("%s\n" % pid)

            return True

        except Exception as e :
            raise


    # --------------------------------------------------------------------------
    #
    def delpid (self) :

        if  self.pidfile :
            os.remove(self.pidfile)


    # --------------------------------------------------------------------------
    #
    def start (self, debug=False) :

      # if  debug :
      #     return self.run ()

        if  self.pidfile :
            
            # Check for a pidfile to see if the daemon already runs
            try :
                pf  = file (self.pidfile,'r')
                pid = int  (pf.read ().strip ())
                pf.close ()
            except IOError :
                pid = None

            if pid :
                raise RuntimeError ("pidfile %s exist - daemon running?\n" \
                                   % self.pidfile)

        # Start the daemon, and in the demon, run the workload
        if self.daemonize (debug=debug) :
            try :
                self.run ()
            except Exception as e :
                # FIXME: we need a logfile :P
                pass


    # --------------------------------------------------------------------------
    #
    def stop (self) :
        """
        Stop the daemon
        """
        # Get the pid from the pidfile
        if  self.pidfile :
            try :
                pf  = file (self.pidfile,'r')
                pid = int  (pf.read ().strip ())
                pf.close ()
            except IOError :
                pid = None

        else :
            pid = self.pid


        if  not pid:
            raise RuntimeError ("no pidfile %s / pid not known.  Daemon not running?\n" \
                               % self.pidfile)

        # Try killing the daemon process    
        try:
            while 1:
                os.kill (pid, signal.SIGTERM)
                time.sleep (0.1)

        except OSError as err :
            err = str(err)

            if  err.find ("No such process") :
                if  os.path.exists (self.pidfile) :
                    os.remove (self.pidfile)

            else :
                raise

    # --------------------------------------------------------------------------
    #
    def restart (self) :
        """
        Restart the daemon
        """
        self.stop  ()
        self.start ()


    # --------------------------------------------------------------------------
    #
    def run (self) :
        """
        You should override this method when you subclass Daemon. 
        It will be called after the process has been
        daemonized by start() or restart().
        """
        raise RuntimeError ("deamon workload undefined")

