
__author__    = "Andre Merzky"
__copyright__ = "Copyright 2012-2013, The SAGA Project"
__license__   = "MIT"


import re
import os
import sys
import pty
import tty
import time
import errno
import shlex
import select
import signal
import termios

import radical.utils         as ru
import radical.utils.logger  as rul

import saga.exceptions       as se

import pty_exceptions        as ptye

# --------------------------------------------------------------------
#
_CHUNKSIZE = 1024*1024  # default size of each read
_POLLDELAY = 0.01       # seconds in between read attempts
_DEBUG_MAX = 600


# --------------------------------------------------------------------
#
class PTYProcess (object) :
    """
    This class spawns a process, providing that child with pty I/O channels --
    it will maintain stdin, stdout and stderr channels to the child.  All
    write-like operations operate on the stdin, all read-like operations operate
    on the stdout stream.  Data from the stderr stream are at this point
    redirected to the stdout channel.

    Example::

        # run an interactive client process
        pty = PTYProcess ("/usr/bin/ssh -t localhost")

        # check client's I/O for one of the following patterns (prompts).
        # Then search again.
        n, match = pty.find (['password\s*:\s*$',
                              'want to continue connecting.*\(yes/no\)\s*$',
                              '[\$#>]\s*$'])

        while True :

            if n == 0 :
                # found password prompt - tell the secret
                pty.write ("secret\\n")
                n, _ = pty.find (['password\s*:\s*$',
                                  'want to continue connecting.*\(yes/no\)\s*$',
                                  '[\$#>]\s*$'])
            elif n == 1 :
                # found request to accept host key - sure we do... (who checks
                # those keys anyways...?).  Then search again.
                pty.write ("yes\\n")
                n, _ = pty.find (['password\s*:\s*$',
                                  'want to continue connecting.*\(yes/no\)\s*$',
                                  '[\$#>]\s*$'])
            elif n == 2 :
                # found shell prompt!  Wohoo!
                break


        while True :
            # go full Dornroeschen (Sleeping Beauty)...
            pty.alive (recover=True) or break      # check / restart process
            pty.find  (['[\$#>]\s*$'])             # find shell prompt
            pty.write ("/bin/sleep "100 years"\\n") # sleep!  SLEEEP!

        # something bad happened
        print pty.autopsy ()

    """

    # ----------------------------------------------------------------
    #
    def __init__ (self, command, logger=None) :
        """
        The class constructor, which runs (execvpe) command in a separately
        forked process.  The bew process will inherit the environment of the
        application process.

        :type  command: string or list of strings
        :param command: The given command is what is run as a child, and
        fed/drained via pty pipes.  If given as string, command is split into an
        array of strings, using :func:`shlex.split`.

        :type  logger:  :class:`radical.utils.logger.Logger` instance
        :param logger:  logger stream to send status messages to.
        """

        self.logger = logger
        if  not  self.logger : self.logger = ru.get_logger ('radical.saga.pty')
        self.logger.debug ("PTYProcess init %s" % self)


        if isinstance (command, basestring) :
            command = shlex.split (command)

        if not isinstance (command, list) :
            raise se.BadParameter ("PTYProcess expects string or list command")

        if len(command) < 1 :
            raise se.BadParameter ("PTYProcess expects non-empty command")

        self.rlock   = ru.RLock ("pty process %s" % command)

        self.command = command # list of strings too run()


        self.cache   = ""      # data cache
        self.tail    = ""      # tail of data data cache for error messages
        self.child   = None    # the process as created by subprocess.Popen
        self.ptyio   = None    # the process' io channel, from pty.fork()

        self.exit_code        = None  # child died with code (may be revived)
        self.exit_signal      = None  # child kill by signal (may be revived)

        self.recover_max      = 3  # TODO: make configure option.  This does not
        self.recover_attempts = 0  # apply for recovers triggered by gc_timeout!


        try :
            self.initialize ()

        except Exception as e :
            raise ptye.translate_exception (e, "pty or process creation failed")

    # --------------------------------------------------------------------
    #
    def __del__ (self) :
        """
        Need to free pty's on destruction, otherwise we might ran out of
        them (see cat /proc/sys/kernel/pty/max)
        """

        self.logger.debug ("PTYProcess del  %s" % self)
        with self.rlock :

            try :
                self.finalize ()
            except :
                pass


    # ----------------------------------------------------------------------
    #
    def flush (self) :

        # clean the *read* data cache -- only use when you know what you are
        # doing!

        try :
             if not self.parent_out :
                 # nothing to read, yet, so there is nothing to cache, so tehre is
                 # nothing to flush...
                 return

             self.logger.debug ("flush: [%5d] [     ] (flush pty read cache)"  \
                     % (self.parent_out))

             # lets see if there are still things to read
             while True :
                 tmp = self.read (timeout=0.1)
                 if tmp :
                     self.logger.warn ("flush: [%5d] [%5d] (discard data: '%s')" \
                             % (self.parent_out, len(tmp), tmp))
                     continue
                 break

             if len(self.cache) :
                  self.logger.warn ("flush: [%5d] [%5d] (discard cache: '%s')" \
                          % (self.parent_out, len(self.cache), self.cache))
             self.cache = ""

        except Exception as e :
            self.logger.exception ('flushing failed')
            raise



    # ----------------------------------------------------------------------
    #
    def _hide_data (self, data, nolog=False) :

        if  nolog :
            import re
            return re.sub (r'([^\n])', 'X', data)

        else :
            return data




    # ----------------------------------------------------------------------
    #
    def initialize (self) :

        with self.rlock :

            # already initialized?
            if  self.child :
                self.logger.warn ("initialization race: %s" % ' '.join (self.command))
                return


            self.logger.info ("running: %s" % ' '.join (self.command))

            # create the child
            try :
                self.child, self.child_fd = pty.fork ()
            except Exception as e:
                raise se.NoSuccess ("Could not run (%s): %s" \
                                 % (' '.join (self.command), e))

            if  not self.child :
                # this is the child

                try :
                    # all I/O set up, have a pty (*fingers crossed*), lift-off!
                    os.execvpe (self.command[0], self.command, os.environ)

                except OSError as e:
                    self.logger.error ("Could not execute (%s): %s" \
                                    % (' '.join (self.command), e))
                    sys.exit (-1)

            else :
                # this is the parent
                new = termios.tcgetattr (self.child_fd)
                new[3] = new[3] & ~termios.ECHO

                termios.tcsetattr (self.child_fd, termios.TCSANOW, new)


                self.parent_in  = self.child_fd
                self.parent_out = self.child_fd


    # --------------------------------------------------------------------
    #
    def finalize (self, wstat=None) :
        """ kill the child, close all I/O channels """

        with self.rlock :

            # now we can safely kill the process -- unless some wait did that before
            if  wstat == None :

                if  self.child :
                    # yes, we have something to kill!
                    try :
                        os.kill (self.child, signal.SIGKILL)
                    except OSError :
                        pass

                    # hey, kiddo, how did that go?
                    max_tries = 10
                    tries     =  0
                    while tries < max_tries :
                        try :
                            wpid, wstat = os.waitpid (self.child, os.WNOHANG)

                        except OSError as e :

                            # this should not have failed -- child disappeared?
                            if e.errno == errno.ECHILD :
                                self.exit_code   = None
                                self.exit_signal = None
                                wstat            = None
                                break
                            else :
                                # other errors are bad, but there is not much to
                                # be done at this point
                                self.logger.warning ("ignore waitpid failure on finalize (%s)" % e)
                                break

                        if  wpid :
                            break

                        time.sleep (0.1)
                        tries += 1


            # at this point, we declare the process to be gone for good
            self.child = None

            # lets see if we can perform some post-mortem analysis
            if  wstat != None :

                if  os.WIFEXITED (wstat) :
                    # child died of natural causes - perform autopsy...
                    self.exit_code   = os.WEXITSTATUS (wstat)
                    self.exit_signal = None

                elif os.WIFSIGNALED (wstat) :
                    # murder!! Child got killed by someone!  recover evidence...
                    self.exit_code   = None
                    self.exit_signal = os.WTERMSIG (wstat)


            try :
                if  self.parent_out :
                    os.close (self.parent_out)
                    self.parent_out = None
            except OSError :
                pass

          # try :
          #     if  self.parent_in :
          #         os.close (self.parent_in)
          #         self.parent_in = None
          # except OSError :
          #     pass

          # try :
          #     os.close (self.parent_err)
          # except OSError :
          #     pass



    # --------------------------------------------------------------------
    #
    def wait (self) :
        """
        blocks forever until the child finishes on its own, or is getting
        killed.

        Actully, we might just as well try to figure out what is going on on the
        remote end of things -- so we read the pipe until the child dies...
        """

        output = ""
        # yes, for ever and ever...
        while True :
            try :
                output += self.read ()
            except :
                break

        # yes, for ever and ever...
        while True :

            if not self.child:
                # this was quick ;-)
              # print "child is gone"
                return output

            # we need to lock, as the SIGCHLD will only arrive once
            with self.rlock :
                # hey, kiddo, whats up?
                try :
                    wpid, wstat = os.waitpid (self.child, 0)
                  # print "wait: %s -- %s" % (wpid, wstat)

                except OSError as e :

                    if e.errno == errno.ECHILD :

                        # child disappeared
                        self.exit_code   = None
                        self.exit_signal = None
                        self.finalize ()
                      # print "no such child"
                        return output

                    # no idea what happened -- it is likely bad
                  # print "waitpid failed"
                    raise se.NoSuccess ("waitpid failed on wait")


                # did we get a note about child termination?
                if 0 == wpid :

                    # nope, all is well - carry on
                    continue


                # Yes, we got a note.
                # Well, maybe the child fooled us and is just playing dead?
                if os.WIFSTOPPED   (wstat) or \
                   os.WIFCONTINUED (wstat)    :
                    # we don't care if someone stopped/resumed the child -- that is up
                    # to higher powers.  For our purposes, the child is alive.  Ha!
                    continue


                # not stopped, poor thing... - soooo, what happened??  But hey,
                # either way, its dead -- make sure it stays dead, to avoid
                # zombie apocalypse...
                self.child = None
                self.finalize (wstat=wstat)

              # print "child is done"
                return output


    # --------------------------------------------------------------------
    #
    def alive (self, recover=False) :
        """
        try to determine if the child process is still active.  If not, mark
        the child as dead and close all IO descriptors etc ("func:`finalize`).

        If `recover` is `True` and the child is indeed dead, we attempt to
        re-initialize it (:func:`initialize`).  We only do that for so many
        times (`self.recover_max`) before giving up -- at that point it seems
        likely that the child exits due to a re-occurring operations condition.

        Note that upstream consumers of the :class:`PTYProcess` should be
        careful to only use `recover=True` when they can indeed handle
        a disconnected/reconnected client at that point, i.e. if there are no
        assumptions on persistent state beyond those in control of the upstream
        consumers themselves.
        """

        with self.rlock :

            # do we have a child which we can check?
            if  self.child :

                wstat = None

                while True :
                  # print 'waitpid %s' % self.child

                    # hey, kiddo, whats up?
                    try :
                        wpid, wstat = os.waitpid (self.child, os.WNOHANG)
                      # print 'waitpid %s : %s - %s' % (self.child, wpid, wstat)

                    except OSError as e :

                        if e.errno == errno.ECHILD :
                            # child disappeared, go to zombie cleanup routine
                            break

                        raise ("waitpid failed on wait (%s)" % e)

                    # did we get a note about child termination?
                    if 0 == wpid :
                      # print 'waitpid %s : %s - %s -- none' % (self.child, wpid, wstat)
                        # nope, all is well - carry on
                        return True


                    # Yes, we got a note.
                    # Well, maybe the child fooled us and is just playing dead?
                    if os.WIFSTOPPED   (wstat) or \
                       os.WIFCONTINUED (wstat)    :
                      # print 'waitpid %s : %s - %s -- stop/cont' % (self.child, wpid, wstat)
                        # we don't care if someone stopped/resumed the child -- that is up
                        # to higher powers.  For our purposes, the child is alive.  Ha!
                        continue

                    break

                # so its dead -- make sure it stays dead, to avoid zombie
                # apocalypse...
              # print "he's dead, honeybunny, jim is dead..."
                self.child = None
                self.finalize (wstat=wstat)


            # check if we can attempt a post-mortem revival though
            if  not recover :
              # print 'not alive, not recover'
                # nope, we are on holy ground - revival not allowed.
                return False

            # we are allowed to revive!  So can we try one more time...  pleeeease??
            # (for cats, allow up to 9 attempts; for Buddhists, always allow to
            # reincarnate, etc.)
            if self.recover_attempts >= self.recover_max :
                # nope, its gone for good - just report the sad news
              # print 'not alive, no recover anymore'
                return False

            # MEDIIIIC!!!!
            self.recover_attempts += 1
            self.initialize ()

            # well, now we don't trust the child anymore, of course!  So we check
            # again.  Yes, this is recursive -- but note that recover_attempts get
            # incremented on every iteration, and this will eventually lead to
            # call termination (tm).
          # print 'alive, or not alive?  Check again!'
            return self.alive (recover=True)



    # --------------------------------------------------------------------
    #
    def autopsy (self) :
        """
        return diagnostics information string for dead child processes
        """

        with self.rlock :

            if  self.child :
                # Boooh!
                return "false alarm, process %s is alive!" % self.child

            # try a last read to grab whatever we can get (from cache)
            data = ''
            try :
                data  = self.tail
                data += self.read (size=0, timeout=-1)
            except :
                pass

            ret  = ""
            ret += "  exit code  : %s\n" % self.exit_code
            ret += "  exit signal: %s\n" % self.exit_signal
            ret += "  last output: %s\n" % data

            return ret


    # --------------------------------------------------------------------
    #
    def read (self, size=0, timeout=0, _force=False) :
        """
        read some data from the child.  By default, the method reads whatever is
        available on the next read, up to _CHUNKSIZE, but other read sizes can
        be specified.

        The method will return whatever data it has at timeout::

          timeout == 0 : return the content of the first successful read, with
                         whatever data up to 'size' have been found.
          timeout <  0 : return after first read attempt, even if no data have
                         been available.

        If no data are found, the method returns an empty string (not None).

        This method will not fill the cache, but will just read whatever data it
        needs (FIXME).

        Note: the returned lines do *not* get '\\\\r' stripped.
        """

        with self.rlock :

            found_eof = False

            try:
                # start the timeout timer right now.  Note that even if timeout is
                # short, and child.poll is slow, we will nevertheless attempt at least
                # one read...
                start = time.time ()
                ret   = ""

                # read until we have enough data, or hit timeout ceiling...
                while True :

                    # first, lets see if we still have data in the cache we can return
                    if len (self.cache) :

                        if  not size :
                            ret        = self.cache
                            self.cache = ""
                            self.tail += ret
                            self.tail  = self.tail[-256:]
                            return ret

                        # we don't even need all of the cache
                        elif size <= len (self.cache) :
                            ret        = self.cache[:size]
                            self.cache = self.cache[size:]
                            self.tail += ret
                            self.tail  = self.tail[-256:]
                            return ret

                    # otherwise we need to read some more data, right?
                    # idle wait 'til the next data chunk arrives, or 'til _POLLDELAY
                    rlist, _, _ = select.select ([self.parent_out], [], [], _POLLDELAY)

                    # got some data?
                    for f in rlist:
                        # read whatever we still need

                        readsize = _CHUNKSIZE
                        if  size:
                            readsize = size-len(ret)

                        buf  = os.read (f, _CHUNKSIZE)

                        if  len(buf) == 0 and sys.platform == 'darwin' :
                            self.logger.debug ("read : MacOS EOF")
                            self.finalize ()
                            found_eof = True
                            raise se.NoSuccess ("unexpected EOF (%s)" % self.tail)


                        self.cache += buf.replace ('\r', '')
                        log         = buf.replace ('\r', '')
                        log         = log.replace ('\n', '\\n')
                      # print "buf: --%s--" % buf
                      # print "log: --%s--" % log
                        if  len(log) > _DEBUG_MAX :
                            self.logger.debug ("read : [%5d] [%5d] (%s ... %s)" \
                                            % (f, len(log), log[:30], log[-30:]))
                        else :
                            self.logger.debug ("read : [%5d] [%5d] (%s)" \
                                            % (f, len(log), log))
                          # for c in log :
                          #     print '%s' % c


                    # lets see if we still got any data in the cache we can return
                    if len (self.cache) :

                        if  not size :
                            ret        = self.cache
                            self.cache = ""
                            self.tail += ret
                            self.tail  = self.tail[-256:]
                            return ret

                        # we don't even need all of the cache
                        elif size <= len (self.cache) :
                            ret        = self.cache[:size]
                            self.cache = self.cache[size:]
                            self.tail += ret
                            self.tail  = self.tail[-256:]
                            return ret

                    # at this point, we do not have sufficient data -- only
                    # return on timeout

                    if  timeout == 0 :
                        # only return if we have data
                        if len (self.cache) :
                            ret        = self.cache
                            self.cache = ""
                            self.tail += ret
                            self.tail  = self.tail[-256:]
                            return ret

                    elif timeout < 0 :
                        # return of we have data or not
                        ret        = self.cache
                        self.cache = ""
                        self.tail += ret
                        self.tail  = self.tail[-256:]
                        return ret

                    else : # timeout > 0
                        # return if timeout is reached
                        now = time.time ()
                        if (now-start) > timeout :
                            ret        = self.cache
                            self.cache = ""
                            self.tail += ret
                            self.tail  = self.tail[-256:]
                            return ret


            except Exception as e :

                if  found_eof :
                    raise e

                raise se.NoSuccess ("read from process failed '%s' : (%s)" \
                                 % (e, self.tail))


    # ----------------------------------------------------------------
    #
    def find (self, patterns, timeout=0) :
        """
        This methods reads bytes from the child process until a string matching
        any of the given patterns is found.  If that is found, all read data are
        returned as a string, up to (and including) the match.  Note that
        pattern can match an empty string, and the call then will return just
        that, an empty string.  If all patterns end with matching a newline,
        this method is effectively matching lines -- but note that '$' will also
        match the end of the (currently available) data stream.

        The call actually returns a tuple, containing the index of the matching
        pattern, and the string up to the match as described above.

        If no pattern is found before timeout, the call returns (None, None).
        Negative timeouts will block until a match is found

        Note that the pattern are interpreted with the re.M (multi-line) and
        re.S (dot matches all) regex flags.

        Performance: the call is doing repeated string regex searches over
        whatever data it finds.  On complex regexes, and large data, and small
        read buffers, this method can be expensive.

        Note: the returned data get '\\\\r' stripped.

        Note: ansi-escape sequences are also stripped before matching, but are
        kept in the returned data.
        """

        def escape (txt) :
            pat = re.compile(r'\x1b[^m]*m')
            return pat.sub ('', txt)

        _debug = False

        with self.rlock :

            try :
                start = time.time ()                       # startup timestamp
                ret   = []                                 # array of read lines
                patts = []                                 # compiled patterns
                data  = self.cache                         # initial data to check
                self.cache = ""

                if  not data : # empty cache?
                    data = self.read (timeout=_POLLDELAY)

                # pre-compile the given pattern, to speed up matching
                for pattern in patterns :
                    patts.append (re.compile (pattern, re.MULTILINE | re.DOTALL))

                # we wait forever -- there are two ways out though: data matches
                # a pattern, or timeout passes
                while True :

                    # skip non-lines
                    if  not data :
                        data += self.read (timeout=_POLLDELAY)

                    if  _debug : print ">>%s<<" % data

                    escaped = escape (data)
                    if _debug : print 'data    ==%s==' % data
                    if _debug : print 'escaped ==%s==' % escaped

                    # check current data for any matching pattern
                    for n in range (0, len(patts)) :

                        escaped = data
                      # escaped = escape (data)
                      # print '-- 1 --%s--' % data
                      # print '-- 2 --%s--' % escaped

                        match = patts[n].search (escaped)
                        if _debug : print "==%s==" % patterns[n]
                        if _debug : print match

                        if match :
                            # a pattern matched the current data: return a tuple of
                            # pattern index and matching data.  The remainder of the
                            # data is cached.
                            ret  = escaped[0:match.end()]
                            self.cache = escaped[match.end():]

                            if _debug : print "~~match!~~ %s" % escaped[match.start():match.end()]
                            if _debug : print "~~match!~~ %s" % (len(escaped))
                            if _debug : print "~~match!~~ %s" % (str(match.span()))
                            if _debug : print "~~match!~~ %s" % (ret)

                            return (n, ret.replace('\r', ''))

                    # if a timeout is given, and actually passed, return
                    # a non-match and a copy of the data we looked at
                    if timeout == 0 :
                        return (None, str(escaped))

                    if timeout > 0 :
                        now = time.time ()
                        if (now-start) > timeout :
                            self.cache = escaped
                            return (None, str(escaped))

                    # no match yet, still time -- read more data
                    data += self.read (timeout=_POLLDELAY)

            except se.NoSuccess as e :
                raise ptye.translate_exception (e, "(%s)" % data)


    # ----------------------------------------------------------------
    #
    def write (self, data, nolog=False) :
        """
        This method will repeatedly attempt to push the given data into the
        child's stdin pipe, until it succeeds to write all data.
        """

        with self.rlock :

            if not self.alive (recover=False) :
                raise ptye.translate_exception (se.NoSuccess ("cannot write to dead process (%s) [%5d]" \
                                                % (self.cache[-256:], self.parent_in)))

            try :

                log = self._hide_data (data, nolog)
                log =  log.replace ('\n', '\\n')
                log =  log.replace ('\r', '')
                if  len(log) > _DEBUG_MAX :
                    self.logger.debug ("write: [%5d] [%5d] (%s ... %s)" \
                                    % (self.parent_in, len(data), log[:30], log[-30:]))
                else :
                    self.logger.debug ("write: [%5d] [%5d] (%s)" \
                                    % (self.parent_in, len(data), log))

                # attempt to write forever -- until we succeeed
                while data :

                    # check if the pty pipe is ready for data
                    _, wlist, _ = select.select ([], [self.parent_in], [], _POLLDELAY)

                    for f in wlist :

                        # write will report the number of written bytes
                        size = os.write (f, data)

                        # otherwise, truncate by written data, and try again
                        data = data[size:]

                        if data :
                            self.logger.info ("write: [%5d] [%5d]" % (f, size))


            except Exception as e :
                raise ptye.translate_exception (e, "write to process failed (%s)" % e)




