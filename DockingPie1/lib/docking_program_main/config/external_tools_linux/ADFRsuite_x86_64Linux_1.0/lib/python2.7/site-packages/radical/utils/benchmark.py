
__author__    = "Radical.Utils Development Team (Andre Merzky)"
__copyright__ = "Copyright 2013, RADICAL@Rutgers"
__license__   = "MIT"


import os
import sys
import math
import time
import socket
import threading
import traceback

import threads as rut
import testing as rutest


# --------------------------------------------------------------------
#
class Benchmark (object) :

    # ----------------------------------------------------------------
    #
    def __init__ (self, config, name, func_pre, func_core, func_post) :

        try:
            import numpy
        except ImportError, e:
            raise RuntimeError ("Benchmark depends on numpy which is not installed")
    
        self.pre    = func_pre
        self.core   = func_core
        self.post   = func_post
        self.name   = name.replace (' ', '_')
        self.lock   = rut.RLock ()
        self.events = dict()
    
        if  isinstance (config, basestring) :
            app_cfg = rutest.TestConfig (config)
        elif isinstance (config, dict) :
            app_cfg = config
        else :
            print 'warning: no valid benchmarl config (neither filename nor dict) -- ignore'
            app_cfg = dict()

        if  'benchmarks' in app_cfg :
            bench_cfg = app_cfg['benchmarks']
        else :
            bench_cfg = dict ()

    
        # RADICAL_BENCHMARK_ environments will overwrite config settings
        if  'RADICAL_BENCHMARK_TAGS' in os.environ :
            bench_cfg['tags'] = os.environ['RADICAL_BENCHMARK_TAGS'].split (' ')
    
        if  'RADICAL_BENCHMARK_CONCURRENCY' in os.environ :
            bench_cfg['concurrency'] = os.environ['RADICAL_BENCHMARK_CONCURRENCY']
    
        if  'RADICAL_BENCHMARK_ITERATIONS' in os.environ :
            bench_cfg['iterations'] = os.environ['RADICAL_BENCHMARK_ITERATIONS']
    
        
        # check benchmark settings for completeness, set some defaults
        if  not 'concurrency' in bench_cfg : 
            bench_cfg['concurrency'] = 1
    
        if  not 'iterations'  in bench_cfg : 
            bench_cfg['iterations'] = 1
    
        if  not 'tags' in bench_cfg : 
            bench_cfg['tags'] = list()
    
        self.app_cfg   = app_cfg
        self.bench_cfg = bench_cfg
    

    # --------------------------------------------------------------------------
    #
    def _thread (self, tid) :
    
        try :
            self.pre (tid, self.app_cfg, self.bench_cfg)

            sys.stdout.write ('-')
            sys.stdout.flush ()
    
            self.events[tid]['event_1'].set  ()  # signal we are done        
            self.events[tid]['event_2'].wait ()  # wait 'til others are done 
    
            iterations = int(self.bench_cfg['iterations']) / int(self.bench_cfg['concurrency'])
    
            # poor-mans ceil()
            if (iterations * int(self.bench_cfg['concurrency'])) < int(self.bench_cfg['iterations']) :
                iterations += 1
    
            for i in range (0, iterations+1) :
                core_timer = self.core (tid, i, self.app_cfg, self.bench_cfg)
                self._tic (tid, core_timer)
    
    
            self.events[tid]['event_3'].set ()   # signal we are done        
            self.events[tid]['event_4'].wait ()  # wait 'til others are done 
    
    
            self.post (tid, self.app_cfg, self.bench_cfg)
            sys.stdout.write ('=')
            sys.stdout.flush ()
    
            self.events[tid]['event_5'].set ()   # signal we are done        

    
        except Exception as e :
    
            sys.stdout.write ("exception in benchmark thread: %s\n\n" % str(e))
            sys.stdout.write (traceback.format_exc())
            sys.stdout.write ("\n")
            sys.stdout.flush ()
    
            # Oops, we are screwed.  Tell main thread that we are done for, and
            # bye-bye...
            self.events[tid]['event_1'].set  ()
            self.events[tid]['event_3'].set  ()
            self.events[tid]['event_5'].set  ()

            raise

        # finish thread
        sys.exit (0)
    
    
    # --------------------------------------------------------------------------
    #
    def run (self) :
        """
        - create 'concurrency' number of threads
        - per thread call pre()
        - sync threads, start timer
        - per thread call core() 'iteration' number of times', tic()
        - stop timer
        - per thread, call post, close threads
        - eval once
        """
    
        threads     = []
        concurrency = int(self.bench_cfg['concurrency'])
    
        self._start ()
    
        for tid in range (0, concurrency) :
    
            self.events[tid] = {}
            self.events[tid]['event_1'] = rut.Event ()
            self.events[tid]['event_2'] = rut.Event ()
            self.events[tid]['event_3'] = rut.Event ()
            self.events[tid]['event_4'] = rut.Event ()
            self.events[tid]['event_5'] = rut.Event ()
            self.start [tid] = time.time ()
            self.times [tid] = list()
    
            t = rut.Thread (self._thread, tid)
            threads.append (t)
    
    
        for t in threads :
            t.start ()
    
        
        # wait for all threads to start up and initialize
        self.t_init = time.time ()
        rut.lout ("\n> " + "="*concurrency)
        rut.lout ("\n> ")
        for tid in range (0, concurrency) :
            self.events[tid]['event_1'].wait ()
    
        # start workload in all threads
        self.t_start = time.time ()
        for tid in range (0, concurrency) :
            self.events[tid]['event_2'].set ()
    
        # wait for all threads to finish core test
        for tid in range (0, concurrency) :
            self.events[tid]['event_3'].wait ()
        self.t_stop = time.time ()
    
        # start shut down
        rut.lout ("\n< " + "-"*concurrency)
        rut.lout ("\n< ")
        for tid in range (0, concurrency) :
            self.events[tid]['event_4'].set ()
    
        # wait for all threads to finish shut down
        for tid in range (0, concurrency) :
            self.events[tid]['event_5'].wait ()
    
    
    
    # ----------------------------------------------------------------
    #
    def _start (self) :
    
        self.start = dict()
        self.times = dict()
        self.idx   = 0
    
        rut.lout ("\n")
        rut.lout ("benchmark   : %s (%s)\n" % (self.name, self.bench_cfg['tags']))
        rut.lout ("concurrency : %s\n"      %  self.bench_cfg['concurrency'])
        rut.lout ("iterations  : %s\n"      %  self.bench_cfg['iterations'])
    
        sys.stdout.flush ()
    
    
    # ----------------------------------------------------------------
    #
    def _tic (self, tid='master_tid', core_timer=None) :
        # if a core_timer is given, we interpret that as time value instead of
        # using 'now-start'.  This assumes that 'core' has its own way of
        # measuring the time, possibly to correct for some systematic error...
    
        with self.lock :

            import numpy

            now   = time.time ()

            if  core_timer :
                timer = core_timer
            else :
                timer = now - self.start[tid]

            self.times[tid].append (timer)
            self.start[tid] = now

            if len(self.times[tid][1:]) :
                numpy_times = numpy.array (self.times[tid][1:])
                vmean = numpy_times.mean ()
            else :
                vmean = timer
    
            if   timer  <  0.75 * vmean : marker = '='
            if   timer  <  0.90 * vmean : marker = '~'
            elif timer  <  0.95 * vmean : marker = '_'
            elif timer  <  1.05 * vmean : marker = '.'
            elif timer  <  1.10 * vmean : marker = '-'
            elif timer  <  1.25 * vmean : marker = '+'
            else                        : marker = '*'
    
    
            if       not ( (self.idx)        ) : sys.stdout.write ('\n* ')
            else :
                if   not ( (self.idx) % 1000 ) : sys.stdout.write (" %7d\n\n# " % self.idx)
                elif not ( (self.idx) %  100 ) : sys.stdout.write (" %7d\n| "   % self.idx)
                elif not ( (self.idx) %   10 ) : sys.stdout.write (' ')
            if  True                                    : sys.stdout.write (marker)
            
            sys.stdout.flush ()
    
            self.idx += 1
    
    # ----------------------------------------------------------------
    #
    def eval (self, error=None) :
    
        import numpy

        times = list()
    
        for tid in self.times :
            times += self.times[tid][1:]

      # import pprint
      # pprint.pprint (times)
    
        if  len(times) < 1 :
            raise ValueError ("min 1 timing value required for benchmark evaluation (%d)" % len(times))
    
        concurrency = int(self.bench_cfg['concurrency'])
        tags        = ' '.join (self.bench_cfg['tags'])
    
        out = "\n"
        top = ""
        tab = ""
        num = ""
    
        out += "Results :\n"

        numpy_times = numpy.array (times)
    
        vtot  = self.t_stop  - self.t_start
        vini  = self.t_start - self.t_init
        vn    = len (times)
        vsum  = sum (times)
        vmin  = min (times)
        vmax  = max (times)

        if len (times) :
            vmean = numpy_times.mean ()
            vsdev = numpy_times.std  ()
        else :
            vmean = 0.0
            vsdev = 0.0
        vrate = vn / vtot
    
        bdat  = "benchmark.%s.dat" % (self.name)
    
        out += "  name    : %s\n"                               % (self.name)
        out += "  tags    : %s\n"                               % (tags)
        out += "  threads : %8d        iters   : %8d\n"        % (concurrency, vn)
        out += "  init    : %8.2fs       min     : %8.2fs\n"    % (vini,        vmin )
        out += "  total   : %8.2fs       max     : %8.2fs\n"    % (vtot,        vmax )
        out += "  rate    : %8.2f/s      mean    : %8.2fs\n"    % (vrate,       vmean)
        out += "                            sdev    : %8.2fs\n" % (             vsdev)
    
        num = "# %7s  %7s  %7s  %7s  %7s  %7s  %8s  %8s  %9s  %10s  %10s" \
            % (1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11)
        top = "# %7s  %7s  %7s  %7s  %7s  %7s  %8s  %8s  %9s  %10s  %10s" \
            % ('n', 'threads', 'init', 'tot', 'min',  'max', 'mean', \
               'std-dev', 'rate', 'name', 'tags')
    
        tab = "  "      \
              "%7d  "   \
              "%7d  "   \
              "%7.2f  " \
              "%7.2f  " \
              "%7.2f  " \
              "%7.2f  " \
              "%8.3f  " \
              "%8.3f  " \
              "%9.3f  " \
              "%10s  "  \
              "%10s  "  \
            % (vn, 
               concurrency, 
               vini,
               vtot,   
               vmin,  
               vmax, 
               vmean, 
               vsdev, 
               vrate, 
               self.name, 
               tags) 
    
        rut.lout ("\n%s" % out)
    
        create_top = True
        try :
            statinfo = os.stat (bdat)
            if  statinfo.st_size > 0 :
                create_top = False
        except :
            pass
    
        f = open (bdat, "a+")
    
        if  create_top :
            f.write ("%s\n" % num)
            f.write ("%s\n" % top)
        f.write ("%s\n" % tab)
    
    
# --------------------------------------------------------------------

