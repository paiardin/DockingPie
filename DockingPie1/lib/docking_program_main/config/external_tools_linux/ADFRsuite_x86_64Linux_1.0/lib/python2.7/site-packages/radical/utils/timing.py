
__author__    = "Radical.Utils Development Team (Andre Merzky)"
__copyright__ = "Copyright 2014, RADICAL@Rutgers"
__license__   = "MIT"


import time
import inspect
import functools


""" 
Provides class decorators to time all public class methods 

The class decorator is provided as decorator method `@timed_class`, and as
decorator class `@TimedClass`.  It thus also serves as documentation on the
different ways to implement class decorators.  
"""


# ------------------------------------------------------------------------------
#
# called on *decorator class* instantiation
def timed_method (func) :
    """ 
    This class decorator will decorate all public class methods with a timing
    function.  That will time the call invocation, and store the respective data
    in a private class dict '__timing'.  Additionally, it will add the class
    method '_timing_last()', which will return the tuple `[method name, method
    timer]` for the last timed invocation of a class method.
    """

    # apply timing decorator to all public methods
    def func_timer(obj, *args, **kwargs):
        
        try  :
            _ = obj.__timing

        except Exception :

            # never seen this one before -- create private timing dict, and add
            # timing_last method
            obj.__timing = dict ()

            def timing_last () :
                last_call = obj.__timing.get ('__last_call', None)
                timer     = obj.__timing.get (last_call, [None])[0]
                return last_call, timer

            def timing_stats () :
                import math
                ret = ""
                for key in sorted (obj.__timing.keys ()) :

                    if  key == '__last_call' :
                        continue

                    tlist = obj.__timing[key]
                    tnum  = len(tlist)
                    tmean = sum(tlist) / tnum
                    tdev  = [x - tmean for x in tlist]
                    tdev2 = [x*x for x in tdev]
                    tstd  = math.sqrt (sum(tdev2) / tnum)

                    ret  += "%-20s : %10.3fs (+/- %5.3fs) [n=%5d]\n" % (key, tmean, tstd, tnum)

                return ret


            def timing_reset () :
                obj.__timing = dict()

            obj._timing_last  = timing_last
            obj._timing_stats = timing_stats
            obj._timing_reset = timing_reset


        # object is set up -- time the call and record timings
        fname = func.__name__
        tdict = obj.__timing
        
        start = time.time()
        ret   = func (obj, *args, **kwargs)
        stop  = time.time()
        timed = stop-start
        
        if not fname in tdict :
            tdict[fname] = list()
        tdict[fname].append (timed)
        tdict['__last_call'] = fname
        
        return ret
        
    return func_timer


