# Copyright 2009 Google Inc.
#
#   Licensed under the Apache License, Version 2.0 (the "License");
#   you may not use this file except in compliance with the License.
#   You may obtain a copy of the License at
#
#       http://www.apache.org/licenses/LICENSE-2.0
#
#   Unless required by applicable law or agreed to in writing, software
#   distributed under the License is distributed on an "AS IS" BASIS,
#   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#   See the License for the specific language governing permissions and
#   limitations under the License.
#
# Licensed to the PSF under a Contributor Agreement.
#
# Author: Gregory P. Smith <greg@krypto.org>

"""
This module implements a pthread_atfork() work-a-like mechanism for all
fork() calls made from the Python os module.  Any time a fork() is called
from Python a set of unique callbacks will be made in each of the following
three states:
    Preparing to fork - Immediately before the fork call is made.
    In the parent after fork - Immediately after the fork (regardless of
                               success or failure) in the parent process.
    In the child after fork - Immediately after the fork in the child process.

To use this module, first import it early on your programs initialization:

    import atfork
    atfork.monkeypatch_os_fork_functions()

That will stub out os.fork and os.forkpty with wrapped versions implementing
the enhanced behavior.

Next, register your atfork actions by calling atfork.atfork:

    atfork.atfork(prepare=my_lock.acquire,
                  parent=my_lock.release,
                  child=my_lock.release)

No API to unregister an atfork call is provided.  If you are concerned
about resource usage by references your callable holds, consider using
weakref's within your callable.
"""

import os
import sys
import threading
import traceback


def monkeypatch_os_fork_functions():
    """
    Replace os.fork* with wrappers that use ForkSafeLock to acquire
    all locks before forking and release them afterwards.
    """
    builtin_function = type(''.join)
    if hasattr(os, 'fork') and isinstance(os.fork, builtin_function):
        global _orig_os_fork
        _orig_os_fork = os.fork
        os.fork = os_fork_wrapper
    if hasattr(os, 'forkpty') and isinstance(os.forkpty, builtin_function):
        global _orig_os_forkpty
        _orig_os_forkpty = os.forkpty
        os.forkpty = os_forkpty_wrapper


# This lock protects all of the lists below.
_fork_lock = threading.Lock()
_prepare_call_list = []
_prepare_call_exceptions = []
_parent_call_list = []
_child_call_list = []


def atfork(prepare=None, parent=None, child=None):
    """A Python work-a-like of pthread_atfork.
    
    Any time a fork() is called from Python, all 'prepare' callables will
    be called in the order they were registered using this function.

    After the fork (successful or not), all 'parent' callables will be called in
    the parent process.  If the fork succeeded, all 'child' callables will be
    called in the child process.

    No exceptions may be raised from any of the registered callables.  If so
    they will be printed to sys.stderr after the fork call once it is safe
    to do so.
    """
    assert not prepare or callable(prepare)
    assert not parent or callable(parent)
    assert not child or callable(child)
    _fork_lock.acquire()
    try:
        if prepare:
            _prepare_call_list.append(prepare)
        if parent:
            _parent_call_list.append(parent)
        if child:
            _child_call_list.append(child)
    finally:
        _fork_lock.release()


def _call_atfork_list(call_list):
    """
    Given a list of callables in call_list, call them all in order and save
    and return a list of sys.exc_info() tuples for each exception raised.
    """
    exception_list = []
    for func in call_list:
        try:
            func()
        except:
            exception_list.append(sys.exc_info())
    return exception_list


def prepare_to_fork_acquire():
    """Acquire our lock and call all prepare callables."""
    _fork_lock.acquire()
    _prepare_call_exceptions.extend(_call_atfork_list(_prepare_call_list))


def parent_after_fork_release():
    """
    Call all parent after fork callables, release the lock and print
    all prepare and parent callback exceptions.
    """
    prepare_exceptions = list(_prepare_call_exceptions)
    del _prepare_call_exceptions[:]
    exceptions = _call_atfork_list(_parent_call_list)
    _fork_lock.release()
   #_print_exception_list(prepare_exceptions, 'before fork')
   #_print_exception_list(exceptions, 'after fork from parent')


def child_after_fork_release():
    """
    Call all child after fork callables, release lock and print all
    all child callback exceptions.
    """
    del _prepare_call_exceptions[:]
    exceptions = _call_atfork_list(_child_call_list)
    _fork_lock.release()
   #_print_exception_list(exceptions, 'after fork from child')


def _print_exception_list(exceptions, message, output_file=None):
    """
    Given a list of sys.exc_info tuples, print them all using the traceback
    module preceeded by a message and separated by a blank line.
    """
    output_file = output_file or sys.stderr
    message = 'Exception %s:\n' % message
    for exc_type, exc_value, exc_traceback in exceptions:
        output_file.write(message)
        traceback.print_exception(exc_type, exc_value, exc_traceback,
                                  file=output_file)
        output_file.write('\n')


def os_fork_wrapper():
    """Wraps os.fork() to run atfork handlers."""
    pid = None
    prepare_to_fork_acquire()
    try:
        pid = _orig_os_fork()
    finally:
        if pid == 0:
            child_after_fork_release()
        else:
            # We call this regardless of fork success in order for
            # the program to be in a sane state afterwards.
            parent_after_fork_release()
    return pid


def os_forkpty_wrapper():
    """Wraps os.forkpty() to run atfork handlers."""
    pid = None
    prepare_to_fork_acquire()
    try:
        pid, fd = _orig_os_forkpty()
    finally:
        if pid == 0:
            child_after_fork_release()
        else:
            parent_after_fork_release()
    return pid, fd


# TODO: Also replace os.fork1() on Solaris.
