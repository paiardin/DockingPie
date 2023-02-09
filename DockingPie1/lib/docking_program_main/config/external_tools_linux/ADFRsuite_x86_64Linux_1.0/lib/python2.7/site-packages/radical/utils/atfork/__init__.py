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

from .atfork       import monkeypatch_os_fork_functions, atfork
from .stdlib_fixer import fix_logging_module


