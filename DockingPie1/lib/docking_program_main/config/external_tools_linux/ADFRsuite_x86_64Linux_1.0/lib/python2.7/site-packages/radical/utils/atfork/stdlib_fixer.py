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
This module provides code to setup appropriate atfork() calls within the
Python standard library in order to make it safe for use in programs that
mix fork() and threads.

Import this and call the appropriate fixers very early in your program's
initialization sequence, preferrably -before- importing the modules you
need to fix.

Provided fixers:

    fix_logging_module()

TODO(gps): Audit more of the stdlib and provide necessary fixers.
In 2.4.5 the following additional stdlib modules use locks:

  threading.Condition and similar classes could use it.
  Queue
  cookielib
  mimetools
  _strptime
"""

import sys
import warnings

from .atfork import atfork


class Error(Exception):
    pass


def fix_logging_module():
    logging = sys.modules.get('logging')
    # Prevent fixing multiple times as that would cause a deadlock.
    if logging and getattr(logging, 'fixed_for_atfork', None):
        return
    if logging:
        warnings.warn('logging module already imported before fixup.')
    import logging
    if logging.getLogger().handlers:
        # We could register each lock with atfork for these handlers but if
        # these exist, other loggers or not yet added handlers could as well.
        # Its safer to insist that this fix is applied before logging has been
        # configured.
        warnings.warn('logging handlers already registered.')
      # raise Error('logging handlers already registered.')

    logging._acquireLock()
    try:
        def fork_safe_createLock(self):
            self._orig_createLock()
            atfork(self.lock.acquire,
                   self.lock.release, self.lock.release)

        # Fix the logging.Handler lock (a major source of deadlocks).
        logging.Handler._orig_createLock = logging.Handler.createLock
        logging.Handler.createLock = fork_safe_createLock

        # Fix the module level lock.
        atfork(logging._acquireLock,
               logging._releaseLock, logging._releaseLock)

        logging.fixed_for_atfork = True
    finally:
        logging._releaseLock()
