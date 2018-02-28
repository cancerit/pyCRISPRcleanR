from __future__ import absolute_import, division, print_function
from builtins import map
from past.builtins import basestring

import contextlib
import logging
import os
import subprocess
import tempfile


# __________________________________________________________________________
# I/O helpers

def call_quiet(*args):
    """Safely run a command and get stdout; print stderr if there's an error.
    Like subprocess.check_output, but silent in the normal case where the
    command logs unimportant stuff to stderr. If there is an error, then the
    full error message(s) is shown in the exception message.
    """
    # args = map(str, args)
    if not len(args):
        raise ValueError("Must supply at least one argument (the command name)")
    try:
        proc = subprocess.Popen(args, stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE)
    except OSError as exc:
        raise RuntimeError("Could not find the executable %r" % args[0]
                           + " -- is it installed correctly?"
                           + "\n(Original error: %s)" % exc)
    out, err = proc.communicate()
    if proc.returncode != 0:
        raise RuntimeError("Subprocess command failed:\n$ %s\n\n%s"
                           % (' '.join(args), err))
    return out

def ensure_path(fname):
    """Create dirs and move an existing file to avoid overwriting, if necessary.
    If a file already exists at the given path, it is renamed with an integer
    suffix to clear the way.
    """
    if '/' in os.path.normpath(fname):
        # Ensure the output directory exists
        dname = os.path.dirname(os.path.abspath(fname))
        if dname and not os.path.isdir(dname):
            try:
                os.makedirs(dname)
            except OSError as exc:
                raise OSError("Output path " + fname +
                              " contains a directory " + dname +
                              " that cannot be created: %s" % exc)
    if os.path.isfile(fname):
        # Add an integer suffix to the existing file name
        cnt = 1
        bak_fname = "%s.%d" % (fname, cnt)
        while os.path.isfile(bak_fname):
            cnt += 1
            bak_fname = "%s.%d" % (fname, cnt)
        os.rename(fname, bak_fname)
        logging.info("Moved existing file %s -> %s", fname, bak_fname)
    return True


@contextlib.contextmanager
def temp_write_text(text, mode="w+b"):
    """Save text to a temporary file.
    NB: This won't work on Windows b/c the file stays open.
    """
    with tempfile.NamedTemporaryFile(mode=mode) as tmp:
        tmp.write(text)
        tmp.flush()
        yield tmp.name
