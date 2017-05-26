# BioPandas
# Author: Sebastian Raschka <mail@sebastianraschka.com>
# License: BSD 3 clause
# Project Website: http://rasbt.github.io/biopandas/
# Code Repository: https://github.com/rasbt/biopandas

import difflib


def assertMultiLineEqual(first, second, preserve_newline=True, msg=None):
    """Assert that two multi-line strings are equal."""
    assert isinstance(first, str), 'First argument is not a string'
    assert isinstance(second, str), 'Second argument is not a string'

    if first != second:
        message = ''.join(difflib.ndiff(first.splitlines(preserve_newline),
                                        second.splitlines(preserve_newline)))
        if msg:
            message += " : " + msg
        raise AssertionError("Multi-line strings are unequal:\n" + message)


def assert_raises(exception_type, message, func, *args, **kwargs):
    """Check that an exception is raised with a specific message

    Parameters
    ----------
    exception_type : exception
        The exception that should be raised
    message : str (default: None)
        The error message that should be raised. Ignored if False or None
    func : callable
        The function that raises the exception
    *args : positional arguments to `func`
    **kwargs : keyword arguments to `func`

    """
    try:
        func(*args, **kwargs)
    except exception_type as e:
        error_message = str(e)
        if message and message not in error_message:
            raise AssertionError("Error message differs from the expected"
                                 " string: %r. Got error message: %r" %
                                 (message, error_message))
    else:
        raise AssertionError('%s not raised.' % exception_type.__name__)
