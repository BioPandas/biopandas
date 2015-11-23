# BioPandas
# Author: Sebastian Raschka <mail@sebastianraschka.com>
# License: BSD 3 clause
# Project Website: http://rasbt.github.io/biopandas/
# Code Repository: https://github.com/rasbt/biopandas

import difflib


def assertMultiLineEqual(first, second, preserve_newline=True, msg=None):
    """Assert that two multi-line strings are equal."""
    assert isinstance(first, str) == True, 'First argument is not a string'
    assert isinstance(second, str) == True, 'Second argument is not a string'

    if first != second:
        message = ''.join(difflib.ndiff(first.splitlines(preserve_newline),
                                        second.splitlines(preserve_newline)))
        if msg:
            message += " : " + msg
        raise AssertionError("Multi-line strings are unequal:\n" + message)
