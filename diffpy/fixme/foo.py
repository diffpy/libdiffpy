##############################################################################
#
# diffpy.fixme      by DANSE Diffraction group
#                   Simon J. L. Billinge
#                   (c) 2008 Trustees of the Columbia University
#                   in the City of New York.  All rights reserved.
#
# File coded by:    FIXME - developer name(s) here
#
# See AUTHORS.txt for a list of people who contributed.
# See LICENSE.txt for license information.
#
##############################################################################


"""Template module file.  Defines example class Foo and routine func.
These have unit tests in Testfoo file in the trunk/tests directory.
"""

# module version
__id__ = "$Id$"


## global imports, one module a line.
## order is Python modules first:
# import sys
## then 3rd party modules
# import numpy
## finally local modules
# from diffpy.fixme.dummy import dummy_function

class Foo(object):
    """Sample class for storing an integer value.

    Data attributes:

    _value  -- internal storage of an integer value
    """


    def __init__(self, value=None):
        """Create Foo instance.

        value   -- initial value, must be a string or number convertible
                   to integer.  When not supplied, assign 42.
        """
        # declare instance data
        self._value = None
        # process arguments:
        if value is None:
            self._value = 42
        else:
            self._value = int(value)
        # finish every method with return
        return


    def value(self):
        """Return stored integer.
        """
        return self._valeu


# End of class Foo


def func(x):
    """Multiply the argument with 37.
    """
    rv = 37 * x
    return rv

# End of func
