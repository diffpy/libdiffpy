#!/usr/bin/env python

"""Unit tests for class Foo from diffpy.fixme.foo
"""

# version
__id__ = '$Id$'

import os
import unittest

# useful variables
thisfile = locals().get('__file__', 'file.py')
tests_dir = os.path.dirname(os.path.abspath(thisfile))
# testdata_dir = os.path.join(tests_dir, 'testdata')

from diffpy.fixme.foo import Foo, func

##############################################################################
class TestRoutines(unittest.TestCase):

    def setUp(self):
        return

    def tearDown(self):
        return

    def test_func(self):
        """check func()
        """
        self.assertEqual(0, func(0))
        self.assertEqual(37, func(1))
        self.assertEqual(-18.5, func(-0.5))
        return

# End of class TestRoutines

##############################################################################
class TestFoo(unittest.TestCase):

    def setUp(self):
        self.foo = Foo()
        return

    def tearDown(self):
        return

    def test___init__(self):
        """check Foo.__init__()
        """
        # self.foo is created before every test in setUp method
        # It should assign default value, which is 42.
        self.assertEqual(42, self.foo.value())
        # check __init__ arguments
        foo1 = Foo(11.3)
        self.assertEqual(11, foo1._value)
        foo1 = Foo('-99')
        self.assertEqual(-99, foo1.value())
        # Foo should raise exception for invalid argument
        self.assertRaises(ValueError, Foo, 'nointeger')
        return

    def test_value(self):
        """check Foo.value()
        """
        foo = Foo(100)
        self.assertEqual(100, foo.value())
        return

# End of class TestFoo

if __name__ == '__main__':
    unittest.main()

# End of file
