#!/usr/bin/env python

# Installation script for diffpy.lib

"""diffpy.lib - template for DiffPy library package
put description of the package here.

Packages:   diffpy.lib
Scripts:    diffpy-config
"""

from setuptools import setup, find_packages
import fix_setuptools_chmod

# define distribution
dist = setup(
        name = 'diffpy.lib',
        version = '0.1a1',
        namespace_packages = ['diffpy'],
        packages = find_packages(),
        entry_points = {
            # define console_scripts here, see setuptools docs for details.
        },
        install_requires = [
            # insert list of DiffPy and DANSE dependencies here
        ],
        dependency_links = [
            # REMOVE dev.danse.us for a public release.
            'http://dev.danse.us/packages/',
            'http://www.diffpy.org/packages/',
        ],

        author = 'Simon J.L. Billinge',
        author_email = 'sb2896@columbia.edu',
        maintainer = 'Pavol Juhas',
        maintainer_email = 'pj2192@columbia.edu',
        url = 'http://www.diffpy.org/',
        download_url = 'http://www.diffpy.org/packages/',
        description = 'DiffPy C++ monolithic library.',
        license = 'BSD',
        keywords = 'diffpy C++ library headers',
        classifiers = [
            # List of possible values at
            # http://pypi.python.org/pypi?:action=list_classifiers
            'Intended Audience :: Science/Research',
            'Programming Language :: C++',
            'Topic :: Scientific/Engineering :: Physics',
        ],
)

# End of file
