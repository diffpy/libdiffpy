#!/usr/bin/env python

# Installation script for diffpy.fixme

"""diffpy.fixme - template for DiffPy library package
put description of the package here.

Packages:   diffpy.fixme
Scripts:    (none yet)
"""

from setuptools import setup, find_packages
import fix_setuptools_chmod

# define distribution
dist = setup(
        name = 'diffpy.fixme',
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
        maintainer = 'FIXME'
        maintainer_email = 'FIXME@columbia.edu',
        url = 'http://www.diffpy.org/',
        download_url = 'http://www.diffpy.org/packages/',
        description = 'FIXME - what this package does.',
        license = 'BSD',
        keywords = 'FIXME',
        classifiers = [
            # List of possible values at
            # http://pypi.python.org/pypi?:action=list_classifiers
            'Intended Audience :: Science/Research',
            'Programming Language :: Python :: 2.5',
            'Topic :: Scientific/Engineering :: Physics',
        ],
)

# End of file
