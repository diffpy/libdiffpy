#!/usr/bin/env python

'''Utility functions used by scons and sphinx for version extraction.
'''

import os

MYDIR = os.path.dirname(os.path.abspath(__file__))

def gitinfo():
    'Extract dictionary of version data from git records.'
    import re
    from subprocess import Popen, PIPE
    kw = dict(stdout=PIPE, cwd=MYDIR)
    proc = Popen(['git', 'describe', '--match=v[[:digit:]]*'], **kw)
    desc = proc.stdout.read()
    proc = Popen(['git', 'log', '-1', '--format=%H %ai'], **kw)
    glog = proc.stdout.read()
    rv = {}
    rv['version'] = '-'.join(desc.strip().split('-')[:2]).lstrip('v')
    rv['commit'], rv['date'] = glog.strip().split(None, 1)
    mx = re.search(r'(?m)^(\d+)\.(\d+)(?:-(\d+))?', rv['version'])
    rv['major'] = int(mx.group(1))
    rv['minor'] = int(mx.group(2))
    rv['number'] = mx.group(3) and int(mx.group(3)) or 0
    return rv
