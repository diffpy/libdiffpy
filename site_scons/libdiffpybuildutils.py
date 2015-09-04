#!/usr/bin/env python

'''Utility functions used by scons and sphinx for version extraction.
'''

import os

MYDIR = os.path.dirname(os.path.abspath(__file__))

def gitinfo():
    'Extract dictionary of version data from git records.'
    import re
    from subprocess import Popen, PIPE
    global _cached_gitinfo
    if _cached_gitinfo is not None:  return _cached_gitinfo
    rv = _cached_gitinfo = {}
    nullfile = open(os.devnull, 'w')
    kw = dict(stdout=PIPE, stderr=nullfile, cwd=MYDIR)
    proc = Popen(['git', 'describe', '--match=v[[:digit:]]*'], **kw)
    desc = proc.stdout.read()
    if proc.wait():  return gitinfo()
    proc = Popen(['git', 'log', '-1', '--format=%H %ai'], **kw)
    glog = proc.stdout.read()
    rv['version'] = '.'.join(desc.strip().split('-')[:2]).lstrip('v')
    rv['commit'], rv['date'] = glog.strip().split(None, 1)
    mx = re.search(r'(?m)^(\d+)\.(\d+)([ab]\d*)?(?:[.](\d+))?', rv['version'])
    rv['major'] = int(mx.group(1))
    rv['minor'] = int(mx.group(2))
    rv['prerelease'] = mx.group(3)
    rv['number'] = mx.group(4) and int(mx.group(4)) or 0
    _cached_gitinfo = rv
    return gitinfo()

_cached_gitinfo = None
