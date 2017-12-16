#!/usr/bin/env python

'''Utility functions used by scons and sphinx for version extraction.
'''

import os
import re

MYDIR = os.path.dirname(os.path.abspath(__file__))

EMSG_GIT_MISSING = """
Cannot determine libdiffpy version.  Compile from a git repository
or use a source archive from

    https://github.com/diffpy/libdiffpy/releases
"""

EMSG_BAD_FALLBACK_VERSION = """
Inconsistent FALLBACK_VERSION {} vs {} from git tag.
"""


def gitinfo():
    'Extract dictionary of version data from git records.'
    from subprocess import Popen, PIPE
    global _cached_gitinfo
    if _cached_gitinfo is not None:
        return _cached_gitinfo
    nullfile = open(os.devnull, 'w')
    kw = dict(stdout=PIPE, stderr=nullfile, cwd=MYDIR,
              universal_newlines=True)
    proc = Popen(['git', 'describe', '--match=v[[:digit:]]*'], **kw)
    desc = proc.stdout.read()
    if proc.wait():
        _cached_gitinfo = {}
        return gitinfo()
    proc = Popen(['git', 'log', '-1', '--format=%H %ai'], **kw)
    glog = proc.stdout.read()
    words = desc.strip().split('-')
    vtag = words[0].lstrip('v')
    vnum = int((words[1:2] or ['0'])[0])
    rv = {}
    rv['commit'], rv['date'] = glog.strip().split(None, 1)
    rv['patchnumber'] = vnum
    rv['version'] = vtag
    if vnum:
        rv['version'] += '.post' + str(vnum)
    _cached_gitinfo = rv
    return gitinfo()

_cached_gitinfo = None


def getversion():
    """Extract version from gitinfo and/or gitarchive.cfg file.

    Process gitinfo first and make sure it matches FALLBACK_VERSION
    in gitarchive.cfg.  Use expanded data from gitarchive.cfg when
    these sources are from git archive bundle.
    """
    from fallback_version import FALLBACK_VERSION
    gitarchivecfgfile = os.path.join(MYDIR, 'gitarchive.cfg')
    assert os.path.isfile(gitarchivecfgfile)
    with open(gitarchivecfgfile) as fp:
        gacontent = fp.read()
    gaitems = re.findall(r'^(\w+) *= *(\S.*?)\s*$', gacontent, re.M)
    ga = dict(gaitems)
    gi = gitinfo()
    rv = {}
    if gi:
        afb = FALLBACK_VERSION
        gfb = gi['version'].split('.post')[0] + '.post0'
        if gfb != afb:
            raise RuntimeError(EMSG_BAD_FALLBACK_VERSION.format(afb, gfb))
        rv.update(gi)
    else:
        # Not a git repository.  Require that gitarchive.cfg is expanded.
        if '$Format:' in ga['commit']:
            raise RuntimeError(EMSG_GIT_MISSING)
        rv['commit'] = ga['commit']
        rv['date'] = ga['date']
        # First assume we have an undetermined post-release.  Keep version
        # suffix as ".post0", but set patchnumber to ensure we are post
        # FALLBACK_VERSION.
        rv['version'] = FALLBACK_VERSION
        rv['patchnumber'] = 1
        # If we are at version tag turn this to regular release.
        mx = re.search(r'\btag: v(\d[^,]*)', ga['refnames'])
        if mx:
            rv['version'] = mx.group(1)
            rv['patchnumber'] = 0
    # rv['version'] is resolved, let's parse its subfields.
    vbase = rv['version'].split('.post')[0]
    mx = re.match(r'^(\d+)\.(\d+)(?:\.(\d+))?((?:[ab]|rc)\d*)?', vbase)
    rv['major'] = int(mx.group(1))
    rv['minor'] = int(mx.group(2))
    rv['micro'] = int(mx.group(3) or 0)
    rv['prerelease'] = mx.group(4)
    return rv
