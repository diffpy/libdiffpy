#!/usr/bin/env python

"""Utility functions used by scons and sphinx for version extraction."""

import os
import re
import subprocess

MYDIR = os.path.dirname(os.path.abspath(__file__))

EMSG_GIT_MISSING = """
Cannot determine libdiffpy version.  Compile from a git repository
or use a source archive from

    https://github.com/diffpy/libdiffpy/releases
"""

EMSG_BAD_FALLBACK_VERSION = """
Inconsistent FALLBACK_VERSION {} vs {} from git tag.
"""

RELEASE_TAG = re.compile(r"([0-9]+\.[0-9]+\.[0-9]+)(?:-rc\.([0-9]+))?")


def release_version(tag):
    """Normalize a release tag to the internal version format."""
    match = RELEASE_TAG.fullmatch(tag)
    if not match:
        raise RuntimeError(
            f"Invalid release tag {tag!r}; expected X.Y.Z or X.Y.Z-rc.N."
        )
    base, candidate = match.groups()
    return base if candidate is None else f"{base}rc{candidate}"


def gitinfo():
    "Extract dictionary of version data from git records."
    global _cached_gitinfo
    if _cached_gitinfo is not None:
        return _cached_gitinfo
    kw = dict(
        stdout=subprocess.PIPE,
        stderr=subprocess.DEVNULL,
        cwd=MYDIR,
        text=True,
    )
    proc = subprocess.run(["git", "log", "-1", "--format=%H %ci"], **kw)
    if proc.returncode:
        _cached_gitinfo = {}
        return _cached_gitinfo
    rv = {}
    rv["commit"], rv["date"] = proc.stdout.strip().split(None, 1)
    # --long makes the commit count explicit even at a tag.  Split from
    # the right so the -rc.N suffix remains part of the release tag.
    proc = subprocess.run(
        ["git", "describe", "--tags", "--long", "--match=[0-9]*"], **kw
    )
    if not proc.returncode:
        tag, count, _ = proc.stdout.strip().rsplit("-", 2)
        rv["patchnumber"] = int(count)
        rv["version"] = release_version(tag)
        if rv["patchnumber"]:
            rv["version"] += ".post" + count
    _cached_gitinfo = rv
    return _cached_gitinfo


_cached_gitinfo = None


def getversion():
    """Extract version from gitinfo and/or gitarchive.cfg file.

    Process gitinfo first and make sure it matches FALLBACK_VERSION in
    fallback_version.py.  Use expanded data from gitarchive.cfg when
    these sources are from git archive bundle.
    """
    from fallback_version import FALLBACK_VERSION

    if not re.fullmatch(
        r"[0-9]+\.[0-9]+\.[0-9]+(?:rc[0-9]+)?\.post0", FALLBACK_VERSION
    ):
        raise RuntimeError(
            "FALLBACK_VERSION must be X.Y.Z.post0 or X.Y.ZrcN.post0."
        )
    gitarchivecfgfile = os.path.join(MYDIR, "gitarchive.cfg")
    assert os.path.isfile(gitarchivecfgfile)
    with open(gitarchivecfgfile) as fp:
        gacontent = fp.read()
    gaitems = re.findall(r"^(\w+) *= *(\S.*?)\s*$", gacontent, re.M)
    ga = dict(gaitems)
    gi = gitinfo()
    rv = {}
    if not gi:
        # Not a git repository.  Require an expanded source archive.
        if "$Format:" in ga["commit"]:
            raise RuntimeError(EMSG_GIT_MISSING)
        gi = dict(commit=ga["commit"], date=ga["date"])
        for ref in ga["refnames"].split(", "):
            if ref.startswith("tag: ") and re.match(r"[0-9]", ref[5:]):
                gi["version"] = release_version(ref[5:])
                gi["patchnumber"] = 0
                break
    rv.update(gi)
    if "version" in gi:
        afb = FALLBACK_VERSION
        gfb = gi["version"].split(".post")[0] + ".post0"
        if gfb != afb:
            raise RuntimeError(EMSG_BAD_FALLBACK_VERSION.format(afb, gfb))
    else:
        # A development checkout without a new-format tag, or an untagged
        # archive, has no known commit count.  Keep the fallback version.
        rv["version"] = FALLBACK_VERSION
        rv["patchnumber"] = 1
    # rv['version'] is resolved, let's parse its subfields.
    vbase = rv["version"].split(".post")[0]
    mx = re.fullmatch(r"([0-9]+)\.([0-9]+)\.([0-9]+)(rc[0-9]+)?", vbase)
    rv["major"] = int(mx.group(1))
    rv["minor"] = int(mx.group(2))
    rv["micro"] = int(mx.group(3))
    rv["prerelease"] = mx.group(4)
    return rv
