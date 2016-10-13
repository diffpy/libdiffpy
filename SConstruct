MY_SCONS_HELP = """\
SCons build rules for the libdiffpy C++ library
Usage: scons [target] [var=value]

Targets:

lib                 build the shared library object [default]
install             install everything under prefix directory
install-lib         install the shared library object
install-include     install the C++ header files
install-data        install data files used by the library
alltests            build the unit test program "alltests"
test                execute unit tests (requires the cxxtest framework)
sdist               create source distribution tarball from git repository

Build configuration variables:
%s
Variables can be also assigned in a user-written script sconsvars.py.
SCons construction environment can be customized in sconscript.local script.
"""

import os
import platform

def subdictionary(d, keyset):
    return dict([kv for kv in d.items() if kv[0] in keyset])

# copy system environment variables related to compilation
DefaultEnvironment(ENV=subdictionary(os.environ, '''
    PATH PYTHONPATH
    CPATH CPLUS_INCLUDE_PATH LIBRARY_PATH LD_RUN_PATH
    LD_LIBRARY_PATH DYLD_LIBRARY_PATH DYLD_FALLBACK_LIBRARY_PATH
    MACOSX_DEPLOYMENT_TARGET
    '''.split())
)

# Create construction environment
env = DefaultEnvironment().Clone()

# Variables definitions below work only with 0.98.1 or later.
env.EnsureSConsVersion(0, 98, 1)

# Customizable compile variables
vars = Variables('sconsvars.py')

vars.Add('tests',
    'Fixed-string patterns for selecting unit test sources.', None)
vars.Add(BoolVariable('test_installed',
    'build unit tests using the installed library.', False))
vars.Add(EnumVariable('build',
    'compiler settings', 'fast',
    allowed_values=('debug', 'fast')))
vars.Add(BoolVariable('profile',
    'build with profiling information', False))
vars.Add(BoolVariable('coverage',
    'build for code coverage analysis', False))
vars.Add(PathVariable('prefix',
    'installation prefix directory', '/usr/local'))
vars.Update(env)
vars.Add(PathVariable('libdir',
    'object code library directory [prefix/lib]',
    env['prefix'] + '/lib',
    PathVariable.PathAccept))
vars.Add(PathVariable('includedir',
    'installation directory for C++ header files [prefix/include]',
    env['prefix'] + '/include',
    PathVariable.PathAccept))
vars.Add(PathVariable('datadir',
    'installation directory for architecture independent data [prefix/share]',
    env['prefix'] + '/share',
    PathVariable.PathAccept))
vars.Add(BoolVariable('enable_objcryst',
    'enable objcryst support, when installed', None))
vars.Update(env)
env.Help(MY_SCONS_HELP % vars.GenerateHelpText(env))

env['has_objcryst'] = None
btags = [env['build'], platform.machine()]
if env['profile']:  btags.append('profile')
if env['coverage']:  btags.append('coverage')
builddir = env.Dir('build/' + '-'.join(btags))

Export('env')

def GlobSources(pattern):
    """Same as Glob but also require that source node is a valid file.
    """
    rv = [f for f in Glob(pattern) if f.srcnode().isfile()]
    return rv

Export('GlobSources')

if os.path.isfile('sconscript.local'):
    env.SConscript('sconscript.local')

env.SConscript('src/SConscript', variant_dir=builddir)

# vim: ft=python
