"""Execute scons before normal setuptools install.
"""

# version
__id__ = "$Id$"


import os
import sys

from distutils.spawn import spawn
from setuptools.command.bdist_egg import bdist_egg
from setuptools.command.easy_install import easy_install


class scons_bdist_egg(bdist_egg):
    '''This class applies when setup.py is sourced from easy_install script.
    We need to walk the trace to recover the original prefix.
    '''
    
    def run(self):
        bdist_egg.run(self)
        if self._root_setup_stack_frame():
            _launch_scons_install(self)
        return


    def _get_scons_prefix(self):
        '''Return installation prefix directory.
        '''
        frame = self._root_setup_stack_frame()
        dist = frame.f_locals['dist']
        cmd = dist.command_obj['easy_install']
        rv = _cut_libdir_prefix(cmd.install_dir)
        return rv


    def _root_setup_stack_frame(self):
        '''Get stack frame of the root setup function in distutils.core.
        Only works when the main executable is easy_install.
        
        Return stack frame when run from easy_install.
        Return None in all other cases.
        '''
        import traceback
        stack = traceback.extract_stack()
        def find_mod_fnc(modname, fncname):
            rv = -1
            modpat = modname.replace('.', os.sep)
            for i, tb in enumerate(stack):
                if modpat in tb[0] and tb[2] == fncname:
                    rv = i
                    break
            return rv
        idxmain = find_mod_fnc('setuptools.command.easy_install', 'main')
        idxsetup = find_mod_fnc('distutils.core', 'setup')
        # recover setup frame only when it exists and was called from
        # easy_install main
        if 0 < idxmain < idxsetup:
            depth = len(stack) - 1 - idxsetup
            rv = sys._getframe(depth)
        else:
            rv = None
        return rv

# class scons_easy_install


class scons_easy_install(easy_install):
    '''This class applies only when setup.py is the main executable script.
    When installed by easy_install script, the distribution object gets
    created before sourcing this file.  Therefore the scons needs to 
    be hooked in the scons_bdist_egg class.
    '''
    
    def run(self):
        easy_install.run(self)
        _launch_scons_install(self)
        return


    def _get_scons_prefix(self):
        '''Return installation prefix directory.
        '''
        rv = _cut_libdir_prefix(self.install_dir)
        return rv

# class scons_easy_install


def _launch_scons_install(cmd):
    '''Execute scons install under a prefix obtained from distutils.

    cmd  -- instance of scons_easy_install or scons_bdist_egg

    No return value.
    '''
    sconsargv = ['scons', 'install',
            'build=fast', 'prefix=%s' % cmd._get_scons_prefix()]
    #print "sconsargv = %r" % sconsargv
    spawn(sconsargv, verbose=True)
    return


def _cut_libdir_prefix(libdirpath):
    '''Get prefix by cutting away lib/pythonX.Y/ from a path.

    libdirpath  -- installation path for python packages

    Return string.
    '''
    pyver = 'python' + sys.version[:3]
    ptail = os.sep.join(['', 'lib', pyver])
    idx = libdirpath.index(ptail)
    rv = libdirpath[0:idx]
    return rv

# End of file
