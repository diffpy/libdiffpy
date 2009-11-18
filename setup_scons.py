"""Execute scons before normal setuptools install.
"""

# version
__id__ = "$Id$"


from distutils.spawn import spawn
from setuptools.command.bdist_egg import bdist_egg

class scons_bdist_egg(bdist_egg):

    def run(self):
        cmdinst = self.get_finalized_command('install')
        sconsargv = ['scons', 'install',
                'build=fast', 'prefix=%s' % cmdinst.prefix]
        spawn(sconsargv, verbose=True)
        bdist_egg.run(self)
        return

# class scons_install

# End of file
