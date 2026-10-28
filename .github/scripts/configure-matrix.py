"""Configure SCons and CxxTest for the CI matrix."""

from pathlib import Path

test_requirements = Path("requirements/tests.txt").read_text().splitlines()
Path("build").mkdir(exist_ok=True)
Path("build/matrix-tests.txt").write_text(
    "\n".join(line for line in test_requirements if line.strip() != "cxxtest")
    + "\n"
)

Path("sconscript.local").write_text("""\
Import('env')
import os
import sys

prefix = os.environ['CONDA_PREFIX']
cxxtest = Dir('#/build/cxxtest').abspath
env.Replace(CXXTEST=os.path.join(cxxtest, 'bin', 'cxxtestgen'),
            CXXTEST_INSTALL_DIR=cxxtest, CXXTEST_PYTHON=sys.executable)

if env['PLATFORM'] != 'win32':
    env.Replace(CXX=os.environ.get('CXX', os.path.join(prefix, 'bin', 'c++')))
    env.Prepend(CCFLAGS=['-isystem' + os.path.join(prefix, 'include')])
    env.Append(LIBPATH=os.path.join(prefix, 'lib'))
    env.Append(LINKFLAGS='-Wl,-rpath,' + os.path.join(prefix, 'lib'))
    if env['PLATFORM'] == 'darwin':
        sdk = os.environ.get('CONDA_BUILD_SYSROOT')
        if sdk:
            env.Append(CXXFLAGS=['-isysroot', sdk])
            env.Append(LINKFLAGS=['-isysroot', sdk])
""")
