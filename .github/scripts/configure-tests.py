"""Configure SCons, CxxTest and the conda environment for CI tests."""

import json
from pathlib import Path

requirements = ["python"]
for filename in ("build.txt", "conda.txt", "tests.txt"):
    for line in (Path("requirements") / filename).read_text().splitlines():
        spec = line.split("#", 1)[0].strip()
        if spec and spec != "cxxtest":
            requirements.append(spec)

Path("build").mkdir(exist_ok=True)
Path("build/test-environment.yml").write_text(
    json.dumps(
        {
            "name": "test",
            "channels": ["conda-forge"],
            "dependencies": requirements,
        },
        indent=2,
    )
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

if env['PLATFORM'] == 'win32':
    from SCons.Platform.win32 import exec_spawn

    # The outer CI shell already activates conda and SCons configures MSVC.
    # /D stops each compiler/test-generator child from repeating cmd AutoRun.
    def spawn_without_autorun(sh, escape, cmd, args, child_env):
        return exec_spawn([sh, '/D', '/C', escape(' '.join(args))], child_env)

    env['SPAWN'] = spawn_without_autorun
else:
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
