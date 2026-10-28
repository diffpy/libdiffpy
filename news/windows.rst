**Added:**

* Support building libdiffpy on Windows with MSVC, clang or MinGW, merged from the upstream ``windows`` branch.
* Export the public API from the Windows DLL with the ``DLL_EXPORT`` macro in ``diffpy/Export.hpp``.
* Install headers, libraries and data under ``$PREFIX`` or ``$CONDA_PREFIX`` by default, using the ``Library/`` layout on Windows.
* Add a Windows build and test job to the pull request CI.
* Add a clang-format pre-commit hook and ``.git-blame-ignore-revs`` for the reformat commits.

**Changed:**

* Format all C++ sources with clang-format (Google style).
* Require ObjCryst by default. Build with ``enable_objcryst=False`` to skip it.

**Deprecated:**

* <news item>

**Removed:**

* <news item>

**Fixed:**

* <news item>

**Security:**

* <news item>
