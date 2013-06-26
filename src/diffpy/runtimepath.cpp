/*****************************************************************************
*
* diffpy.srreal     Complex Modeling Initiative
*                   Pavol Juhas
*                   (c) 2013 Brookhaven National Laboratory,
*                   Upton, New York.  All rights reserved.
*
* File coded by:    Pavol Juhas
*
* See AUTHORS.txt for a list of people who contributed.
* See LICENSE.txt for license information.
*
******************************************************************************
*
* Functions for resolving paths to static data files at runtime.
*
*****************************************************************************/

#include <cstdlib>
#include <stdexcept>
#include <sys/stat.h>
#include <dlfcn.h>
#include <libgen.h>

#include <diffpy/version.hpp>
#include <diffpy/runtimepath.hpp>
#define STRINGIFY(m) STRINGIFY_BRAIN_DAMAGE(m)
#define STRINGIFY_BRAIN_DAMAGE(m) #m

using namespace std;

// Helper Functions ----------------------------------------------------------

namespace {

const char* runtimerelpath =
#ifdef DIFFPYRUNTIMERELPATH
    STRINGIFY(DIFFPYRUNTIMERELPATH)
#else
    "../share/diffpy/diffpy"
    STRINGIFY(DIFFPY_VERSION_MAJOR)
    STRINGIFY(DIFFPY_VERSION_MINOR)
#endif
;

bool isdir(const string& d)
{
    struct stat sd;
    return (stat(d.c_str(), &sd) == 0 && S_ISDIR(sd.st_mode));
}


/// throw runtime error if directory does not exist.
void ensureIsDir(const std::string& d)
{
    if (isdir(d))  return;
    string emsg("Directory '@D' does not exist.");
    emsg = emsg.replace(emsg.find_first_of('@'), 2, d);
    throw runtime_error(emsg);
}


const string& diffpyruntime()
{
    char fpb[PATH_MAX];
    static string rv;
    static bool didruntime = false;
    static string dprt;
    // check the DIFFPYRUNTIME environment variable.
    char* pe = getenv("DIFFPYRUNTIME");
    if (pe)  didruntime = didruntime && pe && dprt == pe;
    if (didruntime)  return rv;
    if (pe)
    {
        dprt = pe;
        rv = dprt;
        size_t pt = rv.find_last_not_of('/');
        if (pt == string::npos)  rv.clear();
        else  rv.erase(pt + 1);
        ensureIsDir(rv);
        rv = realpath(rv.c_str(), fpb);
        didruntime = true;
        return rv;
    }
    // resolve path of the libdiffpy shared library file
    Dl_info i;
    dladdr(reinterpret_cast<void*>(diffpyruntime), &i);
    // first candidate is resolved in relative data path
    strncpy(fpb, i.dli_fname, PATH_MAX);
    string d1 = string(dirname(fpb)) + "/" + runtimerelpath;
    if (isdir(d1))
    {
        rv = realpath(d1.c_str(), fpb);
        didruntime = true;
        return rv;
    }
    // second candidate is resolved with respect to physical libary location
    // according to the source tree layout
    string d2 = dirname(realpath(i.dli_fname, fpb));
    d2 += "/../../src/runtime";
    if (isdir(d2))
    {
        rv = realpath(d2.c_str(), fpb);
        didruntime = true;
        return rv;
    }
    // nothing worked - throw exception about the first candidate path.
    ensureIsDir(d1);
    return rv;
}

}   // namespace

// Implementation ------------------------------------------------------------

namespace diffpy {
namespace runtimepath {

string datapath(const std::string& f)
{
    string rv = diffpyruntime();
    rv += (f.empty() ? "" : "/") + f;
    return rv;
}


}   // namespace runtimepath
}   // namespace diffpy
