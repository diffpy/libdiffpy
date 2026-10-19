/*****************************************************************************
*
* libdiffpy         Complex Modeling Initiative
*                   (c) 2013 Brookhaven Science Associates,
*                   Brookhaven National Laboratory.
*                   All rights reserved.
*
* File coded by:    Pavol Juhas
*
* See AUTHORS.txt for a list of people who contributed.
* See LICENSE.rst for license information.
*
******************************************************************************
*
* Functions for resolving paths to static data files at runtime.
*
*****************************************************************************/

#include <cstdlib>
#include <cassert>
#include <filesystem>
#include <sstream>
#include <stdexcept>

#ifdef _WIN32
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <windows.h>
#include <system_error>
#include <vector>
#else
#include <dlfcn.h>
#endif

#include <diffpy/version.hpp>
#include <diffpy/runtimepath.hpp>
#define STRINGIFY(m) STRINGIFY_BRAIN_DAMAGE(m)
#define STRINGIFY_BRAIN_DAMAGE(m) #m

using namespace std;

// Helper Functions ----------------------------------------------------------

namespace {

namespace fs = std::filesystem;

const char* runtimerelpath =
#ifdef DIFFPYRUNTIMERELPATH
    STRINGIFY(DIFFPYRUNTIMERELPATH)
#else
    "../share/diffpy/libdiffpy-"
    STRINGIFY(DIFFPY_VERSION_MAJOR)
    "."
    STRINGIFY(DIFFPY_VERSION_MINOR)
#endif
;

// Locate the module containing this code, including an embedding extension.
fs::path modulepath()
{
    static const char anchor = 0;
#ifdef _WIN32
    HMODULE module = nullptr;
    if (!GetModuleHandleExW(
            GET_MODULE_HANDLE_EX_FLAG_FROM_ADDRESS |
            GET_MODULE_HANDLE_EX_FLAG_UNCHANGED_REFCOUNT,
            reinterpret_cast<LPCWSTR>(&anchor), &module))
    {
        throw system_error(GetLastError(), system_category(),
                "Cannot locate the libdiffpy module");
    }
    vector<wchar_t> buffer(256);
    for (;;)
    {
        DWORD size = GetModuleFileNameW(module, buffer.data(),
                static_cast<DWORD>(buffer.size()));
        if (size == 0)
        {
            throw system_error(GetLastError(), system_category(),
                    "Cannot resolve the libdiffpy module path");
        }
        if (size < buffer.size())
        {
            return fs::path(wstring(buffer.data(), size));
        }
        buffer.resize(2 * buffer.size());
    }
#else
    Dl_info info = {};
    if (!dladdr(&anchor, &info) || !info.dli_fname)
    {
        throw runtime_error("Cannot resolve the libdiffpy module path");
    }
    return fs::path(info.dli_fname);
#endif
}


// Validate an explicit directory without falling back to packaged data.
string canonicalDirectory(const fs::path& directory)
{
    if (!fs::is_directory(directory))
    {
        throw runtime_error("Directory '" + directory.string() +
                "' does not exist.");
    }
    return fs::canonical(directory).generic_string();
}


string diffpyruntime()
{
    // check the DIFFPYRUNTIME environment variable.
    const char* pe = getenv("DIFFPYRUNTIME");
    if (pe && *pe != '\0')
    {
        return canonicalDirectory(pe);
    }
    static const fs::path module = fs::absolute(modulepath());
    // first candidate is resolved in relative data path
    const fs::path d1 = module.parent_path() / runtimerelpath;
    if (fs::is_directory(d1))
    {
        return fs::canonical(d1).generic_string();
    }
    // second candidate is resolved with respect to physical library location
    // according to the source tree layout
    const fs::path d2 = fs::canonical(module).parent_path() /
        "../../src/runtime";
    if (fs::is_directory(d2))
    {
        return fs::canonical(d2).generic_string();
    }
    // nothing worked - throw exception about the first candidate path.
    return canonicalDirectory(d1);
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

// class LineReader ----------------------------------------------------------

// Constructor

LineReader::LineReader() :
    lineno(0),
    commentmark("#")
{ }

// Methods

bool LineReader::isignored() const
{
    return words.empty() || this->iscomment();
}


bool LineReader::iscomment() const
{
    bool rv = !commentmark.empty() && !words.empty() &&
            (0 == words[0].compare(0, commentmark.size(), commentmark));
    return rv;
}


bool LineReader::isblank() const
{
    return words.empty();
}


size_t LineReader::wcount() const
{
    return words.size();
}


runtime_error LineReader::format_error(
        const string& filename, const string& edetail)
{
    ostringstream emsg;
    emsg << "Invalid data format in " <<
        filename << " line " << this->lineno << '.';
    if (!edetail.empty())  emsg << '\n' << edetail;
    return runtime_error(emsg.str());
}

// Non-member operators

istream& operator>> (istream& fp, LineReader& line)
{
    getline(fp, line.line);
    if (fp)  line.lineno++;
    istringstream fpline(line.line);
    line.words.clear();
    string w;
    if (line.separator.empty()) {
        while (fpline >> w)  line.words.push_back(w);
    }
    else {
        assert(line.separator.size() == 1);
        const char& sep = line.separator[0];
        while (getline(fpline, w, sep))  line.words.push_back(w);
    }
    return fp;
}

}   // namespace runtimepath
}   // namespace diffpy
