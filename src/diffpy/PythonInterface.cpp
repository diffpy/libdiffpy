/*****************************************************************************
*
* diffpy.srreal     by DANSE Diffraction group
*                   Simon J. L. Billinge
*                   (c) 2009 Trustees of the Columbia University
*                   in the City of New York.  All rights reserved.
*
* File coded by:    Pavol Juhas
*
* See AUTHORS.txt for a list of people who contributed.
* See LICENSE.txt for license information.
*
******************************************************************************
*
* Common functions for interfacing with Python interpreter.
*
*****************************************************************************/

#include <diffpy/PythonInterface.hpp>

#include <iostream>
#include <csignal>
#include <map>

using namespace std;
using namespace boost;

namespace diffpy {

// Local Helpers -------------------------------------------------------------

namespace {

// Flag if Python has been launched from C++.  This means there is no
// Python interpreter to handle error_already_set exceptions.

bool python_is_embedded = false;


// Execute import from python without any error handling
python::object do_python_import(const string& modname, const string& item)
{
    typedef std::map<string, python::object> ObjectCache;
    static ObjectCache cacheditems;
    string fullname = modname + ":" + item;
    ObjectCache::iterator ii = cacheditems.find(fullname);
    if (ii == cacheditems.end())
    {
        python::object mod = python::import(modname.c_str());
        python::object obj = mod.attr(item.c_str());
        cacheditems[fullname] = obj;
        ii = cacheditems.find(fullname);
    }
    return ii->second;
}

}   // namespace

// Public Functions ----------------------------------------------------------

void initializePython(int py_argc, char* py_argv[])
{
    if (Py_IsInitialized())  return;
    if (!py_argc && !py_argv)
    {
        static const int initpy_argc = 1;
        static char initpy_arg0[7] = "python";
        static char* initpy_argv[initpy_argc] = {initpy_arg0};
        py_argc = initpy_argc;
        py_argv = initpy_argv;
    }
    Py_Initialize();
    PySys_SetArgv(py_argc, py_argv);
    // Make sure Python does not eat SIGINT.
    python_is_embedded = true;
    signal(SIGINT, SIG_DFL);
}


string getPythonErrorString()
{
    string rv;
    PyObject* ptype = NULL;
    PyObject* pvalue = NULL;
    PyObject* ptraceback = NULL;
    PyErr_Fetch(&ptype, &pvalue, &ptraceback);
    if (ptype == NULL)  return rv;
    PyObject* pemsg = PyObject_Str(pvalue);
    rv = PyString_AsString(pemsg);
    Py_XDECREF(pemsg);
    PyErr_Restore(ptype, pvalue, ptraceback);
    return rv;
}


python::object importFromPyModule(const string& modname, const string& item,
        python::object fallback)
{
    python::object rv;
    try {
        rv = do_python_import(modname, item);
    }
    catch (python::error_already_set e) {
        PyErr_Clear();
        rv = fallback;
    }
    return rv;
}


python::object importFromPyModule(const string& modname, const string& item)
{
    python::object rv;
    try {
        rv = do_python_import(modname, item);
    }
    catch (python::error_already_set e) {
        // display error message when running embedded
        if (python_is_embedded)
        {
            cerr << getPythonErrorString() << endl;
        }
        throw;
    }
    return rv;
}


}   // namespace diffpy

// End of file
