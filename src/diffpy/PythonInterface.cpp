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
* $Id$
*
*****************************************************************************/

#include <diffpy/PythonInterface.hpp>

#include <csignal>
#include <map>

using std::string;
using namespace boost;

namespace diffpy {

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


python::object importFromPyModule(const string& modname, const string& item)
{
    static std::map<string, python::object> cacheditems;
    string fullname = modname + ":" + item;
    // perform import when not in the cache
    if (!cacheditems.count(fullname))
    {
        try {
            cacheditems[fullname] =
                python::import(modname.c_str()).attr(item.c_str());
        }
        // Ignore import errors when item could not be recovered.
        catch (python::error_already_set e) {
            PyErr_Clear();
            cacheditems[fullname] = python::object();
        }
    }
    return cacheditems[fullname];
}


}   // namespace diffpy

// End of file
