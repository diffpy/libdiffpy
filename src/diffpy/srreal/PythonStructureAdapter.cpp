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
* createStructureAdapter - a factory for creating StructureAdapter
*     for recognized Python objects
*
* $Id$
*
*****************************************************************************/

#include <diffpy/srreal/PythonStructureAdapter.hpp>

#include <stdexcept>
#include <sstream>
#include <string>
#include <set>

using namespace boost;
using namespace std;

namespace diffpy {
namespace srreal {

// Local Helpers -------------------------------------------------------------

namespace {

typedef set<PythonStructureAdapterFactory> PSAFactorySet;

PSAFactorySet& getPSAFactorySet()
{
    static PSAFactorySet factories;
    return factories;
}

}   // namespace

// Routines ------------------------------------------------------------------

StructureAdapterPtr createStructureAdapter(const boost::python::object stru)
{
    StructureAdapterPtr rv;
    // First check if stru is already a wrapped StructureAdapterPtr
    python::extract<StructureAdapterPtr> getadapter(stru);
    if (getadapter.check())
    {
        rv = getadapter();
        if (rv.get())  return rv;
    }
    // Loop over all registered factories.  They should return
    // a new StructureAdapter instance or NULL on failure.
    PSAFactorySet::iterator factory = getPSAFactorySet().begin();
    PSAFactorySet::iterator last = getPSAFactorySet().end();
    for (; factory != last; ++factory)
    {
        rv = (*factory)(stru);
        if (rv.get())  return rv;
    }
    // We get here only if nothing worked.
    python::object emsg = "Cannot create structure adapter for %r." %
        python::make_tuple(stru);
    PyErr_SetObject(PyExc_TypeError, emsg.ptr());
    python::throw_error_already_set();
    // This should be never reached
    return rv;
}

// Factory functions registration --------------------------------------------

bool registerPythonStructureAdapterFactory(PythonStructureAdapterFactory fctr)
{
    getPSAFactorySet().insert(fctr);
    return true;
}

}   // namespace srreal
}   // namespace diffpy

// End of file
