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

#include <stdexcept>
#include <sstream>
#include <string>
#include <set>

#include <diffpy/srreal/PythonStructureAdapter.hpp>

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
        return rv;
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
    python::object pytype = python::import("__builtin__").attr("type");
    python::object tp = python::str(pytype(stru));
    ostringstream emsg;
    emsg << "Cannot create structure adapter for Python " <<
        string(python::extract<string>(tp)) << ".";
    throw invalid_argument(emsg.str());
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
