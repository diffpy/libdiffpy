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

#ifndef PYTHONSTRUCTUREADAPTER_HPP_INCLUDED
#define PYTHONSTRUCTUREADAPTER_HPP_INCLUDED

#include <boost/python.hpp>
#include <diffpy/srreal/StructureAdapter.hpp>

namespace diffpy {
namespace srreal {

/// Factory for constructing appropriate StructureAdapter for a Python object.
StructureAdapterPtr createStructureAdapter(const boost::python::object& stru);

/// PythonStructureAdapterFactory type is a function that can
/// create StructureAdapter from boost::python object
typedef StructureAdapterPtr
    (*PythonStructureAdapterFactory)(const boost::python::object&);

/// This function registers a concrete structure adapter factory
bool registerPythonStructureAdapterFactory(PythonStructureAdapterFactory);

}   // namespace srreal
}   // namespace diffpy

#endif  // PYTHONSTRUCTUREADAPTER_HPP_INCLUDED
