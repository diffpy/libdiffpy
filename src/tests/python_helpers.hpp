/*****************************************************************************
*
* diffpy.lib        by DANSE Diffraction group
*                   Simon J. L. Billinge
*                   (c) 2010 Trustees of the Michigan State University.
*                   All rights reserved.
*
* File coded by:    Pavol Juhas
*
* See AUTHORS.txt for a list of people who contributed.
* See LICENSE.txt for license information.
*
******************************************************************************
*
* shared helpers for loading diffpy structure objects.
*
* $Id$
*
*****************************************************************************/

#ifndef PYTHON_HELPERS_HPP_INCLUDED
#define PYTHON_HELPERS_HPP_INCLUDED

#include <boost/python/object.hpp>
#include <string>

boost::python::object loadTestStructure(const std::string& tailname);

#endif  // PYTHON_HELPERS_HPP_INCLUDED
