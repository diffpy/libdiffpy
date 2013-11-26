/*****************************************************************************
*
* diffpy.srreal     by DANSE Diffraction group
*                   Simon J. L. Billinge
*                   (c) 2009 Trustees of the Columbia University
*                   in the City of New York.  All rights reserved.
*
* File coded by:    Chris Farrow
*
* See AUTHORS.txt for a list of people who contributed.
* See LICENSE.txt for license information.
*
******************************************************************************
*
* class ObjCrystStructureAdapter
*   -- adapter to the Crystal class from ObjCryst++.
* class ObjCrystBondGenerator
*   -- Generate bonds from periodic ObjCrystStructureAdapter.
*
*****************************************************************************/

#ifndef OBJCRYSTSTRUCTUREADAPTER_HPP_INCLUDED
#define OBJCRYSTSTRUCTUREADAPTER_HPP_INCLUDED

#include <diffpy/srreal/AtomicStructureAdapter.hpp>
#include <diffpy/srreal/CrystalStructureAdapter.hpp>
#include <ObjCryst/ObjCryst/Crystal.h>
#include <ObjCryst/ObjCryst/Molecule.h>

namespace diffpy {
namespace srreal {

// ObjCryst::Crystal is now adapted with CrystalStructureAdapter

StructureAdapterPtr
createStructureAdapter(const ObjCryst::Crystal& cryst);

// ObjCryst::Molecule can be adapted with AtomicStructureAdapter
//
// Molecules are always considered aperiodic. The anisotropic ADPs are treated
// as if in a cartesian cell. If this is not what is intended, pass the
// molecule as a scattering component within an ObjCryst::Crystal.


StructureAdapterPtr
createStructureAdapter(const ObjCryst::Molecule& molecule);

}   // namespace srreal
}   // namespace diffpy

#endif  // OBJCRYSTSTRUCTUREADAPTER_HPP_INCLUDED
