/*****************************************************************************
*
* diffpy.srreal     by DANSE Diffraction group
*                   Simon J. L. Billinge
*                   (c) 2011 Trustees of the Columbia University
*                   in the City of New York.  All rights reserved.
*
* File coded by:    Pavol Juhas
*
* See AUTHORS.txt for a list of people who contributed.
* See LICENSE.txt for license information.
*
******************************************************************************
*
* class VR3Adapter -- serialization support
*
* $Id$
*
*****************************************************************************/

#include <diffpy/srreal/VR3Structure.hpp>
#include <diffpy/serialization.hpp>

namespace diffpy {
namespace srreal {

//////////////////////////////////////////////////////////////////////////////
// class VR3StructureAdapter
//////////////////////////////////////////////////////////////////////////////

// Constructor ---------------------------------------------------------------

VR3Adapter::VR3Adapter(const VR3Structure& vr3s) : mvr3structure(vr3s)
{ }

// Public Methods ------------------------------------------------------------

int VR3Adapter::countSites() const
{
    return mvr3structure.size();
}


const R3::Vector& VR3Adapter::siteCartesianPosition(int idx) const
{
    return mvr3structure[idx];
}


bool VR3Adapter::siteAnisotropy(int idx) const
{
    return false;
}


const R3::Matrix& VR3Adapter::siteCartesianUij(int idx) const
{
    return R3::zeros();
}


BaseBondGenerator* VR3Adapter::createBondGenerator() const
{
    BaseBondGenerator* bnds = new VR3BondGenerator(this);
    return bnds;
}

//////////////////////////////////////////////////////////////////////////////
// class VR3BondGenerator
//////////////////////////////////////////////////////////////////////////////

// Constructor ---------------------------------------------------------------

VR3BondGenerator::VR3BondGenerator(const VR3Adapter* adpt) :
    BaseBondGenerator(adpt), mvr3structure(adpt->mvr3structure)
{ }

// Public Methods ------------------------------------------------------------

const R3::Vector& VR3BondGenerator::r0() const
{
    return mvr3structure.at(this->site0());
}


const R3::Vector& VR3BondGenerator::r1() const
{
    return mvr3structure.at(this->site1());
}

}   // namespace srreal
}   // namespace diffpy

// Serialization -------------------------------------------------------------

DIFFPY_INSTANTIATE_SERIALIZE(diffpy::srreal::VR3Adapter)
BOOST_CLASS_EXPORT(diffpy::srreal::VR3Adapter)

// End of file
