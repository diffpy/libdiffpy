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
*****************************************************************************/

#include <diffpy/srreal/VR3Structure.hpp>
#include <diffpy/serialization.ipp>

namespace diffpy {
namespace srreal {

//////////////////////////////////////////////////////////////////////////////
// class VR3StructureAdapter
//////////////////////////////////////////////////////////////////////////////

// Constructor ---------------------------------------------------------------

VR3Adapter::VR3Adapter(const VR3Structure& vr3s) : mvr3structure(vr3s)
{ }

// Public Methods ------------------------------------------------------------

StructureAdapterPtr VR3Adapter::clone() const
{
    // VR3Adapter has only constant public methods
    StructureAdapterPtr rv = boost::const_pointer_cast<
        StructureAdapterPtr::element_type>(shared_from_this());
    return rv;
}


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
    return R3::zeromatrix();
}


BaseBondGeneratorPtr VR3Adapter::createBondGenerator() const
{
    BaseBondGeneratorPtr bnds(new VR3BondGenerator(shared_from_this()));
    return bnds;
}

//////////////////////////////////////////////////////////////////////////////
// class VR3BondGenerator
//////////////////////////////////////////////////////////////////////////////

// Constructor ---------------------------------------------------------------

VR3BondGenerator::VR3BondGenerator(StructureAdapterConstPtr adpt) :
    BaseBondGenerator(adpt),
    mvr3structure(dynamic_cast<const VR3Adapter*>(adpt.get())->mvr3structure)
{ }

}   // namespace srreal
}   // namespace diffpy

// Serialization -------------------------------------------------------------

DIFFPY_INSTANTIATE_SERIALIZATION(diffpy::srreal::VR3Adapter)
BOOST_CLASS_EXPORT_IMPLEMENT(diffpy::srreal::VR3Adapter)

// End of file
