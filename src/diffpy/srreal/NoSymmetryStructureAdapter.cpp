/*****************************************************************************
*
* diffpy.srreal     by DANSE Diffraction group
*                   Simon J. L. Billinge
*                   (c) 2010 Trustees of the Columbia University
*                   in the City of New York.  All rights reserved.
*
* File coded by:    Pavol Juhas
*
* See AUTHORS.txt for a list of people who contributed.
* See LICENSE.txt for license information.
*
******************************************************************************
*
* class NoSymmetryStructureAdapter -- StructureAdapter class that removes
*     any symmetry expansions (rotations or periodic translations) from
*     another StructureAdapter instance.  This can be used to use only
*     the asymmetric unit from any adapter to crystal structure.
*
* class NoSymmetryBondGenerator -- bond generator
*
* nosymmetry -- factory function that creates a NoSymmetryStructureAdapter
*     instance inside StructureAdapterPtr
*
* $Id$
*
*****************************************************************************/

#include <diffpy/serialization.hpp>
#include <diffpy/srreal/NoSymmetryStructureAdapter.hpp>

namespace diffpy {
namespace srreal {

//////////////////////////////////////////////////////////////////////////////
// class NoSymmetryStructureAdapter
//////////////////////////////////////////////////////////////////////////////

// Constructor ---------------------------------------------------------------

NoSymmetryStructureAdapter::NoSymmetryStructureAdapter(
        StructureAdapterPtr srcstructure)
{
    msrcstructure = srcstructure;
}

// Public Methods ------------------------------------------------------------

BaseBondGenerator* NoSymmetryStructureAdapter::createBondGenerator() const
{
    BaseBondGenerator* bnds = new BaseBondGenerator(this);
    return bnds;
}


int NoSymmetryStructureAdapter::countSites() const
{
    return msrcstructure->countSites();
}


double NoSymmetryStructureAdapter::numberDensity() const
{
    return 0.0;
}


const R3::Vector& NoSymmetryStructureAdapter::siteCartesianPosition(
        int idx) const
{
    return msrcstructure->siteCartesianPosition(idx);
}


double NoSymmetryStructureAdapter::siteOccupancy(int idx) const
{
    return msrcstructure->siteOccupancy(idx);
}


bool NoSymmetryStructureAdapter::siteAnisotropy(int idx) const
{
    return msrcstructure->siteAnisotropy(idx);
}


const R3::Matrix& NoSymmetryStructureAdapter::siteCartesianUij(int idx) const
{
    return msrcstructure->siteCartesianUij(idx);
}


const std::string& NoSymmetryStructureAdapter::siteAtomType(int idx) const
{
    return msrcstructure->siteAtomType(idx);
}


void NoSymmetryStructureAdapter::customPQConfig(PairQuantity* pq) const
{
    msrcstructure->customPQConfig(pq);
}

// Routines ------------------------------------------------------------------

StructureAdapterPtr nosymmetry(StructureAdapterPtr stru)
{
    StructureAdapterPtr rv =
        boost::dynamic_pointer_cast<NoSymmetryStructureAdapter>(stru) ? stru :
        StructureAdapterPtr(new NoSymmetryStructureAdapter(stru));
    return rv;
}

}   // namespace srreal
}   // namespace diffpy

DIFFPY_INSTANTIATE_SERIALIZE(diffpy::srreal::NoSymmetryStructureAdapter)
BOOST_CLASS_EXPORT(diffpy::srreal::NoSymmetryStructureAdapter)

// End of file
