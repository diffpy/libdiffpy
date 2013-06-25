/*****************************************************************************
*
* diffpy.srreal     by DANSE Diffraction group
*                   Simon J. L. Billinge
*                   (c) 2012 Trustees of the Columbia University
*                   in the City of New York.  All rights reserved.
*
* File coded by:    Pavol Juhas
*
* See AUTHORS.txt for a list of people who contributed.
* See LICENSE.txt for license information.
*
******************************************************************************
*
* class AtomicStructureAdapter -- universal structure adapter for
*     a non-periodic set of atoms.

* class AtomicStructureBondGenerator -- bond generator
*
*****************************************************************************/

#include <cassert>

#include <diffpy/serialization.hpp>
#include <diffpy/srreal/AtomicStructureAdapter.hpp>

using std::string;

namespace diffpy {
namespace srreal {

//////////////////////////////////////////////////////////////////////////////
// class AtomicStructureAdapter
//////////////////////////////////////////////////////////////////////////////

// Public Methods ------------------------------------------------------------

BaseBondGeneratorPtr AtomicStructureAdapter::createBondGenerator() const
{
    // FIXME: hack for handling non-periodic structures
    // should diffpy.Structure get an isPeriodic method?
    BaseBondGeneratorPtr bnds(
            new AtomicStructureBondGenerator(shared_from_this()));
    return bnds;
}


int AtomicStructureAdapter::countSites() const
{
    return matoms.size();
}


const R3::Vector& AtomicStructureAdapter::siteCartesianPosition(int idx) const
{
    assert(0 <= idx && idx < this->countSites());
    return matoms[idx].cartesianposition;
}


double AtomicStructureAdapter::siteOccupancy(int idx) const
{
    assert(0 <= idx && idx < this->countSites());
    return matoms[idx].occupancy;
}


bool AtomicStructureAdapter::siteAnisotropy(int idx) const
{
    assert(0 <= idx && idx < this->countSites());
    return matoms[idx].anisotropy;
}


const R3::Matrix& AtomicStructureAdapter::siteCartesianUij(int idx) const
{
    assert(0 <= idx && idx < this->countSites());
    return matoms[idx].cartesianuij;
}


const string& AtomicStructureAdapter::siteAtomType(int idx) const
{
    assert(0 <= idx && idx < this->countSites());
    return matoms[idx].atomtype;
}


void AtomicStructureAdapter::insert(int idx, const AtomAdapter& atom)
{
    assert(0 <= idx && idx <= this->countSites());
    matoms.insert(matoms.begin() + idx, atom);
}


void AtomicStructureAdapter::append(const AtomAdapter& atom)
{
    matoms.push_back(atom);
}


void AtomicStructureAdapter::remove(int idx)
{
    assert(0 <= idx && idx < this->countSites());
    matoms.erase(matoms.begin() + idx);
}

//////////////////////////////////////////////////////////////////////////////
// class AtomicStructureBondGenerator
//////////////////////////////////////////////////////////////////////////////

// Constructor ---------------------------------------------------------------

AtomicStructureBondGenerator::AtomicStructureBondGenerator(
        StructureAdapterConstPtr adpt) : BaseBondGenerator(adpt)
{
    mastructure = dynamic_cast<const AtomicStructureAdapter*>(adpt.get());
    assert(mastructure);
}

}   // namespace srreal
}   // namespace diffpy

// Serialization -------------------------------------------------------------

DIFFPY_INSTANTIATE_SERIALIZE(diffpy::srreal::AtomicStructureAdapter)
BOOST_CLASS_EXPORT(diffpy::srreal::AtomicStructureAdapter)

// End of file
