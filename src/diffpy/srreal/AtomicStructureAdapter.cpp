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
#include <boost/functional/hash.hpp>

#include <diffpy/serialization.hpp>
#include <diffpy/srreal/AtomicStructureAdapter.hpp>

using std::string;

namespace diffpy {
namespace srreal {

//////////////////////////////////////////////////////////////////////////////
// class Atom
//////////////////////////////////////////////////////////////////////////////

bool operator< (const Atom& a0, const Atom& a1)
{
    if (a0.atomtype < a1.atomtype)  return true;
    if (a0.atomtype > a1.atomtype)  return false;
    R3::Vector::const_iterator p0 = a0.cartesianposition.begin();
    R3::Vector::const_iterator p1 = a1.cartesianposition.begin();
    for (; p0 != a0.cartesianposition.end(); ++p0, ++p1)
    {
        if (*p0 < *p1)  return true;
        if (*p0 > *p1)  return false;
    }
    if (a0.occupancy < a1.occupancy)  return true;
    if (a0.occupancy > a1.occupancy)  return false;
    if (a0.anisotropy < a1.anisotropy)  return true;
    if (a0.anisotropy > a1.anisotropy)  return false;
    const double* u0 = a0.cartesianposition.dataFirst();
    const double* u0last = u0 + R3::Ndim * R3::Ndim;
    const double* u1 = a1.cartesianposition.dataFirst();
    for (; u0 != u0last; ++u0, ++u1)
    {
        if (*u0 < *u1)  return true;
        if (*u0 > *u1)  return false;
    }
    return false;
}


bool operator==(const Atom& a0, const Atom& a1)
{
    bool rv = (&a0 == &a1) || (
            a0.atomtype == a1.atomtype &&
            (a0.cartesianposition == a1.cartesianposition) &&
            a0.occupancy == a1.occupancy &&
            a0.anisotropy == a1.anisotropy &&
            a0.cartesianuij == a1.cartesianuij
            );
    return rv;
}


size_t hash_value(const Atom& a)
{
    size_t seed = 0;
    boost::hash_combine(seed, a.atomtype);
    boost::hash_combine(seed, a.cartesianposition);
    boost::hash_combine(seed, a.occupancy);
    boost::hash_combine(seed, a.anisotropy);
    boost::hash_combine(seed, a.cartesianuij);
    return seed;
}

//////////////////////////////////////////////////////////////////////////////
// class AtomicStructureAdapter
//////////////////////////////////////////////////////////////////////////////

// Public Methods ------------------------------------------------------------

BaseBondGeneratorPtr AtomicStructureAdapter::createBondGenerator() const
{
    BaseBondGeneratorPtr bnds(
            new AtomicStructureBondGenerator(shared_from_this()));
    return bnds;
}


int AtomicStructureAdapter::countSites() const
{
    return matoms.size();
}


const string& AtomicStructureAdapter::siteAtomType(int idx) const
{
    assert(0 <= idx && idx < this->countSites());
    return matoms[idx].atomtype;
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


void AtomicStructureAdapter::insert(int idx, const Atom& atom)
{
    assert(0 <= idx && idx <= this->countSites());
    matoms.insert(matoms.begin() + idx, atom);
}


void AtomicStructureAdapter::append(const Atom& atom)
{
    matoms.push_back(atom);
}


void AtomicStructureAdapter::remove(int idx)
{
    assert(0 <= idx && idx < this->countSites());
    matoms.erase(matoms.begin() + idx);
}


Atom& AtomicStructureAdapter::operator[](int idx)
{
    assert(0 <= idx && idx < this->countSites());
    return matoms[idx];
}


const Atom& AtomicStructureAdapter::operator[](int idx) const
{
    assert(0 <= idx && idx < this->countSites());
    return matoms[idx];
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
