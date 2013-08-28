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
#include <algorithm>
#include <boost/functional/hash.hpp>

#include <diffpy/serialization.hpp>
#include <diffpy/srreal/AtomicStructureAdapter.hpp>
#include <diffpy/srreal/StructureDifference.hpp>

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
    const double* u0 = a0.cartesianuij.dataFirst();
    const double* u0last = u0 + R3::Ndim * R3::Ndim;
    const double* u1 = a1.cartesianuij.dataFirst();
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

// helper for diff
namespace {

typedef std::pair<const Atom*, int> atomindex;

bool cmpatomindex(const atomindex& ai0, const atomindex& ai1)
{
    return (*(ai0.first) < *(ai1.first));
}

}   // namespace

StructureDifference
AtomicStructureAdapter::diff(StructureAdapterConstPtr other) const
{
    using std::min;
    typedef boost::shared_ptr<const class AtomicStructureAdapter> Ptr;
    StructureDifference sd(this->shared_from_this(), other);
    if (sd.stru0 == sd.stru1)  return sd;
    Ptr pstru0 = boost::static_pointer_cast<Ptr::element_type>(sd.stru0);
    Ptr pstru1 = boost::dynamic_pointer_cast<Ptr::element_type>(sd.stru1);
    if (!pstru1)  return this->StructureAdapter::diff(other);
    assert(pstru0);
    assert(pstru1);
    const AtomicStructureAdapter& astru0 = *pstru0;
    const AtomicStructureAdapter& astru1 = *pstru1;
    // build sorted vectors of atoms in stru0 and stru1
    std::vector<atomindex> satoms0, satoms1;
    satoms0.reserve(astru0.countSites());
    std::vector<Atom>::const_iterator ai = astru0.matoms.begin();
    for (int i = 0; ai != astru0.matoms.end(); ++ai, ++i)
    {
        satoms0.push_back(atomindex(&(*ai), i));
    }
    // use negative index for stru1 atoms so we can tell them apart
    // in the output of set_symmetric_difference
    satoms1.reserve(astru1.countSites());
    ai = astru1.matoms.begin();
    for (int i = -1; ai != astru1.matoms.end(); ++ai, --i)
    {
        satoms1.push_back(atomindex(&(*ai), i));
    }
    sort(satoms0.begin(), satoms0.end(), cmpatomindex);
    sort(satoms1.begin(), satoms1.end(), cmpatomindex);
    std::vector<atomindex> symdiffatoms(satoms0.size() + satoms1.size());
    std::vector<atomindex>::iterator ii;
    ii = std::set_symmetric_difference(satoms0.begin(), satoms0.end(),
            satoms1.begin(), satoms1.end(),
            symdiffatoms.begin(), cmpatomindex);
    symdiffatoms.erase(ii, symdiffatoms.end());
    sd.pop0.reserve(min(satoms0.size(), symdiffatoms.size()));
    sd.add1.reserve(min(satoms1.size(), symdiffatoms.size()));
    for (ii = symdiffatoms.begin(); ii != symdiffatoms.end(); ++ii)
    {
        if (ii->second >= 0)  sd.pop0.push_back(ii->second);
        else  sd.add1.push_back(-1 * ii->second - 1);
    }
    assert(sd.pop0.size() <= min(satoms0.size(), symdiffatoms.size()));
    assert(sd.add1.size() <= min(satoms1.size(), symdiffatoms.size()));
    sort(sd.pop0.begin(), sd.pop0.end());
    sort(sd.add1.begin(), sd.add1.end());
    return sd;
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
