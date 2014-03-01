/*****************************************************************************
*
* diffpy.srreal     by DANSE Diffraction group
*                   Simon J. L. Billinge
*                   (c) 2012 The Trustees of Columbia University
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
*
*****************************************************************************/

#include <cassert>
#include <algorithm>
#include <boost/functional/hash.hpp>

#include <diffpy/serialization.ipp>
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
    R3::Vector::const_iterator p0 = a0.xyz_cartn.begin();
    R3::Vector::const_iterator p1 = a1.xyz_cartn.begin();
    for (; p0 != a0.xyz_cartn.end(); ++p0, ++p1)
    {
        if (*p0 < *p1)  return true;
        if (*p0 > *p1)  return false;
    }
    if (a0.occupancy < a1.occupancy)  return true;
    if (a0.occupancy > a1.occupancy)  return false;
    if (a0.anisotropy < a1.anisotropy)  return true;
    if (a0.anisotropy > a1.anisotropy)  return false;
    typedef R3::Matrix::array_type::const_iterator CIter;
    CIter u0 = a0.uij_cartn.data().begin();
    CIter u0last = a0.uij_cartn.data().end();
    CIter u1 = a1.uij_cartn.data().begin();
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
            (a0.xyz_cartn == a1.xyz_cartn) &&
            a0.occupancy == a1.occupancy &&
            a0.anisotropy == a1.anisotropy &&
            a0.uij_cartn == a1.uij_cartn
            );
    return rv;
}


bool operator!=(const Atom& a0, const Atom& a1)
{
    return !(a0 == a1);
}


size_t hash_value(const Atom& a)
{
    size_t seed = 0;
    boost::hash_combine(seed, a.atomtype);
    boost::hash_combine(seed, a.xyz_cartn);
    boost::hash_combine(seed, a.occupancy);
    boost::hash_combine(seed, a.anisotropy);
    boost::hash_combine(seed, a.uij_cartn);
    return seed;
}

//////////////////////////////////////////////////////////////////////////////
// class AtomicStructureAdapter
//////////////////////////////////////////////////////////////////////////////

// Public Methods ------------------------------------------------------------

StructureAdapterPtr AtomicStructureAdapter::clone() const
{
    StructureAdapterPtr rv(new AtomicStructureAdapter(*this));
    return rv;
}


BaseBondGeneratorPtr AtomicStructureAdapter::createBondGenerator() const
{
    BaseBondGeneratorPtr bnds(new BaseBondGenerator(shared_from_this()));
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
    return matoms[idx].xyz_cartn;
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
    return matoms[idx].uij_cartn;
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
    sd.pop0.reserve(astru0.countSites());
    sd.add1.reserve(astru1.countSites());
    // try fast side-by-side comparison for equal length structures
    if (astru0.countSites() == astru1.countSites())
    {
        sd.diffmethod = StructureDifference::Method::SIDEBYSIDE;
        const_iterator ai0 = astru0.matoms.begin();
        const_iterator ai1 = astru1.matoms.begin();
        for (int i = 0; ai0 != astru0.matoms.end(); ++i, ++ai0, ++ai1)
        {
            if (*ai0 != *ai1)
            {
                sd.pop0.push_back(i);
                sd.add1.push_back(i);
            }
            if (!sd.allowsfastupdate())  break;
        }
        if (sd.allowsfastupdate())  return sd;
    }
    // here the structures are either of different length or differ too much
    // when compared side by side.  Let's start from a blank slate.
    sd.pop0.clear();
    sd.add1.clear();
    // let's build sorted vectors of atoms in stru0 and stru1
    sd.diffmethod = StructureDifference::Method::SORTED;
    std::vector<atomindex> satoms0, satoms1;
    satoms0.reserve(astru0.countSites());
    const_iterator ai = astru0.matoms.begin();
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

typedef AtomicStructureAdapter::iterator iterator;

iterator AtomicStructureAdapter::insert(int idx, const Atom& atom)
{
    assert(0 <= idx && idx <= this->countSites());
    return this->insert(matoms.begin() + idx, atom);
}


iterator AtomicStructureAdapter::insert(iterator ii, const Atom& atom)
{
    return matoms.insert(ii, atom);
}


void AtomicStructureAdapter::append(const Atom& atom)
{
    matoms.push_back(atom);
}


void AtomicStructureAdapter::clear()
{
    matoms.clear();
}


iterator AtomicStructureAdapter::erase(int idx)
{
    assert(0 <= idx && idx < this->countSites());
    return matoms.erase(matoms.begin() + idx);
}


iterator AtomicStructureAdapter::erase(iterator pos)
{
    return matoms.erase(pos);
}


iterator AtomicStructureAdapter::erase(iterator first, iterator last)
{
    return matoms.erase(first, last);
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

}   // namespace srreal
}   // namespace diffpy

// Serialization -------------------------------------------------------------

DIFFPY_INSTANTIATE_SERIALIZATION(diffpy::srreal::Atom)
BOOST_CLASS_EXPORT_IMPLEMENT(diffpy::srreal::Atom)
DIFFPY_INSTANTIATE_SERIALIZATION(diffpy::srreal::AtomicStructureAdapter)
BOOST_CLASS_EXPORT_IMPLEMENT(diffpy::srreal::AtomicStructureAdapter)

// End of file
