/*****************************************************************************
*
* libdiffpy         Complex Modeling Initiative
*                   (c) 2013 Brookhaven Science Associates,
*                   Brookhaven National Laboratory.
*                   All rights reserved.
*
* File coded by:    Pavol Juhas
*
* See AUTHORS.txt for a list of people who contributed.
* See LICENSE.txt for license information.
*
******************************************************************************
*
* class CrystalStructureAdapter -- universal adapter for crystal structure
*     composed of asymmetric unit cell and list of symmetry operations in
*     the space group
*
* class CrystalStructureBondGenerator -- bond generator
*
*****************************************************************************/

#include <cassert>

#include <diffpy/serialization.ipp>
#include <diffpy/validators.hpp>
#include <diffpy/srreal/PointsInSphere.hpp>
#include <diffpy/srreal/StructureDifference.hpp>
#include <diffpy/srreal/CrystalStructureAdapter.hpp>

using namespace std;

namespace diffpy {
namespace srreal {

// Constants -----------------------------------------------------------------

const double DEFAULT_SYMMETRY_PRECISION = 5e-5;

//////////////////////////////////////////////////////////////////////////////
// class CrystalStructureAdapter
//////////////////////////////////////////////////////////////////////////////

// Constructor ---------------------------------------------------------------

CrystalStructureAdapter::CrystalStructureAdapter() :
    PeriodicStructureAdapter(),
    msymmetry_cached(false)
{
    this->setSymmetryPrecision(DEFAULT_SYMMETRY_PRECISION);
}

// Public Methods ------------------------------------------------------------

StructureAdapterPtr CrystalStructureAdapter::clone() const
{
    StructureAdapterPtr rv(new CrystalStructureAdapter(*this));
    return rv;
}


BaseBondGeneratorPtr CrystalStructureAdapter::createBondGenerator() const
{
    this->updateSymmetryPositions();
    BaseBondGeneratorPtr bnds(
            new CrystalStructureBondGenerator(shared_from_this()));
    return bnds;
}


int CrystalStructureAdapter::siteMultiplicity(int idx) const
{
    if (!this->isSymmetryCached())  this->updateSymmetryPositions();
    int rv = msymatoms[idx].size();
    return rv;
}


StructureDifference
CrystalStructureAdapter::diff(StructureAdapterConstPtr other) const
{
    StructureDifference sd = this->StructureAdapter::diff(other);
    if (sd.stru0 == sd.stru1)  return sd;
    typedef boost::shared_ptr<const class CrystalStructureAdapter> CPtr;
    CPtr cother = boost::dynamic_pointer_cast<CPtr::element_type>(other);
    if (!cother)  return sd;
    // compare symmetry operations in both adapters
    assert(cother == sd.stru1);
    if (msymops != cother->msymops)  return sd;
    sd = this->PeriodicStructureAdapter::diff(other);
    return sd;
}


void CrystalStructureAdapter::setSymmetryPrecision(double eps)
{
    using namespace diffpy::validators;
    ensureEpsilonPositive("symmetryprecision", eps);
    if (eps != msymmetry_precision)  msymmetry_cached = false;
    msymmetry_precision = eps;
}


const double& CrystalStructureAdapter::getSymmetryPrecision() const
{
    return msymmetry_precision;
}


int CrystalStructureAdapter::countSymOps() const
{
    return msymops.size();
}


void CrystalStructureAdapter::clearSymOps()
{
    msymops.clear();
    msymmetry_cached = false;
}


void CrystalStructureAdapter::addSymOp(const SymOpRotTrans& op)
{
    msymops.push_back(op);
    msymmetry_cached = false;
}


void CrystalStructureAdapter::addSymOp(const R3::Matrix& R, const R3::Vector& t)
{
    SymOpRotTrans op = {R, t};
    this->addSymOp(op);
    assert(!msymmetry_cached);
}


const SymOpRotTrans& CrystalStructureAdapter::getSymOp(int i) const
{
    return msymops[i];
}


const CrystalStructureAdapter::AtomVector&
CrystalStructureAdapter::getEquivalentAtoms(int idx) const
{
    assert(0 <= idx && idx < this->countSites());
    if (!this->isSymmetryCached())  this->updateSymmetryPositions();
    return msymatoms[idx];
}


CrystalStructureAdapter::AtomVector
CrystalStructureAdapter::expandLatticeAtom(const Atom& a0) const
{
    using mathutils::eps_eq;
    AtomVector eqsites;
    vector<R3::Vector> eqsumpos;
    vector<int> eqduplicity;
    eqsumpos.reserve(this->countSymOps());
    eqduplicity.reserve(this->countSymOps());
    const Lattice& L = this->getLattice();
    SymOpVector::const_iterator op = msymops.begin();
    Atom a1 = a0;
    for (; op != msymops.end(); ++op)
    {
        // positions and Uij-s are actually fractional here
        a1.xyz_cartn = R3::mxvecproduct(op->R, a0.xyz_cartn);
        a1.xyz_cartn += op->t;
        // check if a1 is a duplicate of an existing symmetry site
        int ieq = this->findEqualPosition(eqsites, a1);
        if (ieq < 0)
        {
            // a1 is a new symmetry site
            R3::Matrix utmp = R3::prod(a0.uij_cartn, R3::trans(op->R));
            a1.uij_cartn = R3::prod(op->R, utmp);
            eqsites.push_back(a1);
            eqsumpos.push_back(R3::zerovector);
            eqduplicity.push_back(0);
            ieq = eqsites.size() - 1;
        }
        eqsumpos[ieq] += L.ucvFractional(a1.xyz_cartn);
        eqduplicity[ieq] += 1;
    }
    // assume P1 if symmetry operations were not defined
    if (msymops.empty())
    {
        assert(eqsites.empty());
        assert(eqsumpos.empty());
        assert(eqduplicity.empty());
        eqsites.push_back(a0);
        eqsumpos.push_back(a0.xyz_cartn);
        eqduplicity.push_back(1);
    }
    // calculate mean values from equivalent sites and adjust any roundoffs
    assert(eqsites.size() == eqduplicity.size());
    assert(eqsites.size() == eqsumpos.size());
    iterator ai = eqsites.begin();
    vector<R3::Vector>::const_iterator sii = eqsumpos.begin();
    vector<int>::const_iterator dpi = eqduplicity.begin();
    for (; ai != eqsites.end(); ++ai, ++sii, ++dpi)
    {
        ai->xyz_cartn = (*sii) / (*dpi);
    }
    return eqsites;
}


void CrystalStructureAdapter::updateSymmetryPositions() const
{
    // build asymmetric unit in lattice coordinates
    AtomVector lcatoms(this->begin(), this->end());
    AtomVector::iterator lcai = lcatoms.begin();
    for (; lcai != lcatoms.end(); ++lcai)  this->toFractional(*lcai);
    // build symmetry positions for all atoms in the asymmetric unit
    msymatoms.resize(this->countSites());
    assert(lcatoms.size() == msymatoms.size());
    lcai = lcatoms.begin();
    std::vector<AtomVector>::iterator saii = msymatoms.begin();
    for (; lcai != lcatoms.end(); ++lcai, ++saii)
    {
        *saii = this->expandLatticeAtom(*lcai);
        iterator ai = saii->begin();
        for (; ai != saii->end(); ++ai)  this->toCartesian(*ai);
    }
    msymmetry_cached = true;
}

// Private Methods -----------------------------------------------------------

int CrystalStructureAdapter::findEqualPosition(
        const AtomVector& eqsites, const Atom& a0) const
{
    const double symeps = this->getSymmetryPrecision();
    const Lattice& L = this->getLattice();
    R3::Vector dxyz;
    const_iterator ai = eqsites.begin();
    for (; ai != eqsites.end(); ++ai)
    {
        dxyz = ai->xyz_cartn - a0.xyz_cartn;
        dxyz[0] -= round(dxyz[0]);
        dxyz[1] -= round(dxyz[1]);
        dxyz[2] -= round(dxyz[2]);
        if (L.norm(dxyz) <= symeps)  return (ai - eqsites.begin());
    }
    return -1;
}

bool CrystalStructureAdapter::isSymmetryCached() const
{
    msymmetry_cached = msymmetry_cached &&
        (int(msymatoms.size()) == this->countSites());
    return msymmetry_cached;
}

// Comparison functions ------------------------------------------------------

bool operator==(
        const CrystalStructureAdapter& stru0,
        const CrystalStructureAdapter& stru1)
{
    const PeriodicStructureAdapter& pstru0 = stru0;
    const PeriodicStructureAdapter& pstru1 = stru1;
    bool rv = (&stru0 == &stru1) || (
        (pstru0 == pstru1) &&
        (stru0.msymops == stru1.msymops));
    return rv;
}


bool operator!=(
        const CrystalStructureAdapter& stru0,
        const CrystalStructureAdapter& stru1)
{
    return !(stru0 == stru1);
}

//////////////////////////////////////////////////////////////////////////////
// class CrystalStructureBondGenerator
//////////////////////////////////////////////////////////////////////////////

// Constructor ---------------------------------------------------------------

CrystalStructureBondGenerator::CrystalStructureBondGenerator(
        StructureAdapterConstPtr adpt) : PeriodicStructureBondGenerator(adpt)
{
    mcstructure = dynamic_cast<const CrystalStructureAdapter*>(adpt.get());
    assert(mcstructure);
    msymidx = 0;
    mpuc1 = &(R3::zeromatrix());
}

// Public Methods ------------------------------------------------------------

void CrystalStructureBondGenerator::selectAnchorSite(int anchor)
{
    this->BaseBondGenerator::selectAnchorSite(anchor);
    const Atom& a0 = this->symatoms(anchor)[0];
    mr0 = a0.xyz_cartn;
}


const R3::Matrix& CrystalStructureBondGenerator::Ucartesian1() const
{
    return *mpuc1;
}

// Protected Methods ---------------------------------------------------------

bool CrystalStructureBondGenerator::iterateSymmetry()
{
    // Iterate the sphere at a fixed symmetry position.
    if (this->PeriodicStructureBondGenerator::iterateSymmetry())
    {
        this->updater1();
        return true;
    }
    // Advance to the next symmetry position.  We are done if they
    // were all already used.
    const AtomVector& sa = this->symatoms(this->site1());
    ++msymidx;
    if (msymidx >= sa.size())  return false;
    // rewind the sphere for the new symmetry position.
    // this calls CrystalStructureBondGenerator::updater1()
    this->PeriodicStructureBondGenerator::rewindSymmetry();
    return true;
}


void CrystalStructureBondGenerator::rewindSymmetry()
{
    msymidx = 0;
    // this calls CrystalStructureBondGenerator::updater1()
    this->PeriodicStructureBondGenerator::rewindSymmetry();
}


void CrystalStructureBondGenerator::getNextBond()
{
    this->BaseBondGenerator::getNextBond();
}


void CrystalStructureBondGenerator::updater1()
{
    const AtomVector& sa = this->symatoms(this->site1());
    assert(msymidx < sa.size());
    mr1 = mrcsphere + sa[msymidx].xyz_cartn;
    mpuc1 = &(sa[msymidx].uij_cartn);
    this->updateDistance();
}

// Private Methods -----------------------------------------------------------

const CrystalStructureAdapter::AtomVector&
CrystalStructureBondGenerator::symatoms(int idx)
{
    assert(0 <= idx && idx < int(mcstructure->msymatoms.size()));
    return mcstructure->msymatoms[idx];
}

}   // namespace srreal
}   // namespace diffpy

// Serialization -------------------------------------------------------------

DIFFPY_INSTANTIATE_SERIALIZATION(diffpy::srreal::CrystalStructureAdapter)
BOOST_CLASS_EXPORT_IMPLEMENT(diffpy::srreal::CrystalStructureAdapter)

// End of file
