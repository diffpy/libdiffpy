/*****************************************************************************
*
* diffpy.srreal     Complex Modeling Initiative
*                   Pavol Juhas
*                   (c) 2013 Brookhaven National Laboratory,
*                   Upton, New York.  All rights reserved.
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

BaseBondGeneratorPtr CrystalStructureAdapter::createBondGenerator() const
{
    this->updateSymmetryPositions();
    BaseBondGeneratorPtr bnds(
            new CrystalStructureBondGenerator(shared_from_this()));
    return bnds;
}


int CrystalStructureAdapter::siteMultiplicity(int idx) const
{
    if (!msymmetry_cached)  this->updateSymmetryPositions();
    int rv = msymatoms[idx].size();
    return rv;
}


StructureDifference
CrystalStructureAdapter::diff(StructureAdapterConstPtr other) const
{
    StructureDifference alldiffer;
    if (!other || typeid(CrystalStructureAdapter) != typeid(*other))
    {
        return alldiffer;
    }
    StructureDifference sd(this->shared_from_this(), other);
    if (sd.stru0 == sd.stru1)  return sd;
    // compare symmetry operations in both adapters
    typedef boost::shared_ptr<const class CrystalStructureAdapter> CPtr;
    CPtr cstru1 = boost::static_pointer_cast<CPtr::element_type>(sd.stru1);
    assert(cstru1 == sd.stru1);
    if (msymops != cstru1->msymops)  return alldiffer;
    // downcast to allow comparison with PeriodicStructureAdapter::diff
    typedef boost::shared_ptr<const class PeriodicStructureAdapter> PPtr;
    PPtr pstru1 = boost::static_pointer_cast<PPtr::element_type>(sd.stru1);
    assert(typeid(PeriodicStructureAdapter) == typeid(*pstru1));
    sd = this->PeriodicStructureAdapter::diff(pstru1);
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


CrystalStructureAdapter::AtomVector
CrystalStructureAdapter::expandLatticeAtom(const Atom& a0) const
{
    using mathutils::eps_eq;
    AtomVector eqsites;
    vector<int> eqduplicity;
    vector<SymOpRotTrans>::const_iterator op = msymops.begin();
    Atom a1 = a0;
    for (; op != msymops.end(); ++op)
    {
        // positions and Uij-s are actually fractional here
        a1.cartesianposition = R3::mxvecproduct(op->R, a0.cartesianposition);
        a1.cartesianposition += op->t;
        R3::Matrix utmp = R3::prod(a0.cartesianuij, R3::trans(op->R));
        a1.cartesianuij = R3::prod(op->R, utmp);
        // check if a1 is a duplicate of an existing symmetry site
        int ieq = this->findEqualPosition(eqsites, a1);
        if (ieq >= 0)
        {
            Atom& aeq = eqsites[ieq];
            aeq.cartesianposition += a1.cartesianposition;
            aeq.cartesianuij += a1.cartesianuij;
            ++eqduplicity[ieq];
            continue;
        }
        // a1 is a new symmetry site
        eqsites.push_back(a1);
        eqduplicity.push_back(1);
    }
    // calculate mean values from equivalent sites and adjust any roundoffs
    assert(eqsites.size() == eqduplicity.size());
    iterator ai = eqsites.begin();
    vector<int>::const_iterator dpi = eqduplicity.begin();
    for (; ai != eqsites.end(); ++ai, ++dpi)
    {
        R3::Vector& axyz = ai->cartesianposition;
        axyz /= (*dpi);
        axyz -= floor(axyz);
        if (eps_eq(axyz[0], 1.0))  axyz[0] = 0.0;
        if (eps_eq(axyz[1], 1.0))  axyz[1] = 0.0;
        if (eps_eq(axyz[2], 1.0))  axyz[2] = 0.0;
        ai->cartesianuij /= (*dpi);
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
        dxyz = ai->cartesianposition - a0.cartesianposition;
        dxyz[0] -= round(dxyz[0]);
        dxyz[1] -= round(dxyz[1]);
        dxyz[2] -= round(dxyz[2]);
        if (L.norm(dxyz) <= symeps)  return (ai - eqsites.begin());
    }
    return -1;
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
}

// Public Methods ------------------------------------------------------------

void CrystalStructureBondGenerator::rewind()
{
/* FIXME
    // Delay msphere instantiation to here instead of in constructor,
    // so it is possible to use setRmin, setRmax.
    if (!msphere.get())
    {
        const Lattice& L = mpstructure->getLattice();
        double buffzone = L.ucMaxDiagonalLength();
        double rsphmin = this->getRmin() - buffzone;
        double rsphmax = this->getRmax() + buffzone;
        msphere.reset(new PointsInSphere(rsphmin, rsphmax, L));
    }
    // BaseBondGenerator::rewind calls this->rewindSymmetry,
    // which takes care of msphere configuration
    this->BaseBondGenerator::rewind();
*/
}


const R3::Matrix& CrystalStructureBondGenerator::Ucartesian1() const
{
// FIXME
    return R3::zeromatrix();
}


void CrystalStructureBondGenerator::setRmin(double rmin)
{
/* FIXME
    // destroy msphere so it will be created on rewind with new rmin
    if (this->getRmin() != rmin)    msphere.reset();
    this->BaseBondGenerator::setRmin(rmin);
*/
}


void CrystalStructureBondGenerator::setRmax(double rmax)
{
/* FIXME
    // destroy msphere so it will be created on rewind with new rmax
    if (this->getRmax() != rmax)    msphere.reset();
    this->BaseBondGenerator::setRmax(rmax);
*/
}

// Protected Methods ---------------------------------------------------------

bool CrystalStructureBondGenerator::iterateSymmetry()
{
/* FIXME
    msphere->next();
    bool done = msphere->finished();
    mrcsphere = done ? R3::zerovector :
        mpstructure->getLattice().cartesian(msphere->mno());
    return !done;
*/
    return true;
}


void CrystalStructureBondGenerator::rewindSymmetry()
{
/* FIXME
    msphere->rewind();
    mrcsphere = msphere->finished() ? R3::zerovector :
        mpstructure->getLattice().cartesian(msphere->mno());
    this->updater1();
*/
}


// Private Methods -----------------------------------------------------------

/* FIXME
void CrystalStructureBondGenerator::updater1()
{
    mr1 = mrcsphere + mcartesian_positions_uc[this->site1()];
    this->updateDistance();
}
*/

}   // namespace srreal
}   // namespace diffpy

// Serialization -------------------------------------------------------------

DIFFPY_INSTANTIATE_SERIALIZATION(diffpy::srreal::CrystalStructureAdapter)
BOOST_CLASS_EXPORT_IMPLEMENT(diffpy::srreal::CrystalStructureAdapter)

// End of file
