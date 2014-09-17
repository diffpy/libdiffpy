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
* class PeriodicStructureAdapter -- universal adapter for structure with
*     periodic boundary conditions that has no space group symmetry
*
* class PeriodicStructureBondGenerator -- bond generator
*
*****************************************************************************/

#include <cassert>

#include <diffpy/serialization.ipp>
#include <diffpy/srreal/PointsInSphere.hpp>
#include <diffpy/srreal/StructureDifference.hpp>
#include <diffpy/srreal/PeriodicStructureAdapter.hpp>

using namespace std;

namespace diffpy {
namespace srreal {

//////////////////////////////////////////////////////////////////////////////
// class PeriodicStructureAdapter
//////////////////////////////////////////////////////////////////////////////

// Public Methods ------------------------------------------------------------

StructureAdapterPtr PeriodicStructureAdapter::clone() const
{
    StructureAdapterPtr rv(new PeriodicStructureAdapter(*this));
    return rv;
}


BaseBondGeneratorPtr PeriodicStructureAdapter::createBondGenerator() const
{
    BaseBondGeneratorPtr bnds(
            new PeriodicStructureBondGenerator(shared_from_this()));
    return bnds;
}


double PeriodicStructureAdapter::numberDensity() const
{
    const Lattice& L = this->getLattice();
    double rv = this->totalOccupancy() / L.volume();
    return rv;
}


StructureDifference
PeriodicStructureAdapter::diff(StructureAdapterConstPtr other) const
{
    StructureDifference sd = this->StructureAdapter::diff(other);
    if (sd.stru0 == sd.stru1)  return sd;
    typedef boost::shared_ptr<const class PeriodicStructureAdapter> PPtr;
    PPtr pother = boost::dynamic_pointer_cast<PPtr::element_type>(other);
    if (!pother)  return sd;
    assert(pother == sd.stru1);
    if (this->getLattice() != pother->getLattice())  return sd;
    sd = this->AtomicStructureAdapter::diff(other);
    return sd;
}


void PeriodicStructureAdapter::setLatPar(
        double a, double b, double c,
        double alphadeg, double betadeg, double gammadeg)
{
    mlattice.setLatPar(a, b, c, alphadeg, betadeg, gammadeg);
}


const Lattice& PeriodicStructureAdapter::getLattice() const
{
    return mlattice;
}


void PeriodicStructureAdapter::toCartesian(Atom& a) const
{
    const Lattice& L = this->getLattice();
    a.xyz_cartn = L.cartesian(a.xyz_cartn);
    a.uij_cartn = L.cartesianMatrix(a.uij_cartn);
}


void PeriodicStructureAdapter::toFractional(Atom& a) const
{
    const Lattice& L = this->getLattice();
    a.xyz_cartn = L.fractional(a.xyz_cartn);
    a.uij_cartn = L.fractionalMatrix(a.uij_cartn);
}

// Comparison functions ------------------------------------------------------

bool operator==(
        const PeriodicStructureAdapter& stru0,
        const PeriodicStructureAdapter& stru1)
{
    const AtomicStructureAdapter& astru0 = stru0;
    const AtomicStructureAdapter& astru1 = stru1;
    bool rv =
        (astru0 == astru1) &&
        (stru0.getLattice() == stru1.getLattice());
    return rv;
}


bool operator!=(
        const PeriodicStructureAdapter& stru0,
        const PeriodicStructureAdapter& stru1)
{
    return !(stru0 == stru1);
}

//////////////////////////////////////////////////////////////////////////////
// class PeriodicStructureBondGenerator
//////////////////////////////////////////////////////////////////////////////

// Constructor ---------------------------------------------------------------

PeriodicStructureBondGenerator::PeriodicStructureBondGenerator(
        StructureAdapterConstPtr adpt) : BaseBondGenerator(adpt)
{
    mpstructure = dynamic_cast<const PeriodicStructureAdapter*>(adpt.get());
    assert(mpstructure);
    int cntsites = mpstructure->countSites();
    mcartesian_positions_uc.reserve(cntsites);
    const Lattice& L = mpstructure->getLattice();
    PeriodicStructureAdapter::const_iterator ai = mpstructure->begin();
    R3::Vector xyzc;
    for (; ai != mpstructure->end(); ++ai)
    {
        xyzc = L.ucvCartesian(ai->xyz_cartn);
        mcartesian_positions_uc.push_back(xyzc);
    }
}

// Public Methods ------------------------------------------------------------

void PeriodicStructureBondGenerator::rewind()
{
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
}


void PeriodicStructureBondGenerator::selectAnchorSite(int anchor)
{
    this->BaseBondGenerator::selectAnchorSite(anchor);
    mr0 = mcartesian_positions_uc[anchor];
}


void PeriodicStructureBondGenerator::setRmin(double rmin)
{
    // destroy msphere so it will be created on rewind with new rmin
    if (this->getRmin() != rmin)    msphere.reset();
    this->BaseBondGenerator::setRmin(rmin);
}


void PeriodicStructureBondGenerator::setRmax(double rmax)
{
    // destroy msphere so it will be created on rewind with new rmax
    if (this->getRmax() != rmax)    msphere.reset();
    this->BaseBondGenerator::setRmax(rmax);
}

// Protected Methods ---------------------------------------------------------

bool PeriodicStructureBondGenerator::iterateSymmetry()
{
    msphere->next();
    bool done = msphere->finished();
    mrcsphere = done ? R3::zerovector :
        mpstructure->getLattice().cartesian(msphere->mno());
    return !done;
}


void PeriodicStructureBondGenerator::rewindSymmetry()
{
    msphere->rewind();
    mrcsphere = msphere->finished() ? R3::zerovector :
        mpstructure->getLattice().cartesian(msphere->mno());
    this->updater1();
}


void PeriodicStructureBondGenerator::getNextBond()
{
    ++msite_current;
    // go back to the first site if there is next symmetry element
    if (msite_current >= msite_last && this->iterateSymmetry())
    {
        msite_current = msite_first;
    }
    // update values only if not finished
    if (!this->finished())  this->updater1();
}

// Private Methods -----------------------------------------------------------

void PeriodicStructureBondGenerator::updater1()
{
    mr1 = mrcsphere + mcartesian_positions_uc[this->site1()];
    this->updateDistance();
}

}   // namespace srreal
}   // namespace diffpy

// Serialization -------------------------------------------------------------

DIFFPY_INSTANTIATE_SERIALIZATION(diffpy::srreal::PeriodicStructureAdapter)
BOOST_CLASS_EXPORT_IMPLEMENT(diffpy::srreal::PeriodicStructureAdapter)

// End of file
