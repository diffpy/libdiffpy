/*****************************************************************************
*
* diffpy.srreal     by DANSE Diffraction group
*                   Simon J. L. Billinge
*                   (c) 2009 Trustees of the Columbia University
*                   in the City of New York.  All rights reserved.
*
* File coded by:    Pavol Juhas
*
* See AUTHORS.txt for a list of people who contributed.
* See LICENSE.txt for license information.
*
******************************************************************************
*
* class BaseBondGenerator -- semi-abstract class for generation
*     of all atom pairs containing specified anchor atom.
*
*****************************************************************************/

#include <algorithm>
#include <diffpy/srreal/BaseBondGenerator.hpp>
#include <diffpy/srreal/StructureAdapter.hpp>
#include <diffpy/mathutils.hpp>

using diffpy::mathutils::eps_eq;

namespace diffpy {
namespace srreal {

//using namespace std;

// Constructor ---------------------------------------------------------------

BaseBondGenerator::BaseBondGenerator(StructureAdapterConstPtr stru)
{
    msite_anchor = 0;
    msite_first = 0;
    msite_last = 0;
    msite_current = 0;
    mstructure = stru;
    this->setRmin(0.0);
    this->setRmax(DEFAULT_BONDGENERATOR_RMAX);
    msite_selected.resize(mstructure->countSites(), true);
    if (mstructure->countSites())
    {
        this->selectAnchorSite(0);
        this->selectSiteRange(0, mstructure->countSites());
    }
}

// Public Methods ------------------------------------------------------------

// loop control

void BaseBondGenerator::rewind()
{
    msite_current = msite_first;
    this->advanceIfSkippedSite();
    // avoid calling rewindSymmetry at an invalid site
    if (this->finished())   return;
    this->rewindSymmetry();
    this->advanceWhileInvalid();
}


void BaseBondGenerator::next()
{
    this->getNextBond();
    this->advanceWhileInvalid();
}


bool BaseBondGenerator::finished() const
{
    return msite_current >= msite_last;
}

// configuration

void BaseBondGenerator::selectAnchorSite(int anchor)
{
    assert(0 <= anchor && anchor < mstructure->countSites());
    msite_anchor = anchor;
    mr0 = mstructure->siteCartesianPosition(msite_anchor);
    this->setFinishedFlag();
}


void BaseBondGenerator::selectSiteRange(int first, int last)
{
    assert(0 <= first && first < mstructure->countSites());
    assert(last <= mstructure->countSites());
    msite_first = first;
    msite_last = last;
    std::fill(msite_selected.begin() + first,
            msite_selected.begin() + last, true);
    this->setFinishedFlag();
}


void BaseBondGenerator::skipSite(int skipindex)
{
    msite_selected[skipindex] = false;
}

void BaseBondGenerator::setRmin(double rmin)
{
    if (rmin != mrmin)  this->setFinishedFlag();
    mrmin = rmin;
}


void BaseBondGenerator::setRmax(double rmax)
{
    if (rmax != mrmax)  this->setFinishedFlag();
    mrmax = rmax;
}

// data query

const double& BaseBondGenerator::getRmin() const
{
    return mrmin;
}


const double& BaseBondGenerator::getRmax() const
{
    return mrmax;
}


int BaseBondGenerator::site0() const
{
    return msite_anchor;
}


int BaseBondGenerator::site1() const
{
    return msite_current;
}


int BaseBondGenerator::multiplicity() const
{
    return mstructure->siteMultiplicity(this->site0());
}


const R3::Vector& BaseBondGenerator::r0() const
{
    return mr0;
}


const R3::Vector& BaseBondGenerator::r1() const
{
    return mr1;
}


const double& BaseBondGenerator::distance() const
{
    return mdistance;
}


const R3::Vector& BaseBondGenerator::r01() const
{
    return mr01;
}


const R3::Matrix& BaseBondGenerator::Ucartesian0() const
{
    return mstructure->siteCartesianUij(this->site0());
}


const R3::Matrix& BaseBondGenerator::Ucartesian1() const
{
    return mstructure->siteCartesianUij(this->site1());
}


double BaseBondGenerator::msd() const
{
    static R3::Vector s;
    s = this->r01();
    double msd0 = meanSquareDisplacement(this->Ucartesian0(), s,
            mstructure->siteAnisotropy(this->site0()));
    double msd1 = meanSquareDisplacement(this->Ucartesian1(), s,
            mstructure->siteAnisotropy(this->site1()));
    double rv = msd0 + msd1;
    return rv;
}

// Protected Methods ---------------------------------------------------------

bool BaseBondGenerator::iterateSymmetry()
{
    return false;
}


void BaseBondGenerator::rewindSymmetry()
{
    mr1 = mstructure->siteCartesianPosition(this->site1());
    this->updateDistance();
}


void BaseBondGenerator::getNextBond()
{
    if (this->iterateSymmetry())  return;
    msite_current += 1;
    this->advanceIfSkippedSite();
    // avoid calling rewindSymmetry at an invalid site
    if (!(this->finished()))  this->rewindSymmetry();
}


void BaseBondGenerator::updateDistance()
{
    if ((mdistance = fabs(mr01[0] = mr1[0] - mr0[0])) > mrmax)  return;
    if ((mdistance = fabs(mr01[1] = mr1[1] - mr0[1])) > mrmax)  return;
    if ((mdistance = fabs(mr01[2] = mr1[2] - mr0[2])) > mrmax)  return;
    mdistance = R3::norm(mr01);
}

// Private Methods -----------------------------------------------------------

void BaseBondGenerator::advanceIfSkippedSite()
{
    while (!this->finished() && !msite_selected[msite_current])
    {
        ++msite_current;
    }
}


void BaseBondGenerator::advanceWhileInvalid()
{
    while (!this->finished() &&
            (this->bondOutOfRange() || this->atSelfPair()))
    {
        this->getNextBond();
    }
}

// Private Methods -----------------------------------------------------------

bool BaseBondGenerator::bondOutOfRange() const
{
    const double& d = this->distance();
    bool rv = (d < this->getRmin()) || (d > this->getRmax());
    return rv;
}


bool BaseBondGenerator::atSelfPair() const
{
    bool rv = eps_eq(this->distance(), 0.0);
    return rv;
}


void BaseBondGenerator::setFinishedFlag()
{
    msite_current = msite_last;
}

}   // namespace srreal
}   // namespace diffpy

// End of file
