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
* $Id$
*
*****************************************************************************/

#include <diffpy/srreal/BaseBondGenerator.hpp>
#include <diffpy/srreal/StructureAdapter.hpp>
#include <diffpy/mathutils.hpp>

using diffpy::mathutils::eps_eq;

namespace diffpy {
namespace srreal {

//using namespace std;

// Constructor ---------------------------------------------------------------

BaseBondGenerator::BaseBondGenerator(const StructureAdapter* stru)
{
    mstructure = stru;
    this->uncache();
    this->setRmin(0.0);
    this->setRmax(DEFAULT_BONDGENERATOR_RMAX);
    this->selectAnchorSite(0);
    this->selectSiteRange(0, mstructure->countSites());
}

// Public Methods ------------------------------------------------------------

// loop control

void BaseBondGenerator::rewind()
{
    this->uncache();
    msite_current = msite_first;
    // avoid calling rewindSymmetry at an invalid site
    if (this->finished())   return;
    this->rewindSymmetry();
    this->advanceWhileInvalid();
}


void BaseBondGenerator::next()
{
    this->uncache();
    this->getNextBond();
    this->advanceWhileInvalid();
}


void BaseBondGenerator::nextsite()
{
    this->uncache();
    msite_current += 1;
    // avoid calling rewindSymmetry at an invalid site
    if (this->finished())   return;
    this->rewindSymmetry();
}


bool BaseBondGenerator::finished() const
{
    return msite_current >= msite_last;
}

// configuration

void BaseBondGenerator::selectAnchorSite(int anchor)
{
    msite_anchor = anchor;
    this->setFinishedFlag();
}


void BaseBondGenerator::selectSiteRange(int first, int last)
{
    msite_first = first;
    msite_last = last;
    this->setFinishedFlag();
}


void BaseBondGenerator::setRmin(double rmin)
{
    this->uncache();
    if (rmin != mrmin)  this->setFinishedFlag();
    mrmin = rmin;
}


void BaseBondGenerator::setRmax(double rmax)
{
    this->uncache();
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
    const R3::Vector& rv = mstructure->siteCartesianPosition(this->site0());
    return rv;
}


const R3::Vector& BaseBondGenerator::r1() const
{
    const R3::Vector& rv = mstructure->siteCartesianPosition(this->site1());
    return rv;
}


double BaseBondGenerator::distance() const
{
    if (!mdistance_cached)
    {
        mdistance = R3::distance(this->r0(), this->r1());
        mdistance_cached = true;
    }
    return mdistance;
}


const R3::Vector& BaseBondGenerator::r01() const
{
    static R3::Vector rv;
    rv = this->r1() - this->r0();
    return rv;
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
    this->uncache();
    return false;
}


void BaseBondGenerator::rewindSymmetry()
{
    this->uncache();
}


void BaseBondGenerator::uncache()
{
    mdistance_cached = false;
}

// Private Methods -----------------------------------------------------------

void BaseBondGenerator::getNextBond()
{
    this->uncache();
    if (this->iterateSymmetry())  return;
    this->nextsite();
}


void BaseBondGenerator::advanceWhileInvalid()
{
    this->uncache();
    while (!this->finished() &&
            (this->bondOutOfRange() || this->atSelfPair()))
    {
        this->getNextBond();
    }
}


bool BaseBondGenerator::bondOutOfRange() const
{
    double d = this->distance();
    bool rv = (d < this->getRmin()) || (d > this->getRmax());
    return rv;
}


bool BaseBondGenerator::atSelfPair() const
{
    bool rv = (this->site0() == this->site1()) &&
        eps_eq(this->distance(), 0.0);
    return rv;
}


void BaseBondGenerator::setFinishedFlag()
{
    this->uncache();
    msite_current = msite_last;
}

}   // namespace srreal
}   // namespace diffpy

// End of file
