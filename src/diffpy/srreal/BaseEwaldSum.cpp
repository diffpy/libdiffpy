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
* class BaseEwaldSum -- base class for concrete Ewald sum calculators
*
* $Id$
*
*****************************************************************************/

#include <cassert>
#include <string>
#include <stdexcept>
#include <sstream>

#include <diffpy/srreal/BaseEwaldSum.hpp>
#include <diffpy/mathutils.hpp>
#include <diffpy/srreal/StructureAdapter.hpp>
#include <diffpy/srreal/BaseBondGenerator.hpp>

using namespace std;
using diffpy::mathutils::eps_eq;
using diffpy::mathutils::eps_gt;

namespace diffpy {
namespace srreal {

// Declaration of Local Helpers ----------------------------------------------

namespace {

const double DEFAULT_EWALD_PRECISION = 1e-5;
void ensureNonNegative(const string& vname, double value);

}   // namespace

// Constructor ---------------------------------------------------------------

BaseEwaldSum::BaseEwaldSum()
{
    // default configuration
    this->setQmin(0.0);
    this->setQmax(10.0);
    this->setQmin(0.05);
    this->setEwaldPrecision(DEFAULT_EWALD_PRECISION);
}

// Public Methods ------------------------------------------------------------

// results

QuantityType BaseEwaldSum::getF() const
{
    QuantityType rv = this->getExtendedF();
    rv.erase(rv.begin(), rv.begin() + this->qminPoints());
    return rv;
}


QuantityType BaseEwaldSum::getQgrid() const
{
    QuantityType rv = this->getExtendedQgrid();
    rv.erase(rv.begin(), rv.begin() + this->qminPoints());
    return rv;
}


QuantityType BaseEwaldSum::getExtendedF() const
{
    QuantityType rv = mvalue;
    const double& dq = this->getQstep();
    const double totocc = mstructure->totalOccupancy();
    const int npts = this->totalPoints();
    for (int kq = this->qminPoints(); kq < npts; ++kq)
    {
        double q = kq * dq;
        double sfavg = this->sfAverageAtQ(q);
        double fscale = (sfavg * totocc) == 0 ? 0.0 :
            1.0 / (sfavg * sfavg * totocc);
        rv[kq] *= fscale;
    }
    return rv;
}


QuantityType BaseEwaldSum::getExtendedQgrid() const
{
    const int npts = this->totalPoints();
    QuantityType rv;
    rv.reserve(npts);
    for (int kq = 0; kq < npts; ++kq)  rv.push_back(kq * this->getQstep());
    return rv;
}

// Q-range configuration

void BaseEwaldSum::setQmin(double qmin)
{
    ensureNonNegative("Qmin", qmin);
    mqmin = qmin;
    this->cacheQpointsData();
}


const double& BaseEwaldSum::getQmin() const
{
    return mqmin;
}


void BaseEwaldSum::setQmax(double qmax)
{
    ensureNonNegative("Qmax", qmax);
    mqmax = qmax;
    this->cacheQpointsData();
}


const double& BaseEwaldSum::getQmax() const
{
    return mqmax;
}


void BaseEwaldSum::setQstep(double qstep)
{
    if (!eps_gt(qstep, 0))
    {
        const char* emsg = "Qstep must be positive.";
        throw invalid_argument(emsg);
    }
    mqstep = qstep;
    this->cacheQpointsData();
}


const double& BaseEwaldSum::getQstep() const
{
    return mqstep;
}


void BaseEwaldSum::setEwaldPrecision(double precision)
{
    mewaldprecision = precision;
}


const double& BaseEwaldSum::getEwaldPrecision() const
{
    return mewaldprecision;
}

// Protected Methods ---------------------------------------------------------

// PairQuantity overloads

void BaseEwaldSum::resetValue()
{
    this->resizeValue(this->totalPoints());
    this->PairQuantity::resetValue();
}


void BaseEwaldSum::addPairContribution(const BaseBondGenerator& bnds)
{
    const int summationscale = 2;
    const double dist = bnds.distance();
    if (eps_eq(0.0, dist))  return;
    this->setupPairScale(bnds);
    for (int kq = this->qminPoints(); kq < this->totalPoints(); ++kq)
    {
        const double q = kq * this->getQstep();
        double pairscale = this->pairScale(q);
        if (pairscale / dist < this->getEwaldPrecision())   break;
        double scaledsfprod = summationscale * pairscale *
            this->sfSiteAtQ(bnds.site0(), q) *
            this->sfSiteAtQ(bnds.site1(), q);
        mvalue[kq] += scaledsfprod * sin(q * dist) / dist;
    }
}


void BaseEwaldSum::setupPairScale(const BaseBondGenerator& bnds)
{
    return;
}


double BaseEwaldSum::pairScale(const double& q) const
{
    return 1.0;
}


double BaseEwaldSum::sfSiteAtQ(int siteidx, const double& q) const
{
    return 1.0;
}


double BaseEwaldSum::sfAverageAtQ(const double&q) const
{
    return 1.0;
}

// Private Methods -----------------------------------------------------------

int BaseEwaldSum::qminPoints() const
{
    return mqpoints_cache.qminpoints;
}


int BaseEwaldSum::totalPoints() const
{
    return mqpoints_cache.totalpoints;
}


void BaseEwaldSum::cacheQpointsData()
{
    const double& dq = this->getQstep();
    mqpoints_cache.qminpoints = int(this->getQmin() / dq);
    mqpoints_cache.totalpoints = int(ceil(this->getQmax() / dq));
    // include point for qmax when it is a close multiple of dq
    if (eps_eq(this->getQmax(), mqpoints_cache.totalpoints * dq))
    {
        mqpoints_cache.totalpoints += 1;
    }
}

// Local Helpers -------------------------------------------------------------

namespace {

void ensureNonNegative(const string& vname, double value)
{
    if (value < 0.0)
    {
        stringstream emsg;
        emsg << vname << " cannot be negative.";
        throw invalid_argument(emsg.str());
    }
}

}   // namespace

}   // namespace srreal
}   // namespace diffpy

// End of file
