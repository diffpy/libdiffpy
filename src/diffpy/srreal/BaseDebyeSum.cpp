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
* class BaseDebyeSum -- base class for concrete Debye sum calculators
*
*****************************************************************************/

#include <cassert>
#include <string>
#include <stdexcept>
#include <sstream>
#include <functional>

#include <diffpy/srreal/BaseDebyeSum.hpp>
#include <diffpy/mathutils.hpp>
#include <diffpy/validators.hpp>
#include <diffpy/serialization.ipp>

using namespace std;
using namespace diffpy::validators;
using diffpy::mathutils::eps_gt;
using diffpy::mathutils::eps_eq;

namespace diffpy {
namespace srreal {

// Declaration of Local Helpers ----------------------------------------------

namespace {

/// Default cutoff for the Q-decreasing scale of the sine contributions.
const double DEFAULT_DEBYE_PRECISION = 1e-6;

}   // namespace

// Constructor ---------------------------------------------------------------

BaseDebyeSum::BaseDebyeSum()
{
    // default configuration
    this->setPeakWidthModelByType("jeong");
    this->setQmin(0.0);
    this->setQmax(DEFAULT_QGRID_QMAX);
    this->setQstep(DEFAULT_QGRID_QSTEP);
    this->setDebyePrecision(DEFAULT_DEBYE_PRECISION);
    // attributes
    this->registerDoubleAttribute("debyeprecision", this,
            &BaseDebyeSum::getDebyePrecision,
            &BaseDebyeSum::setDebyePrecision);
}

// Public Methods ------------------------------------------------------------

// PairQuantity overloads

eventticker::EventTicker& BaseDebyeSum::ticker() const
{
    const PeakWidthModelPtr& pwm = this->getPeakWidthModel();
    if (pwm)  mticker.updateFrom(pwm->ticker());
    return mticker;
}

// results

QuantityType BaseDebyeSum::getF() const
{
    QuantityType rv = this->value();
    const double& totocc = mstructure_cache.totaloccupancy;
    const int npts = pdfutils_qmaxSteps(this);
    for (int kq = pdfutils_qminSteps(this); kq < npts; ++kq)
    {
        double sfavg = this->sfAverageAtkQ(kq);
        double fscale = (sfavg * totocc) == 0 ? 0.0 :
            1.0 / (sfavg * sfavg * totocc);
        rv[kq] *= fscale;
    }
    return rv;
}

// Q-range methods

QuantityType BaseDebyeSum::getQgrid() const
{
    return pdfutils_getQgrid(this);
}

// Q-range configuration

void BaseDebyeSum::setQmin(double qmin)
{
    ensureNonNegative("Qmin", qmin);
    if (mqmin != qmin)  mticker.click();
    mqmin = qmin;
}


const double& BaseDebyeSum::getQmin() const
{
    return mqmin;
}


void BaseDebyeSum::setQmax(double qmax)
{
    ensureNonNegative("Qmax", qmax);
    if (mqmax != qmax)  mticker.click();
    mqmax = qmax;
}


const double& BaseDebyeSum::getQmax() const
{
    return mqmax;
}


void BaseDebyeSum::setQstep(double qstep)
{
    ensureEpsilonPositive("Qstep", qstep);
    if (mqstep != qstep)  mticker.click();
    mqstep = qstep;
}


const double& BaseDebyeSum::getQstep() const
{
    return mqstep;
}


void BaseDebyeSum::setDebyePrecision(double precision)
{
    if (mdebyeprecision != precision)  mticker.click();
    mdebyeprecision = precision;
}


const double& BaseDebyeSum::getDebyePrecision() const
{
    return mdebyeprecision;
}

// Protected Methods ---------------------------------------------------------

// PairQuantity overloads

void BaseDebyeSum::resetValue()
{
    this->cacheStructureData();
    this->resizeValue(pdfutils_qmaxSteps(this));
    this->PairQuantity::resetValue();
}


void BaseDebyeSum::addPairContribution(const BaseBondGenerator& bnds,
        int summationscale)
{
    const double dist = bnds.distance();
    if (eps_eq(0.0, dist))  return;
    // calculate sigma parameter for the Debye-Waller dampign Gaussian
    const double fwhm = this->getPeakWidthModel()->calculate(bnds);
    const double fwhmtosigma = 1.0 / (2 * sqrt(2 * M_LN2));
    const double dwsigma = fwhmtosigma * fwhm;
    const int nqpts = pdfutils_qmaxSteps(this);
    const int smscale = summationscale * bnds.multiplicity();
    const double& sineprec = this->getDebyePrecision();
    for (int kq = pdfutils_qminSteps(this); kq < nqpts; ++kq)
    {
        const double q = kq * this->getQstep();
        const double dwscale = exp(-0.5 * pow(dwsigma * q, 2));
        const double sinescale = smscale * dwscale *
            this->sfSiteAtkQ(bnds.site0(), kq) *
            this->sfSiteAtkQ(bnds.site1(), kq) / dist;
        if (eps_eq(0.0, sinescale, sineprec))   break;
        mvalue[kq] += sinescale * sin(q * dist);
    }
}


void BaseDebyeSum::stashPartialValue()
{
    mdbsumstash = this->value();
}


void BaseDebyeSum::restorePartialValue()
{
    assert(mdbsumstash.size() == mvalue.size());
    mvalue = mdbsumstash;
    mdbsumstash.clear();
}


double BaseDebyeSum::sfSiteAtQ(int siteidx, const double& q) const
{
    return 1.0;
}

// Private Methods -----------------------------------------------------------

double BaseDebyeSum::sfSiteAtkQ(int siteidx, int kq) const
{
    assert(0 <= siteidx && siteidx < int(mstructure_cache.typeofsite.size()));
    int typeidx = mstructure_cache.typeofsite[siteidx];
    assert(typeidx < int(mstructure_cache.sftypeatkq.size()));
    const QuantityType& sfarray = mstructure_cache.sftypeatkq[typeidx];
    assert(0 <= kq && kq < int(sfarray.size()));
    return sfarray[kq];
}


double BaseDebyeSum::sfAverageAtkQ(int kq) const
{
    assert(kq < int(mstructure_cache.sfaverageatkq.size()));
    return mstructure_cache.sfaverageatkq[kq];
}


void BaseDebyeSum::cacheStructureData()
{
    int cntsites = this->countSites();
    const int nqpts = pdfutils_qmaxSteps(this);
    QuantityType zeros(nqpts, 0.0);
    boost::unordered_map<string,int> atomtypeidx;
    // sftypeatkq
    mstructure_cache.typeofsite.clear();
    mstructure_cache.typeofsite.reserve(cntsites);
    mstructure_cache.sftypeatkq.clear();
    for (int siteidx = 0; siteidx < cntsites; ++siteidx)
    {
        const string& smbl = mstructure->siteAtomType(siteidx);
        if (!atomtypeidx.count(smbl))
        {
            atomtypeidx.insert(make_pair(smbl, int(atomtypeidx.size())));
        }
        int tpidx = atomtypeidx[smbl];
        mstructure_cache.typeofsite.push_back(tpidx);
        // do nothing if the type has been already cached
        if (tpidx < int(mstructure_cache.sftypeatkq.size()))  continue;
        assert(tpidx == int(mstructure_cache.sftypeatkq.size()));
        // here we need to build a new array
        mstructure_cache.sftypeatkq.push_back(zeros);
        QuantityType& sfarray = mstructure_cache.sftypeatkq.back();
        for (int kq = pdfutils_qminSteps(this); kq < nqpts; ++kq)
        {
            double q = this->getQstep() * kq;
            sfarray[kq] = this->sfSiteAtQ(siteidx, q);
        }
    }
    assert(cntsites == int(mstructure_cache.typeofsite.size()));
    assert(atomtypeidx.size() == mstructure_cache.sftypeatkq.size());
    // totaloccupancy
    mstructure_cache.totaloccupancy = mstructure->totalOccupancy();
    // sfaverageatkq
    QuantityType& sfak = mstructure_cache.sfaverageatkq;
    sfak = zeros;
    int ntps = mstructure_cache.sftypeatkq.size();
    vector<int> tpmultipl(ntps, 0);
    for (int siteidx = 0; siteidx < cntsites; ++siteidx)
    {
        int tpidx = mstructure_cache.typeofsite[siteidx];
        assert(tpidx < ntps);
        tpmultipl[tpidx] += mstructure->siteMultiplicity(siteidx);
    }
    for (int tpidx = 0; tpidx < ntps; ++tpidx)
    {
        QuantityType& sfarray = mstructure_cache.sftypeatkq[tpidx];
        const int multipl = tpmultipl[tpidx];
        for (int kq = pdfutils_qminSteps(this); kq < nqpts; ++kq)
        {
            sfak[kq] += sfarray[kq] * multipl;
        }
    }
    const double& totocc = mstructure_cache.totaloccupancy;
    double tosc = eps_gt(totocc, 0.0) ? (1.0 / totocc) : 1.0;
    transform(sfak.begin(), sfak.end(), sfak.begin(),
            bind1st(multiplies<double>(), tosc));
}

}   // namespace srreal
}   // namespace diffpy

// Serialization -------------------------------------------------------------

DIFFPY_INSTANTIATE_SERIALIZATION(diffpy::srreal::BaseDebyeSum)
BOOST_CLASS_EXPORT_IMPLEMENT(diffpy::srreal::BaseDebyeSum)

// End of file
