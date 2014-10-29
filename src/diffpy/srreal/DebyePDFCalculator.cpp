/*****************************************************************************
*
* libdiffpy         by DANSE Diffraction group
*                   Simon J. L. Billinge
*                   (c) 2010 The Trustees of Columbia University
*                   in the City of New York.  All rights reserved.
*
* File coded by:    Pavol Juhas, Chris Farrow
*
* See AUTHORS.txt for a list of people who contributed.
* See LICENSE_DANSE.txt for license information.
*
******************************************************************************
*
* class DebyePDFCalculator -- calculate PDF from the Debye equation.
*
*****************************************************************************/

#include <cassert>
#include <valarray>
#include <stdexcept>

#include <diffpy/srreal/DebyePDFCalculator.hpp>
#include <diffpy/srreal/PDFUtils.hpp>
#include <diffpy/srreal/GaussianProfile.hpp>
#include <diffpy/mathutils.hpp>
#include <diffpy/validators.hpp>
#include <diffpy/serialization.ipp>

using namespace std;
using namespace diffpy::validators;
using diffpy::mathutils::eps_gt;

namespace diffpy {
namespace srreal {

const double DEFAULT_DEBYEPDFCALCULATOR_QMAX = 25.0;

// Constructor ---------------------------------------------------------------

DebyePDFCalculator::DebyePDFCalculator()
{
    // initializations
    mrcalclosteps = 0;
    mrcalchisteps = 0;
    mrlimits_are_cached = false;
    // default configuration
    this->setScatteringFactorTableByType("xray");
    this->setRstep(DEFAULT_PDFCALCULATOR_RSTEP);
    this->setRmax(DEFAULT_PDFCALCULATOR_RMAX);
    this->setMaxExtension(DEFAULT_PDFCALCULATOR_MAXEXTENSION);
    this->setOptimumQstep();
    this->setQmin(0.0);
    this->setQmax(DEFAULT_DEBYEPDFCALCULATOR_QMAX);
    // envelopes
    this->addEnvelopeByType("scale");
    this->addEnvelopeByType("qresolution");
    // cache all internal data according to an empty structure.
    // this rebuilds mstructure_cache and mrlimits_cache.
    this->setStructure(mstructure);
    // attributes
    this->registerDoubleAttribute("qmin", this,
            &DebyePDFCalculator::getQmin, &DebyePDFCalculator::setQmin);
    this->registerDoubleAttribute("qmax", this,
            &DebyePDFCalculator::getQmax, &DebyePDFCalculator::setQmax);
    this->registerDoubleAttribute("qstep", this,
            &DebyePDFCalculator::getQstep, &DebyePDFCalculator::setQstep);
    this->registerDoubleAttribute("rmin", this,
            &DebyePDFCalculator::getRmin, &DebyePDFCalculator::setRmin);
    this->registerDoubleAttribute("rmax", this,
            &DebyePDFCalculator::getRmax, &DebyePDFCalculator::setRmax);
    this->registerDoubleAttribute("rstep", this,
            &DebyePDFCalculator::getRstep, &DebyePDFCalculator::setRstep);
    this->registerDoubleAttribute("maxextension", this,
            &DebyePDFCalculator::getMaxExtension,
            &DebyePDFCalculator::setMaxExtension);
    this->registerDoubleAttribute("extendedrmin", this,
            &DebyePDFCalculator::rcalclo);
    this->registerDoubleAttribute("extendedrmax", this,
            &DebyePDFCalculator::rcalchi);
}

// Public Methods ------------------------------------------------------------

// PairQuantity overloads

eventticker::EventTicker& DebyePDFCalculator::ticker() const
{
    // collect composite tickers for the BaseDebyeSum parent class
    eventticker::EventTicker& tic = this->BaseDebyeSum::ticker();
    assert(&mticker == &tic);
    tic.updateFrom(this->ScatteringFactorTableOwner::ticker());
    return tic;
}

// results

QuantityType DebyePDFCalculator::getPDF() const
{
    QuantityType rgrid = this->getRgrid();
    QuantityType pdf0 = this->getPDFAtQmin(this->getQmin());
    QuantityType pdf1 = this->applyEnvelopes(rgrid, pdf0);
    return pdf1;
}


QuantityType DebyePDFCalculator::getRDF() const
{
    QuantityType rgrid = this->getRgrid();
    QuantityType rv = this->getRDFperR();
    assert(rv.size() == rgrid.size());
    transform(rgrid.begin(), rgrid.end(), rv.begin(), rv.begin(),
            multiplies<double>());
    return rv;
}


QuantityType DebyePDFCalculator::getRDFperR() const
{
    return this->getPDFAtQmin(0.0);
}

// Q-range configuration

void DebyePDFCalculator::setQmin(double qmin)
{
    // NOTE: do not update ticker here, because the actual
    // Debye summation is always conducted from Qmin=0.
    ensureNonNegative("Qmin", qmin);
    this->BaseDebyeSum::setQmin(0.0);
    mqminpdf = qmin;
}


const double& DebyePDFCalculator::getQmin() const
{
    return mqminpdf;
}


void DebyePDFCalculator::setQstep(double qstep)
{
    moptimumqstep = false;
    this->BaseDebyeSum::setQstep(qstep);
    this->updateQstep();
}


void DebyePDFCalculator::setOptimumQstep()
{
    moptimumqstep = true;
    this->updateQstep();
}


bool DebyePDFCalculator::isOptimumQstep() const
{
    return moptimumqstep;
}

// R-range methods

QuantityType DebyePDFCalculator::getRgrid() const
{
    return pdfutils_getRgrid(this);
}

// R-range configuration

void DebyePDFCalculator::setRmin(double rmin)
{
    ensureNonNegative("Rmin", rmin);
    if (mrmin != rmin)  mrlimits_are_cached = false;
    this->PairQuantity::setRmin(rmin);
}


void DebyePDFCalculator::setRmax(double rmax)
{
    ensureNonNegative("Rmax", rmax);
    if (mrmax != rmax)  mrlimits_are_cached = false;
    this->PairQuantity::setRmax(rmax);
    this->updateQstep();
}


void DebyePDFCalculator::setRstep(double rstep)
{
    // NOTE: do not update ticker here, rstep is only used in the FFT
    // and not in the Debye summation.
    ensureEpsilonPositive("Rstep", rstep);
    if (mrstep != rstep)  mrlimits_are_cached = false;
    mrstep = rstep;
}


const double& DebyePDFCalculator::getRstep() const
{
    return mrstep;
}


void DebyePDFCalculator::setMaxExtension(double maxextension)
{
    ensureNonNegative("maxextension", maxextension);
    if (mmaxextension != maxextension)
    {
        mticker.click();
        mrlimits_are_cached = false;
    }
    mmaxextension = maxextension;
}


const double& DebyePDFCalculator::getMaxExtension() const
{
    return mmaxextension;
}

// Protected Methods ---------------------------------------------------------

// attributes overload to direct visitors around data structures

namespace {

// single implementation for constant and non-constant accept
template <class T>
void debyepdfcalc_accept(T* obj, diffpy::BaseAttributesVisitor& v)
{
    obj->getPeakWidthModel()->accept(v);
    // PDF envelopes
    set<string> evnames = obj->usedEnvelopeTypes();
    set<string>::const_iterator nm = evnames.begin();
    for (; nm != evnames.end(); ++nm)
    {
        obj->getEnvelopeByType(*nm)->accept(v);
    }
    // finally call standard accept
    obj->diffpy::Attributes::accept(v);
}

}   // namespace


void DebyePDFCalculator::accept(diffpy::BaseAttributesVisitor& v)
{
    debyepdfcalc_accept(this, v);
}


void DebyePDFCalculator::accept(diffpy::BaseAttributesVisitor& v) const
{
    debyepdfcalc_accept(this, v);
}

// BaseDebyeSum overloads

void DebyePDFCalculator::resetValue()
{
    this->cacheRlimitsData();
    this->updateQstep();
    this->BaseDebyeSum::resetValue();
}


void DebyePDFCalculator::configureBondGenerator(BaseBondGenerator& bnds) const
{
    bnds.setRmin(this->rcalclo());
    bnds.setRmax(this->rcalchi());
}


double DebyePDFCalculator::sfSiteAtQ(int siteidx, const double& Q) const
{
    const ScatteringFactorTablePtr& sftable = this->getScatteringFactorTable();
    const string& smbl = mstructure->siteAtomType(siteidx);
    const double occupancy = mstructure->siteOccupancy(siteidx);
    double rv = sftable->lookup(smbl, Q) * occupancy;
    return rv;
}

// Private Methods -----------------------------------------------------------

QuantityType DebyePDFCalculator::getPDFAtQmin(double qmin) const
{
    // build a zero padded F vector that gives dr <= rstep
    QuantityType fpad = this->getF();
    // zero all F values below qmin
    int nqmin = pdfutils_qminSteps(qmin, this->getQstep());
    if (nqmin > int(fpad.size()))  nqmin = fpad.size();
    fill(fpad.begin(), fpad.begin() + nqmin, 0.0);
    int nfromdr = int(ceil(M_PI / this->getRstep() / this->getQstep()));
    if (nfromdr > int(fpad.size()))  fpad.resize(nfromdr, 0.0);
    QuantityType gpad = fftftog(fpad, this->getQstep());
    const double drpad = M_PI / (gpad.size() * this->getQstep());
    QuantityType rgrid = this->getRgrid();
    QuantityType pdf0(rgrid.size());
    QuantityType::const_iterator ri = rgrid.begin();
    QuantityType::iterator pdfi = pdf0.begin();
    for (; ri != rgrid.end(); ++ri, ++pdfi)
    {
        double xdrp = *ri / drpad;
        int iplo = int(xdrp);
        int iphi = iplo + 1;
        double wphi = xdrp - iplo;
        double wplo = 1.0 - wphi;
        assert(iphi < int(gpad.size()));
        *pdfi = wplo * gpad[iplo] + wphi * gpad[iphi];
    }
    return pdf0;
}


void DebyePDFCalculator::updateQstep()
{
    double rmaxext = this->rcalchi();
    // Use at least 4 steps to Qmax even for tiny rmaxext.
    // Avoid division by zero.
    double oqstep = (this->getQmax() * rmaxext / M_PI > 4) ?
        (M_PI / rmaxext) : (this->getQmax() / 4);
    // if custom qstep is higher than the optimum one,
    // force adjustment to the optimum value.
    if (this->getQstep() > oqstep)  moptimumqstep = true;
    if (moptimumqstep)  this->BaseDebyeSum::setQstep(oqstep);
}


double DebyePDFCalculator::rcalclo() const
{
    if (!mrlimits_are_cached)  this->cacheRlimitsData();
    double rv = mrcalclosteps * this->getRstep();
    return rv;
}


double DebyePDFCalculator::rcalchi() const
{
    if (!mrlimits_are_cached)  this->cacheRlimitsData();
    double rv = mrcalchisteps * this->getRstep();
    return rv;
}


double DebyePDFCalculator::extFromTerminationRipples() const
{
    // number of termination ripples for extending the r-range
    const int nripples = 6;
    // extension due to termination ripples
    double rv = (this->getQmax() > 0.0) ?
        (nripples*2*M_PI / this->getQmax()) : 0.0;
    return rv;
}


double DebyePDFCalculator::extFromPeakTails() const
{
    const PeakWidthModel& pwm = *(this->getPeakWidthModel());
    double maxfwhm = pwm.maxWidth(
            mstructure, this->getRmin(), this->getRmax());
    // assume Gaussian peak profile
    GaussianProfile pkf;
    pkf.setPrecision(DEFAULT_PEAKPRECISION);
    // Gaussian function is symmetric, no need to use xboundlo
    double rv = pkf.xboundhi(maxfwhm);
    return rv;
}


void DebyePDFCalculator::cacheRlimitsData() const
{
    const int minreductionsteps = 50;
    double ext_total = min(this->getMaxExtension(),
        this->extFromTerminationRipples() + this->extFromPeakTails());
    // abbreviations
    const double& rmin = this->getRmin();
    const double& rmax = this->getRmax();
    const double& dr = this->getRstep();
    mrcalclosteps = max(0, pdfutils_rminSteps(rmin - ext_total, dr));
    const int nhi = pdfutils_rmaxSteps(rmax + ext_total, dr);
    if (nhi <= (mrcalchisteps - minreductionsteps) || nhi > mrcalchisteps)
    {
        mrcalchisteps = nhi;
    }
    mrlimits_are_cached = true;
}

}   // namespace srreal
}   // namespace diffpy

// Serialization -------------------------------------------------------------

DIFFPY_INSTANTIATE_SERIALIZATION(diffpy::srreal::DebyePDFCalculator)
BOOST_CLASS_EXPORT_IMPLEMENT(diffpy::srreal::DebyePDFCalculator)

// End of file
