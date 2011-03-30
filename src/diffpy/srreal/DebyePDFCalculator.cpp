/*****************************************************************************
*
* diffpy.srreal     by DANSE Diffraction group
*                   Simon J. L. Billinge
*                   (c) 2010 Trustees of the Columbia University
*                   in the City of New York.  All rights reserved.
*
* File coded by:    Pavol Juhas, Chris Farrow
*
* See AUTHORS.txt for a list of people who contributed.
* See LICENSE.txt for license information.
*
******************************************************************************
*
* class DebyePDFCalculator -- calculate PDF from the Debye equation.
*
* $Id$
*
*****************************************************************************/

#include <cassert>
#include <valarray>
#include <stdexcept>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_fft_complex.h>

#include <diffpy/srreal/DebyePDFCalculator.hpp>
#include <diffpy/srreal/PDFUtils.hpp>
#include <diffpy/srreal/GaussianProfile.hpp>
#include <diffpy/mathutils.hpp>
#include <diffpy/validators.hpp>
#include <diffpy/serialization.hpp>

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
    // default configuration
    this->setScatteringFactorTableByType("periodictablexray");
    this->setRstep(DEFAULT_PDFCALCULATOR_RSTEP);
    this->setRmax(DEFAULT_PDFCALCULATOR_RMAX);
    this->setMaxExtension(DEFAULT_PDFCALCULATOR_MAXEXTENSION);
    this->setOptimumQstep();
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

QuantityType DebyePDFCalculator::getPDF() const
{
    // build a zero padded F vector that gives dr <= rstep
    QuantityType fpad = this->getF();
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
    QuantityType pdf1 = this->applyEnvelopes(rgrid, pdf0);
    return pdf1;
}


QuantityType DebyePDFCalculator::getRDF() const
{
    // FIXME
    return QuantityType();
}


QuantityType DebyePDFCalculator::getRDFperR() const
{
    // FIXME
    return QuantityType();
}

// Q-range configuration

void DebyePDFCalculator::setQstep(double qstep)
{
    moptimumqstep = false;
    this->BaseDebyeSum::setQstep(qstep);
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
    this->PairQuantity::setRmin(rmin);
}


void DebyePDFCalculator::setRmax(double rmax)
{
    ensureNonNegative("Rmax", rmax);
    this->PairQuantity::setRmax(rmax);
    this->updateQstep();
}


void DebyePDFCalculator::setRstep(double rstep)
{
    ensureEpsilonPositive("Rstep", rstep);
    mrstep = rstep;
}


const double& DebyePDFCalculator::getRstep() const
{
    return mrstep;
}


void DebyePDFCalculator::setMaxExtension(double maxextension)
{
    ensureNonNegative("maxextension", maxextension);
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
    const ScatteringFactorTablePtr sftable = this->getScatteringFactorTable();
    const string& smbl = mstructure->siteAtomType(siteidx);
    const double occupancy = mstructure->siteOccupancy(siteidx);
    double rv = sftable->lookupatq(smbl, Q) * occupancy;
    return rv;
}

// Private Methods -----------------------------------------------------------

void DebyePDFCalculator::updateQstep()
{
    if (!moptimumqstep)  return;
    double rmaxext = this->rcalchi();
    // Use at least 4 steps to Qmax even for tiny rmaxext.
    // Avoid division by zero.
    double qstep = (this->getQmax() * rmaxext / M_PI > 4) ?
        (M_PI / rmaxext) : (this->getQmax() / 4);
    this->BaseDebyeSum::setQstep(qstep);
}


double DebyePDFCalculator::rcalclo() const
{
    double rv = mrcalclosteps * this->getRstep();
    return rv;
}


double DebyePDFCalculator::rcalchi() const
{
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
    double maxmsd = 2 * maxUii(mstructure);
    double maxfwhm = this->getPeakWidthModel()->calculateFromMSD(maxmsd);
    // assume Gaussian peak profile
    GaussianProfile pkf;
    pkf.setPrecision(DEFAULT_PEAKPRECISION);
    // Gaussian function is symmetric, no need to use xboundlo
    double rv = pkf.xboundhi(maxfwhm);
    return rv;
}


void DebyePDFCalculator::cacheRlimitsData()
{
    double ext_total = min(this->getMaxExtension(),
        this->extFromTerminationRipples() + this->extFromPeakTails());
    // abbreviations
    const double& rmin = this->getRmin();
    const double& rmax = this->getRmax();
    const double& dr = this->getRstep();
    mrcalclosteps = max(0, pdfutils_rminSteps(rmin - ext_total, dr));
    mrcalchisteps = pdfutils_rmaxSteps(rmax + ext_total, dr);
}

}   // namespace srreal
}   // namespace diffpy

// Serialization -------------------------------------------------------------

DIFFPY_INSTANTIATE_SERIALIZE(diffpy::srreal::DebyePDFCalculator)
BOOST_CLASS_EXPORT(diffpy::srreal::DebyePDFCalculator)

// End of file
