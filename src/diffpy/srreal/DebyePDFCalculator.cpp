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
#include <diffpy/srreal/BaseBondGenerator.hpp>
#include <diffpy/mathutils.hpp>
#include <diffpy/validators.hpp>

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
    mtotalextension = 0.0;
    // default configuration
    this->setScatteringFactorTableByType("SFTperiodictableXray");
    this->setRstep(DEFAULT_PDFCALCULATOR_RSTEP);
    this->setRmax(DEFAULT_PDFCALCULATOR_RMAX);
    this->setMaxExtension(DEFAULT_PDFCALCULATOR_MAXEXTENSION);
    this->setOptimumQstep();
    this->setQmax(DEFAULT_DEBYEPDFCALCULATOR_QMAX);
    // envelopes
    this->addEnvelopeByType("scale");
    this->addEnvelopeByType("qresolution");
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
    // FIXME - check this especially the coefficient 2
    int nfromdr = 2 * M_PI / this->getRstep() / this->getQstep();
    int nrequired = max(pdfutils_qmaxSteps(this), nfromdr);
    int nlog2 = int(floor(log2(nrequired))) + 1;
    int padlen = int(pow(2, nlog2));
    // complex valarray needs to have twice as many elements
    valarray<double> ftog(0.0, 2 * padlen);
    QuantityType F = this->getF();
    assert(F.size() * 2 <= ftog.size());
    QuantityType::const_iterator fe = F.begin();
    double* pfc = &(ftog[0]);
    for (; fe != F.end(); ++fe, pfc += 2)  { *pfc = *fe; }
    // apply inverse fft
    int status;
    status = gsl_fft_complex_radix2_inverse(&(ftog[0]), 1, padlen);
    if (status != GSL_SUCCESS)
    {
        const char* emsgft = "Fourier Transformation failed.";
        throw invalid_argument(emsgft);
    }
    // normalize the complex part
    double qmaxpad = padlen * this->getQstep();
    valarray<double> gpad(0.0, padlen / 2);
    for (int i = 0; i < padlen / 2; ++i)
    {
        gpad[i] = ftog[2 * i + 1] * 2 / M_PI * qmaxpad;
    }
    // interpolate to the output r-grid
    const double drpad = 2 * M_PI / qmaxpad;
    QuantityType rgrid = this->getRgrid();
    QuantityType pdf(rgrid.size());
    QuantityType::const_iterator ri = rgrid.begin();
    QuantityType::iterator pdfi = pdf.begin();
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
    return pdf;
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


void DebyePDFCalculator::configureBondGenerator(BaseBondGenerator& bnds)
{
    bnds.setRmin(this->rcalclo());
    bnds.setRmax(this->rcalchi());
}


double DebyePDFCalculator::sfSiteAtQ(int siteidx, const double& Q) const
{
    const ScatteringFactorTablePtr sftable = this->getScatteringFactorTable();
    const string& smbl = mstructure->siteAtomType(siteidx);
    double rv = sftable->lookup(smbl) * mstructure->siteOccupancy(siteidx);
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
    double rv = this->getRmin() - mtotalextension;
    rv = max(rv, 0.0);
    return rv;
}


double DebyePDFCalculator::rcalchi() const
{
    double rv = this->getRmax() + mtotalextension;
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
    pkf.setPrecision(this->getDebyePrecision());
    // Gaussian function is symmetric, no need to use xboundlo
    double rv = pkf.xboundhi(maxfwhm);
    return rv;
}


void DebyePDFCalculator::cacheRlimitsData()
{
    double totext =
        this->extFromTerminationRipples() + this->extFromPeakTails();
    mtotalextension = min(totext, this->getMaxExtension());
}

}   // namespace srreal
}   // namespace diffpy

// End of file
