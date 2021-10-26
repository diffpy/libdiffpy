/*****************************************************************************
*
* libdiffpy         by DANSE Diffraction group
*                   Simon J. L. Billinge
*                   (c) 2009 The Trustees of Columbia University
*                   in the City of New York.  All rights reserved.
*
* File coded by:    Pavol Juhas
*
* See AUTHORS.txt for a list of people who contributed.
* See LICENSE_DANSE.txt for license information.
*
******************************************************************************
*
* class PDFCalculator -- real space PDF calculator
*
*****************************************************************************/

#include <iostream>
#include <stdexcept>
#include <sstream>
#include <cmath>
#include <cassert>

#include <diffpy/serialization.ipp>
#include <diffpy/srreal/PDFCalculator.hpp>
#include <diffpy/srreal/StructureAdapter.hpp>
#include <diffpy/srreal/R3linalg.hpp>
#include <diffpy/srreal/PDFUtils.hpp>
#include <diffpy/mathutils.hpp>
#include <diffpy/validators.hpp>

namespace diffpy {
namespace srreal {

using namespace std;
using namespace diffpy::validators;
using namespace diffpy::mathutils;

// Local Helpers -------------------------------------------------------------

namespace {

/// Return true if qstep value can be cheaply recomputed in initial setup.
template <class PDFC>
bool _initialQstepUpdate(const PDFC* pc)
{
    // Check if we use the initial empty structure and constructor
    // has completed, i.e., the last attribute is available.
    const bool rv =
        (pc->getStructure() == emptyStructureAdapter()) &&
        pc->hasDoubleAttr("extendedrmax");
    return rv;
}

}   // namespace

// Constructor ---------------------------------------------------------------

PDFCalculator::PDFCalculator() :
    mqmin(0.0),
    mqmax(DOUBLE_MAX),
    mrstep(DEFAULT_PDFCALCULATOR_RSTEP),
    mmaxextension(DEFAULT_PDFCALCULATOR_MAXEXTENSION)
{
    // default configuration
    mrmax = DEFAULT_PDFCALCULATOR_RMAX;
    this->setPeakWidthModelByType("jeong");
    this->setPeakProfileByType("gaussian");
    this->getPeakProfile()->setPrecision(DEFAULT_PEAKPRECISION);
    this->setBaselineByType("linear");
    this->setScatteringFactorTableByType("xray");
    // envelopes
    this->addEnvelopeByType("scale");
    this->addEnvelopeByType("qresolution");
    // cache all internal data according to an empty structure.
    // this rebuilds mstructure_cache and mrlimits_cache.
    this->setStructure(mstructure);
    // setup the OPTIMIZED evaluator last so that its validation passes.
    this->setEvaluatorType(OPTIMIZED);
    // attributes
    this->registerDoubleAttribute("qmin", this,
            &PDFCalculator::getQmin, &PDFCalculator::setQmin);
    this->registerDoubleAttribute("qmax", this,
            &PDFCalculator::getQmax, &PDFCalculator::setQmax);
    this->registerDoubleAttribute("qstep", this,
            &PDFCalculator::getQstep);
    this->registerDoubleAttribute("rmin", this,
            &PDFCalculator::getRmin, &PDFCalculator::setRmin);
    this->registerDoubleAttribute("rmax", this,
            &PDFCalculator::getRmax, &PDFCalculator::setRmax);
    this->registerDoubleAttribute("rstep", this,
            &PDFCalculator::getRstep, &PDFCalculator::setRstep);
    this->registerDoubleAttribute("maxextension", this,
            &PDFCalculator::getMaxExtension, &PDFCalculator::setMaxExtension);
    this->registerDoubleAttribute("extendedrmin", this,
            &PDFCalculator::getExtendedRmin);
    this->registerDoubleAttribute("extendedrmax", this,
            &PDFCalculator::getExtendedRmax);
}

// Public Methods ------------------------------------------------------------

// PairQuantity overloads

eventticker::EventTicker& PDFCalculator::ticker() const
{
    eventticker::EventTicker& tic = this->PairQuantity::ticker();
    assert(&mticker == &tic);
    tic.updateFrom(this->PeakWidthModelOwner::ticker());
    tic.updateFrom(this->ScatteringFactorTableOwner::ticker());
    const PeakProfilePtr& pkpf = this->getPeakProfile();
    if (pkpf)  tic.updateFrom(pkpf->ticker());
    return tic;
}

// results

QuantityType PDFCalculator::getPDF() const
{
    QuantityType pdf = this->getExtendedPDF();
    this->cutRipplePoints(pdf);
    return pdf;
}


QuantityType PDFCalculator::getRDF() const
{
    QuantityType rdf = this->getExtendedRDF();
    this->cutRipplePoints(rdf);
    return rdf;
}


QuantityType PDFCalculator::getRDFperR() const
{
    QuantityType rdfperr = this->getExtendedRDFperR();
    this->cutRipplePoints(rdfperr);
    return rdfperr;
}


QuantityType PDFCalculator::getF() const
{
    QuantityType f_ext = this->getExtendedF();
    assert(pdfutils_qmaxSteps(this) <= int(f_ext.size()));
    QuantityType rv(f_ext.begin(), f_ext.begin() + pdfutils_qmaxSteps(this));
    return rv;
}


QuantityType PDFCalculator::getExtendedPDF() const
{
    QuantityType rgrid_ext = this->getExtendedRgrid();
    // Skip FFT when qmax is not specified and qmin does not exclude the
    // the F(Q=Qstep) point (excluding F(0) == 0 makes no difference to G).
    const bool skipfft =
        !eps_lt(this->getQmax(), M_PI / this->getRstep()) &&
        !(1 < pdfutils_qminSteps(this));
    if (skipfft)
    {
        cout << "WARNING: skip FFT because qmax > pi / rstep.";
        QuantityType rdfpr = this->getExtendedRDFperR();
        QuantityType rdfprb = this->applyBaseline(rgrid_ext, rdfpr);
        QuantityType pdf = this->applyEnvelopes(rgrid_ext, rdfprb);
        return pdf;
    }
    // FFT required here
    // we need a full range PDF to apply termination ripples correctly
    QuantityType f_ext = this->getExtendedF();
    // zero all F points at Q < Qmin
    QuantityType::iterator ii_qmin =
        f_ext.begin() + min(pdfutils_qminSteps(this), int(f_ext.size()));
    fill(f_ext.begin(), ii_qmin, 0.0);
    // zero all F points at Q >= Qmax
    assert(pdfutils_qmaxSteps(this) <= int(f_ext.size()));
    QuantityType::iterator ii_qmax = f_ext.begin() + pdfutils_qmaxSteps(this);
    fill(ii_qmax, f_ext.end(), 0.0);
    QuantityType pdf1 = fftftog(f_ext, this->getQstep());
    // cut away the FFT padded points
    assert(this->extendedRmaxSteps() <= int(pdf1.size()));
    pdf1.erase(pdf1.begin() + this->extendedRmaxSteps(), pdf1.end());
    pdf1.erase(pdf1.begin(), pdf1.begin() + this->extendedRminSteps());
    QuantityType pdf2 = this->applyEnvelopes(rgrid_ext, pdf1);
    return pdf2;
}


QuantityType PDFCalculator::getExtendedRDF() const
{
    QuantityType rdf(this->countExtendedPoints());
    const double& totocc = mstructure_cache.totaloccupancy;
    double sfavg = this->sfAverage();
    double rdf_scale = (totocc * sfavg == 0.0) ? 0.0 :
        1.0 / (totocc * sfavg * sfavg);
    QuantityType::iterator iirdf = rdf.begin();
    QuantityType::const_iterator iival, iival_last;
    iival = this->value().begin() +
        this->extendedRminSteps() - this->rcalcloSteps();
    iival_last = this->value().begin() +
        this->extendedRmaxSteps() - this->rcalcloSteps();
    assert(iival >= this->value().begin());
    assert(iival_last <= this->value().end());
    assert(rdf.size() == size_t(iival_last - iival));
    for (; iirdf != rdf.end(); ++iival, ++iirdf)
    {
        *iirdf = *iival * rdf_scale;
    }
    return rdf;
}


QuantityType PDFCalculator::getExtendedRDFperR() const
{
    QuantityType rdf_ext = this->getExtendedRDF();
    QuantityType rgrid_ext = this->getExtendedRgrid();
    assert(rdf_ext.size() == rgrid_ext.size());
    QuantityType::const_iterator ri = rgrid_ext.begin();
    QuantityType::iterator rdfi = rdf_ext.begin();
    for (; ri != rgrid_ext.end(); ++ri, ++rdfi)
    {
        *rdfi = eps_gt(*ri, 0) ? (*rdfi / *ri) : 0.0;
    }
    return rdf_ext;
}


QuantityType PDFCalculator::getExtendedF() const
{
    QuantityType rdfperr_ext = this->getExtendedRDFperR();
    QuantityType rgrid_ext = this->getExtendedRgrid();
    QuantityType rdfperr_ext1 = this->applyBaseline(rgrid_ext, rdfperr_ext);
    const double rmin_ext = this->getExtendedRmin();
    QuantityType rv = fftgtof(rdfperr_ext1, this->getRstep(), rmin_ext);
    assert(rv.empty() || eps_eq(M_PI,
                this->getQstep() * rv.size() * this->getRstep()));
    return rv;
}


QuantityType PDFCalculator::getExtendedRgrid() const
{
    QuantityType rv;
    rv.reserve(this->countExtendedPoints());
    // make sure exact value of rmin will be in the extended grid
    for (int i = this->extendedRminSteps(); i < this->extendedRmaxSteps(); ++i)
    {
        rv.push_back(i * this->getRstep());
    }
    assert(rv.empty() || !eps_lt(rv.front(), this->getExtendedRmin()));
    assert(rv.empty() || !eps_gt(rv.back(), this->getExtendedRmax()));
    return rv;
}

// Q-range methods

QuantityType PDFCalculator::getQgrid() const
{
    return pdfutils_getQgrid(this);
}

// Q-range configuration

void PDFCalculator::setQmin(double qmin)
{
    ensureNonNegative("Qmin", qmin);
    mqmin = qmin;
}


const double& PDFCalculator::getQmin() const
{
    return mqmin;
}


void PDFCalculator::setQmax(double qmax)
{
    ensureNonNegative("Qmax", qmax);
    double qmax1 = (qmax > 0.0) ? qmax : DOUBLE_MAX;
    if (qmax1 < mqmax)  mticker.click();
    mqmax = qmax1;
    if (_initialQstepUpdate(this))  this->resetValue();
}


const double& PDFCalculator::getQmax() const
{
    static double rv;
    rv = min(mqmax, M_PI / this->getRstep());
    return rv;
}


const double& PDFCalculator::getQstep() const
{
    static double rv;
    // replicate the zero padding as done in fftgtof
    int Npad1 = this->extendedRmaxSteps();
    int Npad2 = (Npad1 > 0) ? (1 << int(ceil(log2(Npad1)))) : 0;
    rv = (Npad2 > 0) ? M_PI / (Npad2 * this->getRstep()) : 0.0;
    return rv;
}

// R-range methods

QuantityType PDFCalculator::getRgrid() const
{
    return pdfutils_getRgrid(this);
}

// R-range configuration

void PDFCalculator::setRmin(double rmin)
{
    ensureNonNegative("Rmin", rmin);
    this->PairQuantity::setRmin(rmin);
}


void PDFCalculator::setRmax(double rmax)
{
    ensureNonNegative("Rmax", rmax);
    this->PairQuantity::setRmax(rmax);
    if (_initialQstepUpdate(this))  this->resetValue();
}


void PDFCalculator::setRstep(double rstep)
{
    ensureEpsilonPositive("Rstep", rstep);
    if (mrstep != rstep)  mticker.click();
    mrstep = rstep;
    if (_initialQstepUpdate(this))  this->resetValue();
}


const double& PDFCalculator::getRstep() const
{
    return mrstep;
}


void PDFCalculator::setMaxExtension(double maxextension)
{
    ensureNonNegative("maxextension", maxextension);
    if (mmaxextension != maxextension)  mticker.click();
    mmaxextension = maxextension;
}


const double& PDFCalculator::getMaxExtension() const
{
    return mmaxextension;
}


double PDFCalculator::getExtendedRmin() const
{
    double rv = this->extendedRminSteps() * this->getRstep();
    return rv;
}


double PDFCalculator::getExtendedRmax() const
{
    double rv = this->extendedRmaxSteps() * this->getRstep();
    return rv;
}

// PDF peak profile configuration

void PDFCalculator::setPeakProfile(PeakProfilePtr pkf)
{
    ensureNonNull("PeakProfile", pkf);
    if (mpeakprofile != pkf)  mticker.click();
    mpeakprofile = pkf;
}


void PDFCalculator::setPeakProfileByType(const string& tp)
{
    PeakProfilePtr pkf = PeakProfile::createByType(tp);
    if (mpeakprofile.get())
    {
        pkf->setPrecision(mpeakprofile->getPrecision());
    }
    this->setPeakProfile(pkf);
}


PeakProfilePtr& PDFCalculator::getPeakProfile()
{
    assert(mpeakprofile.get());
    return mpeakprofile;
}


const PeakProfilePtr& PDFCalculator::getPeakProfile() const
{
    assert(mpeakprofile.get());
    return mpeakprofile;
}

// PDF baseline methods

QuantityType PDFCalculator::applyBaseline(
        const QuantityType& x, const QuantityType& y) const
{
    assert(x.size() == y.size());
    QuantityType z = y;
    const PDFBaseline& baseline = *(this->getBaseline());
    QuantityType::const_iterator xi = x.begin();
    QuantityType::iterator zi = z.begin();
    for (; xi != x.end(); ++xi, ++zi)
    {
        *zi += baseline(*xi);
    }
    return z;
}


void PDFCalculator::setBaseline(PDFBaselinePtr baseline)
{
    ensureNonNull("PDFBaseline", baseline);
    mbaseline = baseline;
}


void PDFCalculator::setBaselineByType(const std::string& tp)
{
    mbaseline = PDFBaseline::createByType(tp);
}


PDFBaselinePtr& PDFCalculator::getBaseline()
{
    assert(mbaseline.get());
    return mbaseline;
}


const PDFBaselinePtr& PDFCalculator::getBaseline() const
{
    assert(mbaseline.get());
    return mbaseline;
}

// Protected Methods ---------------------------------------------------------

// Attributes overloads

namespace {

template <class T>
void pdfcalc_accept(T* obj, diffpy::BaseAttributesVisitor& v)
{
    obj->getPeakWidthModel()->accept(v);
    obj->getPeakProfile()->accept(v);
    obj->getBaseline()->accept(v);
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


void PDFCalculator::accept(diffpy::BaseAttributesVisitor& v)
{
    pdfcalc_accept(this, v);
}


void PDFCalculator::accept(diffpy::BaseAttributesVisitor& v) const
{
    pdfcalc_accept(this, v);
}

// PairQuantity overloads

void PDFCalculator::resetValue()
{
    // calcPoints requires that structure and rlimits data are cached.
    this->cacheStructureData();
    this->cacheRlimitsData();
    // when applicable, configure linear baseline
    if (this->getBaseline()->type() == "linear")
    {
        double partialpdfscale =
            (0.0 == mstructure_cache.totaloccupancy) ? 0.0 :
            mstructure_cache.activeoccupancy / mstructure_cache.totaloccupancy;
        double pnumdensity = partialpdfscale * mstructure->numberDensity();
        PDFBaseline& bl = *(this->getBaseline());
        bl.setDoubleAttr("slope", -4 * M_PI * pnumdensity);
    }
    this->resizeValue(this->countCalcPoints());
    this->PairQuantity::resetValue();
}


void PDFCalculator::configureBondGenerator(BaseBondGenerator& bnds) const
{
    bnds.setRmin(this->rcalclo());
    bnds.setRmax(this->rcalchi());
}


void PDFCalculator::addPairContribution(const BaseBondGenerator& bnds,
        int summationscale)
{
    double sfprod = this->sfSite(bnds.site0()) * this->sfSite(bnds.site1());
    double peakscale = sfprod * bnds.multiplicity() * summationscale;
    double fwhm = this->getPeakWidthModel()->calculate(bnds);
    const PeakProfile& pkf = *(this->getPeakProfile());
    double dist = bnds.distance();
    double xlo = dist + pkf.xboundlo(fwhm);
    double xhi = dist + pkf.xboundhi(fwhm);
    int i = max(0, this->calcIndex(xlo));
    int ilast = min(this->countCalcPoints(), this->calcIndex(xhi) + 1);
    assert(ilast <= int(mvalue.size()));
    assert(eps_gt(dist, 0.0));
    for (; i < ilast; ++i)
    {
        double x = (this->rcalcloSteps() + i) * this->getRstep() - dist;
        double y = pkf(x, fwhm);
        // Contributions in G(r) need to be normalized by pair distance,
        // not by r as done in PDFfit or PDFfit2.  Here we rescale RDF
        // in such way that division by r will give a correct result.
        double yrdf = y * (x / dist + 1);
        mvalue[i] += peakscale * yrdf;
    }
}


void PDFCalculator::stashPartialValue()
{
    mstashedvalue.value = this->value();
    mstashedvalue.rclosteps = this->rcalcloSteps();
}


void PDFCalculator::restorePartialValue()
{
    assert(!mstashedvalue.value.empty());
    assert(!mvalue.empty());
    QuantityType::const_iterator si = mstashedvalue.value.begin();
    QuantityType::const_iterator slast = mstashedvalue.value.end();
    QuantityType::iterator ti = mvalue.begin();
    QuantityType::iterator tlast = mvalue.end();
    int leftshift = this->rcalcloSteps() - mstashedvalue.rclosteps;
    int sz = mstashedvalue.value.size();
    if (leftshift >= 0)  si += min(leftshift, sz);
    else  ti += min(-leftshift, int(mvalue.size()));
    for (; si != slast && ti != tlast; ++si, ++ti)  *ti = *si;
    mstashedvalue.value.clear();
}

// calculation specific

double PDFCalculator::rcalclo() const
{
    double rv = this->rcalcloSteps() * this->getRstep();
    return rv;
}


double PDFCalculator::rcalchi() const
{
    double rv = this->rcalchiSteps() * this->getRstep();
    return rv;
}


double PDFCalculator::extFromTerminationRipples() const
{
    // number of termination ripples for extending the r-range
    const int nripples = 6;
    // extension due to termination ripples.
    // apply only when qmax is below the Nyquist frequency for rstep.
    const double& qmax = this->getQmax();
    const double& dr = this->getRstep();
    double rv = (eps_gt(qmax, 0.0) && eps_lt(qmax, M_PI / dr)) ?
        (nripples * 2 * M_PI / qmax) : 0.0;
    return rv;
}


double PDFCalculator::extFromPeakTails() const
{
    // assume uncorrelated neighbors with maxUii
    const PeakWidthModel& pwm = *(this->getPeakWidthModel());
    double maxfwhm = pwm.maxWidth(
            mstructure, this->getRmin(), this->getRmax());
    const PeakProfile& pkf = *(this->getPeakProfile());
    double xleft = fabs(pkf.xboundlo(maxfwhm));
    double xright = fabs(pkf.xboundhi(maxfwhm));
    double rv = max(xleft, xright);
    return rv;
}


int PDFCalculator::rcalcloSteps() const
{
    return mrlimits_cache.rcalclosteps;
}


int PDFCalculator::rcalchiSteps() const
{
    return mrlimits_cache.rcalchisteps;
}


int PDFCalculator::extendedRminSteps() const
{
    return mrlimits_cache.extendedrminsteps;
}


int PDFCalculator::extendedRmaxSteps() const
{
    return mrlimits_cache.extendedrmaxsteps;
}


int PDFCalculator::countExtendedPoints() const
{
    int rv = this->extendedRmaxSteps() - this->extendedRminSteps();
    assert(rv >= 0);
    return rv;
}


int PDFCalculator::countCalcPoints() const
{
    int rv = this->rcalchiSteps() - this->rcalcloSteps();
    assert(rv >= 0);
    return rv;
}


int PDFCalculator::calcIndex(double r) const
{
    int rv = int(floor(r / this->getRstep())) - this->rcalcloSteps();
    return rv;
}


void PDFCalculator::cutRipplePoints(QuantityType& y) const
{
    if (y.empty())  return;
    assert(int(y.size()) ==
            this->extendedRmaxSteps() - this->extendedRminSteps());
    int ncutlo = pdfutils_rminSteps(this) - this->extendedRminSteps();
    int ncuthi = this->extendedRmaxSteps() - pdfutils_rmaxSteps(this);
    assert(ncutlo + ncuthi <= int(y.size()));
    y.erase(y.end() - ncuthi, y.end());
    y.erase(y.begin(), y.begin() + ncutlo);
}


const double& PDFCalculator::sfSite(int siteidx) const
{
    assert(0 <= siteidx && siteidx < int(mstructure_cache.sfsite.size()));
    return mstructure_cache.sfsite[siteidx];
}


double PDFCalculator::sfAverage() const
{
    return mstructure_cache.sfaverage;
}


void PDFCalculator::cacheStructureData()
{
    int cntsites = this->countSites();
    // sfsite
    unordered_map<string, double> fcache;
    mstructure_cache.sfsite.resize(cntsites);
    const ScatteringFactorTablePtr sftable = this->getScatteringFactorTable();
    for (int i = 0; i < cntsites; ++i)
    {
        const string& smbl = mstructure->siteAtomType(i);
        unordered_map<string, double>::iterator ff = fcache.find(smbl);
        if (ff == fcache.end())
        {
            const double value = sftable->lookup(smbl);
            ff = fcache.insert(make_pair(smbl, value)).first;
        }
        mstructure_cache.sfsite[i] = ff->second * mstructure->siteOccupancy(i);
    }
    // sfaverage
    double totocc = mstructure->totalOccupancy();
    double totsf = 0.0;
    for (int i = 0; i < cntsites; ++i)
    {
        totsf += this->sfSite(i) * mstructure->siteMultiplicity(i);
    }
    mstructure_cache.sfaverage = (totocc == 0.0) ? 0.0 : (totsf / totocc);
    // totaloccupancy
    mstructure_cache.totaloccupancy = totocc;
    // active occupancy
    double invmasktotal = 0.0;
    for (auto&& ij : minvertpairmask)
    {
        const int& i = ij.first;
        const int& j = ij.second;
        bool outofbounds = (i < 0 || i >= cntsites || j < 0 || j >= cntsites);
        if (outofbounds)  continue;
        int sumscale = (i == j) ? 1 : 2;
        double occij = sumscale *
            mstructure->siteOccupancy(i) * mstructure->siteMultiplicity(i) *
            mstructure->siteOccupancy(j) * mstructure->siteMultiplicity(j);
        invmasktotal += occij;
    }
    if (totocc > 0.0)   invmasktotal /= totocc;
    mstructure_cache.activeoccupancy = (mdefaultpairmask) ?
        (totocc - invmasktotal) : invmasktotal;
}


void PDFCalculator::cacheRlimitsData()
{
    mrlimits_cache.extendedrminsteps = 0;
    mrlimits_cache.extendedrmaxsteps = 0;
    mrlimits_cache.rcalclosteps = 0;
    mrlimits_cache.rcalchisteps = 0;
    if (pdfutils_rminSteps(this) >= pdfutils_rmaxSteps(this))   return;
    // obtain extension magnitudes and rescale to fit maximum extension
    double ext_ripples = this->extFromTerminationRipples();
    double ext_pktails = this->extFromPeakTails();
    double ext_total = ext_ripples + ext_pktails;
    if (ext_total > this->getMaxExtension())
    {
        double sc = this->getMaxExtension() / ext_total;
        ext_ripples *= sc;
        ext_pktails *= sc;
        ext_total = this->getMaxExtension();
    }
    // abbreviations
    const double& rmin = this->getRmin();
    const double& rmax = this->getRmax();
    const double& dr = this->getRstep();
    mrlimits_cache.extendedrminsteps = max(0, pdfutils_rminSteps(rmin - ext_ripples, dr));
    mrlimits_cache.extendedrmaxsteps = pdfutils_rmaxSteps(rmax + ext_ripples, dr);
    mrlimits_cache.rcalclosteps = max(0, pdfutils_rminSteps(rmin - ext_total, dr));
    mrlimits_cache.rcalchisteps = pdfutils_rmaxSteps(rmax + ext_total, dr);
}

}   // namespace diffpy
}   // namespace srreal

// Serialization -------------------------------------------------------------

DIFFPY_INSTANTIATE_SERIALIZATION(diffpy::srreal::PDFCalculator)
BOOST_CLASS_EXPORT_IMPLEMENT(diffpy::srreal::PDFCalculator)

// End of file
