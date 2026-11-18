/*****************************************************************************
 *
 * libdiffpy         by DANSE Diffraction group
 *                   Simon J. L. Billinge
 *                   (c) 2010 The Trustees of Columbia University
 *                   in the City of New York.  All rights reserved.
 *
 * File coded by:    Pavol Juhas
 *
 * See AUTHORS.rst for a list of people who contributed.
 * See LICENSE_DANSE.rst for license information.
 *
 ******************************************************************************
 *
 * class BaseDebyeSum -- base class for concrete Debye sum calculators
 *
 *****************************************************************************/

#include <algorithm>
#include <bit>
#include <span>
#include <cassert>
#include <typeinfo>
#include <cstdint>
#include <map>
#include <string>

#include <diffpy/srreal/BaseDebyeSum.hpp>
#include <diffpy/srreal/DebyePDFCalculator.hpp>
#include <diffpy/srreal/kernels/PairTable.hpp>
#include <diffpy/srreal/kernels/PairGeometry.hpp>
#include <diffpy/mathutils.hpp>
#include <diffpy/validators.hpp>
#include <diffpy/serialization.ipp>

using namespace std;
using namespace diffpy::validators;
using diffpy::mathutils::eps_eq;
using diffpy::mathutils::eps_gt;

namespace diffpy {
namespace srreal {

// Declaration of Local Helpers ----------------------------------------------

namespace {

/// Default cutoff for the Q-decreasing scale of the sine contributions.
const double DEFAULT_DEBYE_PRECISION = 1e-6;

std::span<const double> factorRow(const QuantityType& factors, int row,
                                  int nq) {
  return std::span(factors).subspan(std::size_t(row) * nq, nq);
}

}  // namespace

// Constructor ---------------------------------------------------------------

BaseDebyeSum::BaseDebyeSum()
    : mqmin(0.0),
      mqmax(DEFAULT_QGRID_QMAX),
      mqstep(DEFAULT_QGRID_QSTEP),
      mdebyeprecision(DEFAULT_DEBYE_PRECISION) {
  mstructure_cache.totaloccupancy = 0.0;
  // default configuration
  this->setPeakWidthModelByType("jeong");
  this->setEvaluatorType(OPTIMIZED);
  // attributes
  this->registerDoubleAttribute("debyeprecision", this,
                                &BaseDebyeSum::getDebyePrecision,
                                &BaseDebyeSum::setDebyePrecision);
}

// Public Methods ------------------------------------------------------------

// PairQuantity overloads

eventticker::EventTicker& BaseDebyeSum::ticker() const {
  mticker.updateFrom(this->PeakWidthModelOwner::ticker());
  return mticker;
}

// results

QuantityType BaseDebyeSum::getF() const {
  QuantityType rv = this->value();
  const double& totocc = mstructure_cache.totaloccupancy;
  const int npts = pdfutils_qmaxSteps(this);
  for (int kq = pdfutils_qminSteps(this); kq < npts; ++kq) {
    double sfavg = this->sfAverageAtkQ(kq);
    double fscale =
      (sfavg * totocc) == 0 ? 0.0 : 1.0 / (sfavg * sfavg * totocc);
    rv[kq] *= fscale;
  }
  return rv;
}

// Q-range methods

QuantityType BaseDebyeSum::getQgrid() const { return pdfutils_getQgrid(this); }

// Q-range configuration

void BaseDebyeSum::setQmin(double qmin) {
  ensureNonNegative("Qmin", qmin);
  if (mqmin != qmin) mticker.click();
  mqmin = qmin;
}

const double& BaseDebyeSum::getQmin() const { return mqmin; }

void BaseDebyeSum::setQmax(double qmax) {
  ensureNonNegative("Qmax", qmax);
  if (mqmax != qmax) mticker.click();
  mqmax = qmax;
}

const double& BaseDebyeSum::getQmax() const { return mqmax; }

void BaseDebyeSum::setQstep(double qstep) {
  ensureEpsilonPositive("Qstep", qstep);
  if (mqstep != qstep) mticker.click();
  mqstep = qstep;
}

const double& BaseDebyeSum::getQstep() const { return mqstep; }

void BaseDebyeSum::setDebyePrecision(double precision) {
  if (mdebyeprecision != precision) mticker.click();
  mdebyeprecision = precision;
}

const double& BaseDebyeSum::getDebyePrecision() const {
  return mdebyeprecision;
}

// Protected Methods ---------------------------------------------------------

// PairQuantity overloads

void BaseDebyeSum::resetValue() {
  this->cacheStructureData();
  this->resizeValue(pdfutils_qmaxSteps(this));
  this->PairQuantity::resetValue();
}

void BaseDebyeSum::addPairContribution(const BaseBondGenerator& bnds,
                                       int summationscale) {
  const double dist = bnds.distance();
  if (eps_eq(0.0, dist)) return;
  const double fwhm = this->getPeakWidthModel()->calculate(bnds);
  const int nqpts = pdfutils_qmaxSteps(this);
  const auto& cache = mstructure_cache;
  const int t0 = cache.types.typeofsite[bnds.site0()];
  const int t1 = cache.types.typeofsite[bnds.site1()];
  if (!nqpts) return;

  kernels::activeCalculationKernel().debyeSum({
    .out = mvalue.data(),
    .kqlo = pdfutils_qminSteps(this),
    .kqhi = nqpts,
    .qstep = this->getQstep(),
    .dist = dist,
    .dwsigma = kernels::fwhmToSigma(fwhm),
    .scale = double(summationscale * bnds.multiplicity()),
    .prec = this->getDebyePrecision(),
    .sf0 = factorRow(cache.factors, cache.types.curve[t0], nqpts).data(),
    .sf1 = factorRow(cache.factors, cache.types.curve[t1], nqpts).data(),
    .occ0 = cache.types.occupancy[t0],
    .occ1 = cache.types.occupancy[t1],
  });
}

bool BaseDebyeSum::addAllPairContributions() {
  if (typeid(*this) != typeid(BaseDebyeSum)) return false;
  return this->addAllPairsOfDebyeSum();
}

bool BaseDebyeSum::addAllPairsOfDebyeSum() {
  // Other peak width models and structures need per-bond evaluation.
  kernels::DebyePairSumArgs args = {};
  if (!kernels::assignPairWidth(args.w, *this->getPeakWidthModel()) ||
      !kernels::PairTable::supports(*mstructure))
    return false;
  BaseBondGenerator range(mstructure);
  this->configureBondGenerator(range);
  const auto& cache = mstructure_cache;
  const int nq = cache.nq;
  args.typeofsite = cache.types.typeofsite.data();
  args.curve = cache.types.curve.data();
  args.occupancy = cache.types.occupancy.data();
  args.ntypes = int(cache.types.curve.size());
  args.ncurves = nq ? int(cache.factors.size() / nq) : 0;
  args.nq = nq;
  args.sf = cache.factors.data();
  args.kqlo = pdfutils_qminSteps(this);
  args.kqhi = nq;
  args.qstep = this->getQstep();
  args.prec = this->getDebyePrecision();
  const auto& kernel =
    this->preparePairSum(args, range.getRmin(), range.getRmax());
  mallpairs.record(kernel.debyePairSum(args));
  return true;
}

void BaseDebyeSum::stashPartialValue() { mdbsumstash = this->value(); }

void BaseDebyeSum::restorePartialValue() {
  assert(mdbsumstash.size() == mvalue.size());
  mvalue.swap(mdbsumstash);
  mdbsumstash.clear();
}

double BaseDebyeSum::sfSiteAtQ(int siteidx, const double& q) const {
  return 1.0;
}

// Private Methods -----------------------------------------------------------

double BaseDebyeSum::sfAverageAtkQ(int kq) const {
  assert(kq < int(mstructure_cache.sfaverageatkq.size()));
  return mstructure_cache.sfaverageatkq[kq];
}

void BaseDebyeSum::cacheStructureData() {
  auto& c = mstructure_cache;
  const int nq = pdfutils_qmaxSteps(this);
  SiteTypes sites;
  sites.typeofsite.resize(this->countSites());
  if (typeid(*this) == typeid(DebyePDFCalculator) ||
      typeid(*this) == typeid(BaseDebyeSum)) {
    this->cacheElementFactors(sites, nq);
  } else {
    this->cacheCustomFactors(sites, nq);
  }
  c.types = std::move(sites);
  c.nq = nq;
  c.totaloccupancy = mstructure->totalOccupancy();
  c.sfaverageatkq.assign(c.nq, 0.0);
  std::vector<int> multiplicity(c.types.curve.size(), 0);
  for (int i = 0; i < this->countSites(); ++i) {
    multiplicity[c.types.typeofsite[i]] += mstructure->siteMultiplicity(i);
  }
  const int lo = pdfutils_qminSteps(this);
  for (size_t t = 0; t < multiplicity.size(); ++t) {
    const auto row = factorRow(c.factors, c.types.curve[t], c.nq);
    for (int k = lo; k < c.nq; ++k) {
      c.sfaverageatkq[k] += (row[k] * c.types.occupancy[t]) * multiplicity[t];
    }
  }
  if (eps_gt(c.totaloccupancy, 0.0)) {
    const double scale = 1.0 / c.totaloccupancy;
    for (double& f : c.sfaverageatkq) f *= scale;
  }
}

void BaseDebyeSum::cacheElementFactors(SiteTypes& sites, int nq) {
  auto& c = mstructure_cache;
  std::map<SiteKey, int> types;
  std::map<std::string, int> elements;
  const auto* pdf = dynamic_cast<const DebyePDFCalculator*>(this);
  for (int i = 0; i < this->countSites(); ++i) {
    const auto& symbol = mstructure->siteAtomType(i);
    const double occ = pdf ? mstructure->siteOccupancy(i) : 1.0;
    auto [t, added] =
      types.try_emplace({symbol, bit_cast<uint64_t>(occ)}, int(types.size()));
    sites.typeofsite[i] = t->second;
    if (added) {
      auto el = elements.try_emplace(symbol, int(elements.size())).first;
      sites.curve.push_back(el->second);
      sites.occupancy.push_back(occ);
    }
  }
  const auto tick = this->ticker();
  if (elements == melements && tick == msfticker &&
      c.factors.size() == elements.size() * nq)
    return;
  QuantityType factors(elements.size() * nq);
  const int lo = pdfutils_qminSteps(this);
  for (const auto& [symbol, row] : elements) {
    for (int k = lo; k < nq; ++k) {
      factors[size_t(row) * nq + k] =
        pdf ? pdf->getScatteringFactorTable()->lookup(symbol,
                                                      k * this->getQstep())
            : 1.0;
    }
  }
  c.factors.swap(factors);
  melements.swap(elements);
  msfticker = tick;
}

void BaseDebyeSum::cacheCustomFactors(SiteTypes& sites, int nq) {
  auto& c = mstructure_cache;
  std::map<SiteKey, int> types;
  QuantityType factors, row(nq);
  const bool comparePrevious = nq == c.nq && this->ticker() == msfticker;
  const int lo = pdfutils_qminSteps(this);
  bool changed = mcustomAmbiguous;
  bool ambiguous = false;
  for (int i = 0; i < this->countSites(); ++i) {
    SiteKey key{mstructure->siteAtomType(i),
                bit_cast<uint64_t>(mstructure->siteOccupancy(i))};
    for (int k = lo; k < nq; ++k)
      row[k] = this->sfSiteAtQ(i, k * this->getQstep());
    auto old = mcustomTypes.find(key);
    if (comparePrevious && old != mcustomTypes.end() &&
        !std::ranges::equal(row, factorRow(c.factors, old->second, nq)))
      changed = true;
    const int next = int(sites.curve.size());
    auto [t, added] = types.try_emplace(std::move(key), next);
    int id = t->second;
    if (!added && !std::ranges::equal(row, factorRow(factors, id, nq))) {
      id = next;
      ambiguous = true;
    }
    if (id == next) {
      factors.insert(factors.end(), row.begin(), row.end());
      sites.curve.push_back(id);
      sites.occupancy.push_back(1.0);
    }
    sites.typeofsite[i] = id;
  }
  if (changed || ambiguous) mticker.click();
  msfticker = this->ticker();
  mcustomAmbiguous = ambiguous;
  mcustomTypes.swap(types);
  c.factors.swap(factors);
}

std::vector<QuantityType> BaseDebyeSum::archiveFactors() const {
  const auto& c = mstructure_cache;
  std::vector<QuantityType> rows(c.types.curve.size(), QuantityType(c.nq));
  for (size_t t = 0; t < rows.size(); ++t) {
    const auto row = factorRow(c.factors, c.types.curve[t], c.nq);
    std::ranges::transform(row, rows[t].begin(),
                           [&](double f) { return f * c.types.occupancy[t]; });
  }
  return rows;
}

void BaseDebyeSum::loadFactors(const std::vector<QuantityType>& rows) {
  auto& c = mstructure_cache;
  c.nq = rows.empty() ? 0 : int(rows.front().size());
  c.factors.clear();
  c.types.curve.resize(rows.size());
  c.types.occupancy.assign(rows.size(), 1.0);
  for (size_t t = 0; t < rows.size(); ++t) {
    c.types.curve[t] = int(t);
    c.factors.insert(c.factors.end(), rows[t].begin(), rows[t].end());
  }
  melements.clear();
  mcustomTypes.clear();
  mcustomAmbiguous = false;
}

}  // namespace srreal
}  // namespace diffpy

// Serialization -------------------------------------------------------------

DIFFPY_INSTANTIATE_SERIALIZATION(diffpy::srreal::BaseDebyeSum)
BOOST_CLASS_EXPORT_IMPLEMENT(diffpy::srreal::BaseDebyeSum)

// End of file
