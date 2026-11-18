/*****************************************************************************
 *
 * libdiffpy         by DANSE Diffraction group
 *                   Simon J. L. Billinge
 *                   (c) 2009 The Trustees of Columbia University
 *                   in the City of New York.  All rights reserved.
 *
 * File coded by:    Pavol Juhas
 *
 * See AUTHORS.rst for a list of people who contributed.
 * See LICENSE_DANSE.rst for license information.
 *
 ******************************************************************************
 *
 * class PQEvaluatorBasic -- robust PairQuantity evaluator, the result
 *     is always calculated from scratch.
 *
 * class PQEvaluatorOptimized -- optimized PairQuantity evaluator with fast
 *     quantity updates
 *
 *****************************************************************************/

#include <stdexcept>
#include <sstream>

#include <diffpy/serialization.ipp>
#include <diffpy/srreal/PQEvaluator.hpp>
#include <diffpy/srreal/PairQuantity.hpp>
#include <diffpy/srreal/BondCalculator.hpp>
#include <diffpy/srreal/StructureDifference.hpp>
#include <diffpy/srreal/kernels/CalculationKernelImpl.hpp>

using namespace std;

namespace diffpy {
namespace srreal {

// Local Constants and Routines ----------------------------------------------

namespace {

// tolerated load variance for splitting outer loop for parallel evaluation
const double CPU_LOAD_VARIANCE = 0.1;

// Work of a sum of all pairs for each of its pairs and for each of the
// contributions that it merges them to, in the units of the work for a
// pair that is added by itself.
const double ALLPAIRS_PAIR_WORK = 0.2;
const double ALLPAIRS_PEAK_WORK = 0.85;
// The estimate is within about 15 %, so the sum of all pairs must take at
// most this share of the work of the fast update.
const double ALLPAIRS_WORK_MARGIN = 0.85;

SiteIndices complementary_indices(const int sz, const SiteIndices& indices0) {
  SiteIndices rv;
  rv.reserve(sz);
  SiteIndices::const_iterator ii0 = indices0.begin();
  for (int k = 0; k < sz; ++k) {
    if (ii0 == indices0.end() || k < *ii0) {
      rv.push_back(k);
    } else {
      assert(k == *ii0);
      ++ii0;
    }
  }
  return rv;
}

}  // namespace

//////////////////////////////////////////////////////////////////////////////
// class PQEvaluatorBasic
//////////////////////////////////////////////////////////////////////////////

PQEvaluatorBasic::PQEvaluatorBasic()
    : mconfigflags(0), mcpuindex(0), mncpu(1), mtypeused(NONE) {}

PQEvaluatorType PQEvaluatorBasic::typeint() const { return BASIC; }

PQEvaluatorType PQEvaluatorBasic::typeintused() const { return mtypeused; }

void PQEvaluatorBasic::validate(PairQuantity& pq) const { return; }

void PQEvaluatorBasic::updateValue(PairQuantity& pq, StructureAdapterPtr stru) {
  mtypeused = BASIC;
  maddedallpairs = false;
  // The value is not valid until this update completes.
  mvalue_ticker = eventticker::EventTicker();
  pq.setStructure(stru);
  // A batch covers the complete half sum of an unpartitioned evaluation.
  const bool halfsum = !this->getFlag(USEFULLSUM) && !this->isParallel();
  if (halfsum && pq.addAllPairContributions()) {
    maddedallpairs = true;
    mvalue_ticker.click();
    return;
  }
  BaseBondGeneratorPtr bnds = pq.mstructure->createBondGenerator();
  pq.configureBondGenerator(*bnds);
  int cntsites = pq.mstructure->countSites();
  // loop counter
  long n = mcpuindex;
  // split outer loop for many atoms.  The CPUs should have similar load.
  bool chop_outer = (mncpu <= ((cntsites - 1) * CPU_LOAD_VARIANCE + 1));
  bool chop_inner = !chop_outer;
  const bool hasmask = pq.hasMask();
  if (!this->isParallel()) chop_outer = chop_inner = false;
  const bool usefullsum = this->getFlag(USEFULLSUM);
  for (int i0 = 0; i0 < cntsites; ++i0) {
    if (chop_outer && (n++ % mncpu)) continue;
    bnds->selectAnchorSite(i0);
    int i1hi = usefullsum ? cntsites : (i0 + 1);
    bnds->selectSiteRange(0, i1hi);
    for (bnds->rewind(); !bnds->finished(); bnds->next()) {
      if (chop_inner && (n++ % mncpu)) continue;
      int i1 = bnds->site1();
      if (hasmask && !pq.getPairMask(i0, i1)) continue;
      int summationscale = (usefullsum || i0 == i1) ? 1 : 2;
      pq.addPairContribution(*bnds, summationscale);
    }
  }
  mvalue_ticker.click();
}

void PQEvaluatorBasic::setFlag(PQEvaluatorFlag flag, bool value) {
  if (value)
    mconfigflags |= int(flag);
  else
    mconfigflags &= ~int(flag);
}

bool PQEvaluatorBasic::getFlag(PQEvaluatorFlag flag) const {
  return mconfigflags & int(flag);
}

void PQEvaluatorBasic::setupParallelRun(int cpuindex, int ncpu) {
  // make sure ncpu is at least one
  if (ncpu < 1) {
    const char* emsg = "Number of CPU ncpu must be at least 1.";
    throw invalid_argument(emsg);
  }
  mcpuindex = cpuindex;
  mncpu = ncpu;
}

bool PQEvaluatorBasic::isParallel() const { return mncpu > 1; }

//////////////////////////////////////////////////////////////////////////////
// class PQEvaluatorOptimized
//////////////////////////////////////////////////////////////////////////////

PQEvaluatorType PQEvaluatorOptimized::typeint() const { return OPTIMIZED; }

void PQEvaluatorOptimized::validate(PairQuantity& pq) const {
  // Check if PairQuantity provides stashPartialValue.
  pq.stashPartialValue();
  pq.restorePartialValue();
}

void PQEvaluatorOptimized::updateValue(PairQuantity& pq,
                                       StructureAdapterPtr stru) {
  mtypeused = OPTIMIZED;
  if (pq.ticker() >= mvalue_ticker || !mlast_structure ||
      mlastkernel != &kernels::activeCalculationKernel()) {
    return this->updateValueCompletely(pq, stru);
  }
  // do not do fast updates if they take more work
  StructureDifference sd = mlast_structure->diff(stru);
  if (!sd.allowsfastupdate()) {
    return this->updateValueCompletely(pq, stru);
  }
  if ((this->getFlag(FIXEDSITEINDEX) || pq.hasPairMask()) &&
      sd.diffmethod != StructureDifference::Method::SIDEBYSIDE) {
    return this->updateValueCompletely(pq, stru);
  }
  if (this->allPairsTakeLessWork(pq, sd)) {
    return this->updateValueCompletely(pq, stru);
  }
  // Remove contributions from the extra sites in the old structure
  assert(sd.stru0 == mlast_structure);
  // The value belongs to no structure until this update completes.
  const StructureAdapterPtr last_structure = mlast_structure;
  mlast_structure.reset();
  // loop counter
  long n = mcpuindex;
  bool usefullsum = this->getFlag(USEFULLSUM);
  const bool hasmask = pq.hasMask();
  // A structure of the same content has no sites to remove or to add.
  // It needs no bond generators, which copy all positions of a crystal.
  if (!sd.pop0.empty()) {
    int cntsites0 = sd.stru0->countSites();
    BaseBondGeneratorPtr bnds0 = sd.stru0->createBondGenerator();
    pq.configureBondGenerator(*bnds0);
    // the loop is adjusted according to usefullsum and split within
    // the outer loop in case of parallel evaluation.
    SiteIndices anchors = sd.pop0;
    SiteIndices unchanged = complementary_indices(cntsites0, sd.pop0);
    anchors.insert(anchors.end(), unchanged.begin(), unchanged.end());
    bnds0->selectSites(anchors.begin(), anchors.end());
    SiteIndices::const_iterator last_anchor =
      usefullsum ? anchors.end() : (anchors.begin() + sd.pop0.size());
    SiteIndices::const_iterator ii0;
    bool needsreselection = usefullsum;
    for (ii0 = anchors.begin(); ii0 != last_anchor; ++ii0) {
      if (n++ % mncpu) continue;
      const int& i0 = *ii0;
      bnds0->selectAnchorSite(i0);
      // when using half sum, deselect visited popped sites
      if (!usefullsum) bnds0->selectSites(ii0, anchors.end());
      // when using full sum, select only the popped sites when
      // anchored at an unchanged atom
      else if (needsreselection && ii0 >= (anchors.begin() + sd.pop0.size())) {
        bnds0->selectSites(sd.pop0.begin(), sd.pop0.end());
        needsreselection = false;
      }
      for (bnds0->rewind(); !bnds0->finished(); bnds0->next()) {
        int i1 = bnds0->site1();
        if (hasmask && !pq.getPairMask(i0, i1)) continue;
        const int summationscale = (usefullsum || i0 == i1) ? -1 : -2;
        pq.addPairContribution(*bnds0, summationscale);
      }
    }
  }
  // Add contributions from the new atoms in the updated structure
  // save current value to override the resetValue call from setStructure
  assert(sd.stru1);
  pq.stashPartialValue();
  // Fall back to a full calculation if assigning the new structure
  // changes the calculator's configuration ticker.
  assert(pq.ticker() < mvalue_ticker);
  pq.setStructure(sd.stru1);
  if (pq.ticker() >= mvalue_ticker) {
    return this->updateValueCompletely(pq, stru);
  }
  pq.restorePartialValue();
  if (!sd.add1.empty()) {
    int cntsites1 = sd.stru1->countSites();
    BaseBondGeneratorPtr bnds1 = sd.stru1->createBondGenerator();
    pq.configureBondGenerator(*bnds1);
    SiteIndices anchors = sd.add1;
    SiteIndices unchanged = complementary_indices(cntsites1, sd.add1);
    anchors.insert(anchors.begin(), unchanged.begin(), unchanged.end());
    bnds1->selectSites(sd.add1.begin(), sd.add1.end());
    SiteIndices::const_iterator first_anchor =
      usefullsum ? anchors.begin() : (anchors.end() - sd.add1.size());
    SiteIndices::const_iterator ii1;
    bool needsreselection = usefullsum;
    for (ii1 = first_anchor; ii1 != anchors.end(); ++ii1) {
      if (n++ % mncpu) continue;
      const int& i0 = *ii1;
      bnds1->selectAnchorSite(i0);
      // when using half sum, activate the added site
      if (!usefullsum) bnds1->selectSites(anchors.begin(), ii1 + 1);
      // when using full sum select unchanged atoms once anchored
      // at an added atom.
      else if (needsreselection && ii1 >= (anchors.end() - sd.add1.size())) {
        bnds1->selectSites(anchors.begin(), anchors.end());
        needsreselection = false;
      }
      for (bnds1->rewind(); !bnds1->finished(); bnds1->next()) {
        int i1 = bnds1->site1();
        if (hasmask && !pq.getPairMask(i0, i1)) continue;
        const int summationscale = (usefullsum || i0 == i1) ? +1 : +2;
        pq.addPairContribution(*bnds1, summationscale);
      }
    }
  }
  // The saved geometry must use the site order of the new scattering arrays.
  // A sorted comparison can report empty differences after a permutation.
  const bool sameorder =
    sd.diffmethod == StructureDifference::Method::SIDEBYSIDE &&
    sd.pop0.empty() && sd.add1.empty();
  mlast_structure = sameorder ? last_structure : pq.getStructure()->clone();
  mvalue_ticker.click();
}

void PQEvaluatorOptimized::updateValueCompletely(PairQuantity& pq,
                                                 StructureAdapterPtr stru) {
  mlast_structure.reset();
  this->PQEvaluatorBasic::updateValue(pq, stru);
  mlast_structure = pq.getStructure()->clone();
  mlastkernel = &kernels::activeCalculationKernel();
}

// True when a sum of all pairs is expected to take less work than the fast
// update of the changed sites.

bool PQEvaluatorOptimized::allPairsTakeLessWork(
  const PairQuantity& pq, const StructureDifference& sd) const {
  // A change of whether the calculator adds all pairs forces a full
  // calculation anyway.
  if (!maddedallpairs) return false;
  // share of the pairs that have one of the sites
  auto pairshare = [](const SiteIndices& sites, const StructureAdapter& s) {
    if (sites.empty()) return 0.0;
    const double others = 1.0 - double(sites.size()) / s.countSites();
    return 1.0 - others * others;
  };
  const double pairs = pq.mallpairs.pairs;
  const double fastwork =
    (pairshare(sd.pop0, *sd.stru0) + pairshare(sd.add1, *sd.stru1)) * pairs;
  const double fullwork =
    ALLPAIRS_PAIR_WORK * pairs + ALLPAIRS_PEAK_WORK * pq.mallpairs.peaks;
  return fullwork < ALLPAIRS_WORK_MARGIN * fastwork;
}

// Helper classes and functions for PQEvaluatorCheck -------------------------

namespace {

bool cmpfabs(const double& x, const double& y) { return fabs(x) < fabs(y); }

class pqresults {
 public:
  pqresults(const PairQuantity& pq) : msaved_value(pq.value()) {}

  virtual ~pqresults() {}

  virtual bool compare(const PairQuantity& pq) const {
    QuantityType::const_iterator smx =
      max_element(msaved_value.begin(), msaved_value.end(), cmpfabs);
    const double eps =
      R3::SQRT_DOUBLE_EPS * (msaved_value.empty() ? 1.0 : max(*smx, 1.0));
    mathutils::EpsilonEqual allclose(eps);
    const QuantityType& cur_value = pq.value();
    return allclose(msaved_value, cur_value);
  }

 private:
  // data
  QuantityType msaved_value;
};

class pqresults_bondcalculator : public pqresults {
 public:
  pqresults_bondcalculator(const BondCalculator& bc)
      : pqresults(bc),
        msaved_directions(bc.directions()),
        msaved_sites0(bc.sites0()),
        msaved_sites1(bc.sites1()),
        msaved_types0(bc.types0()),
        msaved_types1(bc.types1()) {}

  virtual bool compare(const PairQuantity& pq) const {
    if (!(this->pqresults::compare(pq))) return false;
    const BondCalculator& bc = dynamic_cast<const BondCalculator&>(pq);
    if (msaved_sites0 != bc.sites0()) return false;
    if (msaved_sites1 != bc.sites1()) return false;
    if (msaved_types0 != bc.types0()) return false;
    if (msaved_types1 != bc.types1()) return false;
    const vector<R3::Vector> cur_directions = bc.directions();
    const bool szneq = msaved_directions.size() != cur_directions.size();
    if (szneq) return false;
    vector<R3::Vector>::const_iterator sv, cv;
    sv = msaved_directions.begin();
    cv = cur_directions.begin();
    mathutils::EpsilonEqual allclose;
    for (; sv != msaved_directions.end(); ++sv, ++cv) {
      if (!allclose(*sv, *cv)) return false;
    }
    return true;
  }

 private:
  // data
  vector<R3::Vector> msaved_directions;
  SiteIndices msaved_sites0;
  SiteIndices msaved_sites1;
  vector<string> msaved_types0;
  vector<string> msaved_types1;
};

pqresults* create_pqresults(const PairQuantity& pq) {
  const BondCalculator* pbc = dynamic_cast<const BondCalculator*>(&pq);
  if (pbc) return new pqresults_bondcalculator(*pbc);
  return new pqresults(pq);
}

}  // namespace

//////////////////////////////////////////////////////////////////////////////
// class PQEvaluatorCheck
//////////////////////////////////////////////////////////////////////////////

PQEvaluatorType PQEvaluatorCheck::typeint() const { return CHECK; }

void PQEvaluatorCheck::updateValue(PairQuantity& pq, StructureAdapterPtr stru) {
  this->PQEvaluatorOptimized::updateValue(pq, stru);
  if (mtypeused == BASIC) return;
  unique_ptr<pqresults> results(create_pqresults(pq));
  this->PQEvaluatorBasic::updateValue(pq, stru);
  mtypeused = CHECK;
  if (!results->compare(pq)) {
    const char* emsg = "Inconsistent results from OPTIMIZED evaluation.";
    throw logic_error(emsg);
  }
}

// Factory for PairQuantity evaluators ---------------------------------------

PQEvaluatorPtr createPQEvaluator(PQEvaluatorType pqtp, PQEvaluatorPtr pqevsrc) {
  PQEvaluatorPtr rv;
  switch (pqtp) {
    case BASIC:
      rv.reset(new PQEvaluatorBasic());
      break;

    case OPTIMIZED:
      rv.reset(new PQEvaluatorOptimized());
      break;

    case CHECK:
      rv.reset(new PQEvaluatorCheck());
      break;

    default:
      ostringstream emsg;
      emsg << "Invalid PQEvaluatorType value " << pqtp;
      throw invalid_argument(emsg.str());
  }
  if (pqevsrc) {
    rv->mconfigflags = pqevsrc->mconfigflags;
    rv->mcpuindex = pqevsrc->mcpuindex;
    rv->mncpu = pqevsrc->mncpu;
    rv->mvalue_ticker = pqevsrc->mvalue_ticker;
    rv->mtypeused = pqevsrc->mtypeused;
    rv->maddedallpairs = pqevsrc->maddedallpairs;
  }
  return rv;
}

}  // namespace srreal
}  // namespace diffpy

// Serialization -------------------------------------------------------------

DIFFPY_INSTANTIATE_SERIALIZATION(diffpy::srreal::PQEvaluatorBasic)
BOOST_CLASS_EXPORT_IMPLEMENT(diffpy::srreal::PQEvaluatorBasic)
DIFFPY_INSTANTIATE_SERIALIZATION(diffpy::srreal::PQEvaluatorOptimized)
BOOST_CLASS_EXPORT_IMPLEMENT(diffpy::srreal::PQEvaluatorOptimized)

// End of file
