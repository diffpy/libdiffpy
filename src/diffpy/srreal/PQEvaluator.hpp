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

#ifndef PQEVALUATOR_HPP_INCLUDED
#define PQEVALUATOR_HPP_INCLUDED

#include <memory>
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/assume_abstract.hpp>
#include <boost/serialization/export.hpp>

#include <diffpy/EventTicker.hpp>
#include <diffpy/srreal/QuantityType.hpp>
#include <diffpy/srreal/StructureAdapter.hpp>

#include <diffpy/Export.hpp>

namespace diffpy {
namespace srreal {

class PairQuantity;

namespace kernels {
struct CalculationKernel;
}

/// shared pointer to PQEvaluatorBasic

typedef std::shared_ptr<class PQEvaluatorBasic> PQEvaluatorPtr;

enum PQEvaluatorType { NONE, BASIC, OPTIMIZED, CHECK };

enum PQEvaluatorFlag {
  // sum over full matrix of atom pairs, use pair symmetry otherwise.
  USEFULLSUM = 1,
  // allow fast updates only if unchanged atoms keep their indices.
  FIXEDSITEINDEX = 2,
};

class DLL_EXPORT PQEvaluatorBasic {
 public:
  friend DLL_EXPORT PQEvaluatorPtr createPQEvaluator(PQEvaluatorType,
                                                     PQEvaluatorPtr);
  // constructor
  PQEvaluatorBasic();
  virtual ~PQEvaluatorBasic() {}

  // methods
  virtual PQEvaluatorType typeint() const;
  PQEvaluatorType typeintused() const;
  virtual void updateValue(PairQuantity&, StructureAdapterPtr);
  virtual void validate(PairQuantity&) const;
  void setFlag(PQEvaluatorFlag flag, bool value);
  bool getFlag(PQEvaluatorFlag flag) const;
  void setupParallelRun(int cpuindex, int ncpu);
  bool isParallel() const;

 protected:
  // data
  /// per-bit storage of boolean configuration flags
  int mconfigflags;
  /// zero-based index of this CPU
  int mcpuindex;
  /// total number of the CPU units
  int mncpu;
  /// ticker for recording when was the value updated
  eventticker::EventTicker mvalue_ticker;
  /// type of PQEvaluator that was actually used
  PQEvaluatorType mtypeused;
  /// the calculator added all pairs at once in the last full calculation
  bool maddedallpairs = false;

 private:
  // serialization
  friend class boost::serialization::access;
  template <class Archive>
  void serialize(Archive& ar, const unsigned int version) {
    ar & mconfigflags & mcpuindex & mncpu & mvalue_ticker;
  }
};

class StructureDifference;

class DLL_EXPORT PQEvaluatorOptimized : public PQEvaluatorBasic {
 public:
  // methods
  virtual PQEvaluatorType typeint() const;
  virtual void validate(PairQuantity&) const;
  virtual void updateValue(PairQuantity&, StructureAdapterPtr);

 private:
  // data
  StructureAdapterPtr mlast_structure;
  // Process-local kernel identity. A restored value is recomputed on first use.
  const kernels::CalculationKernel* mlastkernel = nullptr;

  // helper methods
  void updateValueCompletely(PairQuantity&, StructureAdapterPtr);
  bool allPairsTakeLessWork(const PairQuantity&,
                            const StructureDifference&) const;

  // serialization
  friend class boost::serialization::access;
  template <class Archive>
  void serialize(Archive& ar, const unsigned int version) {
    using boost::serialization::base_object;
    ar& base_object<PQEvaluatorBasic>(*this);
    ar & mlast_structure;
    if constexpr (Archive::is_loading::value) mlastkernel = nullptr;
  }
};

class DLL_EXPORT PQEvaluatorCheck : public PQEvaluatorOptimized {
 public:
  // methods
  virtual PQEvaluatorType typeint() const;
  virtual void updateValue(PairQuantity&, StructureAdapterPtr);

 private:
  // serialization
  friend class boost::serialization::access;
  template <class Archive>
  void serialize(Archive& ar, const unsigned int version) {
    using boost::serialization::base_object;
    ar& base_object<PQEvaluatorOptimized>(*this);
  }
};

// Factory function for PairQuantity evaluators ------------------------------

DLL_EXPORT PQEvaluatorPtr createPQEvaluator(
  PQEvaluatorType pqtp, PQEvaluatorPtr pqevsrc = PQEvaluatorPtr());

}  // namespace srreal
}  // namespace diffpy

// Serialization -------------------------------------------------------------

BOOST_SERIALIZATION_ASSUME_ABSTRACT(diffpy::srreal::PQEvaluatorBasic)
BOOST_CLASS_EXPORT_KEY(diffpy::srreal::PQEvaluatorBasic)
BOOST_CLASS_EXPORT_KEY(diffpy::srreal::PQEvaluatorOptimized)

#endif  // PQEVALUATOR_HPP_INCLUDED
