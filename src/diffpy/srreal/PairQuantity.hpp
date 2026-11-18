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
 * class PairQuantity -- general implementation of pair quantity calculator
 *
 *****************************************************************************/

#ifndef PAIRQUANTITY_HPP_INCLUDED
#define PAIRQUANTITY_HPP_INCLUDED

#include <boost/serialization/base_object.hpp>
#include <boost/serialization/assume_abstract.hpp>
#include <boost/serialization/utility.hpp>
#include <boost/serialization/unordered_set.hpp>
#include <boost/serialization/unordered_map.hpp>

#include <diffpy/hash.hpp>
#include <diffpy/srreal/PQEvaluator.hpp>
#include <diffpy/srreal/StructureAdapter.hpp>
#include <diffpy/srreal/QuantityType.hpp>
#include <diffpy/Attributes.hpp>
#include <diffpy/EventTicker.hpp>

#include <diffpy/Export.hpp>

namespace diffpy {
namespace srreal {

class BaseBondGenerator;

namespace kernels {
struct PairSumState;
struct PairSumArgs;
struct PairSumStats;
struct CalculationKernel;
}

class DLL_EXPORT PairQuantity : public diffpy::Attributes {
 public:
  // class constants
  static const int ALLATOMSINT;
  static const std::string ALLATOMSSTR;

  // constructor
  PairQuantity();
  virtual ~PairQuantity() {}

  // methods
  const QuantityType& eval();
  template <class T>
  const QuantityType& eval(const T&);
  const QuantityType& eval(StructureAdapterPtr);
  const QuantityType& value() const;
  void mergeParallelData(const std::string& pdata, int ncpu);
  virtual std::string getParallelData() const;

  // configuration
  template <class T>
  void setStructure(const T&);
  void setStructure(StructureAdapterPtr);
  StructureAdapterPtr& getStructure();
  const StructureAdapterPtr& getStructure() const;
  virtual void setRmin(double);
  const double& getRmin() const;
  virtual void setRmax(double);
  const double& getRmax() const;
  void setEvaluatorType(PQEvaluatorType evtp);
  PQEvaluatorType getEvaluatorType() const;
  PQEvaluatorType getEvaluatorTypeUsed() const;
  void setupParallelRun(int cpuindex, int ncpu);
  void maskAllPairs(bool mask);
  void invertMask();
  void setPairMask(int i, int j, bool mask);
  bool getPairMask(int i, int j) const;
  void setTypeMask(std::string, std::string, bool mask);
  bool getTypeMask(const std::string&, const std::string&) const;

  // ticker for any updates in configuration
  virtual eventticker::EventTicker& ticker() const { return mticker; }

 protected:
  friend class PQEvaluatorBasic;
  friend class PQEvaluatorOptimized;
  friend DLL_EXPORT StructureAdapterPtr
  replacePairQuantityStructure(PairQuantity&, StructureAdapterPtr);

  // methods
  virtual void resizeValue(size_t);
  virtual void resetValue();
  virtual void configureBondGenerator(BaseBondGenerator&) const;
  virtual void addPairContribution(const BaseBondGenerator&, int) {}
  /// For an unpartitioned half sum, add every unmasked pair as
  /// addPairContribution would. Return false before changing the value
  /// to use per-bond evaluation. Derived calculators can opt in when
  /// their batch sum matches their per-bond contributions.
  virtual bool addAllPairContributions() { return false; }
  const kernels::CalculationKernel& preparePairSum(kernels::PairSumArgs&,
                                                   double rmin, double rmax);
  virtual void executeParallelMerge(const std::string& pdata);
  virtual void finishValue() {}
  int countSites() const;
  // support methods for PQEvaluatorOptimized
  bool hasMask() const;
  bool hasPairMask() const;
  bool hasTypeMask() const;
  virtual void stashPartialValue();
  virtual void restorePartialValue();

  // data
  typedef std::unordered_set<std::pair<int, int>, diffpy::PairHash>
    PairMaskStorage;
  typedef std::unordered_map<std::pair<std::string, std::string>, bool,
                             diffpy::PairHash>
    TypeMaskStorage;
  /// Pairs of the last structures for calculators that add all pairs at
  /// once.  Neither copied nor serialized.
  struct AllPairsState {
    AllPairsState() = default;
    void record(kernels::PairSumStats stats, bool usedDevice);
    AllPairsState(const AllPairsState& other)
        : ondevice(other.ondevice), pairs(other.pairs), peaks(other.peaks) {}
    AllPairsState& operator=(const AllPairsState& other) {
      pairsum.reset();
      ondevice = other.ondevice;
      pairs = other.pairs;
      peaks = other.peaks;
      return *this;
    }
    std::shared_ptr<kernels::PairSumState> pairsum;
    /// The last full sum ran on a GPU.
    bool ondevice = false;
    /// Pairs and merged contributions of the last sum of all pairs, from
    /// which the OPTIMIZED evaluator estimates the work of the next one.
    double pairs = 0.0;
    double peaks = 0.0;
  };

  QuantityType mvalue;
  StructureAdapterPtr mstructure;
  double mrmin;
  double mrmax;
  /// The evaluator of the calculator.  A copy gets a fresh evaluator, so
  /// that it does not take the last structure of the original for its own.
  class EvaluatorHolder {
   public:
    EvaluatorHolder() {}
    EvaluatorHolder(const EvaluatorHolder& other) : mptr(other.copy()) {}
    EvaluatorHolder& operator=(const EvaluatorHolder& other) {
      if (this != &other) mptr = other.copy();
      return *this;
    }
    EvaluatorHolder& operator=(PQEvaluatorPtr pqev) {
      mptr = pqev;
      return *this;
    }
    PQEvaluatorBasic* operator->() const { return mptr.get(); }
    PQEvaluatorBasic* get() const { return mptr.get(); }
    PQEvaluatorPtr& ptr() { return mptr; }

   private:
    PQEvaluatorPtr copy() const {
      return mptr ? createPQEvaluator(mptr->typeint(), mptr) : mptr;
    }
    PQEvaluatorPtr mptr;
  };

  EvaluatorHolder mevaluator;
  bool mdefaultpairmask;
  PairMaskStorage minvertpairmask;
  std::unordered_map<int, bool> msiteallmask;
  TypeMaskStorage mtypemask;
  int mmergedvaluescount;
  mutable eventticker::EventTicker mticker;
  AllPairsState mallpairs;

 private:
  // methods
  void updateMaskData();
  bool setPairMaskValue(int i, int j, bool mask);

  // serialization
  friend class boost::serialization::access;
  template <class Archive>
  void serialize(Archive& ar, const unsigned int version) {
    ar & mvalue;
    ar & mstructure;
    ar & mrmin;
    ar & mrmax;
    ar & mevaluator.ptr();
    ar & mdefaultpairmask;
    ar & minvertpairmask;
    ar & msiteallmask;
    ar & mtypemask;
    ar & mmergedvaluescount;
    ar & mticker;
  }
};

// Template Public Methods ---------------------------------------------------

template <class T>
const QuantityType& PairQuantity::eval(const T& stru) {
  StructureAdapterPtr pstru = convertToStructureAdapter(stru);
  return this->eval(pstru);
}

template <class T>
void PairQuantity::setStructure(const T& stru) {
  StructureAdapterPtr pstru = convertToStructureAdapter(stru);
  this->setStructure(pstru);
}

// Other functions -----------------------------------------------------------

/// The purpose of this function is to support Python pickling of
/// PairQuantity objects that hold Python-derived StructureAdapter classes.
/// Use it only if you absolutely have to and you know what you do.
DLL_EXPORT StructureAdapterPtr
replacePairQuantityStructure(PairQuantity& pq, StructureAdapterPtr stru);

}  // namespace srreal
}  // namespace diffpy

// Serialization -------------------------------------------------------------

BOOST_SERIALIZATION_ASSUME_ABSTRACT(diffpy::srreal::PairQuantity)
BOOST_CLASS_EXPORT_KEY(diffpy::srreal::PairQuantity)

#endif  // PAIRQUANTITY_HPP_INCLUDED
