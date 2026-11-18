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
 * Helper functions used in other unit tests.
 *
 *****************************************************************************/

#ifndef TEST_HELPERS_HPP_INCLUDED
#define TEST_HELPERS_HPP_INCLUDED

#include <string>
#include <vector>
#include <catch2/catch_test_macros.hpp>
#include <diffpy/srreal/CalculationKernel.hpp>
#include <diffpy/srreal/QuantityType.hpp>
#include <diffpy/srreal/forwardtypes.hpp>

std::string prepend_tests_dir(const std::string& f);
std::string prepend_testdata_dir(const std::string& f);

diffpy::srreal::StructureAdapterPtr loadTestPeriodicStructure(
  const std::string& tailname);

/// max|actual - reference| / scale, where scale defaults to max|reference|.
/// Infinity for different sizes, a value that is not finite or zero scale.
double maxRelativeError(const diffpy::srreal::QuantityType& actual,
                        const diffpy::srreal::QuantityType& reference,
                        double scale = 0.0);

// Calculation kernels -------------------------------------------------------

/// Restore the kernel that was active when the test started.
class KernelGuard {
 public:
  KernelGuard() : msaved(diffpy::srreal::getCalculationKernel()) {}
  ~KernelGuard() { diffpy::srreal::setCalculationKernel(msaved); }

 private:
  std::string msaved;
};

/// Call fn() with each available CPU kernel active, then restore the kernel
/// that was active before.  Do not use SECTION in fn, it would run only with
/// the first kernel.
template <class Fn>
void testAllKernels(Fn&& fn) {
  KernelGuard guard;
  for (const std::string& kname :
       diffpy::srreal::availableCalculationKernels()) {
    if (kname.starts_with("cuda")) continue;
    diffpy::srreal::setCalculationKernel(kname);
    INFO("kernel " << kname);
    fn();
  }
}

// Sums of all pairs ---------------------------------------------------------

/// A calculator of the class T that adds its pairs one by one, as the
/// calculators of the library do in every class derived from them.
template <class T>
struct PairByPair : T {};

/// A calculator of the class T that adds all pairs at once where T does,
/// and tells whether it did.
template <class T>
class AllPairs : public T {
 public:
  /// True when the last full calculation added all pairs at once.
  bool addedAllPairs() const { return maddedallpairs; }
  /// True when it added them on the GPU.
  bool addedOnDevice() const {
    return maddedallpairs && this->mallpairs.ondevice;
  }

 protected:
  // every calculation starts with the value reset
  void resetValue() override {
    maddedallpairs = false;
    this->T::resetValue();
  }

  bool addAllPairContributions() override {
    if constexpr (requires { this->addAllPairsOfPDF(); }) {
      maddedallpairs = this->addAllPairsOfPDF();
    } else {
      maddedallpairs = this->addAllPairsOfDebyeSum();
    }
    return maddedallpairs;
  }

 private:
  bool maddedallpairs = false;
};

#endif  // TEST_HELPERS_HPP_INCLUDED
