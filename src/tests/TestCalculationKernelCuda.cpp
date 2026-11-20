/*****************************************************************************
 *
 * libdiffpy         Complex Modeling Initiative
 *                   (c) 2026 The Trustees of Columbia University
 *                   in the City of New York.  All rights reserved.
 *
 * See AUTHORS.rst for a list of people who contributed.
 * See LICENSE.rst for license information.
 *
 ******************************************************************************
 *
 * TestCalculationKernelCuda -- the "cuda" and "cuda-mixed" kernels must
 * give the results of the CPU.  Skipped without a CUDA device.
 *
 *****************************************************************************/

#include <catch2/catch_test_macros.hpp>

#ifdef LIBDIFFPY_TEST_CUDA_RUNTIME
#include <cuda_runtime_api.h>
#endif

#include <cmath>
#include <memory>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

#include <diffpy/srreal/AtomicStructureAdapter.hpp>
#include <diffpy/srreal/BVSCalculator.hpp>
#include <diffpy/srreal/CalculationKernel.hpp>
#include <diffpy/srreal/CrystalStructureAdapter.hpp>
#include <diffpy/srreal/DebyePDFCalculator.hpp>
#include <diffpy/srreal/GaussianProfile.hpp>
#include <diffpy/srreal/PDFCalculator.hpp>
#include <diffpy/srreal/PeriodicStructureAdapter.hpp>

// The internal kernel header is missing when testing an installed library.
#if __has_include(<diffpy/srreal/kernels/CalculationKernelImpl.hpp>)
#include <diffpy/srreal/kernels/CalculationKernelImpl.hpp>
#define LIBDIFFPY_TEST_KERNEL_INTERNALS 1
#endif

#include "pair_test_helpers.hpp"

using namespace diffpy::srreal;
using namespace pairtests;

namespace {

/// Agreement with the CPU.  In mixed precision, the results other than the
/// raw PDF peaks cancel larger terms (Debye sums, PDF after the baseline)
/// and get the looser tolerance.
const double TOLERANCE = 1e-12;
const double MIXED_TOLERANCE = 1e-6;
const double MIXED_REST_TOLERANCE = 1e-5;

void skipWithoutDevice() {
  KernelGuard guard;
  try {
    setCalculationKernel("cuda");
  } catch (const std::invalid_argument& e) {
    SKIP(e.what());
  }
}

/// The kernel of the CPU that the calculators use by default.
std::string cpuKernel() { return availableCalculationKernels().front(); }

/// Evaluate with the kernel, then restore the one that was active.
const QuantityType& evalWith(const std::string& kernel, PairQuantity& calc,
                             StructureAdapterPtr stru) {
  KernelGuard guard;
  setCalculationKernel(kernel);
  return calc.eval(stru);
}

/// Ni atoms of the fcc lattice within a sphere.
std::shared_ptr<AtomicStructureAdapter> makeSphere(double radius) {
  auto rv = std::make_shared<AtomicStructureAdapter>();
  Atom a;
  a.atomtype = "Ni";
  a.uij_cartn = 0.005 * R3::identity();
  const double alat = 3.52;
  const int n = int(radius / alat) + 1;
  const R3::Vector basis[] = {R3::Vector(0, 0, 0), R3::Vector(0, 0.5, 0.5),
                              R3::Vector(0.5, 0, 0.5), R3::Vector(0.5, 0.5, 0)};
  for (int i = -n; i <= n; ++i) {
    for (int j = -n; j <= n; ++j) {
      for (int k = -n; k <= n; ++k) {
        for (const R3::Vector& b : basis) {
          a.xyz_cartn = alat * (b + R3::Vector(i, j, k));
          if (R3::norm(a.xyz_cartn) <= radius) rv->append(a);
        }
      }
    }
  }
  return rv;
}

/// Results of a calculator of the class T with the kernel.
template <class T>
std::vector<QuantityType> evaluate(const Setup& setup, StructureAdapterPtr s,
                                   const std::string& kernel) {
  KernelGuard guard;
  setCalculationKernel(kernel);
  AllPairs<T> calc;
  calc.setEvaluatorType(BASIC);
  setup(calc);
  calc.eval(s);
  CHECK(calc.addedAllPairs());
  CHECK((kernel != cpuKernel()) == calc.addedOnDevice());
  std::vector<QuantityType> rv = {calc.value(), calc.getPDF(), calc.getF()};
  // The same again to the last bit, also from the calculator of the
  // library.
  T again;
  again.setEvaluatorType(BASIC);
  setup(again);
  again.eval(s);
  CHECK(rv[0] == again.value());
  return rv;
}

/// Check that the device gives what the CPU does.
template <class T>
void checkAsCPU(const Setup& setup, StructureAdapterPtr s) {
  const auto expected = evaluate<T>(setup, s, cpuKernel());
  REQUIRE(*std::ranges::max_element(expected[0]) > 0.0);
  const auto results = evaluate<T>(setup, s, "cuda");
  const auto mixed = evaluate<T>(setup, s, "cuda-mixed");
  for (size_t i = 0; i < expected.size(); ++i) {
    CHECK(maxRelativeError(results[i], expected[i]) <= TOLERANCE);
    const bool peaks = !i && std::is_same_v<T, PDFCalculator>;
    const double tolerance = peaks ? MIXED_TOLERANCE : MIXED_REST_TOLERANCE;
    CHECK(maxRelativeError(mixed[i], expected[i]) <= tolerance);
  }
}

}  // namespace

TEST_CASE("TestCalculationKernelCuda::test_results_of_calculators", "[cuda]") {
  skipWithoutDevice();
  forEachCalculatorCase([](const Setup& setup, StructureAdapterPtr structure) {
    checkAsCPU<PDFCalculator>(setup, structure);
    checkAsCPU<DebyePDFCalculator>(setup, structure);
  });
}

#ifdef LIBDIFFPY_TEST_CUDA_RUNTIME
TEST_CASE("TestCalculationKernelCuda::test_device_switch", "[cuda]") {
  skipWithoutDevice();
  int count, original;
  REQUIRE(cudaGetDeviceCount(&count) == cudaSuccess);
  if (count < 2) SKIP("Requires two CUDA devices");
  REQUIRE(cudaGetDevice(&original) == cudaSuccess);
  struct RestoreDevice {
    int device;
    ~RestoreDevice() { cudaSetDevice(device); }
  } restore{original};

  const auto stru = loadPeriodic("CaTiO3.stru");
  PDFCalculator pdf;
  DebyePDFCalculator debye;
  PairQuantity* calculators[] = {&pdf, &debye};
  for (PairQuantity* calc : calculators) {
    setupDefault(*calc);
    calc->setEvaluatorType(BASIC);
    const QuantityType expected = evalWith("scalar", *calc, stru);
    for (const int device : {original, (original + 1) % count, original}) {
      REQUIRE(cudaSetDevice(device) == cudaSuccess);
      skipWithoutDevice();
      CHECK(maxRelativeError(evalWith("cuda", *calc, stru), expected) <=
            TOLERANCE);
    }
  }
}
#endif

TEST_CASE("TestCalculationKernelCuda::test_peak_profiles_and_cutoffs",
          "[cuda]") {
  skipWithoutDevice();
  const StructureAdapterPtr stru = loadPeriodic("CaTiO3.stru", true);
  const Setup cropped = [](PairQuantity& c) {
    setupDefault(c);
    dynamic_cast<PDFCalculator&>(c).setPeakProfileByType("croppedgaussian");
  };
  const Setup precise = [](PairQuantity& c) {
    setupDefault(c);
    c.setDoubleAttr("peakprecision", 1e-12);
  };
  // peaks that are cut at one percent of their height
  const Setup cut = [](PairQuantity& c) {
    setupDefault(c);
    dynamic_cast<PDFCalculator&>(c).setPeakProfileByType("croppedgaussian");
    c.setDoubleAttr("peakprecision", 0.01);
  };
  checkAsCPU<PDFCalculator>(cropped, stru);
  checkAsCPU<PDFCalculator>(cut, stru);
  checkAsCPU<PDFCalculator>(precise, stru);
  // Debye sums that stop early, late and never
  for (double prec : {1e-3, 1e-6, 1e-12, 0.0}) {
    INFO("debyeprecision " << prec);
    const Setup debye = [=](PairQuantity& c) {
      setupDefault(c);
      c.setDoubleAttr("debyeprecision", prec);
    };
    checkAsCPU<DebyePDFCalculator>(debye, stru);
    checkAsCPU<DebyePDFCalculator>(debye,
                                   asCluster(*loadPeriodic("Ni.stru", true)));
  }
}

TEST_CASE("TestCalculationKernelCuda::test_narrow_jeong_peak", "[cuda]") {
  skipWithoutDevice();
  auto stru = std::make_shared<AtomicStructureAdapter>();
  Atom atom;
  atom.atomtype = "Ni";
  atom.uij_cartn = 0.05 * R3::identity();
  stru->append(atom);
  atom.xyz_cartn[0] = 1.0;
  stru->append(atom);
  checkAsCPU<PDFCalculator>(
    [](PairQuantity& c) {
      c.setRmax(2.0);
      c.setDoubleAttr("rstep", 0.01);
      c.setDoubleAttr("delta2", 1.0 - 2.5e-8);
    },
    stru);
}

TEST_CASE("TestCalculationKernelCuda::test_sampled_peak_windows", "[cuda]") {
  skipWithoutDevice();
  auto stru = std::make_shared<AtomicStructureAdapter>();
  Atom atom;
  atom.atomtype = "C";
  stru->append(atom);
  stru->append(atom);
  PDFCalculator calc;
  calc.setRmax(4);
  calc.setRstep(0.01);
  calc.setMaxExtension(0);
  calc.setEvaluatorType(BASIC);
  SECTION("an unresolved peak between samples leaves the other peaks intact") {
    (*stru)[1].xyz_cartn[0] = 1.005;
    atom.xyz_cartn[1] = 2;
    atom.uij_cartn = 0.01 * R3::identity();
    stru->append(atom);
    calc.setDoubleAttr("qbroad_seperable", 1e-10);
  }
  SECTION("a 3 angstrom bond at the inclusive upper distance cutoff") {
    (*stru)[1].xyz_cartn = R3::Vector(0.4, 2.0, 2.2);
    (*stru)[1].uij_cartn = 0.02 * R3::identity();
    calc.setRmax(3);
  }
  SECTION("a sample at the cropped Gaussian boundary") {
    (*stru)[1].xyz_cartn[0] = 3;
    calc.setPeakWidthModelByType("constant");
    calc.setDoubleAttr("width", 0.2);
    calc.setPeakProfileByType("croppedgaussian");
    calc.setDoubleAttr("peakprecision", 0.5);
  }
  const QuantityType reference = evalWith("scalar", calc, stru);
  REQUIRE(*std::ranges::max_element(reference) > 0.0);
  CHECK(maxRelativeError(evalWith("cuda", calc, stru), reference) <= TOLERANCE);
  CHECK(maxRelativeError(evalWith("cuda-mixed", calc, stru), reference) <=
        MIXED_TOLERANCE);
}

TEST_CASE("TestCalculationKernelCuda::test_weak_peaks", "[cuda]") {
  skipWithoutDevice();
  auto stru = std::make_shared<AtomicStructureAdapter>();
  Atom atom;
  atom.atomtype = "Ni";
  atom.occupancy = 0.01;
  for (int x = 0; x < 8; ++x)
    for (int y = 0; y < 8; ++y)
      for (int z = 0; z < 8; ++z) {
        atom.xyz_cartn = R3::Vector(x, y, z);
        stru->append(atom);
      }
  atom.occupancy = 1;
  atom.xyz_cartn = R3::Vector(0, 100, 0);
  stru->append(atom);
  atom.xyz_cartn[0] = 2.7;
  stru->append(atom);
  PDFCalculator calc;
  calc.setEvaluatorType(BASIC);
  calc.setRmax(5);
  calc.setMaxExtension(0);
  calc.setPeakWidthModelByType("constant");
  calc.setDoubleAttr("width", 0.1);
  const QuantityType reference = evalWith("scalar", calc, stru);
  const QuantityType actual = evalWith("cuda", calc, stru);
  CHECK(maxRelativeError(actual, reference) <= TOLERANCE);
  CHECK(actual == evalWith("cuda", calc, stru));
}

TEST_CASE("TestCalculationKernelCuda::test_split_site_debye_sum", "[cuda]") {
  skipWithoutDevice();
  auto stru = std::make_shared<AtomicStructureAdapter>();
  Atom atom;
  atom.atomtype = "Ni";
  atom.uij_cartn = 0.01 * R3::identity();
  stru->append(atom);
  atom.xyz_cartn[0] = 2.0;
  stru->append(atom);
  DebyePDFCalculator calc;
  calc.setRmax(4);
  calc.setQstep(0.01);
  calc.setEvaluatorType(BASIC);
  const QuantityType reference = evalWith("scalar", calc, stru);
  // Split the first site into coincident partial occupancies. All its
  // terms for a given Q have the same sign and their sum is unchanged.
  (*stru)[0].occupancy = 1.0 / 2048;
  atom = (*stru)[0];
  for (int i = 1; i < 2048; ++i) stru->append(atom);
  // Put the unsplit site last so that one subchunk sums 1024 equal terms.
  std::swap((*stru)[1], stru->at(2048));
  CHECK(maxRelativeError(evalWith("cuda-mixed", calc, stru), reference) <=
        MIXED_TOLERANCE);
}

TEST_CASE("TestCalculationKernelCuda::test_many_occupancies_on_a_fine_q_grid",
          "[cuda]") {
  skipWithoutDevice();
  auto stru = std::make_shared<AtomicStructureAdapter>();
  Atom atom;
  atom.atomtype = "Ni";
  // Only the first two sites form a bond. All 129 occupancies enter the
  // factor table, which spans more than one Q tile at this Q step.
  for (int i = 0; i < 129; ++i) {
    atom.occupancy = (i + 1) / 129.0;
    atom.xyz_cartn[0] = (i == 1) ? 2.0 : 10.0 * i;
    stru->append(atom);
  }
  DebyePDFCalculator calc;
  calc.setEvaluatorType(BASIC);
  calc.setRmax(4);
  calc.setMaxExtension(0);
  calc.setQmax(9);
  calc.setQstep(0.001);
  const QuantityType expected = evalWith("scalar", calc, stru);
  REQUIRE(*std::ranges::max_element(expected) > 0);
  CHECK(maxRelativeError(evalWith("cuda", calc, stru), expected) <= TOLERANCE);
  CHECK(maxRelativeError(evalWith("cuda-mixed", calc, stru), expected) <=
        MIXED_TOLERANCE);
}

TEST_CASE("TestCalculationKernelCuda::test_change_precision", "[cuda]") {
  skipWithoutDevice();
  const auto stru = loadPeriodic("CaTiO3.stru", true);
  PDFCalculator pdf;
  DebyePDFCalculator debye;
  PairQuantity* calculators[] = {&pdf, &debye};
  for (PairQuantity* calc : calculators) {
    setupDefault(*calc);
    calc->setEvaluatorType(OPTIMIZED);
    const QuantityType mixed = evalWith("cuda-mixed", *calc, stru);
    const QuantityType doubles = evalWith("cuda", *calc, stru);
    CHECK(calc->getEvaluatorTypeUsed() == BASIC);
    calc->setEvaluatorType(BASIC);
    const QuantityType expected = evalWith("cuda", *calc, stru);
    REQUIRE(mixed != expected);
    CHECK(maxRelativeError(doubles, expected) == 0.0);
  }
}

TEST_CASE("TestCalculationKernelCuda::test_many_pairs", "[cuda]") {
  skipWithoutDevice();
  // more pairs than the device takes at once
  const StructureAdapterPtr sphere = makeSphere(19.5);
  const long long n = sphere->countSites();
  INFO(n << " atoms");
  REQUIRE(n * (n - 1) / 2 > (3LL << 20));
  const Setup setup = [](PairQuantity& c) {
    c.setRmax(41.0);
    c.setDoubleAttr("qmax", 20.0);
  };
  checkAsCPU<PDFCalculator>(setup, sphere);
  checkAsCPU<DebyePDFCalculator>(setup, sphere);
  // and of a periodic structure
  const StructureAdapterPtr periodic =
    loadPeriodic("alpha_K2Bi8Se13.stru", true);
  const Setup far = [](PairQuantity& c) {
    c.setRmax(60.0);
    c.setDoubleAttr("qmax", 15.0);
  };
  checkAsCPU<PDFCalculator>(far, periodic);
}

TEST_CASE("TestCalculationKernelCuda::test_evaluations_in_turn", "[cuda]") {
  skipWithoutDevice();
  // The list of pairs stays on the device while the structure changes
  // little.
  const auto periodic = loadPeriodic("CaTiO3.stru", true);
  for (std::shared_ptr<PeriodicStructureAdapter> stru :
       {periodic,
        std::shared_ptr<PeriodicStructureAdapter>(asCrystal(*periodic))}) {
    AllPairs<PDFCalculator> cpu, device;
    AllPairs<DebyePDFCalculator> dcpu, ddevice;
    for (PairQuantity* c : {(PairQuantity*)&cpu, (PairQuantity*)&device,
                            (PairQuantity*)&dcpu, (PairQuantity*)&ddevice}) {
      setupDefault(*c);
    }
    const Lattice L = stru->getLattice();
    for (int k = 0; k < 12; ++k) {
      INFO("evaluation " << k);
      const int n = stru->countSites();
      // every fourth change is beyond the margin of the pairs
      const double step = (k % 4 == 3) ? 0.3 : 0.004;
      stru->at(k % n).xyz_cartn[k % 3] += step * ((k % 5) - 2);
      stru->at((3 * k) % n).uij_cartn *= 1.0 + 0.01 * ((k % 3) - 1);
      for (Atom& a : *stru) stru->toFractional(a);
      const double f = 1.0 + 1e-4 * ((k % 7) - 3);
      stru->setLatPar(f * L.a(), L.b(), L.c() / f, L.alpha(), L.beta(),
                      L.gamma());
      for (Atom& a : *stru) stru->toCartesian(a);
      StructureAdapterPtr copy = stru->clone();
      // every third sum of the device is in mixed precision
      const bool mixed = (k % 3 == 2);
      const std::string kernel = mixed ? "cuda-mixed" : "cuda";
      const double tolerance = mixed ? MIXED_REST_TOLERANCE : TOLERANCE;
      CHECK(maxRelativeError(evalWith(kernel, device, copy), cpu.eval(copy)) <=
            tolerance);
      CHECK(maxRelativeError(evalWith(kernel, ddevice, copy),
                             dcpu.eval(copy)) <= tolerance);
      CHECK(device.addedOnDevice());
      CHECK(ddevice.addedOnDevice());
      CHECK(!cpu.addedOnDevice());
      // The device adds all pairs again rather than those of a few sites.
      if (k > 0) CHECK(BASIC == device.getEvaluatorTypeUsed());
    }
    // a structure of the same content takes no sum
    evalWith("cuda-mixed", device, stru->clone());
    CHECK(OPTIMIZED == device.getEvaluatorTypeUsed());
  }
}

TEST_CASE("TestCalculationKernelCuda::test_sums_on_the_cpu", "[cuda]") {
  skipWithoutDevice();
  struct DerivedProfile : GaussianProfile {};
  const StructureAdapterPtr stru = loadPeriodic("CaTiO3.stru", true);
  // What the device does not take gives the results of the CPU to the
  // last bit.
  auto check = [&](PairQuantity& calc, PairQuantity& cpu) {
    for (PairQuantity* c : {&calc, &cpu}) {
      setupDefault(*c);
      c->setEvaluatorType(BASIC);
    }
    for (const char* kernel : {"cuda", "cuda-mixed"}) {
      INFO("kernel " << kernel);
      CHECK(evalWith(kernel, calc, stru) == cpu.eval(stru));
    }
  };
  // all pairs at once, with the kernel of the CPU
  AllPairs<PDFCalculator> masked, maskedcpu;
  masked.setPairMask(0, 3, false);
  maskedcpu.setPairMask(0, 3, false);
  check(masked, maskedcpu);
  CHECK(masked.addedAllPairs());
  CHECK(!masked.addedOnDevice());
  AllPairs<DebyePDFCalculator> dmasked, dmaskedcpu;
  dmasked.setTypeMask("Ti4+", "O2-", false);
  dmaskedcpu.setTypeMask("Ti4+", "O2-", false);
  check(dmasked, dmaskedcpu);
  CHECK(dmasked.addedAllPairs());
  CHECK(!dmasked.addedOnDevice());
  AllPairs<PDFCalculator> profiled, profiledcpu;
  profiled.setPeakProfile(PeakProfilePtr(new DerivedProfile));
  profiledcpu.setPeakProfile(PeakProfilePtr(new DerivedProfile));
  check(profiled, profiledcpu);
  CHECK_FALSE(profiled.addedAllPairs());
  CHECK(!profiled.addedOnDevice());
  // the pairs one by one
  PairByPair<PDFCalculator> derived, derivedcpu;
  check(derived, derivedcpu);
  BVSCalculator bvsc, bvscpu;
  CHECK(evalWith("cuda", bvsc, stru) == bvscpu.eval(stru));
  // all pairs of a calculator without any of these go to the device
  AllPairs<PDFCalculator> pdfc;
  PDFCalculator pdfcpu;
  for (PairQuantity* c : {(PairQuantity*)&pdfc, (PairQuantity*)&pdfcpu}) {
    setupDefault(*c);
    c->setEvaluatorType(BASIC);
  }
  CHECK(maxRelativeError(evalWith("cuda", pdfc, stru), pdfcpu.eval(stru)) <=
        TOLERANCE);
  CHECK(pdfc.addedOnDevice());
  // and back to the CPU with its kernels
  CHECK(pdfc.eval(stru) == pdfcpu.value());
  CHECK(pdfc.addedAllPairs());
  CHECK(!pdfc.addedOnDevice());
  // then to the device again
  CHECK(maxRelativeError(evalWith("cuda-mixed", pdfc, stru), pdfcpu.value()) <=
        MIXED_TOLERANCE);
  CHECK(pdfc.addedOnDevice());
}

#ifdef LIBDIFFPY_TEST_KERNEL_INTERNALS

TEST_CASE("TestCalculationKernelCuda::test_kernels_of_the_cpu", "[cuda]") {
  skipWithoutDevice();
  // The CUDA kernels are those of the CPU for everything else.
  KernelGuard guard;
  setCalculationKernel(cpuKernel());
  const kernels::CalculationKernel cpu = kernels::activeCalculationKernel();
  CHECK_FALSE(cpu.onDevice);
  CHECK(cpu.pdfPairSum);
  CHECK(cpu.debyePairSum);
  for (const char* name : {"cuda", "cuda-mixed"}) {
    INFO("kernel " << name);
    setCalculationKernel(name);
    const kernels::CalculationKernel& k = kernels::activeCalculationKernel();
    CHECK(std::string(name) == k.name);
    CHECK(k.gaussianWindow == cpu.gaussianWindow);
    CHECK(k.debyeSum == cpu.debyeSum);
    CHECK(k.bondCandidates == cpu.bondCandidates);
    CHECK(k.pdfPairSum);
    CHECK(k.debyePairSum);
    CHECK(k.makeWorkspace);
    CHECK(k.onDevice);
  }
}

#endif  // LIBDIFFPY_TEST_KERNEL_INTERNALS
