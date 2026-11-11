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
 * Tests of the calculation kernels against the scalar kernel
 *
 *****************************************************************************/

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include <catch2/catch_test_macros.hpp>

#include <diffpy/srreal/AtomicStructureAdapter.hpp>
#include <diffpy/srreal/CalculationKernel.hpp>
#include <diffpy/srreal/CroppedGaussianProfile.hpp>
#include <diffpy/srreal/GaussianProfile.hpp>
#include <diffpy/srreal/DebyePDFCalculator.hpp>
#include <diffpy/srreal/PDFCalculator.hpp>

#include "test_helpers.hpp"

// The internal kernel header is missing when testing an installed library.
#if __has_include(<diffpy/srreal/kernels/CalculationKernelImpl.hpp>)
#include <diffpy/srreal/kernels/CalculationKernelImpl.hpp>
#define LIBDIFFPY_TEST_KERNEL_INTERNALS 1
#endif

using namespace diffpy::srreal;

namespace {

const double KERNEL_TOLERANCE = 1e-12;

/// One calculation that is repeated with every kernel.
struct KernelCase {
  std::string name;
  std::string calculator;  // "pdf" or "debye"
  std::string structure;   // .stru file in testdata
  std::string profile;     // peak profile, PDFCalculator only
  std::vector<std::pair<std::string, double>> attributes;
};

const std::vector<KernelCase> KERNEL_CASES = {
  {"pdf_LiTaO3", "pdf", "LiTaO3.stru", "gaussian", {{"rmax", 15.0}}},
  {"pdf_LiTaO3_qmax",
   "pdf",
   "LiTaO3.stru",
   "gaussian",
   {{"rmin", 1.0}, {"qmax", 25.0}, {"delta2", 1.0}}},
  {"pdf_K2Bi8Se13_cropped",
   "pdf",
   "alpha_K2Bi8Se13.stru",
   "croppedgaussian",
   {{"rmax", 12.0}}},
  {"debye_LiTaO3",
   "debye",
   "LiTaO3.stru",
   "",
   {{"qmin", 1.0}, {"qmax", 20.0}, {"debyeprecision", 1e-8}}},
};

template <class Calc>
void setAttributes(Calc& calc, const KernelCase& kc) {
  for (const auto& [name, value] : kc.attributes) {
    calc.setDoubleAttr(name, value);
  }
}

/// PDF and the second result array (RDF or F) for the active kernel.
std::vector<QuantityType> computeCase(const KernelCase& kc) {
  StructureAdapterPtr stru = loadTestPeriodicStructure(kc.structure);
  if (kc.calculator == "pdf") {
    PDFCalculator calc;
    calc.setPeakProfileByType(kc.profile);
    setAttributes(calc, kc);
    calc.eval(stru);
    return {calc.getPDF(), calc.getRDF()};
  }
  DebyePDFCalculator calc;
  setAttributes(calc, kc);
  calc.eval(stru);
  return {calc.getPDF(), calc.getF()};
}

}  // namespace

TEST_CASE("TestCalculationKernel::test_nonfinite_error",
          "[CalculationKernel]") {
  const double inf = std::numeric_limits<double>::infinity();
  const double nan = std::numeric_limits<double>::quiet_NaN();
  const QuantityType reference(2, 1.0);
  for (double invalid : {nan, inf, -inf}) {
    QuantityType actual(reference);
    actual[0] = invalid;
    CHECK(maxRelativeError(actual, reference) == inf);
    CHECK(maxRelativeError(reference, actual) == inf);
    CHECK(maxRelativeError(actual, actual) == inf);
  }
}

TEST_CASE("TestCalculationKernel::test_narrow_peaks", "[CalculationKernel]") {
  KernelGuard guard;
  for (double distance : {2.123, 100.123, 1000.123}) {
    auto stru = std::make_shared<AtomicStructureAdapter>();
    Atom atom;
    atom.atomtype = "C";
    stru->append(atom);
    atom.xyz_cartn[0] = distance;
    stru->append(atom);
    for (const char* profile : {"gaussian", "croppedgaussian"}) {
      for (double width : {0.001, 0.01, 0.1}) {
        PDFCalculator calc;
        calc.setEvaluatorType(BASIC);
        calc.setRmin(distance - 0.25);
        calc.setRmax(distance + 0.25);
        calc.setRstep(0.001);
        calc.setPeakProfileByType(profile);
        calc.setDoubleAttr("peakprecision", 1e-16);
        calc.setPeakWidthModelByType("constant");
        calc.setDoubleAttr("width", width);
        setCalculationKernel("scalar");
        calc.eval(stru);
        const QuantityType reference = calc.getPDF();
        INFO("profile " << profile << ", distance " << distance << ", width "
                        << width);
        testAllKernels([&] {
          calc.eval(stru);
          CHECK(maxRelativeError(calc.getPDF(), reference) <= KERNEL_TOLERANCE);
        });
      }
    }
  }
}

TEST_CASE("TestCalculationKernel::test_wide_windows", "[CalculationKernel]") {
  // addRDFPeak over windows that reach far beyond the peak, where the
  // Gaussian underflows to zero, must match the scalar kernel.
  KernelGuard guard;
  const int npts = 2000;
  const double rstep = 0.01;
  const double dist = 10.0;
  const double fwhm = 0.05;
  GaussianProfile gaussian;
  CroppedGaussianProfile cropped;
  const std::pair<int, int> windows[] = {{0, npts},   {900, npts},  {0, 1010},
                                         {995, 1005}, {1500, npts}, {0, 600}};
  const PeakProfile* profiles[] = {&gaussian, &cropped};
  for (const PeakProfile* pk : profiles) {
    for (const auto& [ilo, ihi] : windows) {
      setCalculationKernel("scalar");
      QuantityType reference(npts, 0.0);
      pk->addRDFPeak(reference.data(), ilo, ihi, 0, rstep, dist, fwhm, 1.0);
      if (ilo <= 1000 && 1000 < ihi) {
        // x is zero at the point 1000, where y equals the peak height
        CHECK(reference[1000] == (*pk)(0.0, fwhm));
      }
      INFO("profile " << pk->type() << ", window [" << ilo << ", " << ihi
                      << ")");
      testAllKernels([&] {
        QuantityType actual(npts, 0.0);
        pk->addRDFPeak(actual.data(), ilo, ihi, 0, rstep, dist, fwhm, 1.0);
        // relative to the peak height, the reference may be all zero
        CHECK(maxRelativeError(actual, reference, (*pk)(0.0, fwhm)) <=
              KERNEL_TOLERANCE);
      });
    }
  }
}

namespace {

/// Gaussian profile with a changed operator(), which must not use the kernel.
class DoubledGaussianProfile : public GaussianProfile {
 public:
  double operator()(double x, double fwhm) const override {
    return 2 * this->GaussianProfile::operator()(x, fwhm);
  }
};

/// Cropped profile with a changed operator(), which must not use the kernel.
class DoubledCroppedGaussianProfile : public CroppedGaussianProfile {
 public:
  double operator()(double x, double fwhm) const override {
    return 2 * this->CroppedGaussianProfile::operator()(x, fwhm);
  }
};

/// Gaussian profile that keeps operator() and asks for the kernel.
class KernelGaussianProfile : public GaussianProfile {
 protected:
  KernelShape kernelShape() const override { return KernelShape::Gaussian; }
};

QuantityType rdfPeak(const PeakProfile& pk) {
  QuantityType rv(300, 0.0);
  pk.addRDFPeak(rv.data(), 50, 250, 900, 0.01, 10.5, 0.37, 1.5);
  return rv;
}

}  // namespace

TEST_CASE("TestCalculationKernel::test_derived_profiles",
          "[CalculationKernel]") {
  testAllKernels([] {
    const QuantityType gaussian = rdfPeak(GaussianProfile());
    const QuantityType cropped = rdfPeak(CroppedGaussianProfile());
    QuantityType doubled(gaussian);
    for (double& y : doubled) y *= 2;
    QuantityType doubledcropped(cropped);
    for (double& y : doubledcropped) y *= 2;
    CHECK(maxRelativeError(rdfPeak(DoubledGaussianProfile()), doubled) <=
          KERNEL_TOLERANCE);
    CHECK(maxRelativeError(rdfPeak(DoubledCroppedGaussianProfile()),
                           doubledcropped) <= KERNEL_TOLERANCE);
    CHECK(rdfPeak(KernelGaussianProfile()) == gaussian);
  });
}

TEST_CASE("TestCalculationKernel::test_available", "[CalculationKernel]") {
  const std::vector<std::string> names = availableCalculationKernels();
  const auto has = [&](const char* name) {
    return std::ranges::find(names, name) != names.end();
  };
  REQUIRE(!names.empty());
  // The scalar kernel is the last resort.
  CHECK(names.back() == "scalar");
  const std::vector<std::string> known = {"avx512", "avx2", "sse2", "neon",
                                          "scalar"};
  for (const std::string& kname : names) {
    INFO("kernel " << kname);
    CHECK(std::ranges::count(names, kname) == 1);
    CHECK(std::ranges::find(known, kname) != known.end());
  }
  // The x86-64 kernels are built together and wider ones need newer CPUs.
  if (has("avx512")) {
    CHECK(has("avx2"));
  }
  if (has("avx2")) {
    CHECK(has("sse2"));
  }
  CHECK(!(has("sse2") && has("neon")));
}

TEST_CASE("TestCalculationKernel::test_selection", "[CalculationKernel]") {
  KernelGuard guard;
  const std::vector<std::string> names = availableCalculationKernels();
  CHECK(std::ranges::find(names, getCalculationKernel()) != names.end());
  for (const std::string& kname : names) {
    setCalculationKernel(kname);
    CHECK(kname == getCalculationKernel());
  }
  // Names are lowercase and kernels of other CPUs are not available.
  const std::string active = getCalculationKernel();
  for (const char* invalid : {"", "sse", "avx", "AVX2", "Scalar", "sve"}) {
    INFO("kernel '" << invalid << "'");
    CHECK_THROWS_AS(setCalculationKernel(invalid), std::invalid_argument);
    CHECK(active == getCalculationKernel());
  }
  if (std::ranges::find(names, "sse2") != names.end()) {
    CHECK_THROWS_AS(setCalculationKernel("neon"), std::invalid_argument);
  } else {
    CHECK_THROWS_AS(setCalculationKernel("sse2"), std::invalid_argument);
  }
}

TEST_CASE("TestCalculationKernel::test_match_scalar", "[CalculationKernel]") {
  KernelGuard guard;
  for (const KernelCase& kc : KERNEL_CASES) {
    setCalculationKernel("scalar");
    const std::vector<QuantityType> reference = computeCase(kc);
    INFO("case " << kc.name);
    testAllKernels([&] {
      const std::vector<QuantityType> actual = computeCase(kc);
      REQUIRE(actual.size() == reference.size());
      for (size_t i = 0; i < actual.size(); ++i) {
        INFO("array " << i);
        CHECK(maxRelativeError(actual[i], reference[i]) <= KERNEL_TOLERANCE);
      }
    });
  }
}

#ifdef LIBDIFFPY_TEST_KERNEL_INTERNALS
namespace {

/// True when both arrays keep the fill value at the same indices.
bool sameUntouched(const QuantityType& actual, const QuantityType& reference,
                   double fill) {
  if (actual.size() != reference.size()) return false;
  for (size_t i = 0; i < actual.size(); ++i) {
    if ((actual[i] == fill) != (reference[i] == fill)) return false;
  }
  return true;
}

}  // namespace

TEST_CASE("TestCalculationKernel::test_window_lengths", "[CalculationKernel]") {
  using namespace diffpy::srreal::kernels;
  KernelGuard guard;
  setCalculationKernel("scalar");
  const CalculationKernel& scalar = activeCalculationKernel();
  // around the vector and the reseed block sizes of all kernels,
  // including the longest last block of 80, 160 or 320 points
  const std::vector<int> lengths = {0,   1,   2,   3,   4,   5,   7,   8,   9,
                                    15,  16,  17,  63,  64,  65,  79,  80,  81,
                                    127, 128, 129, 159, 160, 161, 255, 256, 257,
                                    300, 319, 320, 321, 511, 512, 513, 600};
  const double fill = 0.25;
  const int padding = 20;
  testAllKernels([&] {
    const CalculationKernel& kernel = activeCalculationKernel();
    for (int n : lengths) {
      for (int lo : {0, 1, 3}) {
        const int hi = lo + n;
        // one Gaussian peak across the window
        for (bool cropped : {false, true}) {
          INFO("lo " << lo << ", n " << n << ", cropped " << cropped);
          QuantityType reference(hi + padding, fill);
          QuantityType actual(reference);
          const int rlosteps = 100;
          const double rstep = 0.01;
          GaussianWindowArgs args = {
            .out = reference.data(),
            .ilo = lo,
            .ihi = hi,
            .rlosteps = rlosteps,
            .rstep = rstep,
            .dist = (rlosteps + lo + 0.5 * n + 0.37) * rstep,
            .fwhm = std::max(n, 8) * rstep / 4,
            .scale = -1.5,
            .peakheight = 1.0,
            .cropped = cropped,
            .halfboundrel = 1.5,
          };
          scalar.gaussianWindow(args);
          args.out = actual.data();
          kernel.gaussianWindow(args);
          CHECK(maxRelativeError(actual, reference) <= KERNEL_TOLERANCE);
          CHECK(sameUntouched(actual, reference, fill));
        }
        // Debye sum over the whole range and stopped at kq = 298 or 115.
        // The scale is negative when a calculator removes a pair.
        std::vector<double> sf(hi);
        for (int kq = 0; kq < hi; ++kq) sf[kq] = std::exp(-0.01 * kq);
        for (double scale : {2.0, -2.0}) {
          for (double prec : {0.0, 1e-3, 0.1}) {
            INFO("lo " << lo << ", n " << n << ", scale " << scale << ", prec "
                       << prec);
            QuantityType reference(hi + padding, fill);
            QuantityType actual(reference);
            DebyeSumArgs args = {
              .out = reference.data(),
              .kqlo = lo,
              .kqhi = hi,
              .qstep = 0.05,
              .dist = 1.7,
              .dwsigma = 0.1,
              .scale = scale,
              .prec = prec,
              .sf0 = sf.data(),
              .sf1 = sf.data(),
            };
            scalar.debyeSum(args);
            args.out = actual.data();
            kernel.debyeSum(args);
            CHECK(maxRelativeError(actual, reference) <= KERNEL_TOLERANCE);
            CHECK(sameUntouched(actual, reference, fill));
          }
        }
      }
    }
  });
}

TEST_CASE("TestCalculationKernel::test_bond_candidates",
          "[CalculationKernel]") {
  const double nan = std::numeric_limits<double>::quiet_NaN();
  const double inf = std::numeric_limits<double>::infinity();
  const double d2lo = 0.25;
  const double d2hi = 4.0;
  // squares of 0.5 and 2.0 are the exact bounds
  const std::vector<double> pattern = {0.5,
                                       std::nextafter(0.5, 0.0),
                                       std::nextafter(0.5, 1.0),
                                       2.0,
                                       std::nextafter(2.0, 1.0),
                                       std::nextafter(2.0, 3.0),
                                       nan,
                                       -1.0,
                                       0.0,
                                       3.0,
                                       inf,
                                       -0.5,
                                       nan};
  const std::vector<int> lengths = {0, 1, 2,  3,  4,  5,  7,
                                    8, 9, 15, 16, 17, 31, 100};
  testAllKernels([&] {
    const diffpy::srreal::kernels::CalculationKernel& kernel =
      diffpy::srreal::kernels::activeCalculationKernel();
    for (int n : lengths) {
      for (int shift = 0; shift < 3; ++shift) {
        for (int axis = 0; axis < 3; ++axis) {
          INFO("n = " << n << ", shift " << shift << ", axis " << axis);
          std::vector<double> xyz[3];
          for (std::vector<double>& c : xyz) c.assign(n, 0.0);
          std::vector<int> expected;
          for (int i = 0; i < n; ++i) {
            const double v = pattern[(i + shift) % pattern.size()];
            xyz[axis][i] = v;
            const double d2 = v * v;
            if (std::isnan(d2) || (d2lo <= d2 && d2 <= d2hi)) {
              expected.push_back(i);
            }
          }
          std::vector<int> idx(n + 1, -1);
          diffpy::srreal::kernels::BondCandidatesArgs args;
          args.x = xyz[0].data();
          args.y = xyz[1].data();
          args.z = xyz[2].data();
          args.n = n;
          args.ox = 0.0;
          args.oy = 0.0;
          args.oz = 0.0;
          args.d2lo = d2lo;
          args.d2hi = d2hi;
          args.idx = idx.data();
          const int m = kernel.bondCandidates(args);
          REQUIRE(m == int(expected.size()));
          CHECK(std::vector<int>(idx.begin(), idx.begin() + m) == expected);
          CHECK(-1 == idx[n]);
        }
      }
    }
    // All points are the same, so all or none are candidates.  With the
    // bound between the fused and the separately rounded x^2 + z^2, the
    // count is partial when the last points round unlike the vector lanes.
    int nrounded = 0;
    for (int k = 1; k < 200; ++k) {
      const double vx = 1.0 + k * 0x1p-30;
      const double vz = 3.0 - k * 0x1p-31;
      const volatile double xx = vx * vx;
      const volatile double zz = vz * vz;
      const double fused = std::fma(vz, vz, xx);
      const double separate = zz + xx;
      if (fused == separate) continue;
      ++nrounded;
      for (int n : {9, 15, 17}) {
        INFO("k = " << k << ", n = " << n);
        std::vector<double> x(n, vx), y(n, 0.0), z(n, vz);
        std::vector<int> idx(n);
        const diffpy::srreal::kernels::BondCandidatesArgs args = {
          .x = x.data(),
          .y = y.data(),
          .z = z.data(),
          .n = n,
          .d2hi = std::min(fused, separate),
          .idx = idx.data(),
        };
        const int m = kernel.bondCandidates(args);
        CHECK((m == 0 || m == n));
      }
    }
    CHECK(nrounded > 0);
    std::vector<double> x(9, 1.0), y(9, -2.0), z(9, 2.0);
    std::vector<int> idx(9);
    diffpy::srreal::kernels::BondCandidatesArgs args = {
      x.data(), y.data(), z.data(), 9, -1.0, 2.0, 0.0, d2lo, d2hi, idx.data()};
    CHECK(9 == kernel.bondCandidates(args));
    // z + oz is exactly the next double above 2
    args.oz = std::nextafter(2.0, 3.0) - 2.0;
    CHECK(0 == kernel.bondCandidates(args));
  });
}
#endif  // LIBDIFFPY_TEST_KERNEL_INTERNALS
