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
 * TestPDFUtils -- unit test suite for the Fourier transformations
 *     fftgtof and fftftog
 *
 *****************************************************************************/

#include <cmath>
#include <array>
#include <barrier>
#include <future>

#include <catch2/catch_test_macros.hpp>

#include <diffpy/srreal/PDFUtils.hpp>

#include "test_helpers.hpp"

using namespace std;
using namespace diffpy::srreal;

namespace {

/// Tolerance relative to the largest value, round-off is below 1e-13 here.
const double RELTOL = 1e-12;

/// oscillating signal without any symmetry
QuantityType makeSignal(int npts) {
  QuantityType g(npts);
  for (int i = 0; i < npts; ++i) {
    g[i] = sin(0.37 * i + 0.2) * exp(-0.002 * i) + 0.3 * cos(1.9 * i);
  }
  return g;
}

/// direct summation of the sine transformation evaluated by fftgtof
QuantityType sumgtof(const QuantityType& g, double rstep, int padrmin = 0) {
  int Npad2 = 1;
  while (Npad2 < padrmin + int(g.size())) Npad2 *= 2;
  QuantityType f(Npad2);
  for (int j = 0; j < Npad2; ++j) {
    double fj = 0.0;
    for (int i = 0; i < int(g.size()); ++i) {
      // reduce the sine argument for an accurate result at large j * k
      const int jk = j * (padrmin + i) % (2 * Npad2);
      fj += g[i] * sin(M_PI * jk / Npad2);
    }
    f[j] = fj * rstep;
  }
  return f;
}

}  // namespace

TEST_CASE("TestPDFUtils::test_fftgtof_empty", "[PDFUtils]") {
  QuantityType g;
  CHECK(fftgtof(g, 0.01).empty());
  CHECK(fftftog(g, 0.01).empty());
}

TEST_CASE("TestPDFUtils::test_fftgtof_small", "[PDFUtils]") {
  const double rstep = 0.25;
  QuantityType g(1, 3.0);
  QuantityType f = fftgtof(g, rstep);
  REQUIRE(1u == f.size());
  CHECK(0.0 == f[0]);
  g.push_back(5.0);
  f = fftgtof(g, rstep);
  REQUIRE(2u == f.size());
  CHECK(0.0 == f[0]);
  CHECK(5.0 * rstep == f[1]);
  g.push_back(-2.0);
  f = fftgtof(g, rstep);
  REQUIRE(4u == f.size());
  CHECK(maxRelativeError(f, sumgtof(g, rstep)) < RELTOL);
}

TEST_CASE("TestPDFUtils::test_fftgtof_sizes", "[PDFUtils]") {
  const double rstep = 0.01;
  // one length for each padded size from 4 to 1024
  const int sizes[] = {4, 7, 9, 32, 33, 100, 256, 500, 999};
  for (int npts : sizes) {
    INFO("npts = " << npts);
    QuantityType g = makeSignal(npts);
    QuantityType f = fftgtof(g, rstep);
    QuantityType fsum = sumgtof(g, rstep);
    REQUIRE(fsum.size() == f.size());
    CHECK(0.0 == f[0]);
    CHECK(maxRelativeError(f, fsum) < RELTOL);
  }
}

TEST_CASE("TestPDFUtils::test_fftgtof_rmin", "[PDFUtils]") {
  const double rstep = 0.02;
  QuantityType g = makeSignal(200);
  // rmin is rounded to the nearest multiple of rstep
  const double rmins[] = {0.0, 0.02, 0.029, 0.031, 1.0, 1.12, 5.6};
  const int padrmins[] = {0, 1, 1, 2, 50, 56, 280};
  for (int i = 0; i < 7; ++i) {
    INFO("rmin = " << rmins[i]);
    QuantityType f = fftgtof(g, rstep, rmins[i]);
    QuantityType fsum = sumgtof(g, rstep, padrmins[i]);
    REQUIRE(fsum.size() == f.size());
    CHECK(maxRelativeError(f, fsum) < RELTOL);
  }
  // shifted signal is equivalent to explicit zero padding
  QuantityType gpad(56, 0.0);
  gpad.insert(gpad.end(), g.begin(), g.end());
  CHECK(fftgtof(gpad, rstep) == fftgtof(g, rstep, 1.12));
}

TEST_CASE("TestPDFUtils::test_fftgtof_edge_points", "[PDFUtils]") {
  // Edge-point signals have the largest round-off; it must not grow with N.
  const double rstep = 0.01;
  for (int npts : {2048, 1 << 16}) {
    for (int k : {1, 2, npts - 2, npts - 1}) {
      INFO("npts = " << npts << ", signal at point " << k);
      QuantityType g(npts, 0.0);
      g[k] = 1.0;
      // the exact transformation of a single point
      QuantityType fexact(npts);
      for (int j = 0; j < npts; ++j) {
        const long long jk = (static_cast<long long>(j) * k) % (2 * npts);
        fexact[j] = rstep * sin(M_PI * double(jk) / npts);
      }
      CHECK(maxRelativeError(fftgtof(g, rstep), fexact) < 1e-13);
      // The zero-coordinate sample must not contaminate the sine sum.
      g[0] = 1e100;
      CHECK(maxRelativeError(fftgtof(g, rstep), fexact) < 1e-13);
    }
  }
}

TEST_CASE("TestPDFUtils::test_fftftog", "[PDFUtils]") {
  const double qstep = 0.05;
  QuantityType f = makeSignal(300);
  QuantityType g = fftftog(f, qstep, 0.5);
  QuantityType gsum = sumgtof(f, qstep, 10);
  for (double& gi : gsum) gi *= 2.0 / M_PI;
  REQUIRE(gsum.size() == g.size());
  CHECK(maxRelativeError(g, gsum) < RELTOL);
}

TEST_CASE("TestPDFUtils::test_fft_threads", "[PDFUtils]") {
  // Same length, different data and steps: threads must not share FFT state.
  constexpr int nthreads = 4;
  constexpr int npts = 1024;
  array<QuantityType, nthreads> signals, expected;
  array<double, nthreads> steps;
  for (int i = 0; i < nthreads; ++i) {
    signals[i] = makeSignal(npts);
    for (double& v : signals[i]) v = (i + 1) * v + 0.1 * i;
    signals[i][0] = 0.0;
    steps[i] = 0.01 * (i + 1);
    expected[i] = sumgtof(signals[i], steps[i]);
  }
  barrier ready(nthreads);
  array<future<double>, nthreads> workers;
  for (int i = 0; i < nthreads; ++i) {
    workers[i] = async(launch::async, [&, i] {
      ready.arrive_and_wait();
      double error = 0.0;
      for (int repeat = 0; repeat < 16; ++repeat) {
        const auto f = fftgtof(signals[i], steps[i]);
        error = max(error, maxRelativeError(f, expected[i]));
        const auto g = fftftog(f, M_PI / (npts * steps[i]));
        error = max(error, maxRelativeError(g, signals[i]));
      }
      return error;
    });
  }
  for (auto& worker : workers) CHECK(worker.get() < RELTOL);
}

// End of file
