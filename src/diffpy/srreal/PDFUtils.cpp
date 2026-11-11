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
 * Various common routines useful for PDF calculation.
 *
 *****************************************************************************/

#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <cstdint>
#include <limits>
#include <memory>
#include <optional>
#include <stdexcept>
#include <vector>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_fft_complex.h>

#include <diffpy/srreal/PDFUtils.hpp>

using namespace std;

namespace diffpy {
namespace srreal {

// PDFUtils functions --------------------------------------------------------

QuantityType fftgtof(const QuantityType& g, double rstep, double rmin) {
  if (g.empty()) return g;
  // std::isfinite is always true with -ffast-math.
  const auto isfinitebits = [](double v) {
    const uint64_t expmask = 0x7ff0000000000000;
    return (bit_cast<uint64_t>(v) & expmask) != expmask;
  };
  if (!(isfinitebits(rstep) && rstep > 0.0)) {
    throw invalid_argument("rstep must be finite and positive.");
  }
  if (!isfinitebits(rmin)) {
    throw invalid_argument("rmin must be finite.");
  }
  // rmin within half a step below zero rounds to no padding.
  const double pad = round(rmin / rstep);
  if (pad < 0.0) {
    throw invalid_argument("rmin rounds to a negative offset.");
  }
  // Keep size_t(pad), bit_ceil and the FFT buffer of 2*N doubles in range.
  const size_t maxpoints = bit_floor(vector<double>().max_size() / 2);
  if (!isfinitebits(pad) || pad > maxpoints ||
      g.size() > maxpoints - size_t(pad)) {
    throw length_error("Too many points for the Fourier transformation.");
  }
  const size_t padding = size_t(pad);
  const size_t N = bit_ceil(padding + g.size());
  QuantityType f(N, 0.0);
  copy(g.begin(), g.end(), f.begin() + padding);
  // sin(0) removes f[0] from the sum; zero it to keep it out of round-off.
  f[0] = 0.0;
  if (N < 4) {
    if (N == 2) f[1] *= rstep;
    return f;
  }

  // Per-thread FFT buffers and tables for each power-of-two length.
  struct FFTCache {
    vector<double> data, wr, wi;
    unique_ptr<gsl_fft_complex_wavetable,
               decltype(&gsl_fft_complex_wavetable_free)>
      table;
    unique_ptr<gsl_fft_complex_workspace,
               decltype(&gsl_fft_complex_workspace_free)>
      work;

    explicit FFTCache(size_t npts)
        : data(2 * npts),
          wr(npts / 2 + 1),
          wi(npts / 2 + 1),
          table(gsl_fft_complex_wavetable_alloc(npts),
                gsl_fft_complex_wavetable_free),
          work(gsl_fft_complex_workspace_alloc(npts),
               gsl_fft_complex_workspace_free) {
      if (!table || !work) throw bad_alloc();
      for (size_t j = 0; j <= npts / 2; ++j) {
        wr[j] = cos(M_PI * j / npts);
        wi[j] = sin(M_PI * j / npts);
      }
    }
  };
  thread_local array<optional<FFTCache>, numeric_limits<size_t>::digits>
    transforms;
  auto& plan = transforms[countr_zero(N)];
  if (!plan) plan.emplace(N);
  FFTCache& cache = *plan;

  // Need Im X[m] for X[m] = sum_k f[k] exp(i*pi*m*k/N).  The N-point FFT Z
  // of z[j] = f[2*j] + i*f[2*j+1], padded with zeros, is enough.
  copy(f.begin(), f.end(), cache.data.begin());
  fill(cache.data.begin() + N, cache.data.end(), 0.0);
  const int status = gsl_fft_complex_backward(
    cache.data.data(), 1, N, cache.table.get(), cache.work.get());
  if (status != GSL_SUCCESS)
    throw invalid_argument("Fourier Transformation failed.");

  // Z[m] = E[m] + i*O[m] and conj(Z[N-m]) = E[m] - i*O[m] yield the even
  // and odd sample transforms E, O; X[m] = E[m] + exp(i*pi*m/N) * O[m].
  const double c = 0.5 * rstep;
  const auto& z = cache.data;
  for (size_t m = 1; m <= N / 2; ++m) {
    const size_t mc = N - m;
    const double ei = c * (z[2 * m + 1] - z[2 * mc + 1]);
    const double orr = c * (z[2 * m + 1] + z[2 * mc + 1]);
    const double oi = c * (z[2 * mc] - z[2 * m]);
    const double t = cache.wr[m] * oi + cache.wi[m] * orr;
    f[m] = ei + t;
    f[mc] = t - ei;
  }
  return f;
}

QuantityType fftftog(const QuantityType& f, double qstep, double qmin) {
  // ftog is the same as gtof, just normalized by 2 / PI.
  QuantityType g = fftgtof(f, qstep, qmin);
  for (double& gi : g) gi *= 2.0 / M_PI;
  return g;
}

}  // namespace srreal
}  // namespace diffpy

// End of file
