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
 * CUDA calculation kernels, which add all pairs of a structure on a device
 *
 * The device computes the distance and width of every PairTable pair in
 * double precision.  A periodic pair list stays on the device while its
 * neighbor search remains valid.  Each output point is summed in a fixed
 * order: PDF peaks in per-warp rows, Debye terms in blocks of pairs, then
 * the rows or blocks in order.  The mixed kernel evaluates peaks and Debye
 * terms in single precision from arguments reduced in double precision.
 *
 *****************************************************************************/

#include <cuda_runtime.h>

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

#include <diffpy/srreal/kernels/PairGeometry.hpp>

namespace diffpy {
namespace srreal {
namespace kernels {

namespace {

// Device ----------------------------------------------------------------------

void check(cudaError_t err, const char* what) {
  if (err == cudaSuccess) return;
  throw std::runtime_error(std::string("the CUDA device failed in ") + what +
                           " with '" + cudaGetErrorString(err) + "'");
}

int currentDevice() {
  int device;
  check(cudaGetDevice(&device), "the current device lookup");
  return device;
}

/// Memory on the device that grows as needed.
template <class T>
class Buffer {
 public:
  Buffer() {}
  Buffer(const Buffer&) = delete;
  Buffer& operator=(const Buffer&) = delete;
  ~Buffer() {
    if (mptr) cudaFree(mptr);
  }

  T* reserve(size_t n) {
    if (n > mcapacity) {
      if (mptr) cudaFree(mptr);
      mptr = nullptr;
      mcapacity = 0;
      void* p = nullptr;
      check(cudaMalloc(&p, n * sizeof(T)), "the allocation of memory");
      mptr = static_cast<T*>(p);
      mcapacity = n;
    }
    return mptr;
  }

  T* upload(const T* src, size_t n) {
    this->reserve(n);
    if (n) {
      check(cudaMemcpy(mptr, src, n * sizeof(T), cudaMemcpyHostToDevice),
            "the copy of data to it");
    }
    return mptr;
  }

  T* zero(size_t n) {
    this->reserve(n);
    if (n) check(cudaMemset(mptr, 0, n * sizeof(T)), "the reset of memory");
    return mptr;
  }

  void download(T* dst, size_t n) const {
    if (n) {
      check(cudaMemcpy(dst, mptr, n * sizeof(T), cudaMemcpyDeviceToHost),
            "the copy of data from it");
    }
  }

  T* get() const { return mptr; }

 private:
  T* mptr = nullptr;
  size_t mcapacity = 0;
};

// Pairs -----------------------------------------------------------------------

/// Maximum pairs in a Debye chunk.
const long long CHUNK = 1LL << 20;
/// Threads of a block.
const int THREADS = 128;
/// Pairs that one thread of the Debye sum adds for its Q point.
const int SUBCHUNK = 1024;

const double FOUR_LN2 = 2.7725887222397811;
const double TWO_PI = 6.2831853071795865;

// Peaks of the PDF
// --------------------------------------------------------------

struct PeakGrid {
  int npts;
  int rlosteps;
  double rstep;
  double height;
  int cropped;
  double halfboundrel;
};

// Keep the original lower grid index as the origin for sampling the peak.
struct PeakWindow {
  int origin, first, last;
  double xlo;

  __device__ double x(int i, double step) const {
    return xlo + (i - origin) * step;
  }
};

__device__ PeakWindow peakWindow(const PeakGrid& grid, double dist,
                                 double fwhm) {
  const double half = grid.halfboundrel * fwhm;
  const int lo =
    gridIndex(dist - half, grid.rstep, grid.rlosteps, 0, grid.npts);
  const int hi =
    gridIndex(dist + half, grid.rstep, grid.rlosteps, 1, grid.npts);
  PeakWindow win{lo, lo, hi, (double(grid.rlosteps) + lo) * grid.rstep - dist};
  if (grid.cropped)
    cropGaussianWindow(win.first, win.last, grid.rlosteps, grid.rstep, dist,
                       fwhm, grid.halfboundrel);
  return win;
}

__device__ double samplePeak(double x, double dist, double fwhm,
                             double amplitude) {
  const double xrel = x / fwhm;
  return amplitude * exp(-FOUR_LN2 * xrel * xrel) * (x / dist + 1.0);
}

enum { WARP = 32 };
struct Peak {
  PeakWindow win;
  double dist, fwhm, amplitude, invdist, invwidth;
};

/// Each warp groups bit-identical distances and widths, adds their weights
/// in lane order, and accumulates the resulting peaks in its own output row.
template <bool Mixed>
__global__ void pdfPeaks(PairGeometryArgs g, PairWidthArgs w,
                         const double* sfsite, PeakGrid grid, long long count,
                         double* rows) {
  __shared__ Peak peaks[THREADS];
  __shared__ double scales[THREADS];
  const int lane = threadIdx.x % WARP;
  const int warp = (blockIdx.x * THREADS + threadIdx.x) / WARP;
  double* row = rows + (size_t)warp * grid.npts;
  Peak* batch = peaks + threadIdx.x - lane;
  double* weights = scales + threadIdx.x - lane;
  for (long long first = (long long)warp * WARP; first < count;
       first += (long long)gridDim.x * THREADS) {
    Peak pk{};
    double scale = 0.;
    if (first + lane < count) {
      const Pair p = pairAt(g, first + lane);
      if (p.ok) {
        pk.dist = p.dist;
        pk.fwhm = widthOf(w, p.dist, p.msd);
        scale = sfsite[p.site0] * sfsite[p.site1] * p.weight;
      }
    }
    const unsigned valid = __ballot_sync(0xffffffffu, pk.fwhm > 0.);
    const unsigned distMask =
      __match_any_sync(0xffffffffu, __double_as_longlong(pk.dist));
    const unsigned widthMask =
      __match_any_sync(0xffffffffu, __double_as_longlong(pk.fwhm));
    const unsigned same = valid & distMask & widthMask;
    const bool leader = lane == __ffs(same) - 1;
    unsigned leaders = __ballot_sync(0xffffffffu, leader);
    weights[lane] = scale;
    __syncwarp();
    if (leader) {
      for (unsigned others = same & (same - 1); others; others &= others - 1)
        scale += weights[__ffs(others) - 1];
      pk.win = peakWindow(grid, pk.dist, pk.fwhm);
      pk.invdist = 1.0 / pk.dist;
      pk.invwidth = 1.0 / pk.fwhm;
      pk.amplitude = scale * grid.height / pk.fwhm;
      batch[lane] = pk;
    }
    __syncwarp();
    for (; leaders; leaders &= leaders - 1) {
      const Peak p = batch[__ffs(leaders) - 1];
      const int offset = (lane - p.win.first) & (WARP - 1);
      for (long long i = (long long)p.win.first + offset; i < p.win.last;
           i += WARP) {
        const double x = p.win.x(int(i), grid.rstep);
        // An unresolved peak can have an unrepresentable inverse width.
        if (p.win.last - p.win.first <= 2) {
          row[i] += samplePeak(x, p.dist, p.fwhm, p.amplitude);
        } else {
          using Real = std::conditional_t<Mixed, float, double>;
          const Real xr = Real(x * p.invwidth);
          const Real correction = Real(x) * Real(p.invdist) + Real(1);
          row[i] +=
            Real(p.amplitude) * exp(-Real(FOUR_LN2) * xr * xr) * correction;
        }
      }
    }
    __syncwarp();
  }
}

// Terms of the Debye sum
// --------------------------------------------------------

struct QGrid {
  int nq;
  int kqlo;
  int kqhi;
  double qstep;
  double prec;
};

/// The double kernel stores distance in angstroms. The mixed kernel stores
/// the phase step as a fraction of a turn in units of 2^-64.
template <class Real>
struct DebyePair {
  using Phase = std::conditional_t<std::is_same_v<Real, double>, double,
                                   unsigned long long>;
  Phase phase;
  Real sigma, amplitude;
  int type0, type1, kqstop;
};

template <class Real>
struct DebyeBuffers {
  Buffer<DebyePair<Real>> pairs;
  Buffer<Real> factors;
};

/// Set up the pairs of a chunk. logsf contains log(|sf|) for each type.
template <class Real>
__global__ void debyePairs(PairGeometryArgs g, PairWidthArgs w,
                           const int* typeofsite, const int* curve,
                           const double* occ, const double* sf,
                           const float* logsf, QGrid grid, long long first,
                           long long count, DebyePair<Real>* pairs) {
  const long long j = (long long)blockIdx.x * blockDim.x + threadIdx.x;
  if (j >= count) return;
  const Pair p = pairAt(g, first + j);
  if (!p.ok) {
    pairs[j].kqstop = grid.kqlo;
    return;
  }
  const double sigma = fwhmToSigma(widthOf(w, p.dist, p.msd));
  const int tp0 = typeofsite[p.site0];
  const int tp1 = typeofsite[p.site1];
  const double amplitude = p.weight / p.dist;
  DebyePair<Real> out;
  if constexpr (std::is_same_v<Real, double>) {
    out.phase = p.dist;
  } else {
    const double turns = p.dist * grid.qstep / TWO_PI;
    out.phase =
      (unsigned long long)((turns - floor(turns)) * 18446744073709551616.0);
  }
  out.sigma = Real(sigma);
  out.amplitude = Real(amplitude);
  out.type0 = tp0;
  out.type1 = tp1;
  // The sum of a pair stops at the first kq with
  //     |amplitude * exp(-(sigma q)^2 / 2) * sf0 * sf1| <= prec.
  // Its logarithm in single precision finds where that may be, then the
  // terms themselves decide.
  const double* sf0 = sf + (size_t)curve[tp0] * grid.nq;
  const double* sf1 = sf + (size_t)curve[tp1] * grid.nq;
  const float* lsf0 = logsf + (size_t)curve[tp0] * grid.nq;
  const float* lsf1 = logsf + (size_t)curve[tp1] * grid.nq;
  const float decay = (float)(0.5 * sigma * sigma * grid.qstep * grid.qstep);
  const float limit = (grid.prec > 0.0)
                        ? (float)(log(grid.prec / fabs(amplitude)) -
                                  log(fabs(occ[tp0])) - log(fabs(occ[tp1])))
                        : -800.0f;
  int kq = grid.kqlo;
  for (; kq < grid.kqhi; ++kq) {
    const float lterm = lsf0[kq] + lsf1[kq] - decay * (float)kq * (float)kq;
    if (!(lterm > limit + 0.01f)) break;
  }
  for (; kq < grid.kqhi; ++kq) {
    const double sq = sigma * (kq * grid.qstep);
    const double term = amplitude * exp(-0.5 * sq * sq) * (sf0[kq] * occ[tp0]) *
                        (sf1[kq] * occ[tp1]);
    if (fabs(term) <= grid.prec) break;
  }
  out.kqstop = kq;
  pairs[j] = out;
}

/// Apply each occupancy once per Q tile in the selected precision.
template <class Real>
__global__ void occupiedFactors(const double* sf, const int* curve,
                                const double* occupancy, QGrid tile, int ntypes,
                                Real* out) {
  const int span = tile.kqhi - tile.kqlo;
  const int q = blockIdx.x * blockDim.x + threadIdx.x;
  if (q >= span) return;
  for (int t = blockIdx.y; t < ntypes; t += gridDim.y) {
    out[(size_t)t * span + q] =
      Real(sf[(size_t)curve[t] * tile.nq + tile.kqlo + q] * occupancy[t]);
  }
}

__device__ double debyeSine(double distance, double q, int) {
  return sin(q * distance);
}

__device__ float debyeSine(unsigned long long turns, float, int kq) {
  // Fold the phase to a quarter turn before converting to single precision.
  const long long phase = (long long)((unsigned long long)kq * turns);
  const long long quarter = 1LL << 62;
  const long long folded = (phase > quarter || phase < -quarter)
                             ? (long long)(0x8000000000000000ULL - phase)
                             : phase;
  return sinpif((float)folded * 0x1p-63f);
}

/// Each thread adds one Q point for a contiguous subchunk of pairs.
template <class Real>
__global__ void debyeTerms(const DebyePair<Real>* pairs, const Real* sf,
                           QGrid grid, int count, double* partial) {
  const int i = blockIdx.x * blockDim.x + threadIdx.x;
  const int span = grid.kqhi - grid.kqlo;
  if (i >= span) return;
  const int kq = grid.kqlo + i;
  const Real q = Real(kq * grid.qstep);
  const int first = blockIdx.y * SUBCHUNK;
  const int end = min(first + SUBCHUNK, count);
  double sum = 0.0;
  for (int j = first; j < end; ++j) {
    if (kq >= pairs[j].kqstop) continue;
    const auto& p = pairs[j];
    const Real sq = p.sigma * q;
    sum += p.amplitude * exp(Real(-0.5) * sq * sq) *
           sf[(size_t)p.type0 * span + i] * sf[(size_t)p.type1 * span + i] *
           debyeSine(p.phase, q, kq);
  }
  partial[(size_t)blockIdx.y * span + i] = sum;
}

/// Add the sums of the blocks of pairs in their order.
__global__ void finishRows(const double* partial, int cntsub, int span,
                           double* total) {
  const int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i >= span) return;
  double sum = 0.0;
  for (int sub = 0; sub < cntsub; ++sub) sum += partial[(size_t)sub * span + i];
  total[i] += sum;
}

int blocksFor(long long n) { return (int)((n + THREADS - 1) / THREADS); }

// Sums ------------------------------------------------------------------------

/// What the kernels keep on the device between the sums of a calculator.
struct DeviceState : PairSumWorkspace {
  const int device = currentDevice();
  bool isCurrentDevice() const override { return device == currentDevice(); }

  // structure
  Buffer<int> sitebegin;
  Buffer<int> possite;
  Buffer<double> pos;
  Buffer<double> posuij;
  Buffer<double> siteuij;
  Buffer<unsigned char> anisotropy;
  Buffer<int> multiplicity;
  Buffer<double> trans;
  Buffer<int> pairsite0;
  Buffer<int> pairpos;
  Buffer<int> pairtrans;
  // the list of pairs on the device
  const int* listowner = nullptr;
  int listversion = -1;
  long long listsize = -1;
  // peaks
  Buffer<double> sfsite;
  Buffer<double> total;
  // Debye sum
  Buffer<int> typeofsite, curve;
  Buffer<double> occupancy;
  Buffer<double> sf;
  Buffer<float> logsf;
  DebyeBuffers<double> debye;
  DebyeBuffers<float> debyeMixed;
  Buffer<double> partial;

  /// Copy the structure to the device.  Return the number of its pairs.
  long long upload(const PairGeometryArgs& v, PairGeometryArgs& g) {
    g = v;
    g.sitebegin = sitebegin.upload(v.sitebegin, v.cntsites + 1);
    g.possite = possite.upload(v.possite, v.cntpositions);
    g.pos = pos.upload(v.pos, 3 * (size_t)v.cntpositions);
    g.posuij = posuij.upload(v.posuij, 6 * (size_t)v.cntpositions);
    g.siteuij = siteuij.upload(v.siteuij, 6 * (size_t)v.cntsites);
    g.anisotropy = anisotropy.upload(v.anisotropy, v.cntsites);
    g.multiplicity = multiplicity.upload(v.multiplicity, v.cntsites);
    g.trans = trans.upload(v.trans, 3 * (size_t)v.cnttrans);
    if (!v.translations) {
      listversion = -1;
      g.pairsite0 = g.pairpos = g.pairtrans = nullptr;
      return pairCount(v);
    }
    // The list stays on the device for the structures that keep it.
    const bool same = (listowner == v.pairsite0) &&
                      (listversion == v.pairversion) &&
                      (listsize == v.cntpairs);
    if (!same) {
      pairsite0.upload(v.pairsite0, v.cntpairs);
      pairpos.upload(v.pairpos, v.cntpairs);
      pairtrans.upload(v.pairtrans, v.cntpairs);
      listowner = v.pairsite0;
      listversion = v.pairversion;
      listsize = v.cntpairs;
    }
    g.pairsite0 = pairsite0.get();
    g.pairpos = pairpos.get();
    g.pairtrans = pairtrans.get();
    return v.cntpairs;
  }
};

template <bool Mixed>
void addPDFPairs(const PDFPairSumArgs& job) {
  if (job.npts <= 0 || !pairCount(job.g)) return;
  auto& m = static_cast<DeviceState&>(*job.workspace);
  PairGeometryArgs g;
  const long long cntpairs = m.upload(job.g, g);
  const double* sfsite = m.sfsite.upload(job.sfsite, job.g.cntsites);
  const PeakGrid grid = {job.npts,   job.rlosteps, job.rstep,
                         job.height, job.cropped,  job.halfboundrel};
  // Bound temporary storage to 32 MiB or four rows for a very long grid.
  const int rowLimit = std::max(4, std::min(4096, (4 << 20) / job.npts));
  const int blocks = int(std::min<long long>((cntpairs + THREADS - 1) / THREADS,
                                             rowLimit / (THREADS / WARP)));
  const int rows = blocks * (THREADS / WARP);
  double* partial = m.partial.zero((size_t)rows * job.npts);
  double* total = m.total.zero(job.npts);
  pdfPeaks<Mixed>
    <<<blocks, THREADS>>>(g, job.w, sfsite, grid, cntpairs, partial);
  check(cudaGetLastError(), "the start of a PDF kernel");
  finishRows<<<blocksFor(job.npts), THREADS>>>(partial, rows, job.npts, total);
  check(cudaGetLastError(), "the reduction of the PDF");
  std::vector<double> result(job.npts);
  m.total.download(result.data(), job.npts);
  if (std::any_of(result.begin(), result.end(),
                  [](double x) { return !std::isfinite(x); }))
    throw std::runtime_error("nonfinite PDF sum on the CUDA device");
  for (int i = 0; i < job.npts; ++i) job.out[i] += result[i];
}

template <class Real>
void addDebyePairs(const DebyePairSumArgs& job, DeviceState& m,
                   DebyeBuffers<Real>& work) {
  const int span = job.kqhi - job.kqlo;
  if (span <= 0) return;
  PairGeometryArgs g;
  const long long cntpairs = m.upload(job.g, g);
  if (!cntpairs) return;
  const QGrid grid = {job.nq, job.kqlo, job.kqhi, job.qstep, job.prec};
  const int* typeofsite = m.typeofsite.upload(job.typeofsite, job.g.cntsites);
  const int* curve = m.curve.upload(job.curve, job.ntypes);
  const double* occ = m.occupancy.upload(job.occupancy, job.ntypes);
  const size_t sfsize = (size_t)job.ncurves * job.nq;
  const double* sf = m.sf.upload(job.sf, sfsize);
  std::vector<float> logsf(sfsize);
  for (size_t i = 0; i < sfsize; ++i) {
    logsf[i] = (float)std::log(std::fabs(job.sf[i]));
  }
  const float* dlogsf = m.logsf.upload(logsf.data(), logsf.size());
  const size_t chunk = (size_t)std::min(CHUNK, std::max(cntpairs, 1LL));
  auto* pairs = work.pairs.reserve(chunk);
  const size_t maxsub = (chunk + SUBCHUNK - 1) / SUBCHUNK;
  // Scratch is bounded even when each site has its own occupancy.
  const int tileSize = std::min(span, std::max(1, (1 << 20) / job.ntypes));
  Real* weighted = work.factors.reserve((size_t)tileSize * job.ntypes);
  double* partial = m.partial.reserve(maxsub * tileSize);
  double* total = m.total.zero(span);
  for (long long first = 0; first < cntpairs; first += CHUNK) {
    const long long count = std::min(CHUNK, cntpairs - first);
    const int cntsub = (int)((count + SUBCHUNK - 1) / SUBCHUNK);
    debyePairs<<<blocksFor(count), THREADS>>>(
      g, job.w, typeofsite, curve, occ, sf, dlogsf, grid, first, count, pairs);
    check(cudaGetLastError(), "the start of a kernel");
    for (int lo = grid.kqlo; lo < grid.kqhi;) {
      const int n = std::min(tileSize, grid.kqhi - lo);
      QGrid tile = grid;
      tile.kqlo = lo;
      tile.kqhi = lo + n;
      const dim3 factorBlocks(blocksFor(n), std::min(job.ntypes, 65535));
      occupiedFactors<<<factorBlocks, THREADS>>>(sf, curve, occ, tile,
                                                 job.ntypes, weighted);
      check(cudaGetLastError(), "the occupied scattering factors");
      const dim3 blocks(blocksFor(n), cntsub);
      debyeTerms<<<blocks, THREADS>>>(pairs, weighted, tile, (int)count,
                                      partial);
      check(cudaGetLastError(), "the start of a kernel");
      finishRows<<<blocksFor(n), THREADS>>>(partial, cntsub, n,
                                            total + (lo - grid.kqlo));
      check(cudaGetLastError(), "the reduction of the Debye sum");
      lo += n;
    }
  }
  check(cudaDeviceSynchronize(), "the Debye sum");
  std::vector<double> result(span);
  m.total.download(result.data(), span);
  for (int i = 0; i < span; ++i) job.out[job.kqlo + i] += result[i];
}

}  // namespace

// Public Functions ------------------------------------------------------------

const char* cudaKernelUnavailable() {
  // Loading an actual kernel checks the binary or PTX image against the
  // calling thread's current device and driver without launching any work.
  cudaFuncAttributes attributes;
  const cudaError_t err = cudaFuncGetAttributes(&attributes, pdfPeaks<false>);
  return err == cudaSuccess ? nullptr : cudaGetErrorString(err);
}

PairSumStats cudaPDFPairSum(const PDFPairSumArgs& a) {
  addPDFPairs<false>(a);
  return {};
}

PairSumStats cudaPDFPairSumMixed(const PDFPairSumArgs& a) {
  addPDFPairs<true>(a);
  return {};
}

PairSumStats cudaDebyePairSum(const DebyePairSumArgs& a) {
  auto& state = static_cast<DeviceState&>(*a.workspace);
  addDebyePairs(a, state, state.debye);
  return {};
}

PairSumStats cudaDebyePairSumMixed(const DebyePairSumArgs& a) {
  auto& state = static_cast<DeviceState&>(*a.workspace);
  addDebyePairs(a, state, state.debyeMixed);
  return {};
}

PairSumWorkspace* makeCudaWorkspace() { return new DeviceState; }

}  // namespace kernels
}  // namespace srreal
}  // namespace diffpy

// End of file
