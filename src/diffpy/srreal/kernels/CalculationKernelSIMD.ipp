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
 * SIMD calculation kernel shared by the SSE2, AVX2, AVX-512 and NEON files.
 *
 * On the uniform grid exp and sin are advanced by recurrences instead of
 * being called per point:
 *   g = exp(-b x^2): g(x + h) = g(x) rho(x), rho(x + h) = rho(x) exp(-2 b h^2)
 *   sin, cos of x+d: s' = s cos(d) + c sin(d),  c' = c cos(d) - s sin(d)
 * xsimd::exp and xsimd::sincos reseed them every RESEED_STEPS vector steps
 * to limit the rounding error.  Windows too far from the peak for the
 * recurrence go to the scalar kernel.
 *
 * The linker keeps one copy of each inline function for the whole library,
 * so a copy compiled here could run AVX-512 code on a CPU without it.  To
 * avoid that, each file renames the xsimd namespace to
 * LIBDIFFPY_XSIMD_NAMESPACE and is compiled without link-time optimization.
 * GCC and Clang inline the few stl functions used here, so
 * use only C math functions.  MSVC compiles these files for SSE2,
 * see CalculationKernelMSVC.hpp.
 *
 * Keep scalar floating point code out of the vector loops, there it makes
 * MSVC pass the batches of the loop through the stack.
 *
 *****************************************************************************/

#if defined(__FAST_MATH__) || defined(_M_FP_FAST)
#error "The xsimd math functions give wrong results with fast-math."
#endif

#if (defined(__GNUC__) || defined(__clang__)) && !defined(__OPTIMIZE__)
#error "Compile the SIMD kernels with optimization, see the comment above."
#endif

#ifndef LIBDIFFPY_XSIMD_NAMESPACE
#error "Define LIBDIFFPY_XSIMD_NAMESPACE, see the comment above."
#endif

#include <math.h>
#include <stdint.h>

#include <diffpy/srreal/kernels/CalculationKernelImpl.hpp>

// The xsimd headers and all code below use the renamed namespace.
#define xsimd LIBDIFFPY_XSIMD_NAMESPACE
#if defined(_MSC_VER) && !defined(__clang__) && defined(_M_X64)
#include "CalculationKernelMSVC.hpp"
#endif
#include <xsimd/xsimd.hpp>

// Inlines a helper at every call.
#if defined(_MSC_VER)
#define LIBDIFFPY_KERNEL_INLINE __forceinline
#else
#define LIBDIFFPY_KERNEL_INLINE inline __attribute__((always_inline))
#endif

namespace diffpy {
namespace srreal {
namespace kernels {

namespace {

/// Vector steps in one block of the recurrences.
constexpr int RESEED_STEPS = 32;

/// Extra vector steps allowed in the last block.
constexpr int RESEED_SLACK = 8;

/// Limit of b x^2 at the seeds and of |log(rho)|, which keeps the seeds
/// normal and rho finite.  exp(-SAFE_EXPONENT) is about 1e-261.
constexpr double SAFE_EXPONENT = 600.0;

/// Lane numbers of the widest batch.
constexpr double LANE_INDEX[8] = {0, 1, 2, 3, 4, 5, 6, 7};

template <class Arch>
struct SIMDKernel {
  using B = xsimd::batch<double, Arch>;
  using BB = xsimd::batch_bool<double, Arch>;
  static constexpr int W = static_cast<int>(B::size);

  /// The batch {0, 1, ..., W - 1}.
  static B laneIndex() {
    static_assert(W <= 8, "LANE_INDEX is too short.");
    return B::load_unaligned(LANE_INDEX);
  }

  /// End of the block that starts at i, the last block ends at hi.
  static int blockEnd(int i, int hi) {
    const bool last = (hi - i <= (RESEED_STEPS + RESEED_SLACK) * W);
    return last ? hi : i + RESEED_STEPS * W;
  }

  /// End of the whole vectors in [i, iend), only the last block has a rest.
  static int wholeEnd(int i, int iend) { return iend - (iend - i) % W; }

  // Batches go by value, MSVC keeps them in memory when passed by reference.

  /// Load the first n < W lanes from p, the other lanes are zero.
  static LIBDIFFPY_KERNEL_INLINE B loadPartial(const double* p, int n) {
    double v[W] = {};
    for (int j = 0; j < n; ++j) v[j] = p[j];
    return B::load_unaligned(v);
  }

  /// Add u * v to out.
  static void addProduct(double* out, B u, B v) {
    xsimd::fma(u, v, B::load_unaligned(out)).store_unaligned(out);
  }

  /// Add the first n < W lanes of u * v to out.
  static LIBDIFFPY_KERNEL_INLINE void addProductPartial(double* out, B u, B v,
                                                        int n) {
    double uv[W];
    (u * v).store_unaligned(uv);
    for (int j = 0; j < n; ++j) out[j] += uv[j];
  }

  /// Zero t from the first set lane of small on, small must have one.
  static B truncated(B t, BB small) {
    const uint64_t m = small.mask();
    return xsimd::select(BB::from_mask((m ^ (m - 1)) >> 1), t, B(0.0));
  }

  static void gaussianWindow(const GaussianWindowArgs& args) {
    // Vector stores may alias anything, work on a copy of the arguments.
    const GaussianWindowArgs a = args;
    if (a.fwhm <= 0.0) return;
    // exp(-b x^2) and its ratio over one vector step.
    const double b = 4 * M_LN2 / (a.fwhm * a.fwhm);
    const double wh = W * a.rstep;
    // Use the scalar kernel for NaN and for windows too far from the peak,
    // where b * x^2 or 2 * b * wh * |x| reach SAFE_EXPONENT.
    const double xlo = (a.rlosteps + a.ilo) * a.rstep - a.dist;
    const double xhi = (a.rlosteps + a.ihi) * a.rstep - a.dist;
    const double xfar = (-xlo > xhi) ? -xlo : xhi;
    const double reach = (xfar > 2 * wh) ? xfar : 2 * wh;
    if (!(b * xfar * reach < SAFE_EXPONENT)) {
      scalarCalculationKernel().gaussianWindow(args);
      return;
    }
    const B negb(-b);
    const B twowh(2 * wh);
    const B whwh(wh * wh);
    const B ratiostep(exp(-2 * b * wh * wh));
    const B vw(1.0 * W);
    const B vrstep(a.rstep);
    const B vxlo(xlo);
    const B vinvdist(1.0 / a.dist);
    const B vnorm(a.scale * a.peakheight);
    const B vfwhm(a.fwhm);
    const B vhalfbound(a.halfboundrel);
    const B one(1.0);
    const B zero(0.0);
    // Grid steps of the lanes from ilo, whole numbers are exact in double.
    B steps = laneIndex();
    // x / dist + 1 grows by a constant per vector step.
    const B vrscalestep(wh / a.dist);
    for (int i = a.ilo; i < a.ihi;) {
      const int iend = blockEnd(i, a.ihi);
      const int ivec = wholeEnd(i, iend);
      // Same x as the scalar kernel, so cropping agrees.
      const B xseed = steps * vrstep + vxlo;
      B g = vnorm * xsimd::exp(negb * xseed * xseed);
      B ratio = xsimd::exp(negb * (twowh * xseed + whwh));
      if (a.cropped) {
        for (; i < ivec; i += W) {
          const B x = steps * vrstep + vxlo;
          const B y =
            xsimd::select(xsimd::abs(x / vfwhm) >= vhalfbound, zero, g);
          addProduct(a.out + i, y, xsimd::fma(x, vinvdist, one));
          g *= ratio;
          ratio *= ratiostep;
          steps += vw;
        }
        if (i < iend) {
          const B x = steps * vrstep + vxlo;
          const B y =
            xsimd::select(xsimd::abs(x / vfwhm) >= vhalfbound, zero, g);
          addProductPartial(a.out + i, y, xsimd::fma(x, vinvdist, one),
                            iend - i);
          return;
        }
      } else {
        B rscaled = xsimd::fma(xseed, vinvdist, one);
        // Advance before the loop, the next seeds then do not wait for it.
        steps += B(double(ivec - i));
        for (; i < ivec; i += W) {
          addProduct(a.out + i, g, rscaled);
          g *= ratio;
          ratio *= ratiostep;
          rscaled += vrscalestep;
        }
        if (i < iend) {
          addProductPartial(a.out + i, g, rscaled, iend - i);
          return;
        }
      }
    }
  }

  static void debyeSum(const DebyeSumArgs& a) {
    // Debye-Waller factor exp(-alpha kq^2) and its ratio over W steps.
    const double dwstep = a.dwsigma * a.qstep;
    const double alpha = 0.5 * dwstep * dwstep;
    const B negalpha(-alpha);
    const B twow(2.0 * W);
    const B ww(1.0 * W * W);
    const B ratiostep(exp(-2 * alpha * W * W));
    // Rotation of (sin, cos) of q * dist over W steps.
    const double wangle = W * a.qstep * a.dist;
    const B cosw(cos(wangle));
    const B sinw(sin(wangle));
    const B vqstep(a.qstep);
    const B vdist(a.dist);
    const B vdwsigma(a.dwsigma);
    const B vsinescale(a.scale / a.dist);
    const B vprec(a.prec);
    const B neghalf(-0.5);
    const B lanes = laneIndex();
    for (int kq = a.kqlo; kq < a.kqhi;) {
      const int kqend = blockEnd(kq, a.kqhi);
      const int kqvec = wholeEnd(kq, kqend);
      const B k = B(double(kq)) + lanes;
      const B q = k * vqstep;
      const B dwq = vdwsigma * q;
      // g includes the sine scale.
      B g = vsinescale * xsimd::exp(neghalf * dwq * dwq);
      B ratio = xsimd::exp(negalpha * (twow * k + ww));
      auto [s, c] = xsimd::sincos(q * vdist);
      for (; kq < kqvec; kq += W) {
        const B t =
          g * B::load_unaligned(a.sf0 + kq) * B::load_unaligned(a.sf1 + kq);
        const BB small = xsimd::abs(t) <= vprec;
        if (xsimd::any(small)) {
          // Keep the lanes before the first small term.
          addProduct(a.out + kq, truncated(t, small), s);
          return;
        }
        addProduct(a.out + kq, t, s);
        g *= ratio;
        ratio *= ratiostep;
        const B snext = xsimd::fma(s, cosw, c * sinw);
        c = xsimd::fms(c, cosw, s * sinw);
        s = snext;
      }
      if (kq < kqend) {
        const int n = kqend - kq;
        // Unused lanes have zero t and do not reach the output.
        B t = g * loadPartial(a.sf0 + kq, n) * loadPartial(a.sf1 + kq, n);
        const BB small = xsimd::abs(t) <= vprec;
        if (xsimd::any(small)) t = truncated(t, small);
        addProductPartial(a.out + kq, t, s, n);
        return;
      }
    }
  }

  static int bondCandidates(const BondCandidatesArgs& a) {
    const B vox(a.ox);
    const B voy(a.oy);
    const B voz(a.oz);
    const B vlo(a.d2lo);
    const B vhi(a.d2hi);
    int m = 0;
    int i = 0;
    for (; i + W <= a.n; i += W) {
      const B dx = B::load_unaligned(a.x + i) + vox;
      const B dy = B::load_unaligned(a.y + i) + voy;
      const B dz = B::load_unaligned(a.z + i) + voz;
      const B d2 = xsimd::fma(dz, dz, xsimd::fma(dy, dy, dx * dx));
      // NaN fails both comparisons and stays a candidate.
      const uint64_t bits = (!((d2 < vlo) || (d2 > vhi))).mask();
      if (!bits) continue;
      for (int j = 0; j < W; ++j) {
        a.idx[m] = i + j;
        m += int((bits >> j) & 1);
      }
    }
    // The last points go one by one, a partial vector load would be much
    // slower.
    for (; i < a.n; ++i) {
      const double dx = a.x[i] + a.ox;
      const double dy = a.y[i] + a.oy;
      const double dz = a.z[i] + a.oz;
      const double d2 = dx * dx + dy * dy + dz * dz;
      a.idx[m] = i;
      m += int(!(d2 < a.d2lo || d2 > a.d2hi));
    }
    return m;
  }
};

/// Kernel table for the Arch instruction set.
template <class Arch>
constexpr CalculationKernel makeSIMDKernel(const char* name) {
  using K = SIMDKernel<Arch>;
  return {name, K::gaussianWindow, K::debyeSum, K::bondCandidates};
}

}  // namespace

}  // namespace kernels
}  // namespace srreal
}  // namespace diffpy
