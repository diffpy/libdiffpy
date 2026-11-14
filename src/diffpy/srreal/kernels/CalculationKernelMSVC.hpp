/*****************************************************************************
 *
 * libdiffpy         Complex Modeling Initiative
 *                   (c) 2026 The Trustees of Columbia University
 *                   in the City of New York.  All rights reserved.
 *
 * See AUTHORS.rst for a list of people who contributed.
 * See LICENSE.rst for license information.
 *
 *****************************************************************************/

#ifndef CALCULATIONKERNELMSVC_HPP_INCLUDED
#define CALCULATIONKERNELMSVC_HPP_INCLUDED

// MSVC compiles AVX intrinsics without /arch:AVX*.  All other code of the
// SIMD kernels, such as shared standard library functions, then stays at SSE2.
#if defined(__AVX__) || defined(__AVX2__) || defined(__AVX512F__) ||   \
  defined(__AVX10_VER__) || defined(__SSE3__) || defined(__SSSE3__) || \
  defined(__SSE4_1__) || defined(__SSE4_2__)
#error "Compile the MSVC SIMD kernels for SSE2 without /arch:AVX*."
#endif

// Enable the instruction sets of the kernel in xsimd, which also needs the
// lower ones, before its other headers are included.
#include <xsimd/config/xsimd_config.hpp>

#if defined(LIBDIFFPY_XSIMD_AVX2) || defined(LIBDIFFPY_XSIMD_AVX512)
#undef XSIMD_WITH_SSE3
#define XSIMD_WITH_SSE3 1
#undef XSIMD_WITH_SSSE3
#define XSIMD_WITH_SSSE3 1
#undef XSIMD_WITH_SSE4_1
#define XSIMD_WITH_SSE4_1 1
#undef XSIMD_WITH_SSE4_2
#define XSIMD_WITH_SSE4_2 1
#undef XSIMD_WITH_AVX
#define XSIMD_WITH_AVX 1
#undef XSIMD_WITH_AVX2
#define XSIMD_WITH_AVX2 1
#undef XSIMD_WITH_FMA3_AVX
#define XSIMD_WITH_FMA3_AVX 1
#undef XSIMD_WITH_FMA3_AVX2
#define XSIMD_WITH_FMA3_AVX2 1
#endif

#if defined(LIBDIFFPY_XSIMD_AVX512)
// The subsets that CalculationKernelDispatch.cpp checks.
#undef XSIMD_WITH_AVX512F
#define XSIMD_WITH_AVX512F 1
#undef XSIMD_WITH_AVX512CD
#define XSIMD_WITH_AVX512CD 1
#undef XSIMD_WITH_AVX512DQ
#define XSIMD_WITH_AVX512DQ 1
#endif

#endif  // CALCULATIONKERNELMSVC_HPP_INCLUDED
