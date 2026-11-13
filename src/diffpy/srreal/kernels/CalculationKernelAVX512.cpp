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
 * AVX-512 calculation kernel.  Include nothing else here.
 *
 *****************************************************************************/

// xsimd namespace of this file, see CalculationKernelSIMD.ipp
#define LIBDIFFPY_XSIMD_NAMESPACE xsimd_avx512
#define LIBDIFFPY_XSIMD_AVX512 1
#include "CalculationKernelSIMD.ipp"

namespace diffpy {
namespace srreal {
namespace kernels {

const CalculationKernel& avx512CalculationKernel() {
  static constexpr CalculationKernel kernel =
    makeSIMDKernel<xsimd::avx512dq>("avx512");
  return kernel;
}

}  // namespace kernels
}  // namespace srreal
}  // namespace diffpy
