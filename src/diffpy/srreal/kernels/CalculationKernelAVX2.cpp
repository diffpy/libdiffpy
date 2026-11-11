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
 * AVX2+FMA calculation kernel.  Include nothing else here.
 *
 *****************************************************************************/

// xsimd namespace of this file, see CalculationKernelSIMD.ipp
#define LIBDIFFPY_XSIMD_NAMESPACE xsimd_avx2
#include "CalculationKernelSIMD.ipp"

namespace diffpy {
namespace srreal {
namespace kernels {

const CalculationKernel& avx2CalculationKernel() {
  static constexpr CalculationKernel kernel =
    makeSIMDKernel<xsimd::fma3<xsimd::avx2>>("avx2");
  return kernel;
}

}  // namespace kernels
}  // namespace srreal
}  // namespace diffpy
