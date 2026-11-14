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
 * SSE2 calculation kernel.  Include nothing else here.
 *
 *****************************************************************************/

// xsimd namespace of this file, see CalculationKernelSIMD.ipp
#define LIBDIFFPY_XSIMD_NAMESPACE xsimd_sse2
#include "CalculationKernelSIMD.ipp"

namespace diffpy {
namespace srreal {
namespace kernels {

const CalculationKernel& sse2CalculationKernel() {
  static constexpr CalculationKernel kernel =
    makeSIMDKernel<xsimd::sse2>("sse2");
  return kernel;
}

}  // namespace kernels
}  // namespace srreal
}  // namespace diffpy
