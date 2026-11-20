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
 * Runtime selection of the PDF calculation kernel
 *
 *****************************************************************************/

#ifndef CALCULATIONKERNEL_HPP_INCLUDED
#define CALCULATIONKERNEL_HPP_INCLUDED

#include <string>
#include <vector>

#include <diffpy/Export.hpp>

namespace diffpy {
namespace srreal {

/// Names of the kernels that this build and machine can run, CPU kernels
/// first in order of preference.  A compatible CUDA device adds "cuda" and
/// "cuda-mixed", which only run full PDF and Debye sums.
DLL_EXPORT std::vector<std::string> availableCalculationKernels();

/// Name of the kernel in use, by default the first available one.
DLL_EXPORT std::string getCalculationKernel();

/// Use the named kernel for all calculators.
/// Throw std::invalid_argument when the kernel is not available.
/// Do not call while a calculation runs in another thread.
DLL_EXPORT void setCalculationKernel(const std::string& name);

}  // namespace srreal
}  // namespace diffpy

#endif  // CALCULATIONKERNEL_HPP_INCLUDED
