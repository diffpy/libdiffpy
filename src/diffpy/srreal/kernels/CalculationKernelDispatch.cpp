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
 * Runtime selection of the calculation kernel
 *
 *****************************************************************************/

#include <algorithm>
#include <atomic>
#include <stdexcept>
#include <string>
#include <vector>

#if defined(LIBDIFFPY_SIMD_KERNELS)
#include <xsimd/config/xsimd_cpuid.hpp>
#endif

#include <diffpy/srreal/CalculationKernel.hpp>
#include <diffpy/srreal/kernels/CalculationKernelImpl.hpp>

namespace diffpy {
namespace srreal {

using kernels::CalculationKernel;

namespace {

/// Kernels of this build that the CPU can run, from the most to the least
/// preferred.  CMake defines LIBDIFFPY_KERNEL_<NAME> for each SIMD kernel it
/// compiles.  The scalar kernel is last and always available.
const std::vector<const CalculationKernel*>& availableKernels() {
  static const std::vector<const CalculationKernel*> available = [] {
    std::vector<const CalculationKernel*> rv;
#if defined(LIBDIFFPY_SIMD_KERNELS)
    const auto cpu = xsimd::available_architectures();
#endif
#if defined(LIBDIFFPY_KERNEL_AVX512)
    // The AVX-512 file may also use AVX2 and FMA instructions.
    if (cpu.avx512f && cpu.avx512cd && cpu.avx512dq && cpu.fma3_avx2) {
      rv.push_back(&kernels::avx512CalculationKernel());
    }
#endif
#if defined(LIBDIFFPY_KERNEL_AVX2)
    if (cpu.fma3_avx2) rv.push_back(&kernels::avx2CalculationKernel());
#endif
#if defined(LIBDIFFPY_KERNEL_SSE2)
    if (cpu.sse2) rv.push_back(&kernels::sse2CalculationKernel());
#endif
#if defined(LIBDIFFPY_KERNEL_NEON)
    if (cpu.neon64) rv.push_back(&kernels::neonCalculationKernel());
#endif
    rv.push_back(&kernels::scalarCalculationKernel());
    return rv;
  }();
  return available;
}

/// The active kernel, the most preferred one is selected on first use.
std::atomic<const CalculationKernel*>& activeKernelSlot() {
  static std::atomic<const CalculationKernel*> slot{availableKernels().front()};
  return slot;
}

}  // namespace

// Public functions ----------------------------------------------------------

std::vector<std::string> availableCalculationKernels() {
  std::vector<std::string> rv;
  for (const CalculationKernel* k : availableKernels())
    rv.emplace_back(k->name);
  return rv;
}

std::string getCalculationKernel() {
  return kernels::activeCalculationKernel().name;
}

void setCalculationKernel(const std::string& name) {
  const auto& available = availableKernels();
  const auto it = std::ranges::find_if(
    available, [&](const CalculationKernel* k) { return name == k->name; });
  if (it == available.end()) {
    std::string names;
    for (const CalculationKernel* k : available) {
      if (!names.empty()) names += ", ";
      names += k->name;
    }
    throw std::invalid_argument("Calculation kernel '" + name +
                                "' is not available, use one of: " + names +
                                ".");
  }
  activeKernelSlot().store(*it, std::memory_order_relaxed);
}

namespace kernels {

const CalculationKernel& activeCalculationKernel() {
  return *activeKernelSlot().load(std::memory_order_relaxed);
}

}  // namespace kernels

}  // namespace srreal
}  // namespace diffpy
