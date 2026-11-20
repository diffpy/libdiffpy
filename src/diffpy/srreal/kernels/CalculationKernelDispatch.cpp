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
#include <xsimd/xsimd.hpp>
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
const std::vector<const CalculationKernel*>& cpuKernels() {
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

// Kernel tables have process lifetime. Availability is checked on the
// current device whenever a caller lists or explicitly selects CUDA.
const char* cudaUnavailable() {
#if defined(LIBDIFFPY_KERNEL_CUDA)
  return kernels::cudaKernelUnavailable();
#else
  return "libdiffpy was built without CUDA";
#endif
}

const std::vector<const CalculationKernel*>& cudaKernels() {
  static const std::vector<const CalculationKernel*> cuda = [] {
    std::vector<const CalculationKernel*> rv;
#if defined(LIBDIFFPY_KERNEL_CUDA)
    const CalculationKernel& cpu = *cpuKernels().front();
    static const CalculationKernel doubles = {
      "cuda",
      cpu.gaussianWindow,
      cpu.debyeSum,
      cpu.bondCandidates,
      kernels::cudaPDFPairSum,
      kernels::cudaDebyePairSum,
      kernels::makeCudaWorkspace,
      true,
    };
    static const CalculationKernel mixed = {
      "cuda-mixed",
      cpu.gaussianWindow,
      cpu.debyeSum,
      cpu.bondCandidates,
      kernels::cudaPDFPairSumMixed,
      kernels::cudaDebyePairSumMixed,
      kernels::makeCudaWorkspace,
      true,
    };
    rv = {&doubles, &mixed};
#endif
    return rv;
  }();
  return cuda;
}

bool isCudaName(const std::string& name) {
  return name == "cuda" || name == "cuda-mixed";
}

/// All kernels that setCalculationKernel takes.
std::vector<const CalculationKernel*> availableKernels() {
  std::vector<const CalculationKernel*> rv = cpuKernels();
  if (!cudaUnavailable()) {
    const auto& cuda = cudaKernels();
    rv.insert(rv.end(), cuda.begin(), cuda.end());
  }
  return rv;
}

/// The active kernel, the most preferred one of the CPU is selected on
/// first use.
std::atomic<const CalculationKernel*>& activeKernelSlot() {
  static std::atomic<const CalculationKernel*> slot{cpuKernels().front()};
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
  // No device is looked for when the name is one of the CPU.
  const bool cuda = isCudaName(name);
  if (cuda) {
    if (const char* reason = cudaUnavailable()) {
      throw std::invalid_argument("Calculation kernel '" + name +
                                  "' is not available: " + reason + ".");
    }
  }
  const auto& available = cuda ? cudaKernels() : cpuKernels();
  const auto it = std::ranges::find_if(
    available, [&](const CalculationKernel* k) { return name == k->name; });
  if (it == available.end()) {
    std::string names;
    for (const std::string& kname : availableCalculationKernels()) {
      if (!names.empty()) names += ", ";
      names += kname;
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

const CalculationKernel& pairSumKernel(bool cpuOnly) {
  const auto& kernel = activeCalculationKernel();
  return cpuOnly && kernel.onDevice ? *cpuKernels().front() : kernel;
}

}  // namespace kernels

}  // namespace srreal
}  // namespace diffpy
