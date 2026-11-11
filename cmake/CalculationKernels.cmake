# Calculation kernels.  The scalar kernel is always built.  With LIBDIFFPY_SIMD
# and xsimd, all SIMD kernels of the target architecture are built too and
# the library selects the best one the CPU supports at runtime, see
# src/diffpy/srreal/kernels/CalculationKernelDispatch.cpp.
set(LIBDIFFPY_SIMD "AUTO" CACHE STRING
  "SIMD calculation kernels: AUTO (all kernels of the target architecture \
when xsimd is found) or OFF (scalar kernel only)")
set_property(CACHE LIBDIFFPY_SIMD PROPERTY STRINGS AUTO OFF)

# Target architecture.  A universal macOS build matches no architecture
# below and gets only the scalar kernel.
set(_kernel_arch "${CMAKE_SYSTEM_PROCESSOR}")
if(MSVC AND CMAKE_CXX_COMPILER_ARCHITECTURE_ID)
  set(_kernel_arch "${CMAKE_CXX_COMPILER_ARCHITECTURE_ID}")
endif()
if(APPLE AND CMAKE_OSX_ARCHITECTURES)
  set(_kernel_arch "${CMAKE_OSX_ARCHITECTURES}")
endif()
string(TOLOWER "${_kernel_arch}" _kernel_arch)
if(_kernel_arch MATCHES "^(x86_64|amd64|x64)$")
  set(_kernel_arch "x86_64")
elseif(_kernel_arch MATCHES "^(aarch64|arm64)$")
  set(_kernel_arch "arm64")
else()
  set(_kernel_arch "")
endif()

# The SIMD kernels need xsimd.
set(_kernel_simd FALSE)
if(LIBDIFFPY_SIMD AND _kernel_arch)
  if(NOT TARGET xsimd)
    find_package(xsimd CONFIG QUIET)
  endif()
  if(TARGET xsimd)
    set(_kernel_simd TRUE)
  else()
    message(STATUS "xsimd not found, building only the scalar calculation "
      "kernel.  Install xsimd for the SIMD kernels.")
  endif()
endif()

# The main source glob is not recursive, add the kernels explicitly.
set(_kernel_dir "${PROJECT_SOURCE_DIR}/src/diffpy/srreal/kernels")
target_sources(diffpy PRIVATE
  "${_kernel_dir}/CalculationKernelScalar.cpp"
  "${_kernel_dir}/CalculationKernelDispatch.cpp")
set(_kernel_names scalar)

# All kernels need IEEE arithmetic to give consistent results.
if(MSVC)
  set(_kernel_math_options /fp:precise)
else()
  set(_kernel_math_options -fno-fast-math -ffp-contract=off)
endif()
set_source_files_properties("${_kernel_dir}/CalculationKernelScalar.cpp"
  PROPERTIES COMPILE_OPTIONS "${_kernel_math_options}")

# Inline functions compiled for the SIMD instruction sets must not be shared
# with the rest of the library, see CalculationKernelSIMD.ipp.  -O2 inlines
# them also in unoptimized builds and -fno-lto keeps LTO from sharing them.
set(_simd_isolation_options)
if(_kernel_simd AND NOT MSVC)
  include(CheckCXXCompilerFlag)
  check_cxx_compiler_flag(-fno-lto LIBDIFFPY_HAS_NO_LTO)
  if(LIBDIFFPY_HAS_NO_LTO)
    list(APPEND _simd_isolation_options -fno-lto)
  endif()
  list(APPEND _simd_isolation_options
    "$<$<NOT:$<CONFIG:Release,RelWithDebInfo>>:-O2>")
endif()

# Add the SIMD kernel in CalculationKernel<name>.cpp, compiled with the
# instruction set flags in ARGN.  CalculationKernelDispatch.cpp uses the
# kernel when LIBDIFFPY_KERNEL_<name> is defined.
function(_libdiffpy_add_simd_kernel name)
  set(_source "${_kernel_dir}/CalculationKernel${name}.cpp")
  set(_options ${ARGN} ${_kernel_math_options} ${_simd_isolation_options})
  set_source_files_properties("${_source}"
    PROPERTIES COMPILE_OPTIONS "${_options}")
  target_sources(diffpy PRIVATE "${_source}")
  target_compile_definitions(diffpy PRIVATE LIBDIFFPY_KERNEL_${name}=1)
  string(TOLOWER "${name}" _lower)
  set(_kernel_names ${_kernel_names} ${_lower} PARENT_SCOPE)
endfunction()

if(_kernel_simd AND _kernel_arch STREQUAL "x86_64")
  _libdiffpy_add_simd_kernel(SSE2)
  if(MSVC)
    _libdiffpy_add_simd_kernel(AVX2 /arch:AVX2)
    _libdiffpy_add_simd_kernel(AVX512 /arch:AVX512)
  else()
    _libdiffpy_add_simd_kernel(AVX2 -mavx2 -mfma)
    _libdiffpy_add_simd_kernel(AVX512
      -mavx512f -mavx512cd -mavx512dq -mfma)
  endif()
elseif(_kernel_simd AND _kernel_arch STREQUAL "arm64")
  _libdiffpy_add_simd_kernel(NEON)
endif()
list(JOIN _kernel_names ", " _kernel_names)
message(STATUS "Calculation kernels: ${_kernel_names}")

if(_kernel_simd)
  target_compile_definitions(diffpy PRIVATE LIBDIFFPY_SIMD_KERNELS=1)
  # Include path only, keeps xsimd out of the static library's link interface.
  target_include_directories(diffpy SYSTEM PRIVATE
    $<TARGET_PROPERTY:xsimd,INTERFACE_INCLUDE_DIRECTORIES>)
endif()
