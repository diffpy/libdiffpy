# CUDA calculation kernels, see
# src/diffpy/srreal/kernels/CalculationKernelDispatch.cpp.
#
# The CUDA source makes the static library diffpy_cuda, which the diffpy
# library links.  It is compiled by nvcc and statically linked with its
# own C++ standard and without the options of the C++ sources.

include(CheckLanguage)
check_language(CUDA)
if(NOT CMAKE_CUDA_COMPILER)
  message(FATAL_ERROR
    "LIBDIFFPY_CUDA needs the CUDA compiler nvcc.  Set CMAKE_CUDA_COMPILER "
    "or the environment variable CUDACXX to it.")
endif()
enable_language(CUDA)

add_library(diffpy_cuda STATIC
  "${PROJECT_SOURCE_DIR}/src/diffpy/srreal/kernels/CalculationKernelCuda.cu")
set_target_properties(diffpy_cuda PROPERTIES
  CUDA_STANDARD 17
  CUDA_STANDARD_REQUIRED ON
  CUDA_EXTENSIONS OFF
  CUDA_RUNTIME_LIBRARY Static
  POSITION_INDEPENDENT_CODE ON
  # The objects of nvcc have no data for link-time optimization.
  INTERPROCEDURAL_OPTIMIZATION OFF)
target_include_directories(diffpy_cuda PRIVATE "${PROJECT_SOURCE_DIR}/src")
# Use the CPU kernels' rounding order for distances, widths and grid cutoffs.
target_compile_options(diffpy_cuda PRIVATE
  $<$<CUDA_COMPILER_ID:NVIDIA>:--fmad=false>
  $<$<CUDA_COMPILER_ID:Clang>:-ffp-contract=off>)
# The runtime is named here for the programs that link diffpy_cuda through
# a static or OBJECT diffpy library, whose link step is not one of CUDA.
find_package(CUDAToolkit REQUIRED)
target_link_libraries(diffpy_cuda PRIVATE CUDA::cudart_static)
target_link_libraries(diffpy PRIVATE diffpy_cuda)
target_compile_definitions(diffpy PRIVATE LIBDIFFPY_KERNEL_CUDA=1)
