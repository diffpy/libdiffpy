#pragma once

#if defined(_WIN32)
  #if defined(BUILD_DLL)
    #define DLL_EXPORT __declspec(dllexport)
  #else
    #define DLL_EXPORT __declspec(dllimport)
  #endif
#else
  #define DLL_EXPORT
#endif
