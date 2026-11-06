#pragma once

#if defined(_WIN32) && !defined(DIFFPY_STATIC)
#if defined(BUILD_DLL)
#define DLL_EXPORT __declspec(dllexport)
#else
#define DLL_EXPORT __declspec(dllimport)
#endif
#else
#define DLL_EXPORT
#endif
