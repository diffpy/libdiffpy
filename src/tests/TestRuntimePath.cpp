/*****************************************************************************
 *
 * libdiffpy         Complex Modeling Initiative
 *                   (c) 2013 Brookhaven Science Associates,
 *                   Brookhaven National Laboratory.
 *                   All rights reserved.
 *
 * File coded by:    Pavol Juhas
 *
 * See AUTHORS.rst for a list of people who contributed.
 * See LICENSE.rst for license information.
 *
 ******************************************************************************
 *
 * class TestRuntimePath -- unit tests for the resolution of datafile paths
 *
 *****************************************************************************/

#include <string>
#include <cstdlib>
#include <fstream>
#include <stdexcept>
#include <cxxtest/TestSuite.h>
#include <diffpy/runtimepath.hpp>
#include "test_helpers.hpp"

using diffpy::runtimepath::datapath;

class TestRuntimePath : public CxxTest::TestSuite {
 private:
  bool mhas_diffpyruntime;
  std::string mdiffpyruntime;

  bool filereadable(std::string& fname) {
    std::ifstream fp(fname.c_str());
    return bool(fp);
  }

  void setRuntime(const char* value) {
#ifdef _WIN32
    TS_ASSERT_EQUALS(_putenv_s("DIFFPYRUNTIME", value ? value : ""), 0);
#else
    if (value) {
      TS_ASSERT_EQUALS(setenv("DIFFPYRUNTIME", value, 1), 0);
    } else {
      TS_ASSERT_EQUALS(unsetenv("DIFFPYRUNTIME"), 0);
    }
#endif
  }

 public:
  void setUp() {
    char* pe = getenv("DIFFPYRUNTIME");
    mhas_diffpyruntime = pe;
    if (mhas_diffpyruntime) mdiffpyruntime = pe;
    setRuntime(nullptr);
  }

  void tearDown() {
    setRuntime(mhas_diffpyruntime ? mdiffpyruntime.c_str() : nullptr);
  }

  void test_datapath() {
    std::string f = datapath("f0_WaasKirf.dat");
    TS_ASSERT(filereadable(f));
    std::string d = datapath("");
    TS_ASSERT_EQUALS(f, d + '/' + "f0_WaasKirf.dat");
  }

  void test_envvar_lookup() {
    using namespace std;
    setRuntime("does/not/exist");
    TS_ASSERT_THROWS(datapath("f0_WaasKirf.dat"), runtime_error);
    TS_ASSERT_THROWS(datapath("NaCl.cif"), runtime_error);
    string td = prepend_testdata_dir("");
    setRuntime(td.c_str());
    string fnacl = prepend_testdata_dir("NaCl.cif");
    TS_ASSERT_EQUALS(fnacl, datapath("NaCl.cif"));
  }
};

// End of file
