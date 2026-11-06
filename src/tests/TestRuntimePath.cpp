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
 * TestRuntimePath -- unit tests for the resolution of datafile paths
 *
 *****************************************************************************/

#include <string>
#include <cstdlib>
#include <fstream>
#include <stdexcept>
#include <catch2/catch_test_macros.hpp>
#include <diffpy/runtimepath.hpp>
#include "test_helpers.hpp"

using diffpy::runtimepath::datapath;

class TestRuntimePathFixture {
 protected:
  bool mhas_diffpyruntime;
  std::string mdiffpyruntime;

  bool filereadable(std::string& fname) {
    std::ifstream fp(fname.c_str());
    return bool(fp);
  }

  void setRuntime(const char* value) {
#ifdef _WIN32
    CHECK(_putenv_s("DIFFPYRUNTIME", value ? value : "") == 0);
#else
    if (value) {
      CHECK(setenv("DIFFPYRUNTIME", value, 1) == 0);
    } else {
      CHECK(unsetenv("DIFFPYRUNTIME") == 0);
    }
#endif
  }

 public:
  TestRuntimePathFixture() {
    char* pe = getenv("DIFFPYRUNTIME");
    mhas_diffpyruntime = pe;
    if (mhas_diffpyruntime) mdiffpyruntime = pe;
    setRuntime(nullptr);
  }

  ~TestRuntimePathFixture() {
    setRuntime(mhas_diffpyruntime ? mdiffpyruntime.c_str() : nullptr);
  }
};  // class TestRuntimePathFixture

TEST_CASE_METHOD(TestRuntimePathFixture, "TestRuntimePath::test_datapath",
                 "[RuntimePath]") {
  std::string f = datapath("f0_WaasKirf.dat");
  CHECK(filereadable(f));
  std::string d = datapath("");
  CHECK(f == d + '/' + "f0_WaasKirf.dat");
}

TEST_CASE_METHOD(TestRuntimePathFixture, "TestRuntimePath::test_envvar_lookup",
                 "[RuntimePath]") {
  using namespace std;
  setRuntime("does/not/exist");
  CHECK_THROWS_AS(datapath("f0_WaasKirf.dat"), runtime_error);
  CHECK_THROWS_AS(datapath("NaCl.cif"), runtime_error);
  string td = prepend_testdata_dir("");
  setRuntime(td.c_str());
  string fnacl = prepend_testdata_dir("NaCl.cif");
  CHECK(fnacl == datapath("NaCl.cif"));
}

// End of file
