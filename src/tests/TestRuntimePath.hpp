/*****************************************************************************
*
* libdiffpy         Complex Modeling Initiative
*                   Pavol Juhas
*                   (c) 2013 Brookhaven National Laboratory,
*                   Upton, New York.  All rights reserved.
*
* File coded by:    Pavol Juhas
*
* See AUTHORS.txt for a list of people who contributed.
* See LICENSE.txt for license information.
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


class TestRuntimePath : public CxxTest::TestSuite
{

    private:

        bool mhas_diffpyruntime;
        std::string mdiffpyruntime;

        bool filereadable(std::string& fname)
        {
            std::ifstream fp(fname.c_str());
            return bool(fp);
        }

    public:

        void setUp()
        {
            char* pe = getenv("DIFFPYRUNTIME");
            mhas_diffpyruntime = pe;
            if (mhas_diffpyruntime)  mdiffpyruntime = pe;
            unsetenv("DIFFPYRUNTIME");
        }


        void tearDown()
        {
            unsetenv("DIFFPYRUNTIME");
            if (mhas_diffpyruntime)
            {
                setenv("DIFFPYRUNTIME", mdiffpyruntime.c_str(), 1);
            }
        }


        void test_datapath()
        {
            std::string f = datapath("f0_WaasKirf.dat");
            TS_ASSERT(filereadable(f));
            std::string d = datapath("");
            TS_ASSERT_EQUALS(f, d + '/' + "f0_WaasKirf.dat");
        }


        void test_envvar_lookup()
        {
            using namespace std;
            setenv("DIFFPYRUNTIME", "does/not/exist", 0);
            TS_ASSERT_THROWS(datapath("f0_WaasKirf.dat"), runtime_error);
            TS_ASSERT_THROWS(datapath("NaCl.cif"), runtime_error);
            string td = prepend_testdata_dir("");
            setenv("DIFFPYRUNTIME", td.c_str(), 1);
            string fnacl = prepend_testdata_dir("NaCl.cif");
            TS_ASSERT_EQUALS(fnacl, datapath("NaCl.cif"));
        }

};

// End of file
