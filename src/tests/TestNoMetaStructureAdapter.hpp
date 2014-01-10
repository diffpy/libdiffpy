/*****************************************************************************
*
* diffpy.srreal     by DANSE Diffraction group
*                   Simon J. L. Billinge
*                   (c) 2010 Trustees of the Columbia University
*                   in the City of New York.  All rights reserved.
*
* File coded by:    Pavol Juhas
*
* See AUTHORS.txt for a list of people who contributed.
* See LICENSE.txt for license information.
*
******************************************************************************
*
* class TestNoMetaStructureAdapter -- unit tests for an adapter that
*     avoids metadata forwarding from the diffpy.Structure classes.
*
*****************************************************************************/

#include <typeinfo>
#include <sstream>
#include <cxxtest/TestSuite.h>

#include <diffpy/srreal/PythonStructureAdapter.hpp>
#include <diffpy/srreal/NoMetaStructureAdapter.hpp>
#include <diffpy/srreal/PDFCalculator.hpp>
#include <diffpy/serialization.hpp>
#include "python_helpers.hpp"

using namespace std;
using namespace boost;
using namespace diffpy::srreal;

//////////////////////////////////////////////////////////////////////////////
// class TestNoMetaStructureAdapter
//////////////////////////////////////////////////////////////////////////////

class TestNoMetaStructureAdapter : public CxxTest::TestSuite
{
    private:

        StructureAdapterPtr m_pswt;

    public:

        void setUp()
        {
            CxxTest::setAbortTestOnFail(true);
            if (!m_pswt.get())
            {
                python::object stru;
                stru = loadTestStructure("PbScW25TiO3.stru");
                stru.attr("pdffit")["scale"] = 2.5;
                m_pswt = createStructureAdapter(stru);
            }
            CxxTest::setAbortTestOnFail(false);
        }


        void test_PDFCalcMeta()
        {
            PDFCalculator pdfc;
            pdfc.setRstep(0.1);
            pdfc.setRmax(5.0);
            TS_ASSERT_EQUALS(1.0, pdfc.getDoubleAttr("scale"));
            pdfc.eval(m_pswt);
            TS_ASSERT_EQUALS(2.5, pdfc.getDoubleAttr("scale"));
        }


        void test_PDFCalcNoMeta()
        {
            PDFCalculator pdfc;
            pdfc.setRstep(0.1);
            pdfc.setRmax(5.0);
            TS_ASSERT_EQUALS(1.0, pdfc.getDoubleAttr("scale"));
            pdfc.eval(nometa(m_pswt));
            TS_ASSERT_EQUALS(1.0, pdfc.getDoubleAttr("scale"));
        }


        void test_NoMetaTwice()
        {
            StructureAdapterPtr adpt0 = nometa(m_pswt);
            StructureAdapterPtr adpt1 = nometa(adpt0);
            TS_ASSERT_EQUALS(adpt0, adpt1);
        }


        void test_serialization()
        {
            using namespace diffpy::serialization;
            stringstream storage(ios::in | ios::out | ios::binary);
            oarchive oa(storage, ios::binary);
            StructureAdapterPtr pswtbare0(nometa(m_pswt));
            oa << pswtbare0;
            iarchive ia(storage, ios::binary);
            StructureAdapterPtr pswtbare1;
            TS_ASSERT(!pswtbare1.get());
            ia >> pswtbare1;
            TS_ASSERT_DIFFERS(m_pswt.get(), pswtbare1.get());
            TS_ASSERT_EQUALS(56, pswtbare1->countSites());
            TS_ASSERT_EQUALS(string("Pb"), pswtbare1->siteAtomType(0));
            TS_ASSERT_EQUALS(string("Ti"), pswtbare1->siteAtomType(55));
            TS_ASSERT(typeid(NoMetaStructureAdapter) == typeid(*pswtbare1));
        }

};  // class TestNoMetaStructureAdapter

// End of file
