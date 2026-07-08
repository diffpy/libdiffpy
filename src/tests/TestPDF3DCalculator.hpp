/*****************************************************************************
*
* class TestPDF3DCalculator -- unit tests for PDF3DCalculator class
*
*****************************************************************************/

#include <cxxtest/TestSuite.h>

#include <cmath>
#include <sstream>
#include <stdexcept>
#include <boost/make_shared.hpp>

#include <diffpy/srreal/AtomicStructureAdapter.hpp>
#include <diffpy/srreal/PDF3DCalculator.hpp>
#include <diffpy/serialization.hpp>

using namespace std;
using namespace diffpy::srreal;

class TestPDF3DCalculator : public CxxTest::TestSuite
{
    private:

        boost::shared_ptr<PDF3DCalculator> mpdfc;
        AtomicStructureAdapterPtr mstru2;
        double meps;

    public:

        void setUp()
        {
            meps = diffpy::mathutils::SQRT_DOUBLE_EPS;
            mpdfc.reset(new PDF3DCalculator);

            mstru2 = boost::make_shared<AtomicStructureAdapter>();
            Atom ai;
            ai.atomtype = "Ni";
            ai.xyz_cartn[0] = 0.0;
            ai.xyz_cartn[1] = 0.0;
            ai.xyz_cartn[2] = 0.0;
            ai.uij_cartn = R3::identity();
            ai.uij_cartn(0, 0) = ai.uij_cartn(1, 1) =
                ai.uij_cartn(2, 2) = 0.004;
            mstru2->append(ai);
            ai.xyz_cartn[0] = 1.0;
            mstru2->append(ai);
        }


        void test_defaults_and_attributes()
        {
            TS_ASSERT_DELTA(0.1, mpdfc->getGridStep(), meps);
            TS_ASSERT_EQUALS(32, mpdfc->getAccumBlockSize());
            TS_ASSERT_EQUALS(true, mpdfc->getApplyRho0Background3D());
            TS_ASSERT_EQUALS(true, mpdfc->getUseCQWindow3D());
            TS_ASSERT_EQUALS(false, mpdfc->getEnableNNDelta3D());
            TS_ASSERT_EQUALS(0, mpdfc->getCalculationMode3D());
            TS_ASSERT_EQUALS(0, mpdfc->getHistogramWeightMode3D());

            TS_ASSERT_DELTA(0.0,
                    mpdfc->getDoubleAttr("enable_nn_delta3d"), meps);
            TS_ASSERT_DELTA(0.0, mpdfc->getDoubleAttr("nn_delta3d"), meps);
            TS_ASSERT_DELTA(0.9,
                    mpdfc->getDoubleAttr("nn_delta_positive_eta3d"), meps);
            TS_ASSERT_DELTA(0.0, mpdfc->getDoubleAttr("delta1_3d"), meps);
            TS_ASSERT_DELTA(0.0, mpdfc->getDoubleAttr("delta2_3d"), meps);
            TS_ASSERT_DELTA(1.0,
                    mpdfc->getDoubleAttr("delta_shell_index3d"), meps);
            TS_ASSERT_DELTA(1.0, mpdfc->getDoubleAttr("adp_scale3d"), meps);

            mpdfc->setDoubleAttr("calculation_mode3d", 1.0);
            mpdfc->setDoubleAttr("histogram_weight_mode3d", 1.0);
            mpdfc->setDoubleAttr("enable_nn_delta3d", 1.0);
            TS_ASSERT_EQUALS(1, mpdfc->getCalculationMode3D());
            TS_ASSERT_EQUALS(1, mpdfc->getHistogramWeightMode3D());
            TS_ASSERT_EQUALS(true, mpdfc->getEnableNNDelta3D());
        }


        void test_invalid_configuration()
        {
            TS_ASSERT_THROWS(mpdfc->setGridStep(0.0), invalid_argument);
            TS_ASSERT_THROWS(mpdfc->setAccumBlockSize(0), invalid_argument);
            TS_ASSERT_THROWS(mpdfc->setCalculationMode3D(2), invalid_argument);
            TS_ASSERT_THROWS(mpdfc->setHistogramWeightMode3D(2),
                    invalid_argument);
            TS_ASSERT_THROWS(
                    mpdfc->setDoubleAttr("calculation_mode3d", 0.5),
                    invalid_argument);
            TS_ASSERT_THROWS(mpdfc->setADPScale3D(0.0), invalid_argument);
            TS_ASSERT_THROWS(mpdfc->setRho0BackgroundScale3D(-1.0),
                    invalid_argument);
        }


        void test_vector_histogram_eval()
        {
            mpdfc->setRmax(2.0);
            mpdfc->setGridStep(0.5);
            mpdfc->setCalculationMode3D(1);
            mpdfc->setHistogramWeightMode3D(1);
            mpdfc->eval(mstru2);

            QuantityType pdf3d = mpdfc->get3DPDF();
            TS_ASSERT(!pdf3d.empty());
            TS_ASSERT_EQUALS(0u, pdf3d.size() % 4);

            QuantityType rh = mpdfc->getRadialHistogram3D();
            size_t nr = static_cast<size_t>(
                    std::ceil(mpdfc->getRmax() / mpdfc->getGridStep())) + 1;
            TS_ASSERT_EQUALS(2 * nr, rh.size());

            bool found_distance_one = false;
            for (size_t i = 0; i < rh.size(); i += 2)
            {
                if (std::fabs(rh[i] - 1.0) < meps && rh[i + 1] != 0.0)
                {
                    found_distance_one = true;
                    break;
                }
            }
            TS_ASSERT(found_distance_one);
        }


        void test_serialization()
        {
            mpdfc->setRmax(1.0);
            mpdfc->setGridStep(0.25);
            mpdfc->setAccumBlockSize(7);
            mpdfc->setApplyRho0Background3D(false);
            mpdfc->setUseCQWindow3D(false);
            mpdfc->setCalculationMode3D(1);
            mpdfc->setHistogramWeightMode3D(1);
            mpdfc->setEnableNNDelta3D(true);
            mpdfc->setNNDelta3D(0.02);
            mpdfc->setNNDeltaPositiveEta3D(0.8);
            mpdfc->setDistanceDelta1_3D(0.03);
            mpdfc->setDistanceDelta2_3D(0.04);
            mpdfc->setDeltaPairTypes3D("Ni", "Ni");
            mpdfc->setDeltaShellIndex3D(2);
            mpdfc->setDeltaShellTolerance3D(0.06);
            mpdfc->setDeltaKeyTolerance3D(1.0e-5);
            mpdfc->setUseADPScaleSensitivity3D(true);
            mpdfc->setADPScale3D(1.2);
            mpdfc->setRho0BackgroundScale3D(0.5);

            stringstream storage(ios::in | ios::out | ios::binary);
            diffpy::serialization::oarchive oa(storage, ios::binary);
            oa << mpdfc;
            diffpy::serialization::iarchive ia(storage, ios::binary);
            boost::shared_ptr<PDF3DCalculator> pdfc1;
            ia >> pdfc1;

            TS_ASSERT_DIFFERS(pdfc1.get(), mpdfc.get());
            TS_ASSERT_DELTA(0.25, pdfc1->getGridStep(), meps);
            TS_ASSERT_EQUALS(7, pdfc1->getAccumBlockSize());
            TS_ASSERT_EQUALS(false, pdfc1->getApplyRho0Background3D());
            TS_ASSERT_EQUALS(false, pdfc1->getUseCQWindow3D());
            TS_ASSERT_EQUALS(1, pdfc1->getCalculationMode3D());
            TS_ASSERT_EQUALS(1, pdfc1->getHistogramWeightMode3D());
            TS_ASSERT_EQUALS(true, pdfc1->getEnableNNDelta3D());
            TS_ASSERT_DELTA(0.02, pdfc1->getNNDelta3D(), meps);
            TS_ASSERT_DELTA(0.8, pdfc1->getNNDeltaPositiveEta3D(), meps);
            TS_ASSERT_DELTA(0.03, pdfc1->getDistanceDelta1_3D(), meps);
            TS_ASSERT_DELTA(0.04, pdfc1->getDistanceDelta2_3D(), meps);
            TS_ASSERT_EQUALS(string("Ni"), pdfc1->getDeltaPairA3D());
            TS_ASSERT_EQUALS(string("Ni"), pdfc1->getDeltaPairB3D());
            TS_ASSERT_EQUALS(2, pdfc1->getDeltaShellIndex3D());
            TS_ASSERT_DELTA(0.06, pdfc1->getDeltaShellTolerance3D(), meps);
            TS_ASSERT_DELTA(1.0e-5, pdfc1->getDeltaKeyTolerance3D(), meps);
            TS_ASSERT_EQUALS(true, pdfc1->getUseADPScaleSensitivity3D());
            TS_ASSERT_DELTA(1.2, pdfc1->getADPScale3D(), meps);
            TS_ASSERT_DELTA(0.5, pdfc1->getRho0BackgroundScale3D(), meps);
        }

};  // class TestPDF3DCalculator

// End of file
