/*****************************************************************************
*
* class TestThreeDPDFCalculator -- unit tests for ThreeDPDFCalculator
*
*****************************************************************************/

#include <cxxtest/TestSuite.h>

#include <cmath>
#include <cstdio>
#include <fstream>
#include <limits>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>

#include <diffpy/srreal/AtomicStructureAdapter.hpp>
#include <diffpy/srreal/PeriodicStructureAdapter.hpp>
#include <diffpy/srreal/ThreeDPDFCalculator.hpp>
#include <diffpy/serialization.hpp>

using namespace std;
using namespace diffpy::srreal;

// Local Helpers -------------------------------------------------------------

namespace {

AtomicStructureAdapterPtr
makeDimer(const string& atomtype, double uiso)
{
    AtomicStructureAdapterPtr stru =
        std::make_shared<AtomicStructureAdapter>();
    Atom atom;
    atom.atomtype = atomtype;
    atom.xyz_cartn = R3::Vector(0.0, 0.0, 0.0);
    atom.uij_cartn = R3::zeromatrix();
    atom.uij_cartn(0, 0) = atom.uij_cartn(1, 1) =
        atom.uij_cartn(2, 2) = uiso;
    stru->append(atom);
    atom.xyz_cartn[0] = 1.0;
    stru->append(atom);
    return stru;
}


PeriodicStructureAdapterPtr
makePeriodicDimer()
{
    PeriodicStructureAdapterPtr stru =
        std::make_shared<PeriodicStructureAdapter>();
    stru->setLatPar(10.0, 10.0, 10.0, 90.0, 90.0, 90.0);
    AtomicStructureAdapterPtr atoms = makeDimer("Ni", 0.004);
    stru->append(atoms->at(0));
    stru->append(atoms->at(1));
    return stru;
}


void
configureSmallStandardCalculation(ThreeDPDFCalculator& calc)
{
    calc.setRmax(2.0);
    calc.setGridStep(0.25);
    calc.setUseCQWindow3D(false);
    calc.setApplyRho0Background3D(false);
}


bool
allSparseValuesFinite(const QuantityType& data)
{
    if (data.size() % 4)  return false;
    for (size_t i = 0; i < data.size(); ++i)
    {
        if (!std::isfinite(data[i]))  return false;
    }
    return true;
}


bool
findSparseValue(
        const QuantityType& data,
        double x,
        double y,
        double z,
        double eps,
        double& value)
{
    for (size_t i = 0; i + 3 < data.size(); i += 4)
    {
        if (std::fabs(data[i] - x) <= eps &&
            std::fabs(data[i + 1] - y) <= eps &&
            std::fabs(data[i + 2] - z) <= eps)
        {
            value = data[i + 3];
            return true;
        }
    }
    return false;
}


bool
quantitiesClose(
        const QuantityType& a,
        const QuantityType& b,
        double eps)
{
    if (a.size() != b.size())  return false;
    for (size_t i = 0; i < a.size(); ++i)
    {
        if (std::fabs(a[i] - b[i]) > eps)  return false;
    }
    return true;
}


size_t
binaryFileSize(const string& path)
{
    ifstream stream(path.c_str(), ios::binary | ios::ate);
    if (!stream)  return 0;
    return static_cast<size_t>(stream.tellg());
}

}   // namespace

//////////////////////////////////////////////////////////////////////////////
// class TestThreeDPDFCalculator
//////////////////////////////////////////////////////////////////////////////

class TestThreeDPDFCalculator : public CxxTest::TestSuite
{
    private:

        std::shared_ptr<ThreeDPDFCalculator> mpdfc;
        AtomicStructureAdapterPtr mstru2;
        double meps;
        const string mfloat32path = "threedpdf_test_float32.bin";
        const string mfloat64path = "threedpdf_test_float64.bin";

    public:

        void setUp()
        {
            meps = diffpy::mathutils::SQRT_DOUBLE_EPS;
            mpdfc.reset(new ThreeDPDFCalculator);
            mstru2 = makeDimer("Ni", 0.004);
            std::remove(mfloat32path.c_str());
            std::remove(mfloat64path.c_str());
        }


        void tearDown()
        {
            std::remove(mfloat32path.c_str());
            std::remove(mfloat64path.c_str());
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
            TS_ASSERT_DELTA(0.0, mpdfc->getNNDelta3D(), meps);
            TS_ASSERT_DELTA(0.9, mpdfc->getNNDeltaPositiveEta3D(), meps);
            TS_ASSERT_DELTA(0.0, mpdfc->getDistanceDelta1_3D(), meps);
            TS_ASSERT_DELTA(0.0, mpdfc->getDistanceDelta2_3D(), meps);
            TS_ASSERT_EQUALS(string("*"), mpdfc->getDeltaPairA3D());
            TS_ASSERT_EQUALS(string("*"), mpdfc->getDeltaPairB3D());
            TS_ASSERT_EQUALS(1, mpdfc->getDeltaShellIndex3D());
            TS_ASSERT_DELTA(
                    0.05, mpdfc->getDeltaShellTolerance3D(), meps);
            TS_ASSERT_DELTA(
                    1.0e-6, mpdfc->getDeltaKeyTolerance3D(), meps);
            TS_ASSERT_EQUALS(
                    false, mpdfc->getUseADPScaleSensitivity3D());
            TS_ASSERT_DELTA(1.0, mpdfc->getADPScale3D(), meps);
            TS_ASSERT_DELTA(
                    1.0, mpdfc->getRho0BackgroundScale3D(), meps);

            TS_ASSERT_DELTA(0.0,
                    mpdfc->getDoubleAttr("enable_nn_delta3d"), meps);
            TS_ASSERT_DELTA(0.0,
                    mpdfc->getDoubleAttr("nn_delta3d"), meps);
            TS_ASSERT_DELTA(0.9,
                    mpdfc->getDoubleAttr("nn_delta_positive_eta3d"), meps);
            TS_ASSERT_DELTA(0.0,
                    mpdfc->getDoubleAttr("delta1_3d"), meps);
            TS_ASSERT_DELTA(0.0,
                    mpdfc->getDoubleAttr("delta2_3d"), meps);
            TS_ASSERT_DELTA(1.0,
                    mpdfc->getDoubleAttr("delta_shell_index3d"), meps);
            TS_ASSERT_DELTA(1.0,
                    mpdfc->getDoubleAttr("adp_scale3d"), meps);

            mpdfc->setDoubleAttr("calculation_mode3d", 1.0);
            mpdfc->setDoubleAttr("histogram_weight_mode3d", 1.0);
            mpdfc->setDoubleAttr("enable_nn_delta3d", 1.0);
            mpdfc->setDoubleAttr("use_adp_scale_sensitivity3d", 1.0);
            TS_ASSERT_EQUALS(1, mpdfc->getCalculationMode3D());
            TS_ASSERT_EQUALS(1, mpdfc->getHistogramWeightMode3D());
            TS_ASSERT_EQUALS(true, mpdfc->getEnableNNDelta3D());
            TS_ASSERT_EQUALS(
                    true, mpdfc->getUseADPScaleSensitivity3D());
        }


        void test_invalid_configuration()
        {
            TS_ASSERT_THROWS(mpdfc->setGridStep(0.0), invalid_argument);
            TS_ASSERT_THROWS(mpdfc->setGridStep(-0.1), invalid_argument);
            TS_ASSERT_THROWS(
                    mpdfc->setAccumBlockSize(0), invalid_argument);
            TS_ASSERT_THROWS(
                    mpdfc->setCalculationMode3D(2), invalid_argument);
            TS_ASSERT_THROWS(
                    mpdfc->setHistogramWeightMode3D(-1),
                    invalid_argument);
            TS_ASSERT_THROWS(
                    mpdfc->setDoubleAttr("calculation_mode3d", 0.5),
                    invalid_argument);
            TS_ASSERT_THROWS(
                    mpdfc->setDoubleAttr("histogram_weight_mode3d", 0.5),
                    invalid_argument);
            TS_ASSERT_THROWS(
                    mpdfc->setNNDelta3D(-0.01), invalid_argument);
            TS_ASSERT_THROWS(
                    mpdfc->setNNDeltaPositiveEta3D(0.0),
                    invalid_argument);
            TS_ASSERT_THROWS(
                    mpdfc->setNNDeltaPositiveEta3D(1.0),
                    invalid_argument);
            TS_ASSERT_THROWS(
                    mpdfc->setDistanceDelta1_3D(-0.01),
                    invalid_argument);
            TS_ASSERT_THROWS(
                    mpdfc->setDistanceDelta2_3D(-0.01),
                    invalid_argument);
            TS_ASSERT_THROWS(
                    mpdfc->setDeltaPairTypes3D("", "Ni"),
                    invalid_argument);
            TS_ASSERT_THROWS(
                    mpdfc->setDeltaShellIndex3D(0), invalid_argument);
            TS_ASSERT_THROWS(
                    mpdfc->setDeltaShellTolerance3D(0.0),
                    invalid_argument);
            TS_ASSERT_THROWS(
                    mpdfc->setDeltaKeyTolerance3D(0.0),
                    invalid_argument);
            TS_ASSERT_THROWS(
                    mpdfc->setADPScale3D(0.0), invalid_argument);
            TS_ASSERT_THROWS(
                    mpdfc->setRho0BackgroundScale3D(-1.0),
                    invalid_argument);
        }


        void test_standard_grid_mode()
        {
            configureSmallStandardCalculation(*mpdfc);
            mpdfc->eval(mstru2);

            QuantityType result = mpdfc->getThreeDPDF();
            TS_ASSERT(!result.empty());
            TS_ASSERT_EQUALS(0u, result.size() % 4);
            TS_ASSERT(allSparseValuesFinite(result));

            double positive = 0.0;
            double negative = 0.0;
            TS_ASSERT(findSparseValue(
                    result, 1.0, 0.0, 0.0, meps, positive));
            TS_ASSERT(findSparseValue(
                    result, -1.0, 0.0, 0.0, meps, negative));
            TS_ASSERT(positive > 0.0);
            TS_ASSERT(negative > 0.0);
            TS_ASSERT_DELTA(positive, negative, meps);
        }


        void test_vector_histogram_and_pre_rename_regression()
        {
            mpdfc->setRmax(2.0);
            mpdfc->setGridStep(0.5);
            mpdfc->setCalculationMode3D(1);
            mpdfc->setHistogramWeightMode3D(1);
            mpdfc->eval(mstru2);

            QuantityType result = mpdfc->getThreeDPDF();
            TS_ASSERT_EQUALS(8u, result.size());
            double positive = 0.0;
            double negative = 0.0;
            TS_ASSERT(findSparseValue(
                    result, 1.0, 0.0, 0.0, meps, positive));
            TS_ASSERT(findSparseValue(
                    result, -1.0, 0.0, 0.0, meps, negative));
            TS_ASSERT_DELTA(1.0, positive, meps);
            TS_ASSERT_DELTA(1.0, negative, meps);

            QuantityType radial = mpdfc->getRadialHistogram3D();
            size_t nr = static_cast<size_t>(
                    std::ceil(mpdfc->getRmax() /
                        mpdfc->getGridStep())) + 1;
            TS_ASSERT_EQUALS(2 * nr, radial.size());
            TS_ASSERT_DELTA(2.0, radial[5], meps);
        }


        void test_histogram_mode_does_not_use_adps()
        {
            AtomicStructureAdapterPtr largeadp =
                makeDimer("Ni", 0.4);
            mpdfc->setRmax(2.0);
            mpdfc->setGridStep(0.5);
            mpdfc->setCalculationMode3D(1);
            mpdfc->setHistogramWeightMode3D(1);

            mpdfc->eval(mstru2);
            QuantityType small = mpdfc->getThreeDPDF();
            mpdfc->eval(largeadp);
            QuantityType large = mpdfc->getThreeDPDF();
            TS_ASSERT(quantitiesClose(small, large, 0.0));
        }


        void test_q_window_changes_grid()
        {
            configureSmallStandardCalculation(*mpdfc);
            mpdfc->eval(mstru2);
            QuantityType unwindowed = mpdfc->getThreeDPDF();

            mpdfc->setUseCQWindow3D(true);
            mpdfc->setQmin(0.0);
            mpdfc->setQmax(4.0);
            QuantityType windowed = mpdfc->getThreeDPDF();

            TS_ASSERT(!windowed.empty());
            TS_ASSERT(allSparseValuesFinite(windowed));
            TS_ASSERT(!quantitiesClose(unwindowed, windowed, 1.0e-10));
        }


        void test_rho0_background()
        {
            PeriodicStructureAdapterPtr periodic = makePeriodicDimer();
            configureSmallStandardCalculation(*mpdfc);
            mpdfc->eval(periodic);
            QuantityType raw = mpdfc->getThreeDPDF();
            TS_ASSERT(!raw.empty());

            mpdfc->setApplyRho0Background3D(true);
            mpdfc->setRho0BackgroundScale3D(0.75);
            QuantityType corrected = mpdfc->getThreeDPDF();

            double correctedvalue = 0.0;
            TS_ASSERT(findSparseValue(
                    corrected, raw[0], raw[1], raw[2],
                    meps, correctedvalue));
            const double expected =
                0.75 * periodic->numberDensity();
            TS_ASSERT_DELTA(
                    expected, raw[3] - correctedvalue, 1.0e-10);
        }


        void test_delta_pair_filter_and_shell_covariance()
        {
            configureSmallStandardCalculation(*mpdfc);
            mpdfc->setEnableNNDelta3D(true);
            mpdfc->setDeltaPairTypes3D("O", "O");
            mpdfc->eval(mstru2);
            TS_ASSERT_EQUALS(0, mpdfc->getDeltaEligiblePairCount3D());

            mpdfc->setDeltaPairTypes3D("Ni", "Ni");
            mpdfc->eval(mstru2);
            TS_ASSERT(
                    mpdfc->getDeltaEligiblePairCount3D() > 0);
            const double bound =
                mpdfc->getNNDelta3DUpperBound();
            TS_ASSERT(bound > 0.0);
            QuantityType baseline = mpdfc->getThreeDPDF();

            mpdfc->setNNDelta3D(0.5 * bound);
            mpdfc->eval(mstru2);
            QuantityType correlated = mpdfc->getThreeDPDF();
            TS_ASSERT(!quantitiesClose(
                    baseline, correlated, 1.0e-10));
            TS_ASSERT_THROWS(
                    mpdfc->setNNDelta3D(bound), invalid_argument);
        }


        void test_distance_delta_covariance()
        {
            configureSmallStandardCalculation(*mpdfc);
            mpdfc->setEnableNNDelta3D(true);
            mpdfc->setDeltaPairTypes3D("Ni", "Ni");
            mpdfc->eval(mstru2);
            QuantityType baseline = mpdfc->getThreeDPDF();

            mpdfc->setDistanceDelta1_3D(0.1);
            mpdfc->setDistanceDelta2_3D(0.05);
            mpdfc->eval(mstru2);
            QuantityType correlated = mpdfc->getThreeDPDF();
            TS_ASSERT(!quantitiesClose(
                    baseline, correlated, 1.0e-10));
        }


        void test_rdf_normalization_cancels_atom_scattering_scale()
        {
            AtomicStructureAdapterPtr carbon =
                makeDimer("C", 0.004);
            configureSmallStandardCalculation(*mpdfc);
            mpdfc->eval(mstru2);
            QuantityType nickelresult = mpdfc->getThreeDPDF();
            mpdfc->eval(carbon);
            QuantityType carbonresult = mpdfc->getThreeDPDF();

            TS_ASSERT(quantitiesClose(
                    nickelresult, carbonresult, 1.0e-10));
        }


        void test_binary_export()
        {
            mpdfc->setRmax(1.0);
            mpdfc->setGridStep(0.5);
            mpdfc->setUseCQWindow3D(false);
            mpdfc->setApplyRho0Background3D(false);
            mpdfc->eval(mstru2);

            mpdfc->exportGrid3DBinary(
                    mfloat32path, true, false);
            mpdfc->exportGrid3DBinary(
                    mfloat64path, false, false);
            const size_t gridcells = 5u * 5u * 5u;
            TS_ASSERT_EQUALS(
                    gridcells * sizeof(float),
                    binaryFileSize(mfloat32path));
            TS_ASSERT_EQUALS(
                    gridcells * sizeof(double),
                    binaryFileSize(mfloat64path));
            TS_ASSERT_THROWS(
                    mpdfc->exportGrid3DBinary(
                        "libdiffpy_missing_test_directory/grid.bin"),
                    runtime_error);
        }


        void test_serialization_round_trip()
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
            diffpy::serialization::oarchive oa(
                    storage, ios::binary);
            oa << mpdfc;
            diffpy::serialization::iarchive ia(
                    storage, ios::binary);
            std::shared_ptr<ThreeDPDFCalculator> restored;
            ia >> restored;

            TS_ASSERT_DIFFERS(restored.get(), mpdfc.get());
            TS_ASSERT_DELTA(
                    0.25, restored->getGridStep(), meps);
            TS_ASSERT_EQUALS(7, restored->getAccumBlockSize());
            TS_ASSERT_EQUALS(
                    false, restored->getApplyRho0Background3D());
            TS_ASSERT_EQUALS(
                    false, restored->getUseCQWindow3D());
            TS_ASSERT_EQUALS(
                    1, restored->getCalculationMode3D());
            TS_ASSERT_EQUALS(
                    1, restored->getHistogramWeightMode3D());
            TS_ASSERT_EQUALS(
                    true, restored->getEnableNNDelta3D());
            TS_ASSERT_DELTA(
                    0.02, restored->getNNDelta3D(), meps);
            TS_ASSERT_DELTA(
                    0.8, restored->getNNDeltaPositiveEta3D(), meps);
            TS_ASSERT_DELTA(
                    0.03, restored->getDistanceDelta1_3D(), meps);
            TS_ASSERT_DELTA(
                    0.04, restored->getDistanceDelta2_3D(), meps);
            TS_ASSERT_EQUALS(
                    string("Ni"), restored->getDeltaPairA3D());
            TS_ASSERT_EQUALS(
                    string("Ni"), restored->getDeltaPairB3D());
            TS_ASSERT_EQUALS(
                    2, restored->getDeltaShellIndex3D());
            TS_ASSERT_DELTA(
                    0.06, restored->getDeltaShellTolerance3D(), meps);
            TS_ASSERT_DELTA(
                    1.0e-5, restored->getDeltaKeyTolerance3D(), meps);
            TS_ASSERT_EQUALS(
                    true, restored->getUseADPScaleSensitivity3D());
            TS_ASSERT_DELTA(
                    1.2, restored->getADPScale3D(), meps);
            TS_ASSERT_DELTA(
                    0.5, restored->getRho0BackgroundScale3D(), meps);

            restored->eval(mstru2);
            QuantityType restoredresult =
                restored->getThreeDPDF();
            mpdfc->eval(mstru2);
            QuantityType originalresult =
                mpdfc->getThreeDPDF();
            TS_ASSERT(quantitiesClose(
                    originalresult, restoredresult, 0.0));
        }

};  // class TestThreeDPDFCalculator

// End of file
