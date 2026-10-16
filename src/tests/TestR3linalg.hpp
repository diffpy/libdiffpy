/*****************************************************************************
*
* libdiffpy         by DANSE Diffraction group
*                   Simon J. L. Billinge
*                   (c) 2009 The Trustees of Columbia University
*                   in the City of New York.  All rights reserved.
*
* File coded by:    Pavol Juhas
*
* See AUTHORS.txt for a list of people who contributed.
* See LICENSE_DANSE.txt for license information.
*
******************************************************************************
*
* class TestR3linalg -- unit test suite for the R3linalg module
*
*****************************************************************************/

#include <cxxtest/TestSuite.h>

#include <cmath>
#include <type_traits>

#include <diffpy/srreal/R3linalg.hpp>

using namespace std;
using namespace diffpy::srreal;
using diffpy::mathutils::EpsilonEqual;


class TestR3linalg : public CxxTest::TestSuite
{

private:

    double precision;
    EpsilonEqual allclose;

public:

    void setUp()
    {
        precision = 1.0e-12;
        allclose = EpsilonEqual(precision);
    }


    void test_determinant()
    {
        // default lattice should be cartesian
        R3::Matrix A1(
                9, 1, 9,
                6, 7, 4,
                0, 2, 9);
        R3::Matrix A2(
                9, 1, 9,
                0, 2, 9,
                6, 7, 4);
        double detA = 549;
        TS_ASSERT_EQUALS(detA, R3::determinant(A1));
        TS_ASSERT_EQUALS(-detA, R3::determinant(A2));
    }


    void test_inverse()
    {
        R3::Matrix A(
            0.5359, -0.5904, 0.8670,
           -0.0053,  0.7559, 0.2692,
           -0.8926,  0.9424, 0.9692);
        R3::Matrix invA(
            0.49063197005867, 1.42323870111089, -0.83420736316541,
           -0.24089965988852, 1.32489393619466, -0.15249839300481,
            0.68609361943181, 0.02249568627913,  0.41178393851247);
        TS_ASSERT(allclose(invA, R3::inverse(A)));
    }


    void test_eigen_solve_3x3()
    {
        R3::Matrix A(
                2.0, 1.0, 0.0,
                1.0, 2.0, 0.0,
                0.0, 0.0, 4.0);
        R3::Vector eigenvalues;
        R3::Matrix eigenvectors;
        R3::eigen_solve_3x3(A, eigenvalues, eigenvectors);

        TS_ASSERT_DELTA(1.0, eigenvalues[0], precision);
        TS_ASSERT_DELTA(3.0, eigenvalues[1], precision);
        TS_ASSERT_DELTA(4.0, eigenvalues[2], precision);
        for (int j = 0; j < R3::Ndim; ++j)
        {
            R3::Vector v(
                    eigenvectors(0, j),
                    eigenvectors(1, j),
                    eigenvectors(2, j));
            R3::Vector Av = R3::mxvecproduct(A, v);
            TS_ASSERT(allclose(Av, eigenvalues[j] * v));
            TS_ASSERT_DELTA(1.0, R3::norm(v), precision);
        }
    }


    void test_transpose()
    {
        R3::Matrix A(
            0.5359, -0.5904, 0.8670,
           -0.0053,  0.7559, 0.2692,
           -0.8926,  0.9424, 0.9692);
        R3::Matrix Atrans(
            0.5359, -0.0053, -0.8926,
           -0.5904, 0.7559, 0.9424,
            0.8670, 0.2692, 0.9692);
        TS_ASSERT_EQUALS(Atrans, R3::trans(A));
    }


    void test_norm()
    {
        R3::Vector v1(3.0, 4.0, 0.0);
        R3::Vector v2(0.67538115798129, 0.72108424545413, 0.15458914063315);
        TS_ASSERT_DELTA(5.0, R3::norm(v1), precision);
        TS_ASSERT_DELTA(1.0, R3::norm(v2), precision);
    }


    void test_dot()
    {
        R3::Vector v1(-0.97157650177843, 0.43206192654604, 0.56318686427062);
        R3::Vector v2(-0.04787719419083, 0.55895824010234, -0.34472910285751);
        double dot_v1v2 = 0.09387402846316;
        TS_ASSERT_DELTA(dot_v1v2, R3::dot(v1, v2), precision);
    }


    void test_cross()
    {
        R3::Vector v1(
                -0.55160549932839, -0.58291452407504,  0.12378162306543);
        R3::Vector v2(
                0.60842511285200, -0.97946444006248, -0.02828214306095);
        R3::Vector v1xv2(
                0.13772577008800,  0.05971126233738,  0.89493780662849);
        TS_ASSERT(allclose(v1xv2, R3::cross(v1, v2)));
    }


    void test_product()
    {
        R3::Matrix M(
                0.459631856585519, 0.726448904209060, 0.085844209317482,
                0.806838095807669, 0.240116998848762, 0.305032463662873,
                0.019487235483683, 0.580605953831255, 0.726077578738676);
        R3::Vector v(
                0.608652521912322, 0.519716469261062, 0.842577887601566);
        R3::Vector prod_Mv(
                0.729633980805669, 0.872890409522479, 0.925388343907868);
        R3::Vector prod_vM(
                0.715502648789536, 1.056153454546559, 0.822556602046019);
        TS_ASSERT(allclose(prod_Mv, R3::mxvecproduct(M, v)));
        TS_ASSERT(allclose(prod_vM, R3::mxvecproduct(v, M)));
    }


    void test_result_objects_are_not_shared()
    {
        R3::Matrix M(
                1.0, 2.0, 3.0,
                4.0, 5.0, 6.0,
                7.0, 8.0, 10.0);
        R3::Vector v(1.0, 2.0, 3.0);

        TS_ASSERT(!std::is_reference<decltype(R3::inverse(M))>::value);
        TS_ASSERT(!std::is_reference<decltype(R3::floor(v))>::value);
        TS_ASSERT(!std::is_reference<
                decltype(R3::mxvecproduct(M, v))>::value);
        TS_ASSERT(!std::is_reference<
                decltype(R3::mxvecproduct(v, M))>::value);

        R3::Vector first = R3::mxvecproduct(M, v);
        R3::Vector second = R3::mxvecproduct(v, M);
        TS_ASSERT_EQUALS(14.0, first[0]);
        TS_ASSERT_EQUALS(32.0, first[1]);
        TS_ASSERT_EQUALS(53.0, first[2]);
        TS_ASSERT_EQUALS(30.0, second[0]);
        TS_ASSERT_EQUALS(36.0, second[1]);
        TS_ASSERT_EQUALS(45.0, second[2]);
    }


    void test_helpers_are_thread_safe()
    {
        int failures = 0;
#ifdef _OPENMP
#pragma omp parallel for schedule(static) reduction(+:failures)
#endif
        for (int i = 0; i < 20000; ++i)
        {
            const double a = 2.0 + (i % 7);
            const double b = 11.0 + (i % 11);
            const double c = 23.0 + (i % 13);
            R3::Matrix M(
                    a, 0.0, 0.0,
                    0.0, b, 0.0,
                    0.0, 0.0, c);
            const double x = i * 0.001;
            R3::Vector v(x, x + 0.25, x - 0.5);

            R3::Vector mv = R3::mxvecproduct(M, v);
            R3::Vector vm = R3::mxvecproduct(v, M);
            R3::Vector fv = R3::floor(v);
            R3::Matrix invM = R3::inverse(M);
            const double d = R3::distance(v, fv);
            const double expected_d = std::sqrt(
                    std::pow(v[0] - std::floor(v[0]), 2) +
                    std::pow(v[1] - std::floor(v[1]), 2) +
                    std::pow(v[2] - std::floor(v[2]), 2));

            const bool ok =
                    std::abs(mv[0] - a * v[0]) < precision &&
                    std::abs(mv[1] - b * v[1]) < precision &&
                    std::abs(mv[2] - c * v[2]) < precision &&
                    std::abs(vm[0] - a * v[0]) < precision &&
                    std::abs(vm[1] - b * v[1]) < precision &&
                    std::abs(vm[2] - c * v[2]) < precision &&
                    std::abs(fv[0] - std::floor(v[0])) < precision &&
                    std::abs(fv[1] - std::floor(v[1])) < precision &&
                    std::abs(fv[2] - std::floor(v[2])) < precision &&
                    std::abs(invM(0, 0) - 1.0 / a) < precision &&
                    std::abs(invM(1, 1) - 1.0 / b) < precision &&
                    std::abs(invM(2, 2) - 1.0 / c) < precision &&
                    std::abs(d - expected_d) < precision;
            if (!ok)  ++failures;
        }
        TS_ASSERT_EQUALS(0, failures);
    }


};  // class TestR3linalg

// End of file
