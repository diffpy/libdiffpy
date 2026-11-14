/*****************************************************************************
 *
 * libdiffpy         by DANSE Diffraction group
 *                   Simon J. L. Billinge
 *                   (c) 2009 The Trustees of Columbia University
 *                   in the City of New York.  All rights reserved.
 *
 * File coded by:    Pavol Juhas
 *
 * See AUTHORS.rst for a list of people who contributed.
 * See LICENSE_DANSE.rst for license information.
 *
 ******************************************************************************
 *
 * TestR3linalg -- unit test suite for the R3linalg module
 *
 *****************************************************************************/

#include <bit>
#include <sstream>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>

#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <diffpy/srreal/R3linalg.hpp>

using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

using namespace std;
using namespace diffpy::srreal;
using diffpy::mathutils::EpsilonEqual;

class TestR3linalgFixture {
 protected:
  double precision;
  EpsilonEqual allclose;

 public:
  TestR3linalgFixture() {
    precision = 1.0e-12;
    allclose = EpsilonEqual(precision);
  }
};  // class TestR3linalgFixture

TEST_CASE_METHOD(TestR3linalgFixture, "TestR3linalg::test_determinant",
                 "[R3linalg]") {
  // default lattice should be cartesian
  R3::Matrix A1(9, 1, 9, 6, 7, 4, 0, 2, 9);
  R3::Matrix A2(9, 1, 9, 0, 2, 9, 6, 7, 4);
  double detA = 549;
  CHECK(detA == R3::determinant(A1));
  CHECK(-detA == R3::determinant(A2));
  // singular matrices with dependent rows and with a zero column
  R3::Matrix S1(1, 2, 3, 2, 4, 6, 1, 0, 1);
  R3::Matrix S2(0, 1, 2, 0, 3, 4, 0, 5, 6);
  CHECK(0.0 == R3::determinant(S1));
  CHECK(0.0 == R3::determinant(S2));
}

TEST_CASE_METHOD(TestR3linalgFixture, "TestR3linalg::test_inverse",
                 "[R3linalg]") {
  R3::Matrix A(0.5359, -0.5904, 0.8670, -0.0053, 0.7559, 0.2692, -0.8926,
               0.9424, 0.9692);
  R3::Matrix invA(0.49063197005867, 1.42323870111089, -0.83420736316541,
                  -0.24089965988852, 1.32489393619466, -0.15249839300481,
                  0.68609361943181, 0.02249568627913, 0.41178393851247);
  CHECK(allclose(invA, R3::inverse(A)));
  // the argument may be the matrix returned from the previous call
  CHECK(allclose(A, R3::inverse(R3::inverse(A))));
  R3::Matrix AinvA = R3::prod(A, R3::inverse(A));
  CHECK(allclose(R3::identity(), AinvA));
  // singular matrices with dependent rows and with a zero column
  R3::Matrix S1(1, 2, 3, 2, 4, 6, 1, 0, 1);
  R3::Matrix S2(0, 1, 2, 0, 3, 4, 0, 5, 6);
  CHECK_THROWS_AS(R3::inverse(S1), invalid_argument);
  CHECK_THROWS_AS(R3::inverse(S2), invalid_argument);
}

TEST_CASE_METHOD(TestR3linalgFixture, "TestR3linalg::test_nearly_singular",
                 "[R3linalg]") {
  // Determinant d * d is lost in round-off when evaluated from cofactors.
  const double x = 1 + 1e-8;
  const double d = x - 1;
  R3::Matrix A(1, 1, 1, 1, x, 1, 1, 1, x);
  CHECK_THAT(R3::determinant(A), WithinRel(d * d, 1e-12));
  // exact inverse scaled by d
  R3::Matrix dinvA0(2 + d, -1, -1, -1, 1, 0, -1, 0, 1);
  R3::Matrix dinvA = d * R3::inverse(A);
  CHECK(allclose(dinvA0, dinvA));
}

TEST_CASE_METHOD(TestR3linalgFixture, "TestR3linalg::test_transpose",
                 "[R3linalg]") {
  R3::Matrix A(0.5359, -0.5904, 0.8670, -0.0053, 0.7559, 0.2692, -0.8926,
               0.9424, 0.9692);
  R3::Matrix Atrans(0.5359, -0.0053, -0.8926, -0.5904, 0.7559, 0.9424, 0.8670,
                    0.2692, 0.9692);
  CHECK(Atrans == R3::trans(A));
}

TEST_CASE_METHOD(TestR3linalgFixture, "TestR3linalg::test_norm", "[R3linalg]") {
  R3::Vector v1(3.0, 4.0, 0.0);
  R3::Vector v2(0.67538115798129, 0.72108424545413, 0.15458914063315);
  CHECK_THAT(5.0, WithinAbs(R3::norm(v1), precision));
  CHECK_THAT(1.0, WithinAbs(R3::norm(v2), precision));
}

TEST_CASE_METHOD(TestR3linalgFixture, "TestR3linalg::test_dot", "[R3linalg]") {
  R3::Vector v1(-0.97157650177843, 0.43206192654604, 0.56318686427062);
  R3::Vector v2(-0.04787719419083, 0.55895824010234, -0.34472910285751);
  double dot_v1v2 = 0.09387402846316;
  CHECK_THAT(dot_v1v2, WithinAbs(R3::dot(v1, v2), precision));
}

TEST_CASE_METHOD(TestR3linalgFixture, "TestR3linalg::test_cross",
                 "[R3linalg]") {
  R3::Vector v1(-0.55160549932839, -0.58291452407504, 0.12378162306543);
  R3::Vector v2(0.60842511285200, -0.97946444006248, -0.02828214306095);
  R3::Vector v1xv2(0.13772577008800, 0.05971126233738, 0.89493780662849);
  CHECK(allclose(v1xv2, R3::cross(v1, v2)));
}

TEST_CASE_METHOD(TestR3linalgFixture, "TestR3linalg::test_product",
                 "[R3linalg]") {
  R3::Matrix M(0.459631856585519, 0.726448904209060, 0.085844209317482,
               0.806838095807669, 0.240116998848762, 0.305032463662873,
               0.019487235483683, 0.580605953831255, 0.726077578738676);
  R3::Vector v(0.608652521912322, 0.519716469261062, 0.842577887601566);
  R3::Vector prod_Mv(0.729633980805669, 0.872890409522479, 0.925388343907868);
  R3::Vector prod_vM(0.715502648789536, 1.056153454546559, 0.822556602046019);
  CHECK(allclose(prod_Mv, R3::mxvecproduct(M, v)));
  CHECK(allclose(prod_vM, R3::mxvecproduct(v, M)));
}

TEST_CASE("TestR3linalg::test_default_zero", "[R3linalg]") {
  R3::Vector v;
  R3::Matrix M;
  CHECK(R3::Vector(0, 0, 0) == v);
  CHECK(R3::zeromatrix() == M);
  CHECK(R3::Vector(0, 0, 0) == R3::zerovector);
}

TEST_CASE("TestR3linalg::test_arithmetic", "[R3linalg]") {
  R3::Vector u(1, 2, 3);
  R3::Vector v(4, 5, 7);
  CHECK(R3::Vector(5, 7, 10) == u + v);
  CHECK(R3::Vector(-3, -3, -4) == u - v);
  CHECK(R3::Vector(-1, -2, -3) == -u);
  CHECK(R3::Vector(2, 4, 6) == u * 2.0);
  CHECK(R3::Vector(2, 4, 6) == 2.0 * u);
  CHECK(R3::Vector(0.5, 1, 1.5) == u / 2.0);
  u += v;
  CHECK(R3::Vector(5, 7, 10) == u);
  u -= v;
  CHECK(R3::Vector(1, 2, 3) == u);
  R3::Matrix A(1, 2, 3, 4, 5, 6, 7, 8, 9);
  CHECK(R3::Matrix(2, 4, 6, 8, 10, 12, 14, 16, 18) == A + A);
  CHECK(R3::zeromatrix() == A - A);
  CHECK(A * 2.0 == 2.0 * A);
  CHECK(A == (A * 4.0) / 4.0);
  CHECK(R3::zeromatrix() == A + (-A));
}

TEST_CASE("TestR3linalg::test_rows_columns", "[R3linalg]") {
  R3::Matrix A(1, 2, 3, 4, 5, 6, 7, 8, 9);
  CHECK(R3::Vector(4, 5, 6) == R3::row(A, 1));
  CHECK(R3::Vector(3, 6, 9) == R3::column(A, 2));
  R3::setrow(A, 0, R3::Vector(-1, -2, -3));
  CHECK(R3::Matrix(-1, -2, -3, 4, 5, 6, 7, 8, 9) == A);
  R3::setcolumn(A, 1, R3::Vector(0, 0, 0));
  CHECK(R3::Matrix(-1, 0, -3, 4, 0, 6, 7, 0, 9) == A);
}

namespace {

// T&& is an rvalue when T is not a reference type
template <class T>
concept ScalableInPlace = requires(T&& x) { std::forward<T>(x) *= 2.0; };

template <class T>
concept VectorElementWritable = requires(T&& x) {
  std::forward<T>(x)[1] = 2.0;
  std::forward<T>(x)(1) += 2.0;
  *std::forward<T>(x).begin() = 2.0;
  std::forward<T>(x).data()[1] = 2.0;
};

template <class T>
concept MatrixElementWritable = requires(T&& x) {
  std::forward<T>(x)(0, 1) = 2.0;
  std::forward<T>(x).data()[1] = 2.0;
};

}  // namespace

// Writes to temporaries such as row(A, 0) = v or row(A, 0)[1] = x must not
// compile, because they would silently leave A unchanged.
static_assert(std::is_assignable_v<R3::Vector&, const R3::Vector&>);
static_assert(!std::is_assignable_v<R3::Vector, const R3::Vector&>);
static_assert(ScalableInPlace<R3::Vector&>);
static_assert(!ScalableInPlace<R3::Vector>);
static_assert(VectorElementWritable<R3::Vector&>);
static_assert(!VectorElementWritable<R3::Vector>);
static_assert(std::is_assignable_v<R3::Matrix&, const R3::Matrix&>);
static_assert(!std::is_assignable_v<R3::Matrix, const R3::Matrix&>);
static_assert(ScalableInPlace<R3::Matrix&>);
static_assert(!ScalableInPlace<R3::Matrix>);
static_assert(MatrixElementWritable<R3::Matrix&>);
static_assert(!MatrixElementWritable<R3::Matrix>);
static_assert(std::is_trivially_copyable_v<R3::Vector>);
static_assert(std::is_trivially_copyable_v<R3::Matrix>);

TEST_CASE("TestR3linalg::test_prod", "[R3linalg]") {
  R3::Matrix A(1, 2, 3, 4, 5, 6, 7, 8, 9);
  R3::Matrix B(9, 8, 7, 6, 5, 4, 3, 2, 1);
  R3::Matrix AB(30, 24, 18, 84, 69, 54, 138, 114, 90);
  CHECK(AB == R3::prod(A, B));
  CHECK(A == R3::prod(A, R3::identity()));
}

TEST_CASE_METHOD(TestR3linalgFixture, "TestR3linalg::test_inverse_pivoting",
                 "[R3linalg]") {
  // zero at the first pivot requires row exchanges
  R3::Matrix P(0, 1, 0, 0, 0, 1, 1, 0, 0);
  CHECK(1.0 == R3::determinant(P));
  CHECK(R3::trans(P) == R3::inverse(P));
  R3::Matrix A(0, 2, 1, 3, 0, 4, 5, 6, 0);
  CHECK_THAT(R3::determinant(A), WithinAbs(58.0, precision));
  CHECK(allclose(R3::identity(), R3::prod(A, R3::inverse(A))));
}

namespace {

// Binary archive of a Vector(1, -2.5, 3.25), Matrix(1, 2, ..., 8, 9.5) and
// the same two objects again, as written by the boost::numeric::ublas based
// implementation.  The archive uses no_header, little-endian byte order and
// 8-byte collection sizes.
const char* const LEGACY_ARCHIVE_HEX =
  "0000000000000000000000000000000300000000000000000000000000f03f00"
  "000000000004c00000000000000a400000000000000000000003000000000000"
  "00030000000000000000000000000900000000000000000000000000f03f0000"
  "0000000000400000000000000840000000000000104000000000000014400000"
  "0000000018400000000000001c40000000000000204000000000000023400300"
  "000000000000000000000000f03f00000000000004c00000000000000a400300"
  "00000000000003000000000000000900000000000000000000000000f03f0000"
  "0000000000400000000000000840000000000000104000000000000014400000"
  "0000000018400000000000001c4000000000000020400000000000002340";

std::string fromhex(const std::string& hex) {
  std::string rv;
  for (size_t i = 0; i < hex.size(); i += 2) {
    rv += char(std::stoi(hex.substr(i, 2), nullptr, 16));
  }
  return rv;
}

}  // namespace

TEST_CASE("TestR3linalg::test_legacy_serialization", "[R3linalg]") {
  if (std::endian::native != std::endian::little || sizeof(size_t) != 8) {
    SKIP("legacy archive bytes assume little-endian 64-bit platform");
  }
  const R3::Vector v(1.0, -2.5, 3.25);
  const R3::Matrix M(1, 2, 3, 4, 5, 6, 7, 8, 9.5);
  using boost::archive::no_header;
  const std::string legacy = fromhex(LEGACY_ARCHIVE_HEX);
  // saving produces the same bytes as the ublas based implementation
  std::ostringstream os;
  {
    boost::archive::binary_oarchive oa(os, no_header);
    oa << v << M << v << M;
  }
  CHECK(legacy == os.str());
  // loading of the legacy archive
  std::istringstream is(legacy);
  {
    boost::archive::binary_iarchive ia(is, no_header);
    R3::Vector v1, v2;
    R3::Matrix M1, M2;
    ia >> v1 >> M1 >> v2 >> M2;
    CHECK(v == v1);
    CHECK(M == M1);
    CHECK(v == v2);
    CHECK(M == M2);
  }
  // archive with a wrong Vector size cannot be loaded
  std::string bad = legacy;
  REQUIRE(3 == bad[15]);
  bad[15] = 4;
  std::istringstream isbad(bad);
  boost::archive::binary_iarchive iabad(isbad, no_header);
  R3::Vector vbad;
  CHECK_THROWS_AS(iabad >> vbad, boost::archive::archive_exception);
}

// End of file
