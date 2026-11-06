/***********************************************************************
 *
 * libdiffpy         by DANSE Diffraction group
 *                   Simon J. L. Billinge
 *                   (c) 2006 trustees of the Michigan State University
 *                   All rights reserved.
 *
 * File coded by:    Pavol Juhas
 *
 * See AUTHORS.rst for a list of people who contributed.
 * See LICENSE_DANSE.rst for license information.
 *
 ************************************************************************
 *
 * Unit tests for PointsInSphere module
 *
 * Comments:
 *
 ***********************************************************************/

#include <algorithm>

#include <catch2/catch_test_macros.hpp>

#include <diffpy/srreal/PointsInSphere.hpp>

using namespace std;
using namespace diffpy::srreal;

namespace {

const double eps = 1.0e-8;

struct vidxgroup {
  double vijk[4];
  vidxgroup(double v, const int *ijk) {
    vijk[0] = v;
    for (size_t i = 0; i != 3; ++i) {
      vijk[i + 1] = ijk[i];
    }
  }
  vidxgroup(double v, int i, int j, int k) {
    vijk[0] = v;
    vijk[1] = i;
    vijk[2] = j;
    vijk[3] = k;
  }
};

bool operator<(const vidxgroup &x, const vidxgroup &y) {
  if (x.vijk[0] < y.vijk[0] - eps) return true;
  if (x.vijk[0] > y.vijk[0] + eps) return false;
  return lexicographical_compare(x.vijk, x.vijk + 4, y.vijk, y.vijk + 4);
}

bool operator==(const vidxgroup &x, const vidxgroup &y) {
  bool eq = (fabs(x.vijk[0] - y.vijk[0]) < eps) &&
            equal(x.vijk + 1, x.vijk + 4, y.vijk + 1);
  return eq;
}

/* avoid warnings about unused function.
ostream& operator<<(ostream &s, const vidxgroup &x)
{
    return s << "<" << x.vijk[0] << ";" << int(x.vijk[1])
        << ',' << int(x.vijk[2]) << ',' << int(x.vijk[3]) << '>';
}
*/

}  // namespace

////////////////////////////////////////////////////////////////////////
// TestPointsInSphere
////////////////////////////////////////////////////////////////////////

using diffpy::srreal::pointsinsphere::LatticeParameters;

class TestPointsInSphereFixture {
 protected:
  LatticeParameters *latpar;

 public:
  TestPointsInSphereFixture() {
    latpar = new LatticeParameters(1, 1, 1, 90, 90, 90);
  }

  ~TestPointsInSphereFixture() { delete latpar; }

 protected:
  int count(double Rmin, double Rmax) {
    int c = 0;
    PointsInSphere sph(Rmin, Rmax, *latpar);
    for (; !sph.finished(); sph.next()) ++c;
    return c;
  }

  vector<vidxgroup> sortedPoints(double Rmin, double Rmax) {
    vector<vidxgroup> ridx;
    for (PointsInSphere sph(Rmin, Rmax, *latpar); not sph.finished();
         sph.next()) {
      ridx.push_back(vidxgroup(sph.r(), sph.mno()));
    }
    sort(ridx.begin(), ridx.end());
    return ridx;
  }
};  // class TestPointsInSphereFixture

TEST_CASE_METHOD(TestPointsInSphereFixture, "TestPointsInSphere::test_Cubic",
                 "[PointsInSphere]") {
  latpar->a = latpar->b = latpar->c = 1.0;
  latpar->alpha = latpar->beta = latpar->gamma = 90.0;
  latpar->update();
  CHECK(0 == count(0.0, 0.0));
  CHECK(0 == count(eps, 0.5));
  CHECK(0 == count(1.0 + eps, 1.1));
  CHECK(1 == count(0.0, eps));
  CHECK(7 == count(0.0, 1 + eps));
  CHECK(19 == count(0.0, sqrt(2.0) + eps));
  CHECK(12 == count(1.0 + eps, sqrt(2.0) + eps));
}

TEST_CASE_METHOD(TestPointsInSphereFixture,
                 "TestPointsInSphere::test_Orthorombic", "[PointsInSphere]") {
  latpar->a = 1.0;
  latpar->b = 2.0;
  latpar->c = 3.0;
  latpar->alpha = latpar->beta = latpar->gamma = 90.0;
  latpar->update();
  CHECK(3 == count(0.0, 1.1));
  CHECK(4 == count(1.9, 2.1));
  vidxgroup ep[] = {
    vidxgroup(0, 0, 0, 0),          vidxgroup(1, -1, 0, 0),
    vidxgroup(1, 1, 0, 0),          vidxgroup(2, -2, 0, 0),
    vidxgroup(2, 0, -1, 0),         vidxgroup(2, 0, 1, 0),
    vidxgroup(2, 2, 0, 0),          vidxgroup(sqrt(5.0), -1, -1, 0),
    vidxgroup(sqrt(5.0), -1, 1, 0), vidxgroup(sqrt(5.0), 1, -1, 0),
    vidxgroup(sqrt(5.0), 1, 1, 0),  vidxgroup(sqrt(8.0), -2, -1, 0),
    vidxgroup(sqrt(8.0), -2, 1, 0), vidxgroup(sqrt(8.0), 2, -1, 0),
    vidxgroup(sqrt(8.0), 2, 1, 0),  vidxgroup(3, -3, 0, 0),
    vidxgroup(3, 0, 0, -1),         vidxgroup(3, 0, 0, 1),
    vidxgroup(3, 3, 0, 0),
  };
  vector<vidxgroup> exp_pts(ep, ep + sizeof(ep) / sizeof(vidxgroup));
  vector<vidxgroup> act_pts = sortedPoints(0.0, 3.0 + eps);
  CHECK(exp_pts.size() == act_pts.size());
  for (size_t i = 0; i != exp_pts.size(); ++i) {
    CHECK(exp_pts[i] == act_pts[i]);
  }
}

TEST_CASE_METHOD(TestPointsInSphereFixture,
                 "TestPointsInSphere::test_Hexagonal", "[PointsInSphere]") {
  latpar->a = 1.0;
  latpar->b = 1.0;
  latpar->c = 2.0;
  latpar->alpha = latpar->beta = 90.0;
  latpar->gamma = 120.0;
  latpar->update();
  CHECK(7 == count(0.0, 1 + eps));
  vidxgroup ep[] = {
    vidxgroup(0, 0, 0, 0),           vidxgroup(1, -1, -1, 0),
    vidxgroup(1, -1, 0, 0),          vidxgroup(1, 0, -1, 0),
    vidxgroup(1, 0, 1, 0),           vidxgroup(1, 1, 0, 0),
    vidxgroup(1, 1, 1, 0),           vidxgroup(sqrt(3.0), -2, -1, 0),
    vidxgroup(sqrt(3.0), -1, -2, 0), vidxgroup(sqrt(3.0), -1, 1, 0),
    vidxgroup(sqrt(3.0), 1, -1, 0),  vidxgroup(sqrt(3.0), 1, 2, 0),
    vidxgroup(sqrt(3.0), 2, 1, 0),   vidxgroup(2, -2, -2, 0),
    vidxgroup(2, -2, 0, 0),          vidxgroup(2, 0, -2, 0),
    vidxgroup(2, 0, 0, -1),          vidxgroup(2, 0, 0, 1),
    vidxgroup(2, 0, 2, 0),           vidxgroup(2, 2, 0, 0),
    vidxgroup(2, 2, 2, 0),
  };
  vector<vidxgroup> exp_pts(ep, ep + sizeof(ep) / sizeof(vidxgroup));
  vector<vidxgroup> act_pts = sortedPoints(0.0, 2.0 + eps);
  CHECK(exp_pts.size() == act_pts.size());
  for (size_t i = 0; i != exp_pts.size(); ++i) {
    CHECK(exp_pts[i] == act_pts[i]);
  }
}

TEST_CASE_METHOD(TestPointsInSphereFixture, "TestPointsInSphere::test_FCC",
                 "[PointsInSphere]") {
  latpar->a = latpar->b = latpar->c = sqrt(0.5);
  latpar->alpha = latpar->beta = latpar->gamma = 60.0;
  latpar->update();
  CHECK(13 == count(0.0, sqrt(0.5) + eps));
  CHECK(19 == count(0.0, 1.0 + eps));
}

// End of file
