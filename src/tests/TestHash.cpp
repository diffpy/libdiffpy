/*****************************************************************************
 *
 * libdiffpy         by DANSE Diffraction group
 *                   Simon J. L. Billinge
 *                   (c) 2026 The Trustees of Columbia University
 *                   in the City of New York.  All rights reserved.
 *
 * See AUTHORS.rst for a list of people who contributed.
 * See LICENSE_DANSE.rst for license information.
 *
 ******************************************************************************
 *
 * TestHash -- unit tests for the hash utilities in diffpy/hash.hpp
 *
 *****************************************************************************/

#include <string>
#include <utility>

#include <catch2/catch_test_macros.hpp>

#include <diffpy/hash.hpp>
#include <diffpy/srreal/R3linalg.hpp>
#include <diffpy/srreal/AtomicStructureAdapter.hpp>

using namespace std;
using namespace diffpy::srreal;

TEST_CASE("TestHash::test_vector", "[hash]") {
  R3::Vector v0(1.0, 2.0, 3.0);
  R3::Vector v1(1.0, 2.0, 3.0);
  CHECK(hash_value(v0) == hash_value(v1));
  v1[2] = 4.0;
  CHECK(hash_value(v0) != hash_value(v1));
  // negative zero compares equal, so it must hash equal
  R3::Vector z0(0.0, 0.0, 0.0);
  R3::Vector z1(-0.0, 0.0, -0.0);
  REQUIRE(z0 == z1);
  CHECK(hash_value(z0) == hash_value(z1));
}

TEST_CASE("TestHash::test_atom", "[hash]") {
  Atom a0;
  a0.atomtype = "C";
  a0.xyz_cartn = R3::Vector(1, 2, 3);
  Atom a1 = a0;
  CHECK(hash_value(a0) == hash_value(a1));
  a1.uij_cartn(0, 0) = 0.01;
  CHECK(hash_value(a0) != hash_value(a1));
  // hash_of finds the ADL hash_value overload
  CHECK(diffpy::hash_of(a0) == hash_value(a0));
}

TEST_CASE("TestHash::test_hash_values", "[hash]") {
  const R3::Vector v(1.0, 2.0, 3.0);
  size_t seed = 0;
  diffpy::hash_combine(seed, string("C"));
  diffpy::hash_combine(seed, v);
  diffpy::hash_combine(seed, 7);
  CHECK(diffpy::hash_values(string("C"), v, 7) == seed);
  CHECK(diffpy::hash_values(1, 2) != diffpy::hash_values(2, 1));
  CHECK(diffpy::hash_values() == 0);
}

TEST_CASE("TestHash::test_pair_hash", "[hash]") {
  diffpy::PairHash h;
  CHECK(h(make_pair(1, 2)) != h(make_pair(2, 1)));
}
