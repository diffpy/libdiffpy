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
 * Definitions of helper functions for other unit tests.
 *
 *****************************************************************************/

#include <fstream>
#include <sstream>
#include <cassert>
#include <cctype>
#include <cmath>
#include <algorithm>
#include <limits>
#include <string_view>

#include <diffpy/runtimepath.hpp>
#include <diffpy/srreal/PeriodicStructureAdapter.hpp>

#include "test_helpers.hpp"

namespace {

char lower_char(char c) {
  return static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
}

bool iequals(std::string_view a, std::string_view b) {
  return std::ranges::equal(a, b, {}, lower_char, lower_char);
}

}  // namespace

// Resolve DIFFPYTESTSDIRPATH ------------------------------------------------

#ifndef DIFFPYTESTSDIRPATH
#error Compiler must define the DIFFPYTESTSDIRPATH macro.
#endif
#define STRINGIFY(m) STRINGIFY_BRAIN_DAMAGE(m)
#define STRINGIFY_BRAIN_DAMAGE(m) #m

// Path utilities ------------------------------------------------------------

std::string prepend_tests_dir(const std::string& f) {
  std::string rv = STRINGIFY(DIFFPYTESTSDIRPATH);
  rv = rv + '/' + f;
  return rv;
}

std::string prepend_testdata_dir(const std::string& f) {
  using std::string;
  string ft = "testdata/";
  ft += f;
  string rv = prepend_tests_dir(ft);
  return rv;
}

// Load PeriodicStructureAdapter from a PdfFit formatted file.

diffpy::srreal::StructureAdapterPtr loadTestPeriodicStructure(
  const std::string& tailname) {
  using namespace std;
  using namespace diffpy::srreal;
  using diffpy::mathutils::eps_eq;
  using diffpy::runtimepath::LineReader;
  // instantiate PeriodicStructureAdapter
  StructureAdapterPtr pstru(new PeriodicStructureAdapter);
  PeriodicStructureAdapter& stru =
    static_cast<PeriodicStructureAdapter&>(*pstru);
  // open data file
  string fullpath = prepend_testdata_dir(tailname);
  ifstream fp(fullpath.c_str());
  LineReader line;
  string fileformat;
  // process structure file header
  while (fp >> line) {
    if (line.wcount() >= 2 && iequals(line.words[0], "format")) {
      fileformat = line.words[1];
      assert(iequals(fileformat, "pdffit"));
    }
    if (line.wcount() == 7 && line.words[0] == "cell") {
      std::ranges::replace(line.line, ',', ' ');
      string skip;
      double a, b, c, alpha, beta, gamma;
      istringstream fpwords(line.line);
      fpwords >> skip >> a >> b >> c >> alpha >> beta >> gamma;
      stru.setLatPar(a, b, c, alpha, beta, gamma);
    }
    if (line.wcount() == 1 && line.words[0] == "atoms") break;
  }
  assert(!fileformat.empty());
  // load atoms
  while (true) {
    Atom ai;
    double ignore;
    R3::Vector& xyz = ai.xyz_cartn;
    R3::Matrix& Uc = ai.uij_cartn;
    fp >> ai.atomtype >> xyz[0] >> xyz[1] >> xyz[2] >> ai.occupancy >> ignore >>
      ignore >> ignore >> ignore >> Uc(0, 0) >> Uc(1, 1) >> Uc(2, 2) >>
      ignore >> ignore >> ignore >> Uc(0, 1) >> Uc(0, 2) >> Uc(1, 2) >>
      ignore >> ignore >> ignore;
    Uc(1, 0) = Uc(0, 1);
    Uc(2, 0) = Uc(0, 2);
    Uc(2, 1) = Uc(1, 2);
    if (!fp) break;
    // convert string symbol to standard case
    std::ranges::transform(ai.atomtype, ai.atomtype.begin(), lower_char);
    string::iterator tt = ai.atomtype.begin();
    for (; tt != ai.atomtype.end(); ++tt) {
      if (isalpha(*tt)) {
        *tt = toupper(*tt);
        break;
      }
    }
    // convert coordinates and the ADP matrix and determine anisotropy
    stru.toCartesian(ai);
    ai.anisotropy = !eps_eq(Uc(0, 0), Uc(1, 1)) || !eps_eq(Uc(0, 0), Uc(2, 2));
    for (int i = 0; !ai.anisotropy && i < R3::Ndim; ++i) {
      for (int j = i + 1; j < R3::Ndim; ++j) {
        ai.anisotropy =
          ai.anisotropy || !eps_eq(0.0, Uc(i, j)) || !eps_eq(0.0, Uc(j, i));
      }
    }
    stru.append(ai);
  }
  return pstru;
}

// Array comparison ----------------------------------------------------------

double maxRelativeError(const diffpy::srreal::QuantityType& actual,
                        const diffpy::srreal::QuantityType& reference,
                        double scale) {
  const double inf = std::numeric_limits<double>::infinity();
  if (actual.size() != reference.size()) return inf;
  double maxdiff = 0.0;
  double maxref = 0.0;
  for (size_t i = 0; i < actual.size(); ++i) {
    if (!std::isfinite(actual[i]) || !std::isfinite(reference[i])) return inf;
    maxdiff = std::max(maxdiff, std::fabs(actual[i] - reference[i]));
    maxref = std::max(maxref, std::fabs(reference[i]));
  }
  if (scale == 0.0) scale = maxref;
  return (scale > 0.0) ? maxdiff / scale : inf;
}

// End of test_helpers.cpp
