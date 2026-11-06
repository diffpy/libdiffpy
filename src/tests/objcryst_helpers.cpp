/*****************************************************************************
 *
 * libdiffpy         by DANSE Diffraction group
 *                   Simon J. L. Billinge
 *                   (c) 2010 The Trustees of Columbia University
 *                   in the City of New York.  All rights reserved.
 *
 * File coded by:    Pavol Juhas
 *
 * See AUTHORS.rst for a list of people who contributed.
 * See LICENSE_DANSE.rst for license information.
 *
 ******************************************************************************
 *
 * shared helpers for loading diffpy structure objects.
 *
 *****************************************************************************/

#include <fstream>
#include <map>
#include <memory>
#include <stdexcept>

#include <ObjCryst/ObjCryst/CIF.h>

#include "objcryst_helpers.hpp"
#include "test_helpers.hpp"

using namespace std;

namespace {

struct QuietObjCryst {
  decltype(ObjCryst::fpObjCrystInformUser) saved =
    ObjCryst::fpObjCrystInformUser;
  QuietObjCryst() {
    ObjCryst::fpObjCrystInformUser = [](const string&) {};
  }
  ~QuietObjCryst() { ObjCryst::fpObjCrystInformUser = saved; }
};

}  // namespace

ObjCryst::Crystal* loadTestCrystal(const string& tailname) {
  string fp = prepend_testdata_dir(tailname);
  ifstream in(fp.c_str());
  if (!in) throw runtime_error("Cannot open " + fp);
  QuietObjCryst quiet;
  ObjCryst::CIF cif(in);
  ObjCryst::Crystal* crst = ObjCryst::CreateCrystalFromCIF(cif, false);
  if (!crst) throw runtime_error("No crystal structure in " + fp);
  return crst;
}

const ObjCryst::Crystal& cachedTestCrystal(const string& tailname) {
  static map<string, unique_ptr<ObjCryst::Crystal>> crystals;
  auto& crystal = crystals[tailname];
  if (!crystal) crystal.reset(loadTestCrystal(tailname));
  return *crystal;
}

// End of file
