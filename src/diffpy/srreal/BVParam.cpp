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
* class BVParam -- bond valence parameters for a cation-anion pair
*
* $Id$
*
*****************************************************************************/

#include <cmath>
#include <stdexcept>
#include <sstream>

#include <diffpy/srreal/BVParam.hpp>

using namespace std;

namespace diffpy {
namespace srreal {

// Constructor ---------------------------------------------------------------

BVParam::BVParam()
{
    mvalence0 = 0;
    mvalence1 = 0;
    mRo = 0.0;
    mB = 0.37;
}

BVParam::BVParam(const string& atom0, int valence0,
                const string& atom1, int valence1,
                double Ro, double B, string ref_id)
{
    matom0 = atom0;
    mvalence0 = valence0;
    matom1 = atom1;
    mvalence1 = valence1;
    if (mvalence0 < mvalence1)
    {
        swap(matom0, matom1);
        swap(mvalence0, mvalence1);
    }
    mRo = Ro;
    mB = mB;
    mref_id = ref_id;
}

// Public Methods ------------------------------------------------------------

double BVParam::bondvalence(double distance) const
{
    double rv = (mRo > 0.0) ? exp((mRo - distance) / mB) : 0.0;
    return rv;
}


double BVParam::bonddistance(double bvalence) const
{
    double rv = mRo - mB * log(bvalence);
    return rv;
}


void BVParam::setFromCifLine(const std::string& cifline)
{
    BVParam bp1;
    istringstream linefp(cifline);
    linefp >> bp1.matom0 >> bp1.mvalence0 >>
        bp1.matom1 >> bp1.mvalence1 >> bp1.mRo >> bp1.mB >> bp1.mref_id;
    if (!linefp)
    {
        const char* emsg = "Cannot parse cif line.";
        throw invalid_argument(emsg);
    }
    *this = bp1;
}


// class BVParam::CompareIons
bool BVParam::CompareIons::operator()(
        const BVParam& bp0, const BVParam& bp1) const
{
    if (bp0.matom0 < bp1.matom0)  return true;
    if (bp0.matom0 > bp1.matom0)  return false;
    if (bp0.mvalence0 < bp1.mvalence0)  return true;
    if (bp0.mvalence0 > bp1.mvalence0)  return false;
    if (bp0.matom1 < bp1.matom1)  return true;
    if (bp0.matom1 > bp1.matom1)  return false;
    if (bp0.mvalence1 < bp1.mvalence1)  return true;
    if (bp0.mvalence1 > bp1.mvalence1)  return false;
    return false;
}

}   // namespace srreal
}   // namespace diffpy

// End of file
