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
#include <boost/functional/hash.hpp>
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
    mB = 0.0;
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
    mB = B;
    mref_id = ref_id;
}

// Public Methods ------------------------------------------------------------

bool BVParam::operator==(const BVParam& other) const
{
    return
        mvalence0 == other.mvalence0 && mvalence1 == other.mvalence1 &&
        matom0 == other.matom0 && matom1 == other.matom1;
}


double BVParam::bondvalence(double distance) const
{
    double rv = (mB > 0.0) ? exp((mRo - distance) / mB) : 0.0;
    return rv;
}


double BVParam::bondvalenceToDistance(double bvalence) const
{
    double rv = mRo - mB * log(bvalence);
    return rv;
}


void BVParam::setFromCifLine(const std::string& cifline)
{
    string a0, a1, ref;
    int v0, v1;
    double Ro, B;
    istringstream linefp(cifline);
    linefp >> a0 >> v0 >> a1 >> v1 >> Ro >> B >> ref;
    if (!linefp)
    {
        const char* emsg = "Cannot parse cif line.";
        throw invalid_argument(emsg);
    }
    BVParam bp1(a0, v0, a1, v1, Ro, B, ref);
    *this = bp1;
}

// Functions -----------------------------------------------------------------

size_t hash_value(const BVParam& bp)
{
    size_t seed = 0;
    boost::hash_combine(seed, bp.matom0);
    boost::hash_combine(seed, bp.mvalence0);
    boost::hash_combine(seed, bp.matom1);
    boost::hash_combine(seed, bp.mvalence1);
    return seed;
};

}   // namespace srreal
}   // namespace diffpy

// End of file
