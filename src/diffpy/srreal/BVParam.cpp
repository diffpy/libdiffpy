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

// Public Methods ------------------------------------------------------------

void BVParam::setFromCifLine(const std::string& cifline)
{
    BVParam bv1;
    istringstream linefp(cifline);
    linefp >> bv1.matom0 >> bv1.mvalence0 >>
        bv1.matom1 >> bv1.mvalence1 >> bv1.mRo >> bv1.mB >> bv1.mref_id;
    if (!linefp)
    {
        const char* emsg = "Cannot parse cif line.";
        throw invalid_argument(emsg);
    }
    *this = bv1;
}


double BVParam::atd(double distance)
{
    double rv = exp((mRo - distance) / mB);
    return rv;
}

}   // namespace srreal
}   // namespace diffpy

// End of file
