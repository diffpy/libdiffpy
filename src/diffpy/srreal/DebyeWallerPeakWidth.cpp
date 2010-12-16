/*****************************************************************************
*
* diffpy.srreal     by DANSE Diffraction group
*                   Simon J. L. Billinge
*                   (c) 2009 Trustees of the Columbia University
*                   in the City of New York.  All rights reserved.
*
* File coded by:    Christopher Farrow, Pavol Juhas
*
* See AUTHORS.txt for a list of people who contributed.
* See LICENSE.txt for license information.
*
******************************************************************************
*
* class DebyeWallerPeakWidth -- peak width calculated assuming independent
*     thermal vibrations of atoms forming a pair.
*
* $Id$
*
*****************************************************************************/

#include <diffpy/srreal/DebyeWallerPeakWidth.hpp>
#include <diffpy/serialization.hpp>

namespace diffpy {
namespace srreal {

using namespace std;

// Constructors --------------------------------------------------------------

PeakWidthModelPtr DebyeWallerPeakWidth::create() const
{
    PeakWidthModelPtr rv(new DebyeWallerPeakWidth());
    return rv;
}


PeakWidthModelPtr DebyeWallerPeakWidth::clone() const
{
    PeakWidthModelPtr rv(new DebyeWallerPeakWidth(*this));
    return rv;
}

// Public Methods ------------------------------------------------------------

const string& DebyeWallerPeakWidth::type() const
{
    static const string rv = "debye-waller";
    return rv;
}


double DebyeWallerPeakWidth::calculate(const BaseBondGenerator& bnds) const
{
    double msdval = bnds.msd();
    return this->calculateFromMSD(msdval);
}

// Registration --------------------------------------------------------------

bool reg_DebyeWallerPeakWidth = DebyeWallerPeakWidth().registerThisType();

}   // namespace srreal
}   // namespace diffpy

// End of file
