/*****************************************************************************
*
* libdiffpy         by DANSE Diffraction group
*                   Simon J. L. Billinge
*                   (c) 2009 The Trustees of Columbia University
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
*****************************************************************************/

#include <diffpy/srreal/DebyeWallerPeakWidth.hpp>
#include <diffpy/srreal/StructureAdapter.hpp>
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
    using diffpy::mathutils::GAUSS_SIGMA_TO_FWHM;
    double msdval = bnds.msd();
    double rv = (msdval < 0.0) ? 0.0 : GAUSS_SIGMA_TO_FWHM * sqrt(msdval);
    return rv;
}


double DebyeWallerPeakWidth::maxWidth(StructureAdapterPtr stru,
                double rmin, double rmax) const
{
    using diffpy::mathutils::GAUSS_SIGMA_TO_FWHM;
    double maxmsd = 2 * maxUii(stru);
    double rv = (maxmsd <= 0.0) ? 0.0 : GAUSS_SIGMA_TO_FWHM * sqrt(maxmsd);
    return rv;
}

// Registration --------------------------------------------------------------

bool reg_DebyeWallerPeakWidth = DebyeWallerPeakWidth().registerThisType();

}   // namespace srreal
}   // namespace diffpy

// End of file
