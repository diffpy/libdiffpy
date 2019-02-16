/*****************************************************************************
*
* libdiffpy         by DANSE Diffraction group
*                   Simon J. L. Billinge
*                   (c) 2009 The Trustees of Columbia University
*                   in the City of New York.  All rights reserved.
*
* File coded by:    Pavol Juhas
*
* See AUTHORS.txt for a list of people who contributed.
* See LICENSE_DANSE.txt for license information.
*
******************************************************************************
*
* class ConstantPeakWidth -- constant peak width model for testing
*
*****************************************************************************/

#include <diffpy/srreal/ConstantPeakWidth.hpp>
#include <diffpy/serialization.ipp>

namespace diffpy {
namespace srreal {

using namespace std;

// Local Helpers -------------------------------------------------------------

namespace {

using diffpy::mathutils::GAUSS_SIGMA_TO_FWHM;

double getuisowidth(const ConstantPeakWidth* ppwm)
{
    const double rmsd = ppwm->getWidth() / GAUSS_SIGMA_TO_FWHM;
    const int pm = (rmsd >= 0) ? +1 : -1;
    return pm * 0.5 * rmsd * rmsd;
}

void setuisowidth(ConstantPeakWidth* ppwm, const double& uiso)
{
    const int pm = (uiso >= 0) ? +1 : -1;
    const double uisoplus = pm * uiso;
    const double fwhm = pm * GAUSS_SIGMA_TO_FWHM * sqrt(2.0 * uisoplus);
    ppwm->setWidth(fwhm);
}

}   // namespace

// Constructors --------------------------------------------------------------

ConstantPeakWidth::ConstantPeakWidth() : mwidth(0.0)
{
    this->registerDoubleAttribute("width", this,
            &ConstantPeakWidth::getWidth, &ConstantPeakWidth::setWidth);
    this->registerDoubleAttribute("uisowidth", this,
            &getuisowidth, &setuisowidth);
}


PeakWidthModelPtr ConstantPeakWidth::create() const
{
    PeakWidthModelPtr rv(new ConstantPeakWidth());
    return rv;
}


PeakWidthModelPtr ConstantPeakWidth::clone() const
{
    PeakWidthModelPtr rv(new ConstantPeakWidth(*this));
    return rv;
}

// Public Methods ------------------------------------------------------------

const string& ConstantPeakWidth::type() const
{
    static const string rv = "constant";
    return rv;
}


double ConstantPeakWidth::calculate(const BaseBondGenerator& bnds) const
{
    return this->getWidth();
}


double ConstantPeakWidth::maxWidth(
        StructureAdapterPtr stru, double rmin, double rmax) const
{
    return this->getWidth();
}

// data access

const double& ConstantPeakWidth::getWidth() const
{
    return mwidth;
}


void ConstantPeakWidth::setWidth(double width)
{
    if (mwidth != width)  mticker.click();
    mwidth = width;
}

// Registration --------------------------------------------------------------

bool reg_ConstantPeakWidth = ConstantPeakWidth().registerThisType();

}   // namespace srreal
}   // namespace diffpy

// Serialization -------------------------------------------------------------

DIFFPY_INSTANTIATE_SERIALIZATION(diffpy::srreal::ConstantPeakWidth)
BOOST_CLASS_EXPORT_IMPLEMENT(diffpy::srreal::ConstantPeakWidth)

// End of file
