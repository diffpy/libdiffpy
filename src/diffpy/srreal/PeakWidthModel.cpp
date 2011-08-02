/*****************************************************************************
*
* diffpy.srreal     by DANSE Diffraction group
*                   Simon J. L. Billinge
*                   (c) 2009 Trustees of the Columbia University
*                   in the City of New York.  All rights reserved.
*
* File coded by:    Pavol Juhas
*
* See AUTHORS.txt for a list of people who contributed.
* See LICENSE.txt for license information.
*
******************************************************************************
*
* class PeakWidthModel -- base class for calculation of peak widths.
*
* $Id$
*
*****************************************************************************/

#include <boost/serialization/export.hpp>

#include <diffpy/srreal/PeakWidthModel.hpp>
#include <diffpy/HasClassRegistry.ipp>
#include <diffpy/validators.hpp>

using std::string;
using diffpy::srreal::PeakWidthModel;
using diffpy::validators::ensureNonNull;

// Unique instantiation of the template registry base class.
template class HasClassRegistry<PeakWidthModel>;

namespace diffpy {
namespace srreal {

// class PeakWidthModel ------------------------------------------------------

double PeakWidthModel::calculateFromMSD(double msdval) const
{
    const double tofwhm = 2 * sqrt(2 * M_LN2);
    double fwhm = (msdval < 0.0) ? 0.0 : tofwhm * sqrt(msdval);
    return fwhm;
}

// class PeakWidthModelOwner -------------------------------------------------

void PeakWidthModelOwner::setPeakWidthModel(PeakWidthModelPtr pwm)
{
    ensureNonNull("PeakWidthModel", pwm);
    mpwmodel = pwm;
}


void PeakWidthModelOwner::setPeakWidthModelByType(const string& tp)
{
    mpwmodel = PeakWidthModel::createByType(tp);
}


PeakWidthModelPtr& PeakWidthModelOwner::getPeakWidthModel()
{
    return mpwmodel;
}


const PeakWidthModelPtr& PeakWidthModelOwner::getPeakWidthModel() const
{
    return mpwmodel;
}

}   // namespace srreal
}   // namespace diffpy

// End of file
