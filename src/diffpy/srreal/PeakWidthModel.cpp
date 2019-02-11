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
* class PeakWidthModel -- base class for calculation of peak widths.
*
*****************************************************************************/

#include <diffpy/srreal/PeakWidthModel.hpp>
#include <diffpy/HasClassRegistry.ipp>
#include <diffpy/validators.hpp>
#include <diffpy/serialization.ipp>

using std::string;
using diffpy::validators::ensureNonNull;

namespace diffpy {

// Unique instantiation of the template registry base class.
template class HasClassRegistry<srreal::PeakWidthModel>;

namespace srreal {

// class PeakWidthModelOwner -------------------------------------------------

void PeakWidthModelOwner::setPeakWidthModel(PeakWidthModelPtr pwm)
{
    ensureNonNull("PeakWidthModel", pwm);
    if (mpwmodel != pwm)  mprivateticker.click();
    mpwmodel = pwm;
}


void PeakWidthModelOwner::setPeakWidthModelByType(const string& tp)
{
    mpwmodel = PeakWidthModel::createByType(tp);
    mprivateticker.click();
}


PeakWidthModelPtr& PeakWidthModelOwner::getPeakWidthModel()
{
    return mpwmodel;
}


const PeakWidthModelPtr& PeakWidthModelOwner::getPeakWidthModel() const
{
    return mpwmodel;
}


eventticker::EventTicker& PeakWidthModelOwner::ticker() const
{
    if (mpwmodel)  mprivateticker.updateFrom(mpwmodel->ticker());
    return mprivateticker;
}

}   // namespace srreal
}   // namespace diffpy

// Serialization -------------------------------------------------------------

DIFFPY_INSTANTIATE_SERIALIZATION(diffpy::srreal::PeakWidthModel)
DIFFPY_INSTANTIATE_PTR_SERIALIZATION(diffpy::srreal::PeakWidthModel)

// End of file
