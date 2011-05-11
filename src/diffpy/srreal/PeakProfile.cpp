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
* class PeakProfile -- base class for calculation of peak profiles.
*
* $Id$
*
*****************************************************************************/

#include <boost/serialization/export.hpp>

#include <diffpy/srreal/PeakProfile.hpp>
#include <diffpy/HasClassRegistry.ipp>

// Unique instantiation of the template registry base class.
template class HasClassRegistry<diffpy::srreal::PeakProfile>;

namespace diffpy {
namespace srreal {

//////////////////////////////////////////////////////////////////////////////
// class PeakProfile
//////////////////////////////////////////////////////////////////////////////

// Constructors --------------------------------------------------------------

PeakProfile::PeakProfile() : mprecision(0.0)
{
    this->registerDoubleAttribute("peakprecision",
            this, &PeakProfile::getPrecision, &PeakProfile::setPrecision);
}

// Public Methods ------------------------------------------------------------

void PeakProfile::setPrecision(double eps)
{
    mprecision = eps;
}


const double& PeakProfile::getPrecision() const
{
    return mprecision;
}

}   // namespace srreal
}   // namespace diffpy

// End of file
