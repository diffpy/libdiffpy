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
* Concrete implementations of the abstract PeakProfile class:
*
* class GaussianProfile -- registered as "gaussian"
*
* $Id$
*
*****************************************************************************/

#include <diffpy/srreal/PeakProfile.hpp>
#include <diffpy/HasClassRegistry.ipp>

using diffpy::srreal::PeakProfile;

// Unique instantiation of the template registry base class.
template class HasClassRegistry<PeakProfile>;

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

// End of file
