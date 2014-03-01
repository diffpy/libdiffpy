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
* class SphericalShapeEnvelope -- PDF envelope due to spherical shape factor
*
*****************************************************************************/

#include <cmath>

#include <diffpy/serialization.hpp>
#include <diffpy/srreal/SphericalShapeEnvelope.hpp>

using namespace std;

namespace diffpy {
namespace srreal {

// Constructor ---------------------------------------------------------------

SphericalShapeEnvelope::SphericalShapeEnvelope()
{
    this->setSPDiameter(0.0);
    this->registerDoubleAttribute("spdiameter", this,
            &SphericalShapeEnvelope::getSPDiameter,
            &SphericalShapeEnvelope::setSPDiameter);
}


PDFEnvelopePtr SphericalShapeEnvelope::create() const
{
    PDFEnvelopePtr rv(new SphericalShapeEnvelope());
    return rv;
}


PDFEnvelopePtr SphericalShapeEnvelope::clone() const
{
    PDFEnvelopePtr rv(new SphericalShapeEnvelope(*this));
    return rv;
}

// Public Methods ------------------------------------------------------------

const string& SphericalShapeEnvelope::type() const
{
    static string rv = "sphericalshape";
    return rv;
}


double SphericalShapeEnvelope::operator()(const double& r) const
{
    if (mspdiameter <= 0.0)  return 1.0;
    if (r > mspdiameter)  return 0.0;
    double rdratio = r / mspdiameter;
    double rv = 1.0 - 1.5 * rdratio + 0.5 * pow(rdratio, 3);
    return rv;
}


void SphericalShapeEnvelope::setSPDiameter(double spd)
{
    mspdiameter = spd;
}


const double& SphericalShapeEnvelope::getSPDiameter() const
{
    return mspdiameter;
}

// Registration --------------------------------------------------------------

bool reg_SphericalShapeEnvelope = SphericalShapeEnvelope().registerThisType();

}   // namespace srreal
}   // namespace diffpy

// End of file
