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
* class ScaleEnvelope -- constant scaling factor
*
*****************************************************************************/

#include <diffpy/srreal/ScaleEnvelope.hpp>
#include <diffpy/serialization.ipp>

using namespace std;

namespace diffpy {
namespace srreal {

// Constructors --------------------------------------------------------------

ScaleEnvelope::ScaleEnvelope()
{
    this->setScale(1.0);
    this->registerDoubleAttribute("scale",
            this, &ScaleEnvelope::getScale, &ScaleEnvelope::setScale);
}


PDFEnvelopePtr ScaleEnvelope::create() const
{
    PDFEnvelopePtr rv(new ScaleEnvelope());
    return rv;
}


PDFEnvelopePtr ScaleEnvelope::clone() const
{
    PDFEnvelopePtr rv(new ScaleEnvelope(*this));
    return rv;
}

// Public Methods ------------------------------------------------------------

const string& ScaleEnvelope::type() const
{
    static string rv = "scale";
    return rv;
}


double ScaleEnvelope::operator()(const double& r) const
{
    return this->getScale();
}


void ScaleEnvelope::setScale(double sc)
{
    mscale = sc;
}


const double& ScaleEnvelope::getScale() const
{
    return mscale;
}


// Registration --------------------------------------------------------------

bool reg_ScaleEnvelope = ScaleEnvelope().registerThisType();

}   // namespace srreal
}   // namespace diffpy

// Serialization -------------------------------------------------------------

DIFFPY_INSTANTIATE_SERIALIZATION(diffpy::srreal::ScaleEnvelope)
BOOST_CLASS_EXPORT_IMPLEMENT(diffpy::srreal::ScaleEnvelope)

// End of file
