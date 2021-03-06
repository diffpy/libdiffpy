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
* class QResolutionEnvelope -- Gaussian envelope due to limited Q resolution
*
*****************************************************************************/

#include <cmath>

#include <diffpy/serialization.ipp>
#include <diffpy/srreal/QResolutionEnvelope.hpp>

using namespace std;

namespace diffpy {
namespace srreal {

// Constructor ---------------------------------------------------------------

QResolutionEnvelope::QResolutionEnvelope()
{
    this->setQdamp(0.0);
    this->registerDoubleAttribute("qdamp", this,
            &QResolutionEnvelope::getQdamp, &QResolutionEnvelope::setQdamp);
}


PDFEnvelopePtr QResolutionEnvelope::create() const
{
    PDFEnvelopePtr rv(new QResolutionEnvelope());
    return rv;
}


PDFEnvelopePtr QResolutionEnvelope::clone() const
{
    PDFEnvelopePtr rv(new QResolutionEnvelope(*this));
    return rv;
}

// Public Methods ------------------------------------------------------------

const string& QResolutionEnvelope::type() const
{
    static string rv = "qresolution";
    return rv;
}


double QResolutionEnvelope::operator()(const double& r) const
{
    double rv = (mqdamp > 0.0) ?
        exp(-pow(r * mqdamp, 2) / 2) :
        1.0;
    return rv;
}


void QResolutionEnvelope::setQdamp(double sc)
{
    mqdamp = sc;
}


const double& QResolutionEnvelope::getQdamp() const
{
    return mqdamp;
}

// Registration --------------------------------------------------------------

bool reg_QResolutionEnvelope = QResolutionEnvelope().registerThisType();

}   // namespace srreal
}   // namespace diffpy

// Serialization -------------------------------------------------------------

DIFFPY_INSTANTIATE_SERIALIZATION(diffpy::srreal::QResolutionEnvelope)
BOOST_CLASS_EXPORT_IMPLEMENT(diffpy::srreal::QResolutionEnvelope)

// End of file
