/*****************************************************************************
*
* diffpy.srreal     by DANSE Diffraction group
*                   Simon J. L. Billinge
*                   (c) 2009 The Trustees of Columbia University
*                   in the City of New York.  All rights reserved.
*
* File coded by:    Pavol Juhas
*
* See AUTHORS.txt for a list of people who contributed.
* See LICENSE.txt for license information.
*
******************************************************************************
*
* class LinearBaseline -- linear PDF baseline
*
*****************************************************************************/

#include <diffpy/srreal/LinearBaseline.hpp>

using namespace std;

namespace diffpy {
namespace srreal {

// Constructors --------------------------------------------------------------

LinearBaseline::LinearBaseline()
{
    this->setSlope(0.0);
    this->registerDoubleAttribute("slope", this,
            &LinearBaseline::getSlope,
            &LinearBaseline::setSlope);
}


PDFBaselinePtr LinearBaseline::create() const
{
    PDFBaselinePtr rv(new LinearBaseline());
    return rv;
}


PDFBaselinePtr LinearBaseline::clone() const
{
    PDFBaselinePtr rv(new LinearBaseline(*this));
    return rv;
}

// Public Methods ------------------------------------------------------------

const string& LinearBaseline::type() const
{
    static string rv = "linear";
    return rv;
}


double LinearBaseline::operator()(const double& r) const
{
    return (this->getSlope() * r);
}


void LinearBaseline::setSlope(double sc)
{
    mslope = sc;
}


const double& LinearBaseline::getSlope() const
{
    return mslope;
}


// Registration --------------------------------------------------------------

bool reg_LinearBaseline = LinearBaseline().registerThisType();

}   // namespace srreal
}   // namespace diffpy

// End of file
