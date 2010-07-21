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
* class ConstantPeakWidth -- constant peak width model for testing
*
* $Id$
*
*****************************************************************************/

#include <diffpy/srreal/ConstantPeakWidth.hpp>

using namespace std;
using namespace diffpy::srreal;

// Constructors --------------------------------------------------------------

ConstantPeakWidth::ConstantPeakWidth()
{
    this->setWidth(0.0);
    this->registerDoubleAttribute("width", this,
            &ConstantPeakWidth::getWidth, &ConstantPeakWidth::setWidth);
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


double ConstantPeakWidth::calculateFromMSD(double msdval) const
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
    mwidth = width;
}

// Registration --------------------------------------------------------------

bool reg_ConstantPeakWidth = ConstantPeakWidth().registerThisType();

// End of file
