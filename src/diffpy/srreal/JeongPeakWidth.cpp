/*****************************************************************************
*
* libdiffpy         by DANSE Diffraction group
*                   Simon J. L. Billinge
*                   (c) 2009 The Trustees of Columbia University
*                   in the City of New York.  All rights reserved.
*
* File coded by:    Christopher Farrow, Pavol Juhas
*
* See AUTHORS.txt for a list of people who contributed.
* See LICENSE_DANSE.txt for license information.
*
******************************************************************************
*
* class DebyeWallerPeakWidth -- peak width model based on
*      I.-K. Jeong, et al., Phys. Rev. B 67, 104301 (2003)
*      http://link.aps.org/doi/10.1103/PhysRevB.67.104301
*
*****************************************************************************/

#include <diffpy/srreal/JeongPeakWidth.hpp>
#include <diffpy/mathutils.hpp>
#include <diffpy/serialization.ipp>

namespace diffpy {
namespace srreal {

using namespace std;

// Constructors --------------------------------------------------------------

JeongPeakWidth::JeongPeakWidth() :
    mdelta1(0.0), mdelta2(0.0), mqbroad(0.0)
{
    this->registerDoubleAttribute("delta1",
            this, &JeongPeakWidth::getDelta1, &JeongPeakWidth::setDelta1);
    this->registerDoubleAttribute("delta2",
            this, &JeongPeakWidth::getDelta2, &JeongPeakWidth::setDelta2);
    this->registerDoubleAttribute("qbroad",
            this, &JeongPeakWidth::getQbroad, &JeongPeakWidth::setQbroad);
    this->registerDoubleAttribute("qbroad_new",
            this, &JeongPeakWidth::getQbroad_new, &JeongPeakWidth::setQbroad_new);
}


PeakWidthModelPtr JeongPeakWidth::create() const
{
    PeakWidthModelPtr rv(new JeongPeakWidth());
    return rv;
}


PeakWidthModelPtr JeongPeakWidth::clone() const
{
    PeakWidthModelPtr rv(new JeongPeakWidth(*this));
    return rv;
}

// Public Methods ------------------------------------------------------------

const string& JeongPeakWidth::type() const
{
    static const string rv = "jeong";
    return rv;
}


double JeongPeakWidth::calculate(const BaseBondGenerator& bnds) const
{
    double r = bnds.distance();
    double corr = this->msdSharpeningRatio(r);
    // avoid calculating square root of negative value
    double fwhm = (corr <= 0) ? 0.0 :
        (sqrt(corr) * this->DebyeWallerPeakWidth::calculate(bnds) +
        pow(this->getQbroad_new()*r, 2));
    return fwhm;
}


double JeongPeakWidth::maxWidth(StructureAdapterPtr stru,
        double rmin, double rmax) const
{
    double maxwidth0 = this->DebyeWallerPeakWidth::maxWidth(stru, rmin, rmax);
    double maxmsdsharp = max(
            this->msdSharpeningRatio(rmin),
            this->msdSharpeningRatio(rmax));
    // if the sharpening factor is larger than 1 adjust the maximum width
    double rv = maxwidth0 * sqrt(max(1.0, maxmsdsharp));
    return rv;
}

const double& JeongPeakWidth::getDelta1() const
{
    return mdelta1;
}


void JeongPeakWidth::setDelta1(double delta1)
{
    if (mdelta1 != delta1)  mticker.click();
    mdelta1 = delta1;
}


const double& JeongPeakWidth::getDelta2() const
{
    return mdelta2;
}


void JeongPeakWidth::setDelta2(double delta2)
{
    if (mdelta2 != delta2)  mticker.click();
    mdelta2 = delta2;
}


const double& JeongPeakWidth::getQbroad() const
{
    return mqbroad;
}


const double& JeongPeakWidth::getQbroad_new() const
{
    return mqbroad_new;
}


void JeongPeakWidth::setQbroad_new(double qbroad_new)
{
    if (mqbroad_new != qbroad_new)  mticker.click();
    mqbroad_new = qbroad_new;
}


void JeongPeakWidth::setQbroad(double qbroad)
{
    if (mqbroad != qbroad)  mticker.click();
    mqbroad = qbroad;
}

// Private Methods -----------------------------------------------------------

double JeongPeakWidth::msdSharpeningRatio(const double& r) const
{
    using diffpy::mathutils::DOUBLE_EPS;
    // avoid division by zero
    if (r < DOUBLE_EPS)  return 0.0;
    double rv = 1.0 - this->getDelta1()/r - this->getDelta2()/pow(r, 2) +
         pow(this->getQbroad()*r, 2);
    return rv;
}

// Registration --------------------------------------------------------------

bool reg_JeongPeakWidth = JeongPeakWidth().registerThisType();

}   // namespace srreal
}   // namespace diffpy

// Serialization -------------------------------------------------------------

DIFFPY_INSTANTIATE_SERIALIZATION(diffpy::srreal::JeongPeakWidth)
BOOST_CLASS_EXPORT_IMPLEMENT(diffpy::srreal::JeongPeakWidth)

// End of file
