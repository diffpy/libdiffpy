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
* class CroppedGaussianProfile -- Gaussian profile cropped to zero beyond
*     xboundhi and scaled so that its integrated area equals 1.
*     Registered as "croppedgaussian".
*
* $Id$
*
*****************************************************************************/

#include <cmath>

#include <gsl/gsl_sf_erf.h>

#include <diffpy/srreal/CroppedGaussianProfile.hpp>

using namespace std;

namespace diffpy {
namespace srreal {

// Constructors --------------------------------------------------------------

CroppedGaussianProfile::CroppedGaussianProfile()
{
    mscale = 1.0;
}


PeakProfilePtr CroppedGaussianProfile::create() const
{
    PeakProfilePtr rv(new CroppedGaussianProfile());
    return rv;
}


PeakProfilePtr CroppedGaussianProfile::clone() const
{
    PeakProfilePtr rv(new CroppedGaussianProfile(*this));
    return rv;
}

// Public Methods ------------------------------------------------------------

const string& CroppedGaussianProfile::type() const
{
    static string rv = "croppedgaussian";
    return rv;
}


double CroppedGaussianProfile::yvalue(
        double x, double fwhm, double position) const
{
    double xrel = (x - position) / fwhm;
    double rv = (fabs(xrel) >= mhalfboundrel) ? 0.0 :
        2 * sqrt(M_LN2 / M_PI) / fwhm *
        mscale * exp(-4 * M_LN2 * xrel * xrel);
    // Contributions in G(r) need to be normalized by pair distance,
    // not by r as done in PDFfit or PDFfit2.  Here we rescale RDF
    // in such way that division by x will give a correct result.
    rv *= x / position;
    return rv;
}


void CroppedGaussianProfile::setPrecision(double eps)
{
    this->GaussianProfile::setPrecision(eps);
    double eps1 = this->getPrecision();
    mscale = 1.0;
    if (eps1 < 1)
    {
        double area = gsl_sf_erf(sqrt(-log(eps1)));
        mscale = 1.0 / area;
    }
}

// Registration --------------------------------------------------------------

bool reg_CroppedGaussianProfile = CroppedGaussianProfile().registerThisType();

}   // namespace srreal
}   // namespace diffpy

// End of file
