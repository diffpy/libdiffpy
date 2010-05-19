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
* Various common routines useful for PDF calculation.
*
* $Id$
*
*****************************************************************************/

#include <stdexcept>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_fft_complex.h>

#include <diffpy/srreal/PDFUtils.hpp>
#include <diffpy/mathutils.hpp>
#include <diffpy/validators.hpp>

using namespace std;
using namespace diffpy::validators;
using namespace diffpy::mathutils;

namespace diffpy {
namespace srreal {

// Local Constants -----------------------------------------------------------

namespace {

const char* EMSGFFT = "Fourier Transformation failed.";

}   // namespace

// PDFUtils functions --------------------------------------------------------

double meanSquareDisplacement(const R3::Matrix& Uijcartn,
        const R3::Vector& s, bool anisotropy)
{
    double rv;
    if (anisotropy)
    {
        assert(R3::norm(s) > 0);
        assert(eps_eq(Uijcartn(0,1), Uijcartn(1,0)));
        assert(eps_eq(Uijcartn(0,2), Uijcartn(2,0)));
        assert(eps_eq(Uijcartn(1,2), Uijcartn(2,1)));
        static R3::Vector sn;
        sn = s / R3::norm(s);
        rv = Uijcartn(0,0) * sn(0) * sn(0) +
             Uijcartn(1,1) * sn(1) * sn(1) +
             Uijcartn(2,2) * sn(2) * sn(2) +
             2 * Uijcartn(0,1) * sn(0) * sn(1) +
             2 * Uijcartn(0,2) * sn(0) * sn(2) +
             2 * Uijcartn(1,2) * sn(1) * sn(2);
    }
    else
    {
        assert(eps_eq(Uijcartn(0,0), Uijcartn(1,1)));
        assert(eps_eq(Uijcartn(0,0), Uijcartn(2,2)));
        rv = Uijcartn(0,0);
    }
    return rv;
}


double maxUii(const StructureAdapter* stru)
{
    if (!stru)  return 0.0;
    double rv = 0.0;
    for (int i = 0; i < stru->countSites(); ++i)
    {
        const R3::Matrix& U = stru->siteCartesianUij(i);
        for (int k = 0; k < R3::Ndim; k++)
        {
            if (U(k,k) > rv)   rv = U(k,k);
        }
    }
    return rv;
}


QuantityType fftgtof(const QuantityType& g, double rstep, double rmin)
{
    int padrmin = round(rmin / rstep);
    int Npad1 = padrmin + g.size();
    // pad to the next power of 2 for fast Fourier transformation
    int Npad2 = pow(2, int(ceil(log2(Npad1))));
    // sine transformations needs an odd extension
    // gpadc array has to be doubled for complex coefficients
    int Npad4 = 4 * Npad2;
    valarray<double> gpadc(0.0, Npad4);
    QuantityType::const_iterator gi;
    // copy the original g signal
    int ilo = padrmin;
    for (gi = g.begin(); gi != g.end(); ++gi, ++ilo)
    {
        gpadc[2 * ilo] = *gi;
    }
    // copy the odd part of g skipping the first point,
    // because it is periodic image of gpadc[0]
    int ihi = 2 * Npad2 - 1;
    for (ilo = 1; ilo < Npad2; ++ilo, --ihi)
    {
        gpadc[2 * ihi] = -1 * gpadc[2 * ilo];
    }
    int status;
    status = gsl_fft_complex_radix2_inverse(&(gpadc[0]), 1, 2 * Npad2);
    if (status != GSL_SUCCESS)  throw invalid_argument(EMSGFFT);
    QuantityType f(Npad2);
    for (int i = 0; i < Npad2; ++i)
    {
        f[i] = gpadc[2 * i + 1] * Npad2 * rstep;
    }
    // real components should be all close to zero
    assert(fabs(valarray<double>(gpadc[slice(0, 2 * Npad2, 2)]).min()) <=
            SQRT_DOUBLE_EPS * gpadc.max());
    assert(fabs(valarray<double>(gpadc[slice(0, 2 * Npad2, 2)]).max()) <=
            SQRT_DOUBLE_EPS * gpadc.max());
    return f;
}


QuantityType fftftog(const QuantityType& f, double qstep, double qmin)
{
    // ftog is the same as gtof, just normalized by 2 / PI
    QuantityType g = fftgtof(f, qstep, qmin);
    QuantityType::iterator gi;
    for (gi = g.begin(); gi != g.end(); ++gi)  *gi *= 2.0 / M_PI;
    return g;
}


}   // namespace srreal
}   // namespace diffpy

// End of file
