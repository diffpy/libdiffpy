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
* Various common routines useful for PDF calculation:
*     meanSquareDisplacement
*     bandPassFilter
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


void bandPassFilterCValarray(valarray<double>& ycpa, double dr,
        double qmin, double qmax)
{
    // error message for FT failure
    double* yc = &(ycpa[0]);
    // ycpa is a complex array, its actual length is half the size
    int padlen = ycpa.size() / 2;
    // apply fft
    int status;
    status = gsl_fft_complex_radix2_forward(yc, 1, padlen);
    if (status != GSL_SUCCESS)
    {
        throw invalid_argument(EMSGFFT);
    }
// FIXME: the following 2 lines would force sine FFT, but I would
// rather have a shared FtoG function for all PDF calculators.
//    ycpa[slice(0, padlen, 2)] = 0.0;
//    ycpa *= 2.0;
    // Q step for yc
    double dQ = 2 * M_PI / ((padlen - 1) * dr);
    // Cut away high-Q frequencies -
    // loQmaxidx, hiQmaxidx correspond to Qmax and -Qmax frequencies
    // they need to be integer to catch cases with huge qmax/dQ
    int loQmaxidx = int( ceil(qmax/dQ) );
    int hiQmaxidx = padlen + 1 - loQmaxidx;
    hiQmaxidx = min(padlen, hiQmaxidx);
    // zero high Q components in yc
    for (int i = loQmaxidx; i < hiQmaxidx; ++i)
    {
	yc[2 * i] = yc[2 * i + 1] = 0.0;
    }
    // Cut away low-Q frequencies, while keeping the absolut offset.
    // loQminidx corresponds to the Qmin frequency.
    int loQminidx = (int) min(ceil(qmin / dQ), padlen / 2.0);
    for (int i = 1; i < loQminidx; ++i)
    {
        assert(2 * i + 1 < padlen);
        yc[2 * i] = yc[2 * i + 1] = 0.0;
        assert(padlen - 2 * i >= 0);
        yc[padlen - 2 * i] = yc[padlen - 2 * i + 1] = 0.0;
    }
    // transform back
    status = gsl_fft_complex_radix2_inverse(yc, 1, padlen);
    if (status != GSL_SUCCESS)
    {
        throw invalid_argument(EMSGFFT);
    }
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
    // copy the odd part of g skipping the first point in rmin-padded,
    // because it is periodic image of gpadc[0]
    gi = (padrmin > 0) ? g.begin() : g.empty() ? g.begin() : (g.begin() + 1);
    int ihi = 2 * Npad2 - padrmin - 1;
    for (; gi != g.end(); ++gi, --ihi)
    {
        gpadc[2 * ihi] = -1 * (*gi);
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
    assert(fabs(valarray<double>(gpadc[slice(0, 2 * Npad2, 2)]).min()) <
            SQRT_DOUBLE_EPS * gpadc.max());
    assert(fabs(valarray<double>(gpadc[slice(0, 2 * Npad2, 2)]).max()) <
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
