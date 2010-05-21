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
*     maxUii
*     fftftog  and  fftgtof
*
* $Id$
*
*****************************************************************************/

#ifndef PDFUTILS_HPP_INCLUDED
#define PDFUTILS_HPP_INCLUDED

#include <cmath>
#include <valarray>
#include <diffpy/srreal/R3linalg.hpp>
#include <diffpy/srreal/StructureAdapter.hpp>
#include <diffpy/srreal/PairQuantityUtils.hpp>

namespace diffpy {
namespace srreal {

const double DEFAULT_PDFCALCULATOR_RMAX = 10.0;
const double DEFAULT_PDFCALCULATOR_RSTEP = 0.01;
const double DEFAULT_PDFCALCULATOR_MAXEXTENSION = 10.0;
/// Default peak precision was obtained from the tunePeakPrecision.py script
/// and it was tuned to give average zero slope in the difference curve
/// between pdffit2 and PDFCalculator results.
const double DEFAULT_PEAKPRECISION = 3.33e-6;

const double DEFAULT_QGRID_QMAX = 10.0;
const double DEFAULT_QGRID_QSTEP = 0.05;

/// fast Fourier transformation converting G(r) to F(Q)
QuantityType fftgtof(const QuantityType& g, double rstep, double rmin=0.0);

/// fast Fourier transformation converting F(Q) to G(r)
QuantityType fftftog(const QuantityType& f, double qstep, double qmin=0.0);

/// shared methods for PDFCalculator and DebyePDFCalculator
template <class T> QuantityType pdfutils_getQgrid(const T* pdfc);
template <class T> int pdfutils_qminSteps(const T* pdfc);
template <class T> int pdfutils_qmaxSteps(const T* pdfc);
template <class T> QuantityType pdfutils_getRgrid(const T* pdfc);
template <class T> int pdfutils_rminSteps(const T* pdfc);
template <class T> int pdfutils_rmaxSteps(const T* pdfc);

}   // namespace srreal
}   // namespace diffpy

// Implementation ------------------------------------------------------------

#include <diffpy/srreal/PDFUtils.ipp>

#endif  // PDFUTILS_HPP_INCLUDED
