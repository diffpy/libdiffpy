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
const double DEFAULT_QGRID_QMAX = 10.0;
const double DEFAULT_QGRID_QSTEP = 0.05;

/// Calculate MSD along specified direction in Cartesian space.
double meanSquareDisplacement(const R3::Matrix& Uijcartn, const R3::Vector& s,
        bool anisotropy=true);

/// Maximum diagonal Uii element from all atoms in the structure.
double maxUii(const StructureAdapter* stru);

/// Apply band pass filter to a sequence of doubles
template <class Ti>
void bandPassFilter(Ti first, Ti last, double dr, double qmin, double qmax);

/// Implementation of bandPassFilter using padded complex valarray
void bandPassFilterCValarray(std::valarray<double>& ycpa,
        double dr, double qmin, double qmax);

/// fast Fourier transformation converting G(r) to F(Q)
QuantityType fftgtof(const QuantityType& g, double rstep, double rmin=0.0);
/// fast Fourier transformation converting F(Q) to G(r)
QuantityType fftftog(const QuantityType& f, double qstep, double qmin=0.0);

}   // namespace srreal
}   // namespace diffpy

// Implementation ------------------------------------------------------------

#include <diffpy/srreal/PDFUtils.ipp>

#endif  // PDFUTILS_HPP_INCLUDED
