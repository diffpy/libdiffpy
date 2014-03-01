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
* Various common routines useful for PDF calculation.
* Implementation of template functions.
*
*****************************************************************************/

#ifndef PDFUTILS_IPP_INCLUDED
#define PDFUTILS_IPP_INCLUDED

namespace diffpy {
namespace srreal {

// Shared implementations of DebyePDFCalculator and PDFCalculator methods

template <class T>
QuantityType pdfutils_getQgrid(const T* pdfc)
{
    const int npts = pdfutils_qmaxSteps(pdfc);
    const double dq = pdfc->getQstep();
    QuantityType rv;
    rv.reserve(npts);
    for (int kq = 0; kq < npts; ++kq)  rv.push_back(kq * dq);
    return rv;
}


inline
int pdfutils_qminSteps(const double& qmin, const double& qstep)
{
    using diffpy::mathutils::eps_eq;
    if (qstep <= 0.0)  return 0;
    int rv = int(ceil(qmin / qstep));
    if (eps_eq(qmin, (rv - 1) * qstep))  rv -= 1;
    return rv;
}


template <class T>
int pdfutils_qminSteps(const T* pdfc)
{
    const double qmin = pdfc->getQmin();
    const double qstep = pdfc->getQstep();
    int rv = pdfutils_qminSteps(qmin, qstep);
    return rv;
}


template <class T>
int pdfutils_qmaxSteps(const T* pdfc)
{
    const double dq = pdfc->getQstep();
    if (dq <= 0.0)  return 0;
    int rv = int(ceil(pdfc->getQmax() / dq));
    return rv;
}


template <class T>
QuantityType pdfutils_getRgrid(const T* pdfc)
{
    QuantityType rv;
    int ndrmin = pdfutils_rminSteps(pdfc);
    int ndrmax = pdfutils_rmaxSteps(pdfc);
    const double dr = pdfc->getRstep();
    if (ndrmax > ndrmin)    rv.reserve(ndrmax - ndrmin);
    for (int ndr = ndrmin; ndr < ndrmax; ++ndr)
    {
        rv.push_back(ndr * dr);
    }
    return rv;
}


inline
int pdfutils_rminSteps(const double& rmin, const double& rstep)
{
    using diffpy::mathutils::eps_eq;
    int rv = int(ceil(rmin / rstep));
    if (eps_eq(rmin, (rv - 1) * rstep))  rv -= 1;
    return rv;
}


template <class T>
int pdfutils_rminSteps(const T* pdfc)
{
    const double& rmin = pdfc->getRmin();
    const double& rstep = pdfc->getRstep();
    int rv = pdfutils_rminSteps(rmin, rstep);
    return rv;
}


inline
int pdfutils_rmaxSteps(const double& rmax, const double& rstep)
{
    int rv = int(ceil(rmax / rstep));
    return rv;
}


template <class T>
int pdfutils_rmaxSteps(const T* pdfc)
{
    const double& rmax = pdfc->getRmax();
    const double& rstep = pdfc->getRstep();
    int rv = pdfutils_rmaxSteps(rmax, rstep);
    return rv;
}


}   // namespace srreal
}   // namespace diffpy

// vim:ft=cpp:

#endif  // PDFUTILS_IPP_INCLUDED
