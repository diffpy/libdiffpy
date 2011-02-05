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
* Implementation of template functions.
*
* $Id$
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


template <class T>
int pdfutils_qminSteps(const T* pdfc)
{
    using diffpy::mathutils::eps_eq;
    const double qmin = pdfc->getQmin();
    const double dq = pdfc->getQstep();
    if (dq <= 0.0)  return 0;
    int rv = int(ceil(qmin / dq));
    if (eps_eq(qmin, (rv - 1) * dq))  rv -= 1;
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


template <class T>
int pdfutils_rminSteps(const T* pdfc)
{
    using diffpy::mathutils::eps_eq;
    const double rmin = pdfc->getRmin();
    const double dr = pdfc->getRstep();
    int rv = int(ceil(rmin / dr));
    if (eps_eq(rmin, (rv - 1) * dr))  rv -= 1;
    return rv;
}


template <class T>
int pdfutils_rmaxSteps(const T* pdfc)
{
    int rv = int(ceil(pdfc->getRmax() / pdfc->getRstep()));
    return rv;
}


}   // namespace srreal
}   // namespace diffpy

// vim:ft=cpp:

#endif  // PDFUTILS_IPP_INCLUDED
