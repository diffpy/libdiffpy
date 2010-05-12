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
* Implementation of shared template functions for PDF calculator.
*     bandPassFilter
*
* $Id$
*
*****************************************************************************/

#ifndef PDFUTILS_IPP_INCLUDED
#define PDFUTILS_IPP_INCLUDED

namespace diffpy {
namespace srreal {

template <class Ti>
void bandPassFilter(Ti first, Ti last, double dr, double qmin, double qmax)
{
    if (!(first < last))    return;
    int datalen = last - first;
    // pad data with the same number of zeros up to the next power of 2
    int padlen = (int) pow(2, int(ceil(log2(datalen) + 1)));
    // ycpad is complex, so it needs to be twice as long
    std::valarray<double> ycpa(0.0, 2 * padlen);
    double* ycfirst = &(ycpa[0]);
    double* yci = ycfirst;
    for (Ti p = first; p != last; ++p, yci += 2)  { *yci = *p; }
    // perform the filtering
    bandPassFilterCValarray(ycpa, dr, qmin, qmax);
    // copy real components back to the input sequence
    yci = &(ycpa[0]);
    for (Ti p = first; p != last; ++p, yci += 2)  { *p = *yci; }
}

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
    int rv = int(ceil(qmin / dq));
    if (eps_eq(qmin, (rv - 1) * dq))  rv -= 1;
    return rv;
}


template <class T>
int pdfutils_qmaxSteps(const T* pdfc)
{
    int rv = int(ceil(pdfc->getQmax() / pdfc->getQstep()));
    return rv;
}


template <class T>
QuantityType pdfutils_getRgrid(const T* pdfc)
{
    QuantityType rv;
    int ndrmin = pdfutils_rminSteps(pdfc);
    int ndrmax = pdfutils_rmaxSteps(pdfc);
    const double dr = pdfc->getRstep();
    rv.reserve(ndrmax - ndrmin);
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
