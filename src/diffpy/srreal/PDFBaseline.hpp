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
* class PDFBaseline -- abstract base class for PDF baseline functions
*     A concrete instance of PDFBaseline is a functor, that calculates
*     baseline value at a given pair distance r.  The baseline is added to
*     (R(r) * r) before multiplication by any envelope functions.
*
* $Id$
*
*****************************************************************************/

#ifndef PDFBASELINE_HPP_INCLUDED
#define PDFBASELINE_HPP_INCLUDED

#include <string>
#include <set>

#include <diffpy/Attributes.hpp>
#include <diffpy/HasClassRegistry.hpp>

namespace diffpy {
namespace srreal {

/// @class PDFBaseline
/// @brief abstract base class for PDF baseline function

class PDFBaseline :
    public Attributes,
    public HasClassRegistry<PDFBaseline>
{
    public:
        virtual double operator()(const double& r) const = 0;
};

typedef PDFBaseline::SharedPtr PDFBaselinePtr;

}   // namespace srreal
}   // namespace diffpy

#endif  // PDFBASELINE_HPP_INCLUDED
