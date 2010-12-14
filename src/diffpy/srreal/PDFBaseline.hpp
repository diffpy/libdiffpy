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
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/assume_abstract.hpp>
#include <boost/serialization/shared_ptr.hpp>
#include <boost/serialization/export.hpp>

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

    private:
        // serialization
        friend class boost::serialization::access;
        template<class Archive>
            void serialize(Archive& ar, const unsigned int version)  { }
};

typedef PDFBaseline::SharedPtr PDFBaselinePtr;

}   // namespace srreal
}   // namespace diffpy

// Serialization -------------------------------------------------------------

BOOST_SERIALIZATION_ASSUME_ABSTRACT(diffpy::srreal::PDFBaseline)

#endif  // PDFBASELINE_HPP_INCLUDED
