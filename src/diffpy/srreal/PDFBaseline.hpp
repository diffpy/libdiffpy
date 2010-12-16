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
#include <boost/serialization/map.hpp>
#include <boost/serialization/split_free.hpp>

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

// Serialization -------------------------------------------------------------

namespace boost {
namespace serialization {

template<class Archive>
void save(Archive& ar,
        const diffpy::srreal::PDFBaselinePtr& ptr, unsigned int version)
{
    using namespace diffpy::attributes;
    std::string tp;
    AttributesDataMap dt;
    if (ptr.get())
    {
        tp = ptr->type();
        dt = saveAttributesData(*ptr);
    }
    ar & tp & dt;
}


template<class Archive>
void load(Archive& ar,
        diffpy::srreal::PDFBaselinePtr& ptr, unsigned int version)
{
    using namespace diffpy::attributes;
    using namespace diffpy::srreal;
    std::string tp;
    AttributesDataMap dt;
    ar & tp & dt;
    if (!tp.empty())
    {
        ptr = PDFBaseline::createByType(tp);
        loadAttributesData(*ptr, dt);
    }
    else  ptr.reset();
}

}   // namespace serialization
}   // namespace boost

BOOST_SERIALIZATION_SPLIT_FREE(diffpy::srreal::PDFBaselinePtr)

#endif  // PDFBASELINE_HPP_INCLUDED
