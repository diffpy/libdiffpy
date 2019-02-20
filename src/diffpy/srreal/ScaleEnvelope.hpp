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
* class ScaleEnvelope -- constant scaling factor
*
*****************************************************************************/

#ifndef SCALEENVELOPE_HPP_INCLUDED
#define SCALEENVELOPE_HPP_INCLUDED

#include <diffpy/srreal/PDFEnvelope.hpp>

namespace diffpy {
namespace srreal {

/// @class ScaleEnvelope
/// @brief constant PDF scaling factor

class ScaleEnvelope : public PDFEnvelope
{
    public:

        // constructors
        ScaleEnvelope();
        PDFEnvelopePtr create() const;
        PDFEnvelopePtr clone() const;

        // methods
        const std::string& type() const;
        double operator()(const double& r) const;
        void setScale(double sc);
        const double& getScale() const;

    private:

        // data
        double mscale;

        // serialization
        friend class boost::serialization::access;

        template<class Archive>
            void serialize(Archive& ar, const unsigned int version)
        {
            using boost::serialization::base_object;
            ar & base_object<PDFEnvelope>(*this);
            ar & mscale;
        }

};  // class ScaleEnvelope

}   // namespace srreal
}   // namespace diffpy

// Serialization -------------------------------------------------------------

BOOST_CLASS_EXPORT_KEY(diffpy::srreal::ScaleEnvelope)

#endif  // SCALEENVELOPE_HPP_INCLUDED
