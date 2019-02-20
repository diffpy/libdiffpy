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
* class LinearBaseline -- linear PDF baseline
*
*****************************************************************************/

#ifndef LINEARBASELINE_HPP_INCLUDED
#define LINEARBASELINE_HPP_INCLUDED

#include <diffpy/srreal/PDFBaseline.hpp>

namespace diffpy {
namespace srreal {

/// @class LinearBaseline
/// @brief linear PDF baseline

class LinearBaseline : public PDFBaseline
{
    public:

        // constructors
        LinearBaseline();
        PDFBaselinePtr create() const;
        PDFBaselinePtr clone() const;

        // methods
        const std::string& type() const;
        double operator()(const double& r) const;
        void setSlope(double sc);
        const double& getSlope() const;

    private:

        // data
        double mslope;

        // serialization
        friend class boost::serialization::access;

        template<class Archive>
            void serialize(Archive& ar, const unsigned int version)
        {
            using boost::serialization::base_object;
            ar & base_object<PDFBaseline>(*this);
            ar & mslope;
        }

};  // class LinearBaseline

}   // namespace srreal
}   // namespace diffpy

// Serialization -------------------------------------------------------------

BOOST_CLASS_EXPORT_KEY(diffpy::srreal::LinearBaseline)

#endif  // LINEARBASELINE_HPP_INCLUDED
