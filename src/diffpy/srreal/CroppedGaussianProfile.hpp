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
* class CroppedGaussianProfile -- Gaussian profile cropped to zero beyond
*     xboundhi and scaled so that its integrated area equals 1.
*     Registered as "croppedgaussian".
*
*****************************************************************************/

#ifndef CROPPEDGAUSSIANPROFILE_HPP_INCLUDED
#define CROPPEDGAUSSIANPROFILE_HPP_INCLUDED

#include <diffpy/srreal/GaussianProfile.hpp>

namespace diffpy {
namespace srreal {

class CroppedGaussianProfile : public GaussianProfile
{
    public:

        // constructors
        CroppedGaussianProfile();
        PeakProfilePtr create() const;
        PeakProfilePtr clone() const;

        // methods
        const std::string& type() const;
        double operator()(double x, double fwhm) const;
        void setPrecision(double eps);

    private:

        // data
        double mscale;

        // serialization
        friend class boost::serialization::access;

        template<class Archive>
            void serialize(Archive& ar, const unsigned int version)
        {
            using boost::serialization::base_object;
            ar & base_object<GaussianProfile>(*this);
            ar & mscale;
        }

};

}   // namespace srreal
}   // namespace diffpy

// Serialization -------------------------------------------------------------

BOOST_CLASS_EXPORT_KEY(diffpy::srreal::CroppedGaussianProfile)

#endif  // CROPPEDGAUSSIANPROFILE_HPP_INCLUDED
