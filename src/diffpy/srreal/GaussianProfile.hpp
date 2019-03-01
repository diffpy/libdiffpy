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
* class GaussianProfile -- concrete implementation of the PeakProfile class.
*     GaussianProfile is registered as "gaussian".
*
*****************************************************************************/

#ifndef GAUSSIANPROFILE_HPP_INCLUDED
#define GAUSSIANPROFILE_HPP_INCLUDED

#include <diffpy/srreal/PeakProfile.hpp>

namespace diffpy {
namespace srreal {

class GaussianProfile : public PeakProfile
{
    public:

        // constructors
        GaussianProfile();
        PeakProfilePtr create() const;
        PeakProfilePtr clone() const;

        // methods
        const std::string& type() const;
        double operator()(double x, double fwhm) const;
        double xboundlo(double fwhm) const;
        double xboundhi(double fwhm) const;
        void setPrecision(double eps);

    protected:

        // data
        double mhalfboundrel;

    private:

        // serialization
        friend class boost::serialization::access;

        template<class Archive>
            void serialize(Archive& ar, const unsigned int version)
        {
            using boost::serialization::base_object;
            ar & base_object<PeakProfile>(*this);
            ar & mhalfboundrel;
        }

};

}   // namespace srreal
}   // namespace diffpy

// Serialization -------------------------------------------------------------

BOOST_CLASS_EXPORT_KEY(diffpy::srreal::GaussianProfile)

#endif  // GAUSSIANPROFILE_HPP_INCLUDED
