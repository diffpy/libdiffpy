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
* class PeakProfile -- base class for calculation of peak profiles.
*     When possible total integrated area of the profile should eqaul 1.
*     The operator()(x, fwhm) returns amplitude of a zero-centered profile.
*     Methods xboundlo(fwhm), xboundhi(fwhm) return low and high x-boundaries,
*     where amplitude relative to the maximum becomes smaller than precision
*     set by setPrecision().
*
*****************************************************************************/

#ifndef PEAKPROFILE_HPP_INCLUDED
#define PEAKPROFILE_HPP_INCLUDED

#include <boost/serialization/base_object.hpp>
#include <boost/serialization/assume_abstract.hpp>
#include <boost/serialization/export.hpp>
#include <boost/serialization/shared_ptr.hpp>

#include <diffpy/Attributes.hpp>
#include <diffpy/HasClassRegistry.hpp>
#include <diffpy/EventTicker.hpp>

namespace diffpy {
namespace srreal {

class PeakProfile :
    public diffpy::Attributes,
    public diffpy::HasClassRegistry<PeakProfile>
{
    public:

        // constructors
        PeakProfile();
        // methods
        virtual double operator()(double x, double fwhm) const = 0;
        virtual double xboundlo(double fwhm) const = 0;
        virtual double xboundhi(double fwhm) const = 0;
        virtual void setPrecision(double eps);
        const double& getPrecision() const;
        virtual eventticker::EventTicker& ticker() const  { return mticker; }

    protected:

        // data
        mutable eventticker::EventTicker mticker;

    private:

        // data
        double mprecision;

        // serialization
        friend class boost::serialization::access;
        template<class Archive>
            void serialize(Archive& ar, const unsigned int version)
        {
            ar & mticker & mprecision;
        }

};

typedef PeakProfile::SharedPtr PeakProfilePtr;

}   // namespace srreal
}   // namespace diffpy

// Serialization -------------------------------------------------------------

BOOST_SERIALIZATION_ASSUME_ABSTRACT(diffpy::srreal::PeakProfile)

#endif  // PEAKPROFILE_HPP_INCLUDED
