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
* See LICENSE.txt for license information.
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
#include <boost/serialization/map.hpp>
#include <boost/serialization/split_free.hpp>

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
};

typedef PeakProfile::SharedPtr PeakProfilePtr;

}   // namespace srreal
}   // namespace diffpy

// Serialization -------------------------------------------------------------

BOOST_SERIALIZATION_SPLIT_FREE(diffpy::srreal::PeakProfilePtr)

namespace boost {
namespace serialization {

template<class Archive>
void save(Archive& ar,
        const diffpy::srreal::PeakProfilePtr& ptr,
        const unsigned int version)
{
    using namespace diffpy::attributes;
    std::string tp;
    AttributesDataMap dt;
    diffpy::eventticker::EventTicker tc;
    if (ptr.get())
    {
        tp = ptr->type();
        dt = saveAttributesData(*ptr);
        tc = ptr->ticker();
    }
    ar & tp & dt & tc;
}


template<class Archive>
void load(Archive& ar,
        diffpy::srreal::PeakProfilePtr& ptr,
        const unsigned int version)
{
    using namespace diffpy::attributes;
    using namespace diffpy::srreal;
    std::string tp;
    AttributesDataMap dt;
    diffpy::eventticker::EventTicker tc;
    ar & tp & dt & tc;
    if (!tp.empty())
    {
        ptr = PeakProfile::createByType(tp);
        loadAttributesData(*ptr, dt);
        ptr->ticker() = tc;
    }
    else  ptr.reset();
}

}   // namespace serialization
}   // namespace boost

#endif  // PEAKPROFILE_HPP_INCLUDED
