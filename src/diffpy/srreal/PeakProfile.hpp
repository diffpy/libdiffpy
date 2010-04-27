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
* class PeakProfile -- base class for calculation of peak profiles.
*     When possible total integrated area of the profile should eqaul 1.
*     The method yvalue(x, fwhm) returns amplitude of a zero-centered profile.
*     Methods xboundlo(fwhm), xboundhi(fwhm) return low and high x-boundaries,
*     where amplitude relative to the maximum becomes smaller than precision
*     set by setPrecision().
*
* $Id$
*
*****************************************************************************/

#ifndef PEAKPROFILE_HPP_INCLUDED
#define PEAKPROFILE_HPP_INCLUDED

#include <diffpy/Attributes.hpp>
#include <diffpy/HasClassRegistry.hpp>

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
        virtual double yvalue(double x, double fwhm) const = 0;
        virtual double xboundlo(double fwhm) const = 0;
        virtual double xboundhi(double fwhm) const = 0;
        virtual void setPrecision(double eps);
        const double& getPrecision() const;

    private:

        // data
        double mprecision;
};

typedef PeakProfile::SharedPtr PeakProfilePtr;

}   // namespace srreal
}   // namespace diffpy

#endif  // PEAKPROFILE_HPP_INCLUDED
