/*****************************************************************************
*
* libdiffpy         by DANSE Diffraction group
*                   Simon J. L. Billinge
*                   (c) 2009 The Trustees of Columbia University
*                   in the City of New York.  All rights reserved.
*
* File coded by:    Christopher Farrow, Pavol Juhas
*
* See AUTHORS.txt for a list of people who contributed.
* See LICENSE_DANSE.txt for license information.
*
******************************************************************************
*
* class DebyeWallerPeakWidth -- peak width model based on
*      I.-K. Jeong, et al., Phys. Rev. B 67, 104301 (2003)
*      http://link.aps.org/doi/10.1103/PhysRevB.67.104301
*
*****************************************************************************/

#ifndef JEONGPEAKWIDTH_HPP_INCLUDED
#define JEONGPEAKWIDTH_HPP_INCLUDED

#include <diffpy/srreal/DebyeWallerPeakWidth.hpp>

namespace diffpy {
namespace srreal {


class JeongPeakWidth : public DebyeWallerPeakWidth
{
    public:

        // constructors
        JeongPeakWidth();
        virtual PeakWidthModelPtr create() const;
        virtual PeakWidthModelPtr clone() const;

        // methods
        virtual const std::string& type() const;
        virtual double calculate(const BaseBondGenerator&) const;
        virtual double maxWidth(StructureAdapterPtr,
                double rmin, double rmax) const;

        // data access
        const double& getDelta1() const;
        void setDelta1(double);
        const double& getDelta2() const;
        void setDelta2(double);
        const double& getQbroad() const;
        void setQbroad(double);
        const double& getQbroad_seperable() const;
        void setQbroad_seperable(double);

    private:

        // data
        double mdelta1;
        double mdelta2;
        double mqbroad;
        double mqbroad_seperable;

        // methods
        double msdSharpeningRatio(const double& r) const;

        // serialization
        friend class boost::serialization::access;

        template<class Archive>
            void serialize(Archive& ar, const unsigned int version)
        {
            using boost::serialization::base_object;
            ar & base_object<DebyeWallerPeakWidth>(*this);
            ar & mdelta1;
            ar & mdelta2;
            ar & mqbroad;
            ar & mqbroad_seperable;
        }

};


}   // namespace srreal
}   // namespace diffpy

// Serialization -------------------------------------------------------------

BOOST_CLASS_EXPORT_KEY(diffpy::srreal::JeongPeakWidth)

#endif  // JEONGPEAKWIDTH_HPP_INCLUDED
