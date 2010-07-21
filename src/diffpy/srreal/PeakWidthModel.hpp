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
* class PeakWidthModel -- base class for calculation of peak widths.
*     The calculate function takes a BondGenerator instance and
*     returns full width at half maximum, based on peak model parameters
*     and anisotropic displacement parameters of atoms in the pair.
*
* class PeakWidthModelOwner -- to be used as a base class for classes
*     that own PeakWidthModel
*
* $Id$
*
*****************************************************************************/

#ifndef PEAKWIDTHMODEL_HPP_INCLUDED
#define PEAKWIDTHMODEL_HPP_INCLUDED

#include <diffpy/Attributes.hpp>
#include <diffpy/HasClassRegistry.hpp>
#include <diffpy/srreal/BaseBondGenerator.hpp>

namespace diffpy {
namespace srreal {

class PeakWidthModel :
    public diffpy::Attributes,
    public diffpy::HasClassRegistry<PeakWidthModel>
{
    public:

        // methods
        virtual double calculate(const BaseBondGenerator&) const = 0;
        virtual double calculateFromMSD(double msdval) const = 0;

};


typedef PeakWidthModel::SharedPtr PeakWidthModelPtr;


class PeakWidthModelOwner
{
    public:

        // PDF peak width configuration
        void setPeakWidthModel(PeakWidthModelPtr);
        void setPeakWidthModelByType(const std::string& tp);
        PeakWidthModelPtr getPeakWidthModel();
        const PeakWidthModelPtr getPeakWidthModel() const;

    private:

        // data
        PeakWidthModelPtr mpwmodel;
};

}   // namespace srreal
}   // namespace diffpy

#endif  // PEAKWIDTHMODEL_HPP_INCLUDED
