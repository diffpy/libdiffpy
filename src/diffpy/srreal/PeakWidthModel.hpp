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

#include <string>
#include <set>
#include <boost/shared_ptr.hpp>

#include <diffpy/Attributes.hpp>

namespace diffpy {
namespace srreal {

class BaseBondGenerator;

class PeakWidthModel : public diffpy::Attributes
{
    public:

        // constructors
        virtual boost::shared_ptr<PeakWidthModel> create() const = 0;
        virtual boost::shared_ptr<PeakWidthModel> clone() const = 0;
        virtual ~PeakWidthModel()  { }

        // methods
        virtual const std::string& type() const = 0;
        virtual double calculate(const BaseBondGenerator&) const = 0;
        virtual double calculateFromMSD(double msdval) const = 0;

        // comparison with derived classes
        virtual bool operator==(const PeakWidthModel&) const = 0;
};


class PeakWidthModelOwner
{
    public:

        // PDF peak width configuration
        void setPeakWidthModel(const PeakWidthModel&);
        void setPeakWidthModel(const std::string& tp);
        PeakWidthModel& getPeakWidthModel();
        const PeakWidthModel& getPeakWidthModel() const;

    private:

        // data
        boost::shared_ptr<PeakWidthModel> mpwmodel;
};

// Factory functions for Peak Width Models -----------------------------------

boost::shared_ptr<PeakWidthModel> createPeakWidthModel(const std::string& tp);
bool registerPeakWidthModel(const PeakWidthModel&);
bool aliasPeakWidthModel(const std::string& tp, const std::string& al);
std::set<std::string> getPeakWidthModelTypes();

}   // namespace srreal
}   // namespace diffpy

#endif  // PEAKWIDTHMODEL_HPP_INCLUDED
