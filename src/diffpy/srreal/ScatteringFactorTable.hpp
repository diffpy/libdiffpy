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
* class ScatteringFactorTable -- base class for looking up scattering factors

* class ScatteringFactorTableOwner -- to be used as a base class for classes
*     that own ScatteringFactorTable
*
* $Id$
*
*****************************************************************************/

#ifndef SCATTERINGFACTORTABLE_HPP_INCLUDED
#define SCATTERINGFACTORTABLE_HPP_INCLUDED

#include <diffpy/HasClassRegistry.hpp>

namespace diffpy {
namespace srreal {

class ScatteringFactorTable :
    public diffpy::HasClassRegistry<ScatteringFactorTable>
{
    public:

        // methods
        virtual const std::string& radiationType() const = 0;
        const double& lookup(const std::string& smbl) const;
        virtual const double& lookupatq(const std::string&, double) const;
        void setCustom(const std::string& smbl, double value);
        void resetCustom(const std::string& smbl);
        void resetAll();
        virtual double fetch(const std::string& smbl) const = 0;

    protected:

        mutable std::map<std::string,double> mtable;
};

typedef ScatteringFactorTable::SharedPtr ScatteringFactorTablePtr;

class ScatteringFactorTableOwner
{
    public:

        // access and configuration of scattering factors
        void setScatteringFactorTable(ScatteringFactorTablePtr);
        void setScatteringFactorTableByType(const std::string& tp);
        ScatteringFactorTablePtr getScatteringFactorTable();
        const ScatteringFactorTablePtr getScatteringFactorTable() const;
        const std::string& getRadiationType() const;


    private:

        // data
        ScatteringFactorTablePtr msftable;
};

}   // namespace srreal
}   // namespace diffpy

#endif  // SCATTERINGFACTORTABLE_HPP_INCLUDED
