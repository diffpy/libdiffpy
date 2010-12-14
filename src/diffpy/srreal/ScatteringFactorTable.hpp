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

#include <boost/serialization/base_object.hpp>
#include <boost/serialization/assume_abstract.hpp>
#include <boost/serialization/shared_ptr.hpp>
#include <boost/serialization/export.hpp>
#include <boost/serialization/map.hpp>
#include <boost/serialization/set.hpp>

#include <diffpy/HasClassRegistry.hpp>

namespace diffpy {
namespace srreal {

class ScatteringFactorTable :
    public diffpy::HasClassRegistry<ScatteringFactorTable>
{
    public:

        // methods
        virtual const std::string& radiationType() const = 0;
        double lookup(const std::string& smbl) const;
        virtual double lookupatq(const std::string&, double) const = 0;
        void setCustom(const std::string& smbl, double value);
        void resetCustom(const std::string& smbl);
        std::map<std::string,double> getAllCustom() const;
        void resetAll();

    protected:

        mutable std::map<std::string,double> mtable;
        std::set<std::string> mcustomsymbols;

    private:

        // serialization
        friend class boost::serialization::access;
        template<class Archive>
            void serialize(Archive& ar, const unsigned int version)
        {
            ar & mtable;
            ar & mcustomsymbols;
        }
};

typedef ScatteringFactorTable::SharedPtr ScatteringFactorTablePtr;

class ScatteringFactorTableOwner
{
    public:

        // access and configuration of scattering factors
        void setScatteringFactorTable(ScatteringFactorTablePtr);
        void setScatteringFactorTableByType(const std::string& tp);
        ScatteringFactorTablePtr& getScatteringFactorTable();
        const ScatteringFactorTablePtr& getScatteringFactorTable() const;
        const std::string& getRadiationType() const;


    private:

        // data
        ScatteringFactorTablePtr msftable;

        // serialization
        friend class boost::serialization::access;
        template<class Archive>
            void serialize(Archive& ar, const unsigned int version)
        {
            ar & msftable;
        }
};

}   // namespace srreal
}   // namespace diffpy

// Serialization -------------------------------------------------------------

BOOST_SERIALIZATION_ASSUME_ABSTRACT(diffpy::srreal::ScatteringFactorTable)
BOOST_SERIALIZATION_ASSUME_ABSTRACT(diffpy::srreal::ScatteringFactorTableOwner)

#endif  // SCATTERINGFACTORTABLE_HPP_INCLUDED
