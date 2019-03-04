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
* class ScatteringFactorTable -- base class for looking up scattering factors
*
* class ScatteringFactorTableOwner -- to be used as a base class for classes
*     that own ScatteringFactorTable
*
*****************************************************************************/

#ifndef SCATTERINGFACTORTABLE_HPP_INCLUDED
#define SCATTERINGFACTORTABLE_HPP_INCLUDED

#include <unordered_set>

#include <boost/serialization/base_object.hpp>
#include <boost/serialization/assume_abstract.hpp>
#include <boost/serialization/shared_ptr.hpp>
#include <boost/serialization/unordered_map.hpp>
#include <boost/serialization/export.hpp>

#include <diffpy/HasClassRegistry.hpp>
#include <diffpy/EventTicker.hpp>

namespace diffpy {
namespace srreal {

class ScatteringFactorTable :
    public diffpy::HasClassRegistry<ScatteringFactorTable>
{
    public:

        // HasClassRegistry override
        bool registerThisType() const;

        // own methods
        virtual const std::string& radiationType() const = 0;
        double lookup(const std::string& smbl, double q=0.0) const;
        virtual double standardLookup(const std::string&, double) const = 0;
        void setCustomAs(const std::string& smbl, const std::string& srcsmbl);
        void setCustomAs(const std::string& smbl, const std::string& srcsmbl,
                double value, double q=0.0);
        void resetCustom(const std::string& smbl);
        void resetAll();
        std::unordered_set<std::string> getCustomSymbols() const;

        typedef std::unordered_map<std::string,
                std::pair<std::string, double> > CustomDataStorage;
        virtual eventticker::EventTicker& ticker() const  { return mticker; }

    protected:

        // data
        CustomDataStorage mcustom;
        mutable eventticker::EventTicker mticker;

    private:

        // serialization
        friend class boost::serialization::access;
        template<class Archive>
            void serialize(Archive& ar, const unsigned int version)
        {
            ar & mcustom & mticker;
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
        eventticker::EventTicker& ticker() const;

    private:

        // data
        ScatteringFactorTablePtr msftable;
        mutable eventticker::EventTicker mprivateticker;

        // serialization
        friend class boost::serialization::access;
        template<class Archive>
            void serialize(Archive& ar, const unsigned int version)
        {
            ar & msftable & mprivateticker;
        }
};

}   // namespace srreal
}   // namespace diffpy

// Serialization -------------------------------------------------------------

BOOST_SERIALIZATION_ASSUME_ABSTRACT(diffpy::srreal::ScatteringFactorTable)

#endif  // SCATTERINGFACTORTABLE_HPP_INCLUDED
