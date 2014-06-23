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

#include <boost/unordered_set.hpp>
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/split_free.hpp>

#include <diffpy/boostextensions/serialize_unordered_map.hpp>
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
        boost::unordered_set<std::string> getCustomSymbols() const;

        typedef boost::unordered_map<std::string,
                std::pair<std::string, double> > CustomDataStorage;
        virtual eventticker::EventTicker& ticker() const  { return mticker; }

    protected:

        // data
        CustomDataStorage mcustom;
        mutable eventticker::EventTicker mticker;

        // serialization helpers for accessing mcustom
        friend const CustomDataStorage& getsftcustomdata(const SharedPtr&);
        friend void setsftcustomdata(SharedPtr&, const CustomDataStorage&);
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

// serialization helpers for accessing ScatteringFactorTable data

inline
const ScatteringFactorTable::CustomDataStorage&
getsftcustomdata(const ScatteringFactorTablePtr& ptr)
{
    return ptr->mcustom;
}


inline
void setsftcustomdata(ScatteringFactorTablePtr& ptr,
        const ScatteringFactorTable::CustomDataStorage& dt)
{
    ptr->mcustom = dt;
}

}   // namespace srreal
}   // namespace diffpy

// Serialization -------------------------------------------------------------

namespace boost {
namespace serialization {

template<class Archive>
void save(Archive& ar,
        const diffpy::srreal::ScatteringFactorTablePtr& ptr,
        const unsigned int version)
{
    using namespace diffpy::srreal;
    std::string tp;
    ScatteringFactorTable::CustomDataStorage dt;
    diffpy::eventticker::EventTicker tc;
    if (ptr.get())
    {
        tp = ptr->type();
        dt = getsftcustomdata(ptr);
        tc = ptr->ticker();
    }
    ar & tp & dt & tc;
}


template<class Archive>
void load(Archive& ar,
        diffpy::srreal::ScatteringFactorTablePtr& ptr,
        const unsigned int version)
{
    using namespace diffpy::srreal;
    std::string tp;
    ScatteringFactorTable::CustomDataStorage dt;
    diffpy::eventticker::EventTicker tc;
    ar & tp & dt & tc;
    if (!tp.empty())
    {
        ptr = ScatteringFactorTable::createByType(tp);
        setsftcustomdata(ptr, dt);
        ptr->ticker() = tc;
    }
    else  ptr.reset();
}

}   // namespace serialization
}   // namespace boost

BOOST_SERIALIZATION_SPLIT_FREE(diffpy::srreal::ScatteringFactorTablePtr)

#endif  // SCATTERINGFACTORTABLE_HPP_INCLUDED
