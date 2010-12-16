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
#include <boost/serialization/map.hpp>
#include <boost/serialization/split_free.hpp>

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

namespace boost {
namespace serialization {

template<class Archive>
void save(Archive& ar,
        const diffpy::srreal::ScatteringFactorTablePtr& ptr,
        unsigned int version)
{
    std::string tp;
    std::map<std::string,double> dt;
    if (ptr.get())
    {
        tp = ptr->type();
        dt = ptr->getAllCustom();
    }
    ar & tp & dt;
}


template<class Archive>
void load(Archive& ar,
        diffpy::srreal::ScatteringFactorTablePtr& ptr,
        unsigned int version)
{
    using namespace diffpy::srreal;
    std::string tp;
    std::map<std::string,double> dt;
    ar & tp & dt;
    if (!tp.empty())
    {
        ptr = ScatteringFactorTable::createByType(tp);
        std::map<std::string,double>::const_iterator kv;
        for (kv = dt.begin(); kv != dt.end(); ++kv)
        {
            ptr->setCustom(kv->first, kv->second);
        }
    }
    else  ptr.reset();
}

}   // namespace serialization
}   // namespace boost

BOOST_SERIALIZATION_SPLIT_FREE(diffpy::srreal::ScatteringFactorTablePtr)

#endif  // SCATTERINGFACTORTABLE_HPP_INCLUDED
