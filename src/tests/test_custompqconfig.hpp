/*****************************************************************************
*
* libdiffpy         Complex Modeling Initiative
*                   Pavol Juhas
*                   (c) 2014 Brookhaven National Laboratory,
*                   Upton, New York.  All rights reserved.
*
* File coded by:    Pavol Juhas
*
* See AUTHORS.txt for a list of people who contributed.
* See LICENSE.txt for license information.
*
******************************************************************************
*
* Helper class for testing the effect of customPQConfig.
*
*****************************************************************************/

#ifndef TEST_CUSTOMPQCONFIG_HPP_INCLUDED
#define TEST_CUSTOMPQCONFIG_HPP_INCLUDED

#include <boost/serialization/map.hpp>
#include <diffpy/srreal/PeriodicStructureAdapter.hpp>

using diffpy::srreal::StructureAdapterPtr;
using diffpy::srreal::PeriodicStructureAdapter;
using diffpy::srreal::PeriodicStructureAdapterPtr;

// ---------------------------------------------------------------------------

class PeriodicAdapterWithPQConfig :
    public PeriodicStructureAdapter
{
    public:

        // constructors
        PeriodicAdapterWithPQConfig()  { }
        PeriodicAdapterWithPQConfig(const PeriodicStructureAdapter& src)
            : PeriodicStructureAdapter(src)
        { }

        // overloads
        void customPQConfig(diffpy::srreal::PairQuantity* pq) const
        {
            std::map<std::string, double>::const_iterator it;
            for (it = cfg.begin(); it != cfg.end(); ++it)
            {
                pq->setDoubleAttr(it->first, it->second);
            }
        }

        // data
        std::map<std::string, double> cfg;

    private:

        // serialization
        friend class boost::serialization::access;
        template<class Archive>
            void serialize(Archive& ar, const unsigned int version)
        {
            using boost::serialization::base_object;
            ar & base_object<PeriodicStructureAdapter>(*this);
            ar & cfg;
        }

};

BOOST_CLASS_EXPORT(PeriodicAdapterWithPQConfig)

// Typedefs and Functions ----------------------------------------------------

typedef boost::shared_ptr<PeriodicAdapterWithPQConfig>
    PeriodicAdapterWithPQConfigPtr;

inline
PeriodicAdapterWithPQConfigPtr
addCustomPQConfig(StructureAdapterPtr stru)
{
    PeriodicStructureAdapterPtr pstru =
        boost::dynamic_pointer_cast<PeriodicStructureAdapter>(stru);
    return PeriodicAdapterWithPQConfigPtr(
            new PeriodicAdapterWithPQConfig(*pstru));
}

#endif  // TEST_CUSTOMPQCONFIG_HPP_INCLUDED
