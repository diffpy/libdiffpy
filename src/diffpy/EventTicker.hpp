/*****************************************************************************
*
* diffpy.srreal     Complex Modeling Initiative
*                   Pavol Juhas
*                   (c) 2013 Brookhaven National Laboratory,
*                   Upton, New York.  All rights reserved.
*
* File coded by:    Pavol Juhas
*
* See AUTHORS.txt for a list of people who contributed.
* See LICENSE.txt for license information.
*
******************************************************************************
*
* class EventTicker for managing modification times of dependent objects.
*
* This is used to determine if cached values need update, due to a change
* of their dependencies.  The EventTicker class is inspired by the
* ObjCryst::RefinableObjClock class.
*
*****************************************************************************/

#ifndef EVENTTICKER_HPP_INCLUDED
#define EVENTTICKER_HPP_INCLUDED

#include <boost/serialization/utility.hpp>
#include <boost/serialization/split_member.hpp>

namespace diffpy {
namespace eventticker {

class EventTicker
{
    public:

        // constructor
        EventTicker();

        // methods
        /// record an event and advance the global ticker.
        void click();
        /// advance internal ticker if the other ticker is newer.
        void updateFrom(const EventTicker& other);

        // comparison
        bool operator<(const EventTicker&) const;
        bool operator<=(const EventTicker&) const;
        bool operator>(const EventTicker&) const;
        bool operator>=(const EventTicker&) const;
        bool operator==(const EventTicker&) const;

    private:

        // global counter
        static std::pair<long, long> gtick;

        // data
        std::pair<long, long> mtick;

        // serialization
        friend class boost::serialization::access;
        template<class Archive>
            void save(Archive& ar, const unsigned int version) const
        {
            ar << mtick << gtick;
        }

        template<class Archive>
            void load(Archive& ar, const unsigned int version)
        {
            std::pair<long, long> ga;
            ar >> mtick >> ga;
            if (ga > gtick)
            {
                if (ga.first != gtick.first)  mtick.first = mtick.second = 0;
                else  mtick.second += gtick.second - ga.second;
            }
        }

        BOOST_SERIALIZATION_SPLIT_MEMBER()

};

}   // namespace eventticker
}   // namespace diffpy

#endif  // EVENTTICKER_HPP_INCLUDED
