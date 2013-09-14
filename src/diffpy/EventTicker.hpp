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

#include <utility>

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

    private:

        // global counter
        static std::pair<unsigned long,unsigned long> gtick;

        // data
        std::pair<unsigned long,unsigned long> mtick;
};

}   // namespace eventticker
}   // namespace diffpy

#endif  // EVENTTICKER_HPP_INCLUDED
