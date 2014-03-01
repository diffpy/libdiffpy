/*****************************************************************************
*
* libdiffpy         Complex Modeling Initiative
*                   (c) 2013 Brookhaven Science Associates,
*                   Brookhaven National Laboratory.
*                   All rights reserved.
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

#include <diffpy/EventTicker.hpp>
#include <diffpy/serialization.ipp>

namespace diffpy {
namespace eventticker {

//////////////////////////////////////////////////////////////////////////////
// class EventTicker
//////////////////////////////////////////////////////////////////////////////

// Constructor ---------------------------------------------------------------

EventTicker::EventTicker() : mtick(0, 0)
{ }

// Public Methods ------------------------------------------------------------

void EventTicker::click()
{
    ++gtick.second;
    if (0 >= gtick.second)
    {
        ++gtick.first;
        gtick.second = 0;
    }
    mtick = gtick;
}


void EventTicker::updateFrom(const EventTicker& other)
{
    if (other > *this)  mtick = other.mtick;
}


EventTicker::value_type EventTicker::value() const
{
    return mtick;
}


bool EventTicker::operator<(const EventTicker& other) const
{
    return mtick < other.mtick;
}


bool EventTicker::operator<=(const EventTicker& other) const
{
    return mtick <= other.mtick;
}


bool EventTicker::operator>(const EventTicker& other) const
{
    return mtick > other.mtick;
}


bool EventTicker::operator>=(const EventTicker& other) const
{
    return mtick >= other.mtick;
}


bool EventTicker::operator==(const EventTicker& other) const
{
    return mtick == other.mtick;
}


bool EventTicker::operator!=(const EventTicker& other) const
{
    return mtick != other.mtick;
}


// Static Global Data --------------------------------------------------------

EventTicker::value_type EventTicker::gtick(0, 0);

}   // namespace eventticker
}   // namespace diffpy

// Serialization -------------------------------------------------------------

DIFFPY_INSTANTIATE_SERIALIZATION(diffpy::eventticker::EventTicker)

// End of file
