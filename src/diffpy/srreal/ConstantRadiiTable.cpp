/*****************************************************************************
*
* diffpy.srreal     by DANSE Diffraction group
*                   Simon J. L. Billinge
*                   (c) 2011 Trustees of the Columbia University
*                   in the City of New York.  All rights reserved.
*
* File coded by:    Pavol Juhas
*
* See AUTHORS.txt for a list of people who contributed.
* See LICENSE.txt for license information.
*
******************************************************************************
*
* class ConstantRadiiTable -- this AtomRadiiTable returns same value
*   for all elements unless redefined by the setCustom method.
*
*****************************************************************************/

#include <algorithm>
#include <cctype>
#include <sstream>
#include <stdexcept>
#include <vector>
#include <boost/smart_ptr/make_shared.hpp>
#include <boost/serialization/export.hpp>

#include <diffpy/srreal/ConstantRadiiTable.hpp>
#include <diffpy/serialization.hpp>

namespace diffpy {
namespace srreal {

using namespace std;

// Constructors --------------------------------------------------------------

ConstantRadiiTable::ConstantRadiiTable() : mdefaultradius(0.0)
{ }

// Public Methods ------------------------------------------------------------

// HasClassRegistry methods

AtomRadiiTablePtr ConstantRadiiTable::create() const
{
    return boost::make_shared<ConstantRadiiTable>();
}


AtomRadiiTablePtr ConstantRadiiTable::clone() const
{
    return boost::make_shared<ConstantRadiiTable>(*this);
}


const string& ConstantRadiiTable::type() const
{
    static string rv = "constant";
    return rv;
}

// own methods

double ConstantRadiiTable::standardLookup(const string& smbl) const
{
    return mdefaultradius;
}


void ConstantRadiiTable::setDefault(double value)
{
    mdefaultradius = value;
}


double ConstantRadiiTable::getDefault() const
{
    return mdefaultradius;
}

// Registration --------------------------------------------------------------

bool reg_ConstantRadiiTable = ConstantRadiiTable().registerThisType();

}   // srreal
}   // diffpy

// Serialization -------------------------------------------------------------

DIFFPY_INSTANTIATE_SERIALIZE(diffpy::srreal::ConstantRadiiTable)
BOOST_CLASS_EXPORT(diffpy::srreal::ConstantRadiiTable)

// End of file
