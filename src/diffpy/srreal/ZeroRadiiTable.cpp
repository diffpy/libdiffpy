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
* class ZeroRadiiTable -- this AtomRadiiTable returns zero for all elements
*
* $Id$
*
*****************************************************************************/

#include <algorithm>
#include <cctype>
#include <sstream>
#include <stdexcept>
#include <vector>
#include <boost/smart_ptr/make_shared.hpp>

#include <diffpy/srreal/ZeroRadiiTable.hpp>

namespace diffpy {
namespace srreal {

using namespace std;

// Public Methods ------------------------------------------------------------

// HasClassRegistry methods

AtomRadiiTablePtr ZeroRadiiTable::create() const
{
    return boost::make_shared<ZeroRadiiTable>();
}


AtomRadiiTablePtr ZeroRadiiTable::clone() const
{
    return boost::make_shared<ZeroRadiiTable>(*this);
}


const string& ZeroRadiiTable::type() const
{
    static string rv = "zeroradii";
    return rv;
}

// own methods

double ZeroRadiiTable::tableLookup(const string& smbl) const
{
    return 0.0;
}

// Registration --------------------------------------------------------------

bool reg_ZeroRadiiTable = ZeroRadiiTable().registerThisType();

}   // srreal
}   // diffpy

// End of file
