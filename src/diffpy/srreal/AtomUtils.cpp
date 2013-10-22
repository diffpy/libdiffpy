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
* Small routines related to atom properties.
*
*****************************************************************************/

#include <diffpy/srreal/AtomUtils.hpp>

namespace diffpy {
namespace srreal {

std::string atomBareSymbol(const std::string& atomtype)
{
    std::string::size_type pb, pe;
    pb = atomtype.find_first_not_of("0123456789- \t");
    pe = atomtype.find_last_not_of("+-012345678 \t");
    std::string rv = atomtype.substr(pb, pe - pb + 1);
    return rv;
}


/// Return valence of possibly ionic symbol such as "S2-" or "Cl-".
int atomValence(const std::string& atomtype)
{
    std::string::const_reverse_iterator ci;
    int rv = 0;
    for (ci = atomtype.rbegin(); ci != atomtype.rend(); ++ci)
    {
        if (isspace(*ci))   continue;
        // read trailing sign
        if (rv == 0)
        {
            if (*ci == '+')
            {
                rv = 1;
                continue;
            }
            if (*ci == '-')
            {
                rv = -1;
                continue;
            }
            break;
        }
        // allow one [1-8] digit before the +/- sign
        if (rv && '0' <= *ci && *ci <= '8')
        {
            rv *= (*ci - '0');
        }
        break;
    }
    return rv;
}

}   // namespace srreal
}   // namespace diffpy

// End of file
