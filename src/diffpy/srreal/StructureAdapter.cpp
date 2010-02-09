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
* class StructureAdapter -- abstract base class for interfacing general
*     structure objects with srreal classes such as PairQuantity
*
* $Id$
*
*****************************************************************************/

#include <cstdlib>
#include <cctype>

#include <diffpy/srreal/StructureAdapter.hpp>
#include <diffpy/srreal/BaseBondGenerator.hpp>

using namespace std;

namespace diffpy {
namespace srreal {

// Public Methods ------------------------------------------------------------

double StructureAdapter::totalOccupancy() const
{
    double total_occupancy = 0.0;
    int cntsites = this->countSites();
    for (int i = 0; i < cntsites; ++i)
    {
        total_occupancy += this->siteOccupancy(i) * this->siteMultiplicity(i);
    }
    return total_occupancy;
}


double StructureAdapter::numberDensity() const
{
    return 0.0;
}


const string& StructureAdapter::siteAtomType(int idx) const
{
    static string rv = "";
    return rv;
}


double StructureAdapter::siteMultiplicity(int idx) const
{
    return 1.0;
}


double StructureAdapter::siteOccupancy(int idx) const
{
    return 1.0;
}

// Routines ------------------------------------------------------------------

std::string atomBareSymbol(const std::string& atomtype)
{
    string::size_type pb, pe;
    pb = atomtype.find_first_not_of("0123456789- \t");
    pe = atomtype.find_last_not_of("+-12345678 \t");
    string rv = atomtype.substr(pb, pe - pb + 1);
    return rv;
}


/// Return valence of possibly ionic symbol such as "S2-" or "Cl-".
int atomValence(const std::string& atomtype)
{
    string::const_reverse_iterator ci;
    int rv = 0;
    for (ci = atomtype.rbegin(); ci != atomtype.rend(); ++ci)
    {
        if (isblank(*ci))   continue;
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
        if (rv && '1' <= *ci && *ci <= '8')
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
