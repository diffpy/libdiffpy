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
*****************************************************************************/

#include <cassert>
#include <cctype>

#include <diffpy/mathutils.hpp>
#include <diffpy/srreal/StructureAdapter.hpp>
#include <diffpy/srreal/VR3Structure.hpp>

using namespace std;
using diffpy::mathutils::eps_eq;

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


int StructureAdapter::siteMultiplicity(int idx) const
{
    return 1;
}


double StructureAdapter::siteOccupancy(int idx) const
{
    return 1.0;
}

// Routines ------------------------------------------------------------------

StructureAdapterPtr emptyStructureAdapter()
{
    static StructureAdapterPtr stru(new VR3Adapter);
    assert(stru.get() && stru->countSites() == 0);
    return stru;
}


double meanSquareDisplacement(const R3::Matrix& Uijcartn,
        const R3::Vector& s, bool anisotropy)
{
    double rv;
    if (anisotropy)
    {
        assert(R3::norm(s) > 0);
        assert(eps_eq(Uijcartn(0,1), Uijcartn(1,0)));
        assert(eps_eq(Uijcartn(0,2), Uijcartn(2,0)));
        assert(eps_eq(Uijcartn(1,2), Uijcartn(2,1)));
        static R3::Vector sn;
        sn = s / R3::norm(s);
        rv = Uijcartn(0,0) * sn(0) * sn(0) +
             Uijcartn(1,1) * sn(1) * sn(1) +
             Uijcartn(2,2) * sn(2) * sn(2) +
             2 * Uijcartn(0,1) * sn(0) * sn(1) +
             2 * Uijcartn(0,2) * sn(0) * sn(2) +
             2 * Uijcartn(1,2) * sn(1) * sn(2);
    }
    else
    {
        assert(eps_eq(Uijcartn(0,0), Uijcartn(1,1)));
        assert(eps_eq(Uijcartn(0,0), Uijcartn(2,2)));
        rv = Uijcartn(0,0);
    }
    return rv;
}


double maxUii(StructureAdapterPtr stru)
{
    if (!stru)  return 0.0;
    double rv = 0.0;
    for (int i = 0; i < stru->countSites(); ++i)
    {
        const R3::Matrix& U = stru->siteCartesianUij(i);
        for (int k = 0; k < R3::Ndim; k++)
        {
            if (U(k,k) > rv)   rv = U(k,k);
        }
    }
    return rv;
}


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
