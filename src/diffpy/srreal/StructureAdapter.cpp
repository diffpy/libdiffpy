/*****************************************************************************
*
* diffpy.srreal     by DANSE Diffraction group
*                   Simon J. L. Billinge
*                   (c) 2009 The Trustees of Columbia University
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

#include <diffpy/serialization.ipp>
#include <diffpy/mathutils.hpp>
#include <diffpy/srreal/StructureAdapter.hpp>
#include <diffpy/srreal/AtomicStructureAdapter.hpp>
#include <diffpy/srreal/StructureDifference.hpp>

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


StructureDifference
StructureAdapter::diff(StructureAdapterConstPtr other) const
{
    StructureDifference sd;
    return sd;
}

// Routines ------------------------------------------------------------------

StructureAdapterPtr emptyStructureAdapter()
{
    static StructureAdapterPtr stru(new AtomicStructureAdapter);
    assert(stru && stru->countSites() == 0);
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

}   // namespace srreal
}   // namespace diffpy

DIFFPY_INSTANTIATE_PTR_SERIALIZATION(diffpy::srreal::StructureAdapter)

// End of file
