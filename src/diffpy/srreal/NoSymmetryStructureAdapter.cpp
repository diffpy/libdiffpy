/*****************************************************************************
*
* libdiffpy         by DANSE Diffraction group
*                   Simon J. L. Billinge
*                   (c) 2010 The Trustees of Columbia University
*                   in the City of New York.  All rights reserved.
*
* File coded by:    Pavol Juhas
*
* See AUTHORS.txt for a list of people who contributed.
* See LICENSE_DANSE.txt for license information.
*
******************************************************************************
*
* class NoSymmetryStructureAdapter -- StructureAdapter class that removes
*     any symmetry expansions (rotations or periodic translations) from
*     another StructureAdapter instance.  This can be used to use only
*     the asymmetric unit from any adapter to crystal structure.
*
* class NoSymmetryBondGenerator -- bond generator
*
* nosymmetry -- factory function that creates a NoSymmetryStructureAdapter
*     instance inside StructureAdapterPtr
*
*****************************************************************************/

#include <diffpy/serialization.ipp>
#include <diffpy/srreal/StructureDifference.hpp>
#include <diffpy/srreal/NoSymmetryStructureAdapter.hpp>

namespace diffpy {
namespace srreal {

//////////////////////////////////////////////////////////////////////////////
// class NoSymmetryStructureAdapter
//////////////////////////////////////////////////////////////////////////////

// Constructor ---------------------------------------------------------------

NoSymmetryStructureAdapter::NoSymmetryStructureAdapter(
        StructureAdapterPtr srcstructure)
{
    boost::shared_ptr<NoSymmetryStructureAdapter> nmptr =
        boost::dynamic_pointer_cast<NoSymmetryStructureAdapter>(srcstructure);
    if (nmptr)  srcstructure = nmptr->getSourceStructure();
    msrcstructure = srcstructure.get() ?
        srcstructure : emptyStructureAdapter();
}

// Public Methods ------------------------------------------------------------

StructureAdapterPtr NoSymmetryStructureAdapter::clone() const
{
    boost::shared_ptr<NoSymmetryStructureAdapter>
        rv(new NoSymmetryStructureAdapter);
    if (msrcstructure)  rv->msrcstructure = msrcstructure->clone();
    return rv;
}


BaseBondGeneratorPtr NoSymmetryStructureAdapter::createBondGenerator() const
{
    BaseBondGeneratorPtr bnds(new BaseBondGenerator(shared_from_this()));
    return bnds;
}


int NoSymmetryStructureAdapter::countSites() const
{
    return msrcstructure->countSites();
}


double NoSymmetryStructureAdapter::numberDensity() const
{
    return 0.0;
}


const std::string& NoSymmetryStructureAdapter::siteAtomType(int idx) const
{
    return msrcstructure->siteAtomType(idx);
}


const R3::Vector& NoSymmetryStructureAdapter::siteCartesianPosition(
        int idx) const
{
    return msrcstructure->siteCartesianPosition(idx);
}


double NoSymmetryStructureAdapter::siteOccupancy(int idx) const
{
    return msrcstructure->siteOccupancy(idx);
}


bool NoSymmetryStructureAdapter::siteAnisotropy(int idx) const
{
    return msrcstructure->siteAnisotropy(idx);
}


const R3::Matrix& NoSymmetryStructureAdapter::siteCartesianUij(int idx) const
{
    return msrcstructure->siteCartesianUij(idx);
}


void NoSymmetryStructureAdapter::customPQConfig(PairQuantity* pq) const
{
    msrcstructure->customPQConfig(pq);
}


StructureDifference
NoSymmetryStructureAdapter::diff(StructureAdapterConstPtr other) const
{
    StructureDifference alldiffer;
    typedef boost::shared_ptr<const class NoSymmetryStructureAdapter> PPtr;
    PPtr pother = boost::dynamic_pointer_cast<PPtr::element_type>(other);
    if (!pother)  return alldiffer;
    StructureDifference sd =
        this->getSourceStructure()->diff(pother->getSourceStructure());
    assert(sd.stru0 == this->getSourceStructure());
    assert(sd.stru1 == pother->getSourceStructure());
    sd.stru0 = this->shared_from_this();
    sd.stru1 = other;
    return sd;
}


StructureAdapterPtr
NoSymmetryStructureAdapter::getSourceStructure()
{
    return msrcstructure;
}


StructureAdapterConstPtr
NoSymmetryStructureAdapter::getSourceStructure() const
{
    return msrcstructure;
}

// Routines ------------------------------------------------------------------

StructureAdapterPtr nosymmetry(StructureAdapterPtr stru)
{
    StructureAdapterPtr rv =
        boost::dynamic_pointer_cast<NoSymmetryStructureAdapter>(stru) ? stru :
        StructureAdapterPtr(new NoSymmetryStructureAdapter(stru));
    return rv;
}

}   // namespace srreal
}   // namespace diffpy

// Serialization -------------------------------------------------------------

DIFFPY_INSTANTIATE_SERIALIZATION(diffpy::srreal::NoSymmetryStructureAdapter)
BOOST_CLASS_EXPORT_IMPLEMENT(diffpy::srreal::NoSymmetryStructureAdapter)

// End of file
