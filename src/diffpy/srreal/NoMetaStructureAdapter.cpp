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
* class NoMetaStructureAdapter -- deprecated StructureAdapter proxy retained
*     for source and serialization compatibility.  Structures no longer
*     configure calculators through metadata.
*
* nometa -- compatibility function returning the source StructureAdapter
*
*****************************************************************************/

#include <diffpy/serialization.ipp>
#include <diffpy/srreal/StructureDifference.hpp>
#include <diffpy/srreal/NoMetaStructureAdapter.hpp>

namespace diffpy {
namespace srreal {

//////////////////////////////////////////////////////////////////////////////
// class NoMetaStructureAdapter
//////////////////////////////////////////////////////////////////////////////

// Constructor ---------------------------------------------------------------

NoMetaStructureAdapter::NoMetaStructureAdapter(
        StructureAdapterPtr srcstructure)
{
    std::shared_ptr<NoMetaStructureAdapter> nmptr =
        std::dynamic_pointer_cast<NoMetaStructureAdapter>(srcstructure);
    if (nmptr)  srcstructure = nmptr->getSourceStructure();
    msrcstructure = srcstructure.get() ?
        srcstructure : emptyStructureAdapter();
}

// Public Methods ------------------------------------------------------------

StructureAdapterPtr NoMetaStructureAdapter::clone() const
{
    std::shared_ptr<NoMetaStructureAdapter> rv(new NoMetaStructureAdapter);
    if (msrcstructure)  rv->msrcstructure = msrcstructure->clone();
    return rv;
}


BaseBondGeneratorPtr NoMetaStructureAdapter::createBondGenerator() const
{
    BaseBondGeneratorPtr bnds = msrcstructure->createBondGenerator();
    return bnds;
}


int NoMetaStructureAdapter::countSites() const
{
    return msrcstructure->countSites();
}


double NoMetaStructureAdapter::numberDensity() const
{
    return msrcstructure->numberDensity();
}


const std::string& NoMetaStructureAdapter::siteAtomType(int idx) const
{
    return msrcstructure->siteAtomType(idx);
}


const R3::Vector& NoMetaStructureAdapter::siteCartesianPosition(
        int idx) const
{
    return msrcstructure->siteCartesianPosition(idx);
}


int NoMetaStructureAdapter::siteMultiplicity(int idx) const
{
    return msrcstructure->siteMultiplicity(idx);
}


double NoMetaStructureAdapter::siteOccupancy(int idx) const
{
    return msrcstructure->siteOccupancy(idx);
}


bool NoMetaStructureAdapter::siteAnisotropy(int idx) const
{
    return msrcstructure->siteAnisotropy(idx);
}


const R3::Matrix& NoMetaStructureAdapter::siteCartesianUij(int idx) const
{
    return msrcstructure->siteCartesianUij(idx);
}


StructureDifference
NoMetaStructureAdapter::diff(StructureAdapterConstPtr other) const
{
    StructureDifference sd = this->StructureAdapter::diff(other);
    if (sd.stru0 == sd.stru1)  return sd;
    typedef std::shared_ptr<const class NoMetaStructureAdapter> PPtr;
    PPtr pother = std::dynamic_pointer_cast<PPtr::element_type>(other);
    if (!pother)  return sd;
    sd = this->getSourceStructure()->diff(pother->getSourceStructure());
    assert(sd.stru0 == this->getSourceStructure());
    assert(sd.stru1 == pother->getSourceStructure());
    sd.stru0 = this->shared_from_this();
    sd.stru1 = other;
    return sd;
}


StructureAdapterPtr
NoMetaStructureAdapter::getSourceStructure()
{
    return msrcstructure;
}


StructureAdapterConstPtr
NoMetaStructureAdapter::getSourceStructure() const
{
    return msrcstructure;
}

// Routines ------------------------------------------------------------------

StructureAdapterPtr nometa(StructureAdapterPtr stru)
{
    return stru ? stru : emptyStructureAdapter();
}

}   // namespace srreal
}   // namespace diffpy

// Serialization -------------------------------------------------------------

DIFFPY_INSTANTIATE_SERIALIZATION(diffpy::srreal::NoMetaStructureAdapter)
BOOST_CLASS_EXPORT_IMPLEMENT(diffpy::srreal::NoMetaStructureAdapter)

// End of file
