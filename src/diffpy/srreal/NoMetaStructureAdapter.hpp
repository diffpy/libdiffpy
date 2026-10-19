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

#ifndef NOMETASTRUCTUREADAPTER_HPP_INCLUDED
#define NOMETASTRUCTUREADAPTER_HPP_INCLUDED

#include <diffpy/srreal/StructureAdapter.hpp>

namespace diffpy {
namespace srreal {

class NoMetaStructureAdapter : public StructureAdapter
{
    public:

        // constructors
        NoMetaStructureAdapter()  { }
        NoMetaStructureAdapter(StructureAdapterPtr);

        // methods - overloaded
        virtual StructureAdapterPtr clone() const;
        virtual BaseBondGeneratorPtr createBondGenerator() const;
        virtual int countSites() const;
        virtual double numberDensity() const;
        virtual const std::string& siteAtomType(int idx) const;
        virtual const R3::Vector& siteCartesianPosition(int idx) const;
        virtual int siteMultiplicity(int idx) const;
        virtual double siteOccupancy(int idx) const;
        virtual bool siteAnisotropy(int idx) const;
        virtual const R3::Matrix& siteCartesianUij(int idx) const;
        virtual StructureDifference diff(StructureAdapterConstPtr) const;

        // methods - own
        StructureAdapterPtr getSourceStructure();
        StructureAdapterConstPtr getSourceStructure() const;

    private:

        // data
        StructureAdapterPtr msrcstructure;

        // serialization
        friend class boost::serialization::access;
        template<class Archive>
            void serialize(Archive& ar, const unsigned int version)
        {
            ar & boost::serialization::base_object<StructureAdapter>(*this);
            ar & msrcstructure;
        }

};

// Routines ------------------------------------------------------------------

/// Return the input adapter unchanged, or the empty adapter for a null input.
StructureAdapterPtr nometa(StructureAdapterPtr stru);


/// Convert an adaptable structure without creating a proxy.
template <class T>
StructureAdapterPtr nometa(const T& stru)
{
    StructureAdapterPtr bstru = convertToStructureAdapter(stru);
    return nometa(bstru);
}

}   // namespace srreal
}   // namespace diffpy

// Serialization -------------------------------------------------------------

BOOST_CLASS_EXPORT_KEY(diffpy::srreal::NoMetaStructureAdapter)

#endif  // NOMETASTRUCTUREADAPTER_HPP_INCLUDED
