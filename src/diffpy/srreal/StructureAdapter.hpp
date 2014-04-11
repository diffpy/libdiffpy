/*****************************************************************************
*
* libdiffpy         by DANSE Diffraction group
*                   Simon J. L. Billinge
*                   (c) 2009 The Trustees of Columbia University
*                   in the City of New York.  All rights reserved.
*
* File coded by:    Pavol Juhas
*
* See AUTHORS.txt for a list of people who contributed.
* See LICENSE_DANSE.txt for license information.
*
******************************************************************************
*
* class StructureAdapter -- abstract base class for interfacing general
*     structure objects with srreal classes such as PairQuantity
*
*****************************************************************************/

#ifndef STRUCTUREADAPTER_HPP_INCLUDED
#define STRUCTUREADAPTER_HPP_INCLUDED

#include <vector>
#include <boost/shared_ptr.hpp>
#include <boost/enable_shared_from_this.hpp>
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/assume_abstract.hpp>
#include <boost/serialization/shared_ptr.hpp>
#include <boost/serialization/export.hpp>

#include <diffpy/srreal/R3linalg.hpp>
#include <diffpy/srreal/BaseBondGenerator.hpp>

namespace diffpy {
namespace srreal {

class PairQuantity;
class StructureDifference;


/// @class StructureAdapter
/// @brief abstract adaptor to structure data needed by
/// PairQuantity calculator

class StructureAdapter :
    public boost::enable_shared_from_this<StructureAdapter>
{
    public:

        virtual ~StructureAdapter()  { }

        // methods
        virtual StructureAdapterPtr clone() const = 0;

        /// factory for creating compatible BondGenerator instance.
        virtual BaseBondGeneratorPtr createBondGenerator() const = 0;

        /// number of independent sites in the structure, before
        /// any symmetry expansion.
        virtual int countSites() const = 0;

        /// total number of atoms in the structure unit accounting
        /// for possibly fractional occupancies.
        double totalOccupancy() const;

        /// number density in the structure model or 0 when not defined.
        virtual double numberDensity() const;

        /// symbol for element or ion at the independent site @param idx
        virtual const std::string& siteAtomType(int idx) const;

        /// Cartesian coordinates of the independent site @param idx
        virtual const R3::Vector& siteCartesianPosition(int idx) const = 0;

        /// multiplicity of the independent site @param idx in the structure
        virtual int siteMultiplicity(int idx) const;

        /// site occupancy at the independent site @param idx
        virtual double siteOccupancy(int idx) const;

        /// flag for anisotropic atom displacements at independent site
        /// @param idx
        virtual bool siteAnisotropy(int idx) const = 0;

        /// tensor of atom displacement parameters converted in
        /// Cartesian coordinate system at independent site @param idx
        virtual const R3::Matrix& siteCartesianUij(int idx) const = 0;

        /// this method allows custom special configuration for a concrete
        /// pair of StructureAdapter and PairQuantity objects.
        virtual void customPQConfig(PairQuantity* pq) const  { }

        /// Return difference from the other StructureAdapter
        virtual StructureDifference diff(StructureAdapterConstPtr) const;

    private:

        // serialization
        friend class boost::serialization::access;
        template<class Archive>
            void serialize(Archive& ar, const unsigned int version)  { }

};

// Routines ------------------------------------------------------------------

/// Return a singleton instance of an empty StructureAdapter
StructureAdapterPtr emptyStructureAdapter();

/// Calculate MSD along specified direction in Cartesian space.
double meanSquareDisplacement(const R3::Matrix& Uijcartn, const R3::Vector& s,
        bool anisotropy=true);

/// Maximum diagonal Uii element from all atoms in the structure.
double maxUii(StructureAdapterPtr stru);

/// Translate an index container to a vector of string symbols
template <class T>
std::vector<std::string>
siteIndicesToTypes(const StructureAdapterPtr& stru, const T& indices)
{
    std::vector<std::string> rv;
    rv.reserve(stru->countSites());
    typename T::const_iterator i;
    for (i = indices.begin(); i != indices.end(); ++i)
    {
        const std::string& smbl = stru->siteAtomType(*i);
        rv.push_back(smbl);
    }
    return rv;
}

// Conversion helpers --------------------------------------------------------

/// createStructureAdapter a last-resort conversion of the argument type
/// to StructureAdapterPtr.  Classes convertible to StructureAdapter should
/// define a corresponding createStructureAdapter function for seamless
/// integration with PairQuantity calculators.

/// This converter provides by-value support for StructureAdapter instances.

inline
StructureAdapterPtr createStructureAdapter(const StructureAdapter& stru)
{
    return stru.clone();
}

/// convertToStructureAdapter is an interface function that should be used
/// by the clients that need to convert to StructureAdapterPtr.

template <class T>
StructureAdapterPtr convertToStructureAdapter(const T& stru)
{
    StructureAdapterPtr rv = createStructureAdapter(stru);
    return rv;
}


template <class T>
StructureAdapterPtr convertToStructureAdapter(const boost::shared_ptr<T>& stru)
{
    StructureAdapterPtr rv =
        boost::dynamic_pointer_cast<StructureAdapterPtr::element_type>(stru);
    assert(rv);
    return rv;
}


inline
StructureAdapterPtr convertToStructureAdapter(StructureAdapterConstPtr stru)
{
    StructureAdapterPtr rv =
        boost::const_pointer_cast<StructureAdapterPtr::element_type>(stru);
    return rv;
}

}   // namespace srreal
}   // namespace diffpy

// Serialization -------------------------------------------------------------

BOOST_SERIALIZATION_ASSUME_ABSTRACT(diffpy::srreal::StructureAdapter)

#endif  // STRUCTUREADAPTER_HPP_INCLUDED
