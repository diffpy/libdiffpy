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

#ifndef STRUCTUREADAPTER_HPP_INCLUDED
#define STRUCTUREADAPTER_HPP_INCLUDED

#include <vector>
#include <boost/shared_ptr.hpp>
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/assume_abstract.hpp>
#include <boost/serialization/shared_ptr.hpp>
#include <boost/serialization/export.hpp>

#include <diffpy/srreal/R3linalg.hpp>
#include <diffpy/srreal/BaseBondGenerator.hpp>

namespace diffpy {
namespace srreal {

class PairQuantity;

/// shared pointer to StructureAdapter

typedef boost::shared_ptr<class StructureAdapter> StructureAdapterPtr;

/// @class StructureAdapter
/// @brief abstract adaptor to structure data needed by
/// PairQuantity calculator

class StructureAdapter
{
    public:

        virtual ~StructureAdapter()  { }

        // methods

        /// factory for creating compatible BondGenerator instance.
        virtual BaseBondGenerator* createBondGenerator() const = 0;

        /// number of independent sites in the structure, before
        /// any symmetry expansion.
        virtual int countSites() const = 0;

        /// total number of atoms in the structure unit accounting
        /// for possibly fractional occupancies.
        virtual double totalOccupancy() const;

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

    private:

        // serialization
        friend class boost::serialization::access;
        template<class Archive>
            void serialize(Archive& ar, const unsigned int version)  { }

};

// Routines ------------------------------------------------------------------

/// Calculate MSD along specified direction in Cartesian space.
double meanSquareDisplacement(const R3::Matrix& Uijcartn, const R3::Vector& s,
        bool anisotropy=true);

/// Maximum diagonal Uii element from all atoms in the structure.
double maxUii(StructureAdapterPtr stru);

/// Return bare element symbol from atom symbol that may be isotope or ion,
/// for example "Ca2+" or "12-C".
std::string atomBareSymbol(const std::string& atomtype);

/// Return valence of possibly ionic symbol such as "S2-" or "Cl-".
int atomValence(const std::string& atomtype);

/// Translate an index container to a vector of string symbols
template <class T>
std::vector<std::string>
siteIndicesToTypes(const StructureAdapterPtr stru, const T& indices)
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

}   // namespace srreal
}   // namespace diffpy

// Serialization -------------------------------------------------------------

BOOST_SERIALIZATION_ASSUME_ABSTRACT(diffpy::srreal::StructureAdapter)

#endif  // STRUCTUREADAPTER_HPP_INCLUDED
