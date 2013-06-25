/*****************************************************************************
*
* diffpy.srreal     by DANSE Diffraction group
*                   Simon J. L. Billinge
*                   (c) 2012 Trustees of the Columbia University
*                   in the City of New York.  All rights reserved.
*
* File coded by:    Pavol Juhas
*
* See AUTHORS.txt for a list of people who contributed.
* See LICENSE.txt for license information.
*
******************************************************************************
*
* class AtomicStructureAdapter -- universal structure adapter for
*     a non-periodic set of atoms.
*
* class AtomicStructureBondGenerator -- bond generator
*
*****************************************************************************/

#ifndef ATOMICSTRUCTUREADAPTER_HPP_INCLUDED
#define ATOMICSTRUCTUREADAPTER_HPP_INCLUDED

#include <boost/serialization/vector.hpp>

#include <diffpy/srreal/StructureAdapter.hpp>

namespace diffpy {
namespace srreal {

class AtomAdapter
{
    public:
        
        // constructor
        AtomAdapter() :
            cartesianposition(0.0, 0.0, 0.0),
            occupancy(1.0),
            anisotropy(false),
            cartesianuij(R3::zeros())
        { };

        // data
        std::string atomtype;
        R3::Vector cartesianposition;
        double occupancy;
        bool anisotropy;
        R3::Matrix cartesianuij;

    private:

        // serialization
        friend class boost::serialization::access;
        template<class Archive>
            void serialize(Archive& ar, const unsigned int version)
        {
            ar & atomtype;
            ar & cartesianposition;
            ar & occupancy;
            ar & anisotropy;
            ar & cartesianuij;
        }

};


class AtomicStructureAdapter : public StructureAdapter
{
    public:

        // methods - overloaded
        virtual BaseBondGeneratorPtr createBondGenerator() const;
        virtual int countSites() const;
        virtual const R3::Vector& siteCartesianPosition(int idx) const;
        virtual double siteOccupancy(int idx) const;
        virtual bool siteAnisotropy(int idx) const;
        virtual const R3::Matrix& siteCartesianUij(int idx) const;
        virtual const std::string& siteAtomType(int idx) const;

        // methods - own
        void insert(int, const AtomAdapter&);
        void append(const AtomAdapter&);
        void remove(int);

    private:

        // data
        std::vector<AtomAdapter> matoms;

        // serialization
        friend class boost::serialization::access;
        template<class Archive>
            void serialize(Archive& ar, const unsigned int version)
        {
            ar & boost::serialization::base_object<StructureAdapter>(*this);
            ar & matoms;
        }

};


class AtomicStructureBondGenerator : public BaseBondGenerator
{
    public:

        // constructors
        AtomicStructureBondGenerator(StructureAdapterConstPtr);

    protected:

        // data
        const AtomicStructureAdapter* mastructure;
};


}   // namespace srreal
}   // namespace diffpy

#endif  // ATOMICSTRUCTUREADAPTER_HPP_INCLUDED
