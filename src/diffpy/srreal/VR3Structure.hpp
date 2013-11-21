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
* class VR3Structure -- trivial structure representation.
* class VR3Adapter -- concrete StructureAdapter for VR3Structure
* class VR3BondGenerator -- concrete BaseBondGenerator for VR3Structure
*
*****************************************************************************/

#ifndef VR3STRUCTURE_HPP_INCLUDED
#define VR3STRUCTURE_HPP_INCLUDED

#include <boost/serialization/vector.hpp>

#include <diffpy/srreal/R3linalg.hpp>
#include <diffpy/srreal/StructureAdapter.hpp>

namespace diffpy {
namespace srreal {

typedef std::vector<R3::Vector> VR3Structure;


class VR3Adapter : public StructureAdapter
{
    friend class VR3BondGenerator;

    public:

        // constructors
        VR3Adapter()  { }
        VR3Adapter(const VR3Structure& vr3s);

        // methods
        virtual int countSites() const;
        virtual const R3::Vector& siteCartesianPosition(int idx) const;
        virtual bool siteAnisotropy(int idx) const;
        virtual const R3::Matrix& siteCartesianUij(int idx) const;

        virtual BaseBondGeneratorPtr createBondGenerator() const;

    private:

        // data
        VR3Structure mvr3structure;

        // serialization
        friend class boost::serialization::access;
        template<class Archive>
            void serialize(Archive& ar, const unsigned int version)
        {
            ar & boost::serialization::base_object<StructureAdapter>(*this);
            ar & mvr3structure;
        }

};


class VR3BondGenerator : public BaseBondGenerator
{
    public:

        // constructors
        VR3BondGenerator(StructureAdapterConstPtr);

    private:
        // data
        const VR3Structure& mvr3structure;
};

/// Factory for constructing StructureAdapter for VR3Structure
inline
StructureAdapterPtr createStructureAdapter(const VR3Structure& vr3stru)
{
    StructureAdapterPtr adapter(new VR3Adapter(vr3stru));
    return adapter;
}

}   // namespace srreal
}   // namespace diffpy

// Serialization -------------------------------------------------------------

BOOST_CLASS_EXPORT_KEY(diffpy::srreal::VR3Adapter)

#endif  // VR3STRUCTURE_HPP_INCLUDED
