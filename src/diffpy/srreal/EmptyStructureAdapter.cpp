/*****************************************************************************
*
* libdiffpy         Complex Modeling Initiative
*                   (c) 2016 Brookhaven Science Associates,
*                   Brookhaven National Laboratory.
*                   All rights reserved.
*
* File coded by:    Pavol Juhas
*
* See AUTHORS.txt for a list of people who contributed.
* See LICENSE.txt for license information.
*
******************************************************************************
*
* class EmptyStructureAdapter -- immutable class that is always empty.
*
*****************************************************************************/

#include <diffpy/serialization.ipp>
#include <diffpy/srreal/StructureAdapter.hpp>

namespace diffpy {
namespace srreal {

// Local constants -----------------------------------------------------------

namespace {

std::out_of_range outofrange("Index out of range.");

}   // namespace

//////////////////////////////////////////////////////////////////////////////
// class EmptyStructureAdapter
//////////////////////////////////////////////////////////////////////////////

class EmptyStructureAdapter : public StructureAdapter
{
    public:

        // virtual methods from the base class

        virtual StructureAdapterPtr clone() const
        {
            StructureAdapterPtr rv =
                boost::const_pointer_cast<StructureAdapter>(shared_from_this());
            return rv;
        }


        virtual BaseBondGeneratorPtr createBondGenerator() const
        {
            return boost::make_shared<BaseBondGenerator>(shared_from_this());
        }


        virtual int countSites() const
        {
            return 0;
        }

        // reusing StructureAdapter::totalOccupancy()
        // reusing StructureAdapter::numberDensity()

        virtual const std::string& siteAtomType(int idx) const
        {
            throw outofrange;
        }


        virtual const R3::Vector& siteCartesianPosition(int idx) const
        {
            throw outofrange;
        }


        virtual int siteMultiplicity(int idx) const
        {
            throw outofrange;
        }


        virtual double siteOccupancy(int idx) const
        {
            throw outofrange;
        }


        virtual bool siteAnisotropy(int idx) const
        {
            throw outofrange;
        }


        virtual const R3::Matrix& siteCartesianUij(int idx) const
        {
            throw outofrange;
        }

        // reusing StructureAdapter::customPQConfig()
        // reusing StructureAdapter::diff()

    private:

        // serialization
        friend class boost::serialization::access;
        template<class Archive>
            void serialize(Archive& ar, const unsigned int version)
        {
            ar & boost::serialization::base_object<StructureAdapter>(*this);
        }

};

// Routines ------------------------------------------------------------------

StructureAdapterPtr emptyStructureAdapter()
{
    static StructureAdapterPtr stru(new EmptyStructureAdapter);
    assert(stru && stru->countSites() == 0);
    return stru;
}

}   // namespace srreal
}   // namespace diffpy

BOOST_CLASS_EXPORT(diffpy::srreal::EmptyStructureAdapter)
DIFFPY_INSTANTIATE_PTR_SERIALIZATION(diffpy::srreal::EmptyStructureAdapter)

// End of file
