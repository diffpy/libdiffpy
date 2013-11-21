/*****************************************************************************
*
* diffpy.srreal     by DANSE Diffraction group
*                   Simon J. L. Billinge
*                   (c) 2011 Trustees of the Columbia University
*                   in the City of New York.  All rights reserved.
*
* File coded by:    Pavol Juhas
*
* See AUTHORS.txt for a list of people who contributed.
* See LICENSE.txt for license information.
*
******************************************************************************
*
* class AtomRadiiTable -- storage of empirical atomic radii
*   This is an abstract base class that provides a createByType factory
*   function for creating a registered concrete instances.
*   The derived classes has to overload the standardLookup method and also
*   create, clone and type methods from the HasClassRegistry base.
*
*****************************************************************************/

#ifndef ATOMRADIITABLE_HPP_INCLUDED
#define ATOMRADIITABLE_HPP_INCLUDED

#include <string>
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/shared_ptr.hpp>
#include <boost/serialization/assume_abstract.hpp>
#include <boost/serialization/export.hpp>

#include <diffpy/boostextensions/serialize_unordered_map.hpp>
#include <diffpy/HasClassRegistry.hpp>

namespace diffpy {
namespace srreal {

class AtomRadiiTable :
    public diffpy::HasClassRegistry<AtomRadiiTable>
{
    public:

        // types
        typedef boost::unordered_map<std::string,double> CustomRadiiStorage;

        // methods
        /// fast value lookup, which does not change the table.
        double lookup(const std::string& smbl) const;
        /// overloadable lookup function that retrieved standard values
        virtual double standardLookup(const std::string& smbl) const = 0;
        /// set custom radius for a specified atom symbol
        void setCustom(const std::string& smbl, double radius);
        /// set custom radii from a string in (A1:r1, A2:r2, ...) format
        void fromString(const std::string& s);
        /// reset custom value for the specified atom type
        void resetCustom(const std::string& smbl);
        /// reset all custom values
        void resetAll();
        /// return all custom radii defined in this table
        const CustomRadiiStorage& getAllCustom() const;
        /// convert all custom radii to a string in (A1:r1,A2:r2,...) format
        std::string toString(std::string separator=",") const;

    private:

        // data
        CustomRadiiStorage mcustomradius;

        // serialization
        friend class boost::serialization::access;
        template<class Archive>
            void serialize(Archive& ar, const unsigned int version)
        {
            ar & mcustomradius;
        }

};

typedef AtomRadiiTable::SharedPtr AtomRadiiTablePtr;

}   // namespace srreal
}   // namespace diffpy

// Serialization -------------------------------------------------------------

BOOST_SERIALIZATION_ASSUME_ABSTRACT(diffpy::srreal::AtomRadiiTable)

#endif  // ATOMRADIITABLE_HPP_INCLUDED
