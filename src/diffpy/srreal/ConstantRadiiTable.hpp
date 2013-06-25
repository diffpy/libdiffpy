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
* class ConstantRadiiTable -- concrete AtomRadiiTable for a uniform radius
*
*****************************************************************************/

#ifndef CONSTANTRADIITABLE_HPP_INCLUDED
#define CONSTANTRADIITABLE_HPP_INCLUDED

#include <diffpy/srreal/AtomRadiiTable.hpp>

namespace diffpy {
namespace srreal {

class ConstantRadiiTable : public AtomRadiiTable
{
    public:

        // constructor
        ConstantRadiiTable();

        // methods - HasClassRegistry
        SharedPtr create() const;
        SharedPtr clone() const;
        const std::string& type() const;
        // method overloads
        double tableLookup(const std::string& smbl) const;
        // methods specific for this class
        void setDefault(double);
        double getDefault() const;

    private:

        // data
        double mdefaultradius;

        // serialization
        friend class boost::serialization::access;
        template<class Archive>
            void serialize(Archive& ar, const unsigned int version)
        {
            using boost::serialization::base_object;
            ar & base_object<AtomRadiiTable>(*this);
            ar & mdefaultradius;
        }

};

}   // namespace srreal
}   // namespace diffpy

#endif  // CONSTANTRADIITABLE_HPP_INCLUDED
