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
* Shared types and small functions used by PairQuantity calculators.
*
* QuantityType -- type that stores PairQuantity results, an array of doubles
*   It is a unique derived class from vector<double> to avoid conflicts with
*   boost_python convertors in cctbx.
*
*****************************************************************************/

#ifndef QUANTITYTYPE_HPP_INCLUDED
#define QUANTITYTYPE_HPP_INCLUDED

#include <boost/serialization/base_object.hpp>
#include <boost/serialization/vector.hpp>

namespace diffpy {
namespace srreal {

class QuantityType : public std::vector<double>
{
    private:

        typedef std::vector<double> Base;

    public:

        // Constructors from std::vector
        QuantityType() : Base()  { }
        QuantityType(const Base& src) : Base(src)  { }
        explicit QuantityType(Base::size_type n) : Base(n)  { }
        explicit QuantityType(Base::size_type n, double x) : Base(n, x)  { }
        template <class Iter>
            QuantityType(Iter first, Iter last) : Base(first, last)  { }

    private:

        // serialization
        friend class boost::serialization::access;
        template<class Archive>
            void serialize(Archive& ar, const unsigned int version)
        {
            ar & boost::serialization::base_object<Base>(*this);
        }

};

}   // namespace srreal
}   // namespace diffpy

#endif  // QUANTITYTYPE_HPP_INCLUDED
