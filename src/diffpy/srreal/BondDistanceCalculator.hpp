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
* class BondDistanceCalculator -- bond distance calculator
*
* $Id$
*
*****************************************************************************/

#ifndef BONDDISTANCECALCULATOR_HPP_INCLUDED
#define BONDDISTANCECALCULATOR_HPP_INCLUDED

#include <diffpy/srreal/PairQuantity.hpp>

namespace diffpy {
namespace srreal {

class BondDistanceCalculator : public PairQuantity
{
    public:

        // constructor
        BondDistanceCalculator();

        // methods
        template <class T> QuantityType operator()(const T&);
        virtual void mergeParallelValue(const QuantityType&);
        QuantityType distances() const;
        std::vector<R3::Vector> directions() const;
        std::vector<int> sites0() const;
        std::vector<int> sites1() const;

    protected:

        // PairQuantity overloads
        virtual void resetValue();
        virtual void addPairContribution(const BaseBondGenerator&, int);
        virtual void finishValue();
        int count() const;

    private:

        // serialization
        friend class boost::serialization::access;
        template<class Archive>
            void serialize(Archive& ar, const unsigned int version)
        {
            using boost::serialization::base_object;
            ar & base_object<PairQuantity>(*this);
        }

};

// Public Template Methods ---------------------------------------------------

template <class T>
QuantityType BondDistanceCalculator::operator()(const T& stru)
{
    this->eval(stru);
    return this->distances();
}


}   // namespace srreal
}   // namespace diffpy

#endif  // BONDDISTANCECALCULATOR_HPP_INCLUDED
