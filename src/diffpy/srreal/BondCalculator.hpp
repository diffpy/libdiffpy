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
* class BondCalculator -- bond distance calculator
*
* $Id$
*
*****************************************************************************/

#ifndef BONDCALCULATOR_HPP_INCLUDED
#define BONDCALCULATOR_HPP_INCLUDED

#include <diffpy/srreal/PairQuantity.hpp>

namespace diffpy {
namespace srreal {

class BondCalculator : public PairQuantity
{
    public:

        // constructor
        BondCalculator();

        // methods
        template <class T> QuantityType operator()(const T&);
        QuantityType distances() const;
        std::vector<R3::Vector> directions() const;
        std::vector<int> sites0() const;
        std::vector<int> sites1() const;
        std::vector<std::string> types0() const;
        std::vector<std::string> types1() const;
        void filterCone(R3::Vector coneaxis, double degrees);
        void filterOff();

    protected:

        // PairQuantity overloads
        virtual void resetValue();
        virtual void addPairContribution(const BaseBondGenerator&, int);
        virtual void executeParallelMerge(const std::string& pdata);
        virtual void finishValue();

    private:

        // serialization
        friend class boost::serialization::access;
        template<class Archive>
            void serialize(Archive& ar, const unsigned int version)
        {
            using boost::serialization::base_object;
            ar & base_object<PairQuantity>(*this);
            ar & mfilter_directions;
            ar & mfilter_degrees;
        }

        // methods
        int count() const;
        bool checkConeFilters(const R3::Vector& ru01) const;

        // data
        std::vector<R3::Vector> mfilter_directions;
        std::vector<double> mfilter_degrees;

};

// Public Template Methods ---------------------------------------------------

template <class T>
QuantityType BondCalculator::operator()(const T& stru)
{
    this->eval(stru);
    return this->distances();
}


}   // namespace srreal
}   // namespace diffpy

#endif  // BONDCALCULATOR_HPP_INCLUDED
