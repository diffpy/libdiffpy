/*****************************************************************************
*
* libdiffpy         by DANSE Diffraction group
*                   Simon J. L. Billinge
*                   (c) 2011 The Trustees of Columbia University
*                   in the City of New York.  All rights reserved.
*
* File coded by:    Pavol Juhas
*
* See AUTHORS.txt for a list of people who contributed.
* See LICENSE_DANSE.txt for license information.
*
******************************************************************************
*
* class BondCalculator -- bond distance calculator
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

        // PairQuantity overloads
        virtual std::string getParallelData() const;

    protected:

        // PairQuantity overloads
        virtual void resetValue();
        virtual void addPairContribution(const BaseBondGenerator&, int);
        virtual void executeParallelMerge(const std::string& pdata);
        virtual void finishValue();

        // support for PQEvaluatorOptimized
        virtual void stashPartialValue();
        virtual void restorePartialValue();

        friend class BondOp;
        class BondEntry {

            public:

                double distance;
                int site0;
                int site1;
                double direction0;
                double direction1;
                double direction2;

            private:

                friend class boost::serialization::access;
                template<class Archive>
                void serialize(Archive& ar, const unsigned int version)
                {
                    ar & distance & site0 & site1;
                    ar & direction0 & direction1 & direction2;
                }

        };

        typedef std::vector<BondEntry> BondDataStorage;

    private:

        // serialization
        friend class boost::serialization::access;
        template<class Archive>
            void serialize(Archive& ar, const unsigned int version)
        {
            using boost::serialization::base_object;
            ar & base_object<PairQuantity>(*this);
            ar & mbonds;
            ar & mfilter_directions;
            ar & mfilter_degrees;
        }

        // methods
        int count() const;
        bool checkConeFilters(const R3::Vector& ru01) const;

        // data
        std::vector<R3::Vector> mfilter_directions;
        std::vector<double> mfilter_degrees;
        BondDataStorage mbonds;
        BondDataStorage mstashedbonds;
        BondDataStorage mpopbonds;
        BondDataStorage maddbonds;

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

// Serialization -------------------------------------------------------------

BOOST_CLASS_EXPORT_KEY(diffpy::srreal::BondCalculator)

#endif  // BONDCALCULATOR_HPP_INCLUDED
