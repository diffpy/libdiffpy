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
* class PairQuantity -- general implementation of pair quantity calculator
*
* $Id$
*
*****************************************************************************/

#ifndef PAIRQUANTITY_HPP_INCLUDED
#define PAIRQUANTITY_HPP_INCLUDED

#include <boost/unordered_set.hpp>
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/assume_abstract.hpp>
#include <boost/serialization/utility.hpp>
#include <boost/serialization/vector.hpp>

#include <diffpy/boostextensions/serialize_unordered_set.hpp>
#include <diffpy/srreal/PQEvaluator.hpp>
#include <diffpy/srreal/StructureAdapter.hpp>
#include <diffpy/srreal/PairQuantityUtils.hpp>
#include <diffpy/Attributes.hpp>

namespace diffpy {
namespace srreal {

class BaseBondGenerator;

class PairQuantity : public diffpy::Attributes
{
    public:

        // constructor
        PairQuantity();
        virtual ~PairQuantity()  { }

        // methods
        template <class T> const QuantityType& eval(const T&);
        const QuantityType& value() const;
        void mergeParallelValue(const QuantityType&);

        // configuration
        void setStructure(StructureAdapterPtr);
        template <class T> void setStructure(const T&);
        virtual void setRmin(double);
        const double& getRmin() const;
        virtual void setRmax(double);
        const double& getRmax() const;
        void setEvaluator(PQEvaluatorType evtp);
        void setupParallelRun(int cpuindex, int ncpu);
        int countSites() const;
        void maskAllPairs(bool mask);
        void invertMask();
        void setPairMask(int i, int j, bool mask);
        bool getPairMask(int i, int j);
        // FIXME: remove this when serialization works
        const boost::unordered_set< std::pair<int,int> >&
            getMaskData() const  { return minvertpairmask; }

    protected:

        friend class PQEvaluatorBasic;
        friend class PQEvaluatorOptimized;

        // methods
        virtual void resizeValue(size_t);
        virtual void resetValue();
        virtual void configureBondGenerator(BaseBondGenerator&) const;
        virtual void addPairContribution(const BaseBondGenerator&, int) { }

        // data
        QuantityType mvalue;
        StructureAdapterPtr mstructure;
        double mrmin;
        double mrmax;
        PQEvaluatorPtr mevaluator;
        int mcountsites;
        boost::unordered_set< std::pair<int,int> > minvertpairmask;
        bool mdefaultpairmask;

    private:

        // serialization
        friend class boost::serialization::access;
        template<class Archive>
            void serialize(Archive& ar, const unsigned int version)
        {
            ar & mvalue;
            ar & mstructure;
            ar & mrmin;
            ar & mrmax;
            ar & mevaluator;
            ar & mcountsites;
            ar & minvertpairmask;
            ar & mdefaultpairmask;
        }

};

// Template Public Methods ---------------------------------------------------

template <class T>
const QuantityType& PairQuantity::eval(const T& stru)
{
    this->setStructure(stru);
    mevaluator->updateValue(*this);
    return this->value();
}


template <class T>
void PairQuantity::setStructure(const T& stru)
{
    StructureAdapterPtr bstru = createStructureAdapter(stru);
    this->setStructure(bstru);
}


}   // namespace srreal
}   // namespace diffpy

// Serialization -------------------------------------------------------------

BOOST_SERIALIZATION_ASSUME_ABSTRACT(diffpy::srreal::PairQuantity)

#endif  // PAIRQUANTITY_HPP_INCLUDED
