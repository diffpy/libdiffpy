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

#include <boost/serialization/base_object.hpp>
#include <boost/serialization/assume_abstract.hpp>
#include <boost/serialization/utility.hpp>
#include <boost/serialization/vector.hpp>

#include <diffpy/boostextensions/serialize_unordered_set.hpp>
#include <diffpy/boostextensions/serialize_unordered_map.hpp>
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

        // class constants
        static const int ALLATOMSINT;
        static const std::string ALLATOMSSTR;

        // constructor
        PairQuantity();
        virtual ~PairQuantity()  { }

        // methods
        const QuantityType& eval();
        template <class T> const QuantityType& eval(const T&);
        const QuantityType& value() const;
        void mergeParallelData(const std::string& pdata, int ncpu);
        virtual std::string getParallelData() const;

        // configuration
        template <class T> void setStructure(const T&);
        void setStructure(StructureAdapterPtr);
        const StructureAdapterPtr& getStructure() const;
        virtual void setRmin(double);
        const double& getRmin() const;
        virtual void setRmax(double);
        const double& getRmax() const;
        void setEvaluator(PQEvaluatorType evtp);
        void setupParallelRun(int cpuindex, int ncpu);
        void maskAllPairs(bool mask);
        void invertMask();
        void setPairMask(int i, int j, bool mask);
        bool getPairMask(int i, int j) const;
        void setTypeMask(std::string, std::string, bool mask);
        bool getTypeMask(const std::string&, const std::string&) const;

    protected:

        friend class PQEvaluatorBasic;
        friend class PQEvaluatorOptimized;

        // methods
        virtual void resizeValue(size_t);
        virtual void resetValue();
        virtual void configureBondGenerator(BaseBondGenerator&) const;
        virtual void addPairContribution(const BaseBondGenerator&, int) { }
        virtual void executeParallelMerge(const std::string& pdata);
        virtual void finishValue() { }
        int countSites() const;

        // data
        typedef boost::unordered_map<
            std::pair<std::string,std::string>, bool> TypeMaskStorage;
        QuantityType mvalue;
        StructureAdapterPtr mstructure;
        double mrmin;
        double mrmax;
        PQEvaluatorPtr mevaluator;
        bool mdefaultpairmask;
        boost::unordered_set< std::pair<int,int> > minvertpairmask;
        boost::unordered_map<int, bool> msiteallmask;
        TypeMaskStorage mtypemask;
        int mmergedvaluescount;

    private:

        // methods
        void updateMaskData();
        void setPairMaskValue(int i, int j, bool mask);

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
            ar & mdefaultpairmask;
            ar & minvertpairmask;
            ar & msiteallmask;
            ar & mtypemask;
            ar & mmergedvaluescount;
        }

};

// Template Public Methods ---------------------------------------------------

template <class T>
const QuantityType& PairQuantity::eval(const T& stru)
{
    this->setStructure(stru);
    mevaluator->updateValue(*this);
    this->finishValue();
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
