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

#include <boost/scoped_ptr.hpp>
#include <boost/shared_ptr.hpp>

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
        const QuantityType& eval(const StructureAdapter&);
        template <class T> const QuantityType& eval(const T&);
        const QuantityType& value() const;

        // configuration
        virtual void setRmin(double);
        const double& getRmin() const;
        virtual void setRmax(double);
        const double& getRmax() const;
        void setEvaluator(PQEvaluatorType evtp);
        int countSites() const;
        void maskAllPairs(bool mask);
        void maskSitePair(int i, int j, bool mask);
        bool getPairMask(int i, int j);

    protected:

        friend class PQEvaluatorBasic;
        friend class PQEvaluatorOptimized;

        // methods
        virtual void resizeValue(size_t);
        virtual void resetValue();
        virtual void configureBondGenerator(BaseBondGenerator&);
        virtual void addPairContribution(const BaseBondGenerator&, int) { }

        // data
        QuantityType mvalue;
        const StructureAdapter* mstructure;
        double mrmin;
        double mrmax;
        boost::shared_ptr<PQEvaluatorBasic> mevaluator;
        int mcountsites;
        std::set< std::pair<int,int> > minvertpairmask;
        bool mdefaultpairmask;

};

// Template Public Methods ---------------------------------------------------

template <class T>
const QuantityType& PairQuantity::eval(const T& stru)
{
    boost::scoped_ptr<StructureAdapter> bstru(createPQAdapter(stru));
    return this->eval(*bstru);
}


}   // namespace srreal
}   // namespace diffpy

#endif  // PAIRQUANTITY_HPP_INCLUDED
