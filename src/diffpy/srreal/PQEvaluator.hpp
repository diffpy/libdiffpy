/*****************************************************************************
*
* libdiffpy         by DANSE Diffraction group
*                   Simon J. L. Billinge
*                   (c) 2009 The Trustees of Columbia University
*                   in the City of New York.  All rights reserved.
*
* File coded by:    Pavol Juhas
*
* See AUTHORS.txt for a list of people who contributed.
* See LICENSE_DANSE.txt for license information.
*
******************************************************************************
*
* class PQEvaluatorBasic -- robust PairQuantity evaluator, the result
*     is always calculated from scratch.
*
* class PQEvaluatorOptimized -- optimized PairQuantity evaluator with fast
*     quantity updates
*
*****************************************************************************/


#ifndef PQEVALUATOR_HPP_INCLUDED
#define PQEVALUATOR_HPP_INCLUDED

#include <boost/shared_ptr.hpp>
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/assume_abstract.hpp>
#include <boost/serialization/export.hpp>

#include <diffpy/EventTicker.hpp>
#include <diffpy/srreal/QuantityType.hpp>
#include <diffpy/srreal/StructureAdapter.hpp>

namespace diffpy {
namespace srreal {

class PairQuantity;

/// shared pointer to PQEvaluatorBasic

typedef boost::shared_ptr<class PQEvaluatorBasic> PQEvaluatorPtr;

enum PQEvaluatorType {NONE, BASIC, OPTIMIZED};

class PQEvaluatorBasic
{
    public:

        friend
            PQEvaluatorPtr createPQEvaluator(PQEvaluatorType, PQEvaluatorPtr);
        // constructor
        PQEvaluatorBasic();
        virtual ~PQEvaluatorBasic()  { }

        // methods
        virtual PQEvaluatorType typeint() const;
        PQEvaluatorType typeintused() const;
        virtual void updateValue(PairQuantity&, StructureAdapterPtr);
        void useFullSum(bool flag);
        void setupParallelRun(int cpuindex, int ncpu);

    protected:

        // data
        /// flag for performing full double summation in updateValue
        bool musefullsum;
        /// zero-based index of this CPU
        int mcpuindex;
        /// total number of the CPU units
        int mncpu;
        /// ticker for recording when was the value updated
        eventticker::EventTicker mvalue_ticker;
        /// type of PQEvaluator that was actually used
        PQEvaluatorType mtypeused;

    private:

        // serialization
        friend class boost::serialization::access;
        template<class Archive>
            void serialize(Archive& ar, const unsigned int version)
        {
            ar & musefullsum & mcpuindex & mncpu & mvalue_ticker;
        }
};


class PQEvaluatorOptimized : public PQEvaluatorBasic
{
    public:

        // methods
        virtual PQEvaluatorType typeint() const;
        virtual void updateValue(PairQuantity&, StructureAdapterPtr);

    private:

        // data
        StructureAdapterPtr mlast_structure;

        // helper method
        void updateValueCompletely(PairQuantity&, StructureAdapterPtr);

        // serialization
        friend class boost::serialization::access;
        template<class Archive>
            void serialize(Archive& ar, const unsigned int version)
        {
            using boost::serialization::base_object;
            ar & base_object<PQEvaluatorBasic>(*this);
            ar & mlast_structure;
        }
};

// Factory function for PairQuantity evaluators ------------------------------

PQEvaluatorPtr createPQEvaluator(
        PQEvaluatorType pqtp, PQEvaluatorPtr pqevsrc=PQEvaluatorPtr());

}   // namespace srreal
}   // namespace diffpy

// Serialization -------------------------------------------------------------

BOOST_SERIALIZATION_ASSUME_ABSTRACT(diffpy::srreal::PQEvaluatorBasic)
BOOST_CLASS_EXPORT_KEY(diffpy::srreal::PQEvaluatorBasic)
BOOST_CLASS_EXPORT_KEY(diffpy::srreal::PQEvaluatorOptimized)

#endif  // PQEVALUATOR_HPP_INCLUDED
