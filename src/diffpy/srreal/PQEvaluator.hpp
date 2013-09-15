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

#include <diffpy/EventTicker.hpp>
#include <diffpy/srreal/QuantityType.hpp>
#include <diffpy/srreal/forwardtypes.hpp>

namespace diffpy {
namespace srreal {

class PairQuantity;

/// shared pointer to PQEvaluatorBasic

typedef boost::shared_ptr<class PQEvaluatorBasic> PQEvaluatorPtr;

enum PQEvaluatorType {BASIC, OPTIMIZED};

class PQEvaluatorBasic
{
    public:

        // constructor
        PQEvaluatorBasic();
        virtual ~PQEvaluatorBasic()  { }

        // methods
        virtual PQEvaluatorType typeint() const;
        virtual void updateValue(PairQuantity&, StructureAdapterConstPtr);
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

        // testing
        friend class TestPQEvaluator;
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
        virtual void updateValue(PairQuantity&, StructureAdapterConstPtr);

    private:

        // serialization
        friend class boost::serialization::access;
        template<class Archive>
            void serialize(Archive& ar, const unsigned int version)
        {
            using boost::serialization::base_object;
            ar & base_object<PQEvaluatorBasic>(*this);
        }
};

// Factory function for PairQuantity evaluators ------------------------------

PQEvaluatorPtr createPQEvaluator(PQEvaluatorType pqtp);

}   // namespace srreal
}   // namespace diffpy

// Serialization -------------------------------------------------------------

BOOST_SERIALIZATION_ASSUME_ABSTRACT(diffpy::srreal::PQEvaluatorBasic)

#endif  // PQEVALUATOR_HPP_INCLUDED
