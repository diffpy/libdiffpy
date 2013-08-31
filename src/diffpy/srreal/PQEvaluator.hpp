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
#include <boost/serialization/shared_ptr.hpp>
#include <boost/serialization/export.hpp>

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
        virtual void reset()  { }
        virtual void updateValue(PairQuantity& pq);
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

    private:

        // serialization
        friend class boost::serialization::access;
        template<class Archive>
            void serialize(Archive& ar, const unsigned int version)
        {
            ar & musefullsum & mcpuindex & mncpu;
        }
};


class PQEvaluatorOptimized : public PQEvaluatorBasic
{
    public:

        // methods
        virtual PQEvaluatorType typeint() const;
        virtual void reset();
        virtual void updateValue(PairQuantity& pq);

    private:

        // data
        /// last structure used for evaluation
        StructureAdapterConstPtr mstructure0;
        /// last evaluated value
        QuantityType mvalue0;
};

// Factory function for PairQuantity evaluators ------------------------------

PQEvaluatorPtr createPQEvaluator(PQEvaluatorType pqtp);

}   // namespace srreal
}   // namespace diffpy

// Serialization -------------------------------------------------------------

BOOST_SERIALIZATION_ASSUME_ABSTRACT(diffpy::srreal::PQEvaluatorBasic)

#endif  // PQEVALUATOR_HPP_INCLUDED
