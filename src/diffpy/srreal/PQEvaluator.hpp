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

* class PQEvaluatorOptimized -- optimized PairQuantity evaluator with fast
*     quantity updates
*
* $Id$
*
*****************************************************************************/


#ifndef PQEVALUATOR_HPP_INCLUDED
#define PQEVALUATOR_HPP_INCLUDED

#include <boost/shared_ptr.hpp>
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/assume_abstract.hpp>
#include <boost/serialization/shared_ptr.hpp>
#include <boost/serialization/export.hpp>

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
        PQEvaluatorBasic() : mcpuindex(0), mncpu(1)  { }
        virtual ~PQEvaluatorBasic()  { }

        // methods
        virtual PQEvaluatorType typeint() const;
        virtual void updateValue(PairQuantity& pq);
        void setupParallelRun(int cpuindex, int ncpu);

    protected:

        // data
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
            ar & mcpuindex & mncpu;
        }
};


class PQEvaluatorOptimized : public PQEvaluatorBasic
{
    public:

        // methods
        virtual PQEvaluatorType typeint() const;
        virtual void updateValue(PairQuantity& pq);

};


// Factory function for PairQuantity evaluators ------------------------------

PQEvaluatorPtr createPQEvaluator(PQEvaluatorType pqtp);


}   // namespace srreal
}   // namespace diffpy

// Serialization -------------------------------------------------------------

BOOST_SERIALIZATION_ASSUME_ABSTRACT(diffpy::srreal::PQEvaluatorBasic)

#endif  // PQEVALUATOR_HPP_INCLUDED
