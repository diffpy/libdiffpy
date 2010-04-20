/*****************************************************************************
*
* diffpy.srreal     by DANSE Diffraction group
*                   Simon J. L. Billinge
*                   (c) 2010 Trustees of the Columbia University
*                   in the City of New York.  All rights reserved.
*
* File coded by:    Pavol Juhas
*
* See AUTHORS.txt for a list of people who contributed.
* See LICENSE.txt for license information.
*
******************************************************************************
*
* class BVSCalculator -- bond valence sums calculator
*
* $Id$
*
*****************************************************************************/

#ifndef BVSCALCULATOR_HPP_INCLUDED
#define BVSCALCULATOR_HPP_INCLUDED

#include <boost/shared_ptr.hpp>

#include <diffpy/srreal/PairQuantity.hpp>
#include <diffpy/srreal/BVParametersTable.hpp>

namespace diffpy {
namespace srreal {

class BVSCalculator : public PairQuantity
{
    public:

        // constructor
        BVSCalculator();

        // results
        /// expected valence per each site
        QuantityType valences() const;
        /// difference between expected and calculated absolute valence at
        /// each site.  Positive for underbonding, negative for overbonding.
        QuantityType bvdiff() const;
        /// mean square difference of BVS from the expected values
        double bvmsdiff() const;
        /// root mean square difference of BVS from the expected values
        double bvrmsdiff() const;

        // access and configuration of BVS parameters
        void setBVParamTable(const BVParametersTable&);
        const BVParametersTable& getBVParamTable() const;

        // R-range configuration using the valence precision
        /// set cutoff value for bond valence contributions
        void setValencePrecision(double);
        /// return cutoff value for bond valence contributions
        double getValencePrecision() const;

    protected:

        // PairQuantity overloads
        virtual void resetValue();
        virtual void configureBondGenerator(BaseBondGenerator&);
        virtual void addPairContribution(const BaseBondGenerator&, int);

    private:

        // methods
        void cacheStructureData();
        /// rmax necessary for achieving the specified valence precision
        double rmaxFromPrecision(double) const;

        // data
        // configuration
        boost::shared_ptr<BVParametersTable> mbvptable;
        double mvalenceprecision;
        // cache
        struct {
            std::vector<std::string> baresymbols;
            std::vector<int> valences;
            QuantityType multiplicities;
            QuantityType occupancies;
            double total_occupancy;
        } mstructure_cache;

};  // class BVSCalculator

}   // namespace srreal
}   // namespace diffpy

#endif  // BVSCALCULATOR_HPP_INCLUDED
