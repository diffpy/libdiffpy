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

namespace diffpy {
namespace srreal {

class BVParametersTable { 
    public:
        BVParametersTable* copy() const { return new BVParametersTable; }
};

class BVPar { };

class BVSCalculator : public PairQuantity
{
    public:

        // constructor
        BVSCalculator();

        // results
        /// mean square difference of BVS from the expected values
        double bvmsd() const;
        /// root mean square difference of BVS from the expected values
        double bvrmsd() const;

        // access and configuration of BVS parameters
        void setBVParamTable(const BVParametersTable&);
        const BVParametersTable& getBVParamTable() const;
        // scattering factors lookup
        const BVPar& bvpar(int idx0, int idx1) const;

    protected:

        // PairQuantity overloads
        virtual void resetValue();
        virtual void configureBondGenerator(BaseBondGenerator&);
        virtual void addPairContribution(const BaseBondGenerator&);

    private:

        // data
        // configuration
        boost::shared_ptr<BVParametersTable> mbvptable;

};  // class BVSCalculator

}   // namespace srreal
}   // namespace diffpy

#endif  // BVSCALCULATOR_HPP_INCLUDED
