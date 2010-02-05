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
        /// difference between expected and calculated valence per each site
        QuantityType bvdiff() const;
        /// mean square difference of BVS from the expected values
        double bvmsdiff() const;
        /// root mean square difference of BVS from the expected values
        double bvrmsdiff() const;

        // access and configuration of BVS parameters
        void setBVParamTable(const BVParametersTable&);
        const BVParametersTable& getBVParamTable() const;
        // scattering factors lookup
        const BVParam& bvpar(int idx0, int idx1) const;

        // R-range configuration using the valence precision
        void setValencePrecision(double);
        double getValencePrecision() const;

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
