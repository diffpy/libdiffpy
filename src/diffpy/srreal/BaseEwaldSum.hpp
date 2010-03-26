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
* class BaseEwaldSum -- base class for concrete Ewald sum calculators
*
* $Id$
*
*****************************************************************************/

#ifndef BASEEWALDSUM_HPP_INCLUDED
#define BASEEWALDSUM_HPP_INCLUDED

#include <diffpy/srreal/PairQuantity.hpp>

namespace diffpy {
namespace srreal {

class BaseEwaldSum : public PairQuantity
{
    public:

        // constructor
        BaseEwaldSum();

        // results
        /// F(Q) values on a Q-grid cut off at Qmin
        QuantityType getF() const;
        /// output Q-grid that is cut off at Qmin
        QuantityType getQgrid() const;
        /// F values on a full Q-grid starting at 0
        QuantityType getExtendedF() const;
        /// Full Q-grid starting at 0
        QuantityType getExtendedQgrid() const;

        // Q-range configuration
        void setQmin(double);
        const double& getQmin() const;
        void setQmax(double);
        const double& getQmax() const;
        void setQstep(double);
        const double& getQstep() const;

        // Summation cutoff due to Q-dependent pair scaling or large distance
        /// set relative cutoff value for Ewald sum contribution
        void setEwaldPrecision(double);
        /// return relative cutoff value for Ewald sum contribution
        const double& getEwaldPrecision() const;

    protected:

        // PairQuantity overloads
        virtual void resetValue();
        virtual void addPairContribution(const BaseBondGenerator&);

        // own methods
        virtual double pairScale(const BaseBondGenerator&, const double&) const;
        virtual double sfSiteAtQ(int, const double&) const;
        virtual double sfAverageAtQ(const double&) const;

    private:

        // methods
        /// number of ignored points in the extended r-grid below qmin
        int qminPoints() const;
        /// total number of points in mvalue
        int totalPoints() const;
        /// cache integer Q-grid limits from Qmin, Qmax, Qstep
        void cacheQpointsData();

        // data
        // configuration
        double mqmin;
        double mqmax;
        double mqstep;
        double mewaldprecision;
        struct {
            int qminpoints;
            int totalpoints;
        } mqpoints_cache;

};  // class BaseEwaldSum

}   // namespace srreal
}   // namespace diffpy

#endif  // BASEEWALDSUM_HPP_INCLUDED
