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
* class BaseDebyeSum -- base class for concrete Debye sum calculators
*
* $Id$
*
*****************************************************************************/

#ifndef BASEDEBYESUM_HPP_INCLUDED
#define BASEDEBYESUM_HPP_INCLUDED

#include <diffpy/srreal/PairQuantity.hpp>

namespace diffpy {
namespace srreal {

class BaseDebyeSum : public PairQuantity
{
    public:

        // constructor
        BaseDebyeSum();

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
        /// set relative cutoff value for Debye sum contribution
        void setDebyePrecision(double);
        /// return relative cutoff value for Debye sum contribution
        const double& getDebyePrecision() const;

    protected:

        // PairQuantity overloads
        virtual void resetValue();
        virtual void addPairContribution(const BaseBondGenerator&);

        // own methods
        virtual void setupPairScale(const BaseBondGenerator&);
        virtual double pairScale(const double& Q) const;
        virtual double sfSiteAtQ(int, const double& Q) const;
        /// number of ignored points in the extended r-grid below qmin
        int qminPoints() const;
        /// total number of points in mvalue
        int totalPoints() const;

    private:

        // methods
        /// cache integer Q-grid limits from Qmin, Qmax, Qstep
        void cacheQpointsData();
        /// cache structure factors data for a quick access during summation
        double sfSiteAtkQ(int siteidx, int kq) const;
        double sfAverageAtkQ(int kq) const;
        void cacheStructureData();

        // data
        // configuration
        double mqmin;
        double mqmax;
        double mqstep;
        double mdebyeprecision;
        struct {
            int qminpoints;
            int totalpoints;
        } mqpoints_cache;
        struct {
            std::vector< boost::shared_ptr<QuantityType> > sfsiteatkq;
            QuantityType sfaverageatkq;
        } mstructure_cache;

};  // class BaseDebyeSum

}   // namespace srreal
}   // namespace diffpy

#endif  // BASEDEBYESUM_HPP_INCLUDED
