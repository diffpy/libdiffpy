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
#include <diffpy/srreal/PeakWidthModel.hpp>
#include <diffpy/srreal/PDFUtils.hpp>

namespace diffpy {
namespace srreal {

class BaseDebyeSum :
    public virtual PairQuantity,
    public virtual PeakWidthModelOwner
{
    public:

        // constructor
        BaseDebyeSum();

        // results
        /// F values on a full Q-grid starting at 0
        QuantityType getF() const;

        // Q-range methods
        /// Full Q-grid starting at 0
        QuantityType getQgrid() const;
        // Q-range configuration
        void setQmin(double);
        const double& getQmin() const;
        void setQmax(double);
        const double& getQmax() const;
        virtual void setQstep(double);
        const double& getQstep() const;

        // Summation cutoff due to Q-dependent pair scaling or large distance
        /// set relative cutoff value for Debye sum contribution
        void setDebyePrecision(double);
        /// return relative cutoff value for Debye sum contribution
        const double& getDebyePrecision() const;

    protected:

        // PairQuantity overloads
        virtual void resetValue();
        virtual void addPairContribution(const BaseBondGenerator&, int);

        // own methods
        virtual double sfSiteAtQ(int, const double& Q) const;

    private:

        // methods
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
            std::vector< boost::shared_ptr<QuantityType> > sfsiteatkq;
            QuantityType sfaverageatkq;
            double totaloccupancy;
        } mstructure_cache;

};  // class BaseDebyeSum

}   // namespace srreal
}   // namespace diffpy

#endif  // BASEDEBYESUM_HPP_INCLUDED
