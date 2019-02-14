/*****************************************************************************
*
* libdiffpy         by DANSE Diffraction group
*                   Simon J. L. Billinge
*                   (c) 2010 The Trustees of Columbia University
*                   in the City of New York.  All rights reserved.
*
* File coded by:    Pavol Juhas
*
* See AUTHORS.txt for a list of people who contributed.
* See LICENSE_DANSE.txt for license information.
*
******************************************************************************
*
* class BaseDebyeSum -- base class for concrete Debye sum calculators
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
    public PairQuantity,
    public PeakWidthModelOwner
{
    public:

        // constructor
        BaseDebyeSum();

        // PairQuantity overloads
        virtual eventticker::EventTicker& ticker() const;

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
        virtual void addPairContribution(const BaseBondGenerator&, int);
        // support for PQEvaluatorOptimized
        virtual void stashPartialValue();
        virtual void restorePartialValue();

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
            std::vector<int> typeofsite;
            std::vector<QuantityType> sftypeatkq;
            QuantityType sfaverageatkq;
            double totaloccupancy;
        } mstructure_cache;
        QuantityType mdbsumstash;

        // serialization
        friend class boost::serialization::access;
        template<class Archive>
            void serialize(Archive& ar, const unsigned int version)
        {
            using boost::serialization::base_object;
            ar & base_object<PairQuantity>(*this);
            ar & base_object<PeakWidthModelOwner>(*this);
            ar & mqmin;
            ar & mqmax;
            ar & mqstep;
            ar & mdebyeprecision;
            ar & mstructure_cache.typeofsite;
            ar & mstructure_cache.sftypeatkq;
            ar & mstructure_cache.sfaverageatkq;
            ar & mstructure_cache.totaloccupancy;
        }

};  // class BaseDebyeSum

}   // namespace srreal
}   // namespace diffpy

// Serialization -------------------------------------------------------------

BOOST_CLASS_EXPORT_KEY(diffpy::srreal::BaseDebyeSum)

#endif  // BASEDEBYESUM_HPP_INCLUDED
