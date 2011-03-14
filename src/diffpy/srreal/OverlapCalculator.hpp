/*****************************************************************************
*
* diffpy.srreal     by DANSE Diffraction group
*                   Simon J. L. Billinge
*                   (c) 2011 Trustees of the Columbia University
*                   in the City of New York.  All rights reserved.
*
* File coded by:    Pavol Juhas
*
* See AUTHORS.txt for a list of people who contributed.
* See LICENSE.txt for license information.
*
******************************************************************************
*
* class OverlapCalculator -- calculator of atom radii overlaps.
*
* $Id$
*
*****************************************************************************/

#ifndef OVERLAPCALCULATOR_HPP_INCLUDED
#define OVERLAPCALCULATOR_HPP_INCLUDED

#include <diffpy/srreal/PairQuantity.hpp>

namespace diffpy {
namespace srreal {

// FIXME: temporary types
class AtomRadiiTable
{
    public:
        double lookup(const std::string& smbl)  { return 1.0; }
};
typedef boost::shared_ptr<AtomRadiiTable> AtomRadiiTablePtr;


class OverlapCalculator : public PairQuantity
{
    public:

        // constructor
        OverlapCalculator();

        // methods
        /// return siteSquareOverlaps for the specified structure
        template <class T> QuantityType operator()(const T&);
        /// overlap values for all overlapping sites in the structure
        QuantityType overlaps() const;
        /// distances all overlapping sites in the structure
        QuantityType distances() const;
        /// directions between the centers of overlapping sites
        std::vector<R3::Vector> directions() const;
        /// indices of the first site for all overlapping pairs
        std::vector<int> sites0() const;
        /// indices of the second site for all overlapping pairs
        std::vector<int> sites1() const;
        /// sum of squared overlaps at each structure site
        QuantityType siteSquareOverlaps() const;
        /// sum of squared overlaps for all atoms in the structure
        double totalSquareOverlap() const;
        /// difference in the totalSquareOverlap for a flip of two sites
        double totalFlipDiff(int i, int j) const;
        /// gradients of totalFlipDiff at each site in the structure
        std::vector<R3::Vector> gradients() const;
        /// mean square per one atom in the structure
        double msoverlap() const;
        /// root mean square atom-radii overlap per one atom
        double rmsoverlap() const;

        // access and configuration of the atom radii
        void setAtomRadiiTable(AtomRadiiTablePtr);
        AtomRadiiTablePtr& getAtomRadiiTable();
        const AtomRadiiTablePtr& getAtomRadiiTable() const;

        /// effective rmax value, usually a double of the maximum atom radius.
        double getRmaxUsed() const;


    protected:

        // PairQuantity overloads
        virtual void resetValue();
        virtual void addPairContribution(const BaseBondGenerator&, int);
        virtual void executeParallelMerge(const QuantityType&);
        virtual void finishValue();

    private:

        // methods
        void cacheStructureData();

        // data
        AtomRadiiTablePtr matomradiitable;
        // cache
        struct {
            std::vector<double> siteradii;
        } mstructure_cache;

        // serialization
        friend class boost::serialization::access;
        template<class Archive>
            void serialize(Archive& ar, const unsigned int version)
        {
            using boost::serialization::base_object;
            ar & base_object<PairQuantity>(*this);
            // FIXME: ar & matomradiitable;
            ar & mstructure_cache.siteradii;
        }

};

// Public Template Methods ---------------------------------------------------

template <class T>
QuantityType OverlapCalculator::operator()(const T& stru)
{
    this->eval(stru);
    return this->siteSquareOverlaps();
}


}   // namespace srreal
}   // namespace diffpy

#endif  // OVERLAPCALCULATOR_HPP_INCLUDED
