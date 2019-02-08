/*****************************************************************************
*
* libdiffpy         by DANSE Diffraction group
*                   Simon J. L. Billinge
*                   (c) 2011 The Trustees of Columbia University
*                   in the City of New York.  All rights reserved.
*
* File coded by:    Pavol Juhas
*
* See AUTHORS.txt for a list of people who contributed.
* See LICENSE_DANSE.txt for license information.
*
******************************************************************************
*
* class OverlapCalculator -- calculator of atom radii overlaps
*
*****************************************************************************/

#ifndef OVERLAPCALCULATOR_HPP_INCLUDED
#define OVERLAPCALCULATOR_HPP_INCLUDED

#include <boost/serialization/list.hpp>

#include <diffpy/srreal/PairQuantity.hpp>
#include <diffpy/srreal/AtomRadiiTable.hpp>

namespace diffpy {
namespace srreal {

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
        SiteIndices sites0() const;
        /// indices of the second site for all overlapping pairs
        SiteIndices sites1() const;
        /// atom symbols for the first site of all overlapping pairs
        std::vector<std::string> types0() const;
        /// atom symbols for the second site of all overlapping pairs
        std::vector<std::string> types1() const;
        /// sum of squared overlaps at each structure site
        QuantityType siteSquareOverlaps() const;
        /// sum of squared overlaps for all atoms in the structure
        double totalSquareOverlap() const;
        /// mean square per one atom in the structure
        double meanSquareOverlap() const;
        /// difference in the totalSquareOverlap for a flip of two sites
        double flipDiffTotal(int i, int j) const;
        /// difference in the meanSquareOverlap for a flip of two sites
        double flipDiffMean(int i, int j) const;
        /// gradients of totalSquareOverlap at each site in the structure
        std::vector<R3::Vector> gradients() const;
        /// indices of the neighboring sites
        std::unordered_set<int> getNeighborSites(int i) const;
        /// coordination number at each site of the structure
        QuantityType coordinations() const;
        /// coordination number split per each type of neighboring atoms
        std::unordered_map<std::string,double>
            coordinationByTypes(int i) const;
        /// sets of site inidices per each sites neighborhood in the structure
        std::vector< std::unordered_set<int> > neighborhoods() const;

        // access and configuration of the atom radii
        void setAtomRadiiTable(AtomRadiiTablePtr);
        void setAtomRadiiTableByType(const std::string& tp);
        AtomRadiiTablePtr& getAtomRadiiTable();
        const AtomRadiiTablePtr& getAtomRadiiTable() const;

        /// effective rmax value, usually a double of the maximum atom radius.
        double getRmaxUsed() const;


    protected:

        // PairQuantity overloads
        virtual void resetValue();
        virtual void configureBondGenerator(BaseBondGenerator&) const;
        virtual void addPairContribution(const BaseBondGenerator&, int);
        virtual void executeParallelMerge(const std::string&);

    private:

        // types
        enum OverlapFlag {ALLVALUES, OVERLAPPING};
        typedef std::unordered_map<int, std::list<int> > NeighborIdsStorage;

        // methods
        int count() const;
        QuantityType subvector(int offset, OverlapFlag flag) const;
        const double& subvalue(int offset, int index) const;
        const R3::Vector& subdirection(int index) const;
        double suboverlap(int index, int iflip=0, int jflip=0) const;
        void cacheStructureData();
        const std::list<int>& getNeighborIds(int i) const;

        // data
        AtomRadiiTablePtr matomradiitable;
        mutable NeighborIdsStorage mneighborids;
        mutable bool mneighborids_cached;
        // cache
        struct {
            QuantityType siteradii;
            double maxseparation;
        } mstructure_cache;

        // serialization
        friend class boost::serialization::access;
        template<class Archive>
            void serialize(Archive& ar, const unsigned int version)
        {
            using boost::serialization::base_object;
            ar & base_object<PairQuantity>(*this);
            ar & matomradiitable;
            ar & mneighborids;
            ar & mstructure_cache.siteradii;
            ar & mstructure_cache.maxseparation;
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

// Serialization -------------------------------------------------------------

BOOST_CLASS_EXPORT_KEY(diffpy::srreal::OverlapCalculator)

#endif  // OVERLAPCALCULATOR_HPP_INCLUDED
