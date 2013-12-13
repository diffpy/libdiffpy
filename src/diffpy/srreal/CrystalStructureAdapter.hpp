/*****************************************************************************
*
* diffpy.srreal     Complex Modeling Initiative
*                   Pavol Juhas
*                   (c) 2013 Brookhaven National Laboratory,
*                   Upton, New York.  All rights reserved.
*
* File coded by:    Pavol Juhas
*
* See AUTHORS.txt for a list of people who contributed.
* See LICENSE.txt for license information.
*
******************************************************************************
*
* class CrystalStructureAdapter -- universal adapter for crystal structure
*     composed of asymmetric unit cell and list of symmetry operations in
*     the space group
*
* class CrystalStructureBondGenerator -- bond generator
*
*****************************************************************************/

#ifndef CRYSTALSTRUCTUREADAPTER_HPP_INCLUDED
#define CRYSTALSTRUCTUREADAPTER_HPP_INCLUDED

#include <diffpy/srreal/PeriodicStructureAdapter.hpp>

namespace diffpy {
namespace srreal {

/// store rotation matrix and translation vector for a symmetry operation
class SymOpRotTrans
{
    public:

        R3::Matrix R;
        R3::Vector t;

    private:

        // serialization
        friend class boost::serialization::access;
        template<class Archive>
            void serialize(Archive& ar, const unsigned int version)
        {
            ar & R & t;
        }

};

// Comparison functions

inline
bool operator==(const SymOpRotTrans& op0, const SymOpRotTrans& op1)
{
    return (op0.R == op1.R) && (op0.t == op1.t);
}

inline
bool operator!=(const SymOpRotTrans& op0, const SymOpRotTrans& op1)
{
    return !(op0 == op1);
}


class CrystalStructureAdapter : public PeriodicStructureAdapter
{
    friend class CrystalStructureBondGenerator;

    public:

        typedef std::vector<SymOpRotTrans> SymOpVector;
        // constructor
        CrystalStructureAdapter();

        // methods - overloaded
        virtual StructureAdapterPtr clone() const;
        virtual BaseBondGeneratorPtr createBondGenerator() const;
        virtual int siteMultiplicity(int idx) const;
        virtual StructureDifference diff(StructureAdapterConstPtr other) const;

        // methods - own
        void setSymmetryPrecision(double eps);
        const double& getSymmetryPrecision() const;
        int countSymOps() const;
        void clearSymOps();
        void addSymOp(const SymOpRotTrans&);
        void addSymOp(const R3::Matrix& R, const R3::Vector& t);
        const SymOpRotTrans& getSymOp(int i) const;
        /// return all symmetry equivalent atoms in the unit cell for site i
        const AtomVector& getEquivalentAtoms(int idx) const;
        /// return all symmetry related atoms in fractional coordinates
        AtomVector expandLatticeAtom(const Atom&) const;
        void updateSymmetryPositions() const;

    private:

        // data
        /// array of symmetry operations
        SymOpVector msymops;
        double msymmetry_precision;
        mutable std::vector<AtomVector> msymatoms;
        mutable bool msymmetry_cached;

        // symmetry helpers
        /// return index of AtomVector atom at an equal position or -1
        int findEqualPosition(const AtomVector&, const Atom&) const;
        /// fuzzy check if symmetry positions are up to date
        /// this only detects addition or removal of atom in the asymmetric
        /// unit, but does not check for changes in atom positions.
        bool isSymmetryCached() const;

        // comparison
        friend bool operator==(
                const CrystalStructureAdapter&, const CrystalStructureAdapter&);

        // serialization
        friend class boost::serialization::access;
        template<class Archive>
            void serialize(Archive& ar, const unsigned int version)
        {
            ar & boost::serialization::base_object<PeriodicStructureAdapter>(*this);
            ar & msymops;
            ar & msymmetry_precision;
            ar & msymatoms;
            ar & msymmetry_cached;
        }

};


typedef boost::shared_ptr<CrystalStructureAdapter> CrystalStructureAdapterPtr;

// Comparison functions

bool operator==(const CrystalStructureAdapter&, const CrystalStructureAdapter&);
bool operator!=(const CrystalStructureAdapter&, const CrystalStructureAdapter&);


class CrystalStructureBondGenerator : public PeriodicStructureBondGenerator
{
    public:

        // constructors
        CrystalStructureBondGenerator(StructureAdapterConstPtr);

        // configuration
        virtual void selectAnchorSite(int);

        // data access
        virtual const R3::Matrix& Ucartesian1() const;

    protected:

        // methods
        virtual bool iterateSymmetry();
        virtual void rewindSymmetry();
        virtual void getNextBond();
        virtual void updater1();

        // data
        const CrystalStructureAdapter* mcstructure;
        size_t msymidx;
        const R3::Matrix* mpuc1;

    private:

        typedef CrystalStructureAdapter::AtomVector AtomVector;

        // methods
        const AtomVector& symatoms(int idx);

};

}   // namespace srreal
}   // namespace diffpy

BOOST_CLASS_EXPORT_KEY(diffpy::srreal::CrystalStructureAdapter)

#endif  // CRYSTALSTRUCTUREADAPTER_HPP_INCLUDED
