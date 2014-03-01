/*****************************************************************************
*
* libdiffpy         Complex Modeling Initiative
*                   (c) 2013 Brookhaven Science Associates,
*                   Brookhaven National Laboratory.
*                   All rights reserved.
*
* File coded by:    Pavol Juhas
*
* See AUTHORS.txt for a list of people who contributed.
* See LICENSE.txt for license information.
*
******************************************************************************
*
* class PeriodicStructureAdapter -- universal adapter for structure with
*     periodic boundary conditions that has no space group symmetry
*
* class PeriodicStructureBondGenerator -- bond generator
*
*****************************************************************************/

#ifndef PERIODICSTRUCTUREADAPTER_HPP_INCLUDED
#define PERIODICSTRUCTUREADAPTER_HPP_INCLUDED

#include <boost/scoped_ptr.hpp>

#include <diffpy/srreal/forwardtypes.hpp>
#include <diffpy/srreal/AtomicStructureAdapter.hpp>
#include <diffpy/srreal/Lattice.hpp>

namespace diffpy {
namespace srreal {

class PointsInSphere;

class PeriodicStructureAdapter : public AtomicStructureAdapter
{
    public:

        // methods - overloaded
        virtual StructureAdapterPtr clone() const;
        virtual BaseBondGeneratorPtr createBondGenerator() const;
        virtual double numberDensity() const;
        virtual StructureDifference diff(StructureAdapterConstPtr other) const;

        // methods - own
        void setLatPar(
                double a, double b, double c,
                double alphadeg, double betadeg, double gammadeg);
        const Lattice& getLattice() const;
        void toCartesian(Atom&) const;
        void toFractional(Atom&) const;

    private:

        // data
        Lattice mlattice;

        // serialization
        friend class boost::serialization::access;
        template<class Archive>
            void serialize(Archive& ar, const unsigned int version)
        {
            ar & boost::serialization::base_object<AtomicStructureAdapter>(*this);
            ar & mlattice;
        }

};

typedef boost::shared_ptr<PeriodicStructureAdapter> PeriodicStructureAdapterPtr;

// Comparison functions

bool operator==(const PeriodicStructureAdapter&, const PeriodicStructureAdapter&);
bool operator!=(const PeriodicStructureAdapter&, const PeriodicStructureAdapter&);


class PeriodicStructureBondGenerator : public BaseBondGenerator
{
    public:

        // constructors
        PeriodicStructureBondGenerator(StructureAdapterConstPtr);

        // methods
        // loop control
        virtual void rewind();

        // configuration
        virtual void selectAnchorSite(int);
        virtual void setRmin(double);
        virtual void setRmax(double);

    protected:

        // data
        const PeriodicStructureAdapter* mpstructure;
        boost::scoped_ptr<PointsInSphere> msphere;
        R3::Vector mrcsphere;

        // methods
        virtual bool iterateSymmetry();
        virtual void rewindSymmetry();
        virtual void getNextBond();
        virtual void updater1();

    private:

        // data
        std::vector<R3::Vector> mcartesian_positions_uc;
};

}   // namespace srreal
}   // namespace diffpy

// Serialization -------------------------------------------------------------

BOOST_CLASS_EXPORT_KEY(diffpy::srreal::PeriodicStructureAdapter)

#endif  // PERIODICSTRUCTUREADAPTER_HPP_INCLUDED
