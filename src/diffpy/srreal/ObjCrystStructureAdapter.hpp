/*****************************************************************************
*
* diffpy.srreal     by DANSE Diffraction group
*                   Simon J. L. Billinge
*                   (c) 2009 Trustees of the Columbia University
*                   in the City of New York.  All rights reserved.
*
* File coded by:    Chris Farrow
*
* See AUTHORS.txt for a list of people who contributed.
* See LICENSE.txt for license information.
*
******************************************************************************
*
* class ObjCrystStructureAdapter   
*   -- adapter to the Crystal class from ObjCryst++.
* class ObjCrystBondGenerator     
*   -- Generate bonds from periodic ObjCrystStructureAdapter.
* class ObjCrystMoleculeAdapter
*   -- adapter class for Molecule class from ObjCryst++.
* class ObjCrystMoleculeBondGenerator
*   -- Generate bonds from ObjCrystMoleculeAdapter
*
* $Id$
*
*****************************************************************************/

#ifndef OBJCRYSTSTRUCTUREADAPTER_HPP_INCLUDED
#define OBJCRYSTSTRUCTUREADAPTER_HPP_INCLUDED

#include <memory>
#include <set>
#include <vector>

#include <diffpy/srreal/R3linalg.hpp>
#include <diffpy/srreal/StructureAdapter.hpp>
#include <diffpy/srreal/Lattice.hpp>

#include <ObjCryst/ObjCryst/Crystal.h>
#include <ObjCryst/ObjCryst/Molecule.h>
#include <ObjCryst/ObjCryst/ScatteringPower.h>
#include <ObjCryst/ObjCryst/SpaceGroup.h>

namespace diffpy {
namespace srreal {


class PointsInSphere;

class ObjCrystStructureAdapter : public StructureAdapter
{
    friend class ObjCrystBondGenerator;

    public:

        // constructors
        ObjCrystStructureAdapter(const ObjCryst::Crystal&);

        // methods - overloaded
        virtual BaseBondGenerator* createBondGenerator() const;
        virtual int countSites() const;
        virtual double numberDensity() const;
        virtual const R3::Vector& siteCartesianPosition(int idx) const;
        virtual double siteOccupancy(int idx) const;
        virtual bool siteAnisotropy(int idx) const;
        virtual int siteMultiplicity(int idx) const;
        virtual const R3::Matrix& siteCartesianUij(int idx) const;
        virtual const std::string& siteAtomType(int idx) const;
        const Lattice& getLattice() const;


    private:

        typedef std::vector<R3::Vector> SymPosVec;
        typedef std::vector<R3::Matrix> SymUijVec;

        // methods - own
        void getUnitCell(const ObjCryst::Crystal&);

        // Tolerance on distance measurements.  Two sites are the same if
        // their fractional coordinates are within this tolerance
        static const double mtoler;
        // The asymmetric unit cell of ScatteringComponent instances
        std::vector< ObjCryst::ScatteringComponent > mvsc;
        // The symmetry-related positions of the asymmetric unit cell
        std::vector<SymPosVec> mvsym;
        // The Uij for scatterers. Stored in same order as mvsym.
        std::vector<SymUijVec> mvuij;
        // The Lattice instance needed by the bond generator
        Lattice mlattice;

};


class ObjCrystBondGenerator : public BaseBondGenerator
{

    public:

        // constructors
        ObjCrystBondGenerator(const ObjCrystStructureAdapter*);

        // loop control
        virtual void rewind();

        // data access
        virtual const R3::Matrix& Ucartesian1() const;

        // configuration
        virtual void setRmin(double);
        virtual void setRmax(double);

    protected:

        // methods
        virtual bool iterateSymmetry();
        virtual void rewindSymmetry();

        // data
        const ObjCrystStructureAdapter* mpstructure;
        size_t msymidx;

    private:

        // methods
        void updater1();

        // data
        std::auto_ptr<PointsInSphere> msphere;
};

// Adapter for ObjCryst::Molecule
//
// Molecules are always considered aperiodic. The anisotropic ADPs are treated
// as if in a cartesian cell. If this is not what is intended, pass the
// molecule as a scattering component within an ObjCryst::Crystal. 

class ObjCrystMoleculeAdapter : public StructureAdapter
{

    friend class ObjCrystMoleculeBondGenerator;

    public:

        // constructors
        ObjCrystMoleculeAdapter(const ObjCryst::Molecule&);

        // methods - overloaded
        virtual BaseBondGenerator* createBondGenerator() const;
        virtual int countSites() const;
        virtual double numberDensity() const;
        virtual const R3::Vector& siteCartesianPosition(int idx) const;
        virtual double siteOccupancy(int idx) const;
        virtual bool siteAnisotropy(int idx) const;
        virtual const R3::Matrix& siteCartesianUij(int idx) const;
        virtual const std::string& siteAtomType(int idx) const;
        const Lattice& getLattice() const;


    private:

        // The MolAtom instances
        std::vector< ObjCryst::MolAtom > mvatoms;
        // The positions of the scatterers. Same order as mvatoms.
        std::vector<R3::Vector> mvpos;
        // The Uij for scatterers. Same order as mvatoms.
        std::vector<R3::Matrix> mvuij;
        // The Lattice instance needed by the bond generator
        Lattice mlattice;

};


class ObjCrystMoleculeBondGenerator : public BaseBondGenerator
{
    public:

        // constructors
        ObjCrystMoleculeBondGenerator(const ObjCrystMoleculeAdapter*);

    protected:

        // The adapted structure
        const ObjCrystMoleculeAdapter* mpstructure;

};


namespace objcrystutil
{
    R3::Matrix getUij(const ObjCryst::ScatteringPower* sp);
} // namespace objcrystutil


inline
StructureAdapterPtr 
createStructureAdapter(const ObjCryst::Crystal& cryst)
{
    StructureAdapterPtr adapter(new ObjCrystStructureAdapter(cryst));
    return adapter;
}


inline
StructureAdapterPtr 
createStructureAdapter(const ObjCryst::Molecule& molecule)
{
    StructureAdapterPtr adapter(new ObjCrystMoleculeAdapter(molecule));
    return adapter;
}


}   // namespace srreal
}   // namespace diffpy

#endif  // OBJCRYSTSTRUCTUREADAPTER_HPP_INCLUDED
