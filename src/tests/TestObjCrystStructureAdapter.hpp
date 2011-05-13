/*****************************************************************************
*
* diffpy.srreal     by DANSE Diffraction group
*                   Simon J. L. Billinge
*                   (c) 2009 Trustees of the Columbia University
*                   in the City of New York.  All rights reserved.
*
* File coded by:    Pavol Juhas
*
* See AUTHORS.txt for a list of people who contributed.
* See LICENSE.txt for license information.
*
******************************************************************************
*
* class TestObjCrystStructureAdapter -- unit tests for an adapter
*     to the ObjCryst::Crystal class
*
* class TestObjCrystStructureBondGenerator -- unit tests for bond generator
*
*****************************************************************************/

#include <typeinfo>
#include <cxxtest/TestSuite.h>

#include <diffpy/srreal/ObjCrystStructureAdapter.hpp>
#include <diffpy/srreal/PointsInSphere.hpp>
#include <diffpy/serialization.hpp>
#include "objcryst_helpers.hpp"
#include <ObjCryst/ObjCryst/Crystal.h>
#include <ObjCryst/ObjCryst/Molecule.h>
#include <ObjCryst/ObjCryst/ScatteringPower.h>

using namespace std;
using namespace diffpy::srreal;
using ObjCryst::Crystal;
using ObjCryst::Molecule;
using ObjCryst::ScatteringPowerAtom;

// Local Helpers -------------------------------------------------------------

namespace {

int countBonds(BaseBondGenerator& bnds)
{
    int rv = 0;
    for (bnds.rewind(); !bnds.finished(); bnds.next())  ++rv;
    return rv;
}


template <class Tstru, class Tbnds>
double testmsd0(const Tstru& stru, const Tbnds& bnds)
{
    bool anisotropy0 = stru->siteAnisotropy(bnds->site0());
    double rv = meanSquareDisplacement(
            bnds->Ucartesian0(), bnds->r01(), anisotropy0);
    return rv;
}


template <class Tstru, class Tbnds>
double testmsd1(const Tstru& stru, const Tbnds& bnds)
{
    bool anisotropy1 = stru->siteAnisotropy(bnds->site1());
    double rv = meanSquareDisplacement(
            bnds->Ucartesian1(), bnds->r01(), anisotropy1);
    return rv;
}

// Defined below
Molecule* makeC60Molecule();

}   // namespace

//////////////////////////////////////////////////////////////////////////////
// class TestObjCrystStructureAdapter
//////////////////////////////////////////////////////////////////////////////

class TestObjCrystStructureAdapter : public CxxTest::TestSuite
{
    private:

        auto_ptr<Crystal> mcryst_ni;
        StructureAdapterPtr m_ni;
        auto_ptr<Crystal> mcryst_kbise;
        StructureAdapterPtr m_kbise;
        auto_ptr<Crystal> mcryst_catio3;
        StructureAdapterPtr m_catio3;

    public:

        void setUp()
        {
            if (!m_ni.get())
            {
                mcryst_ni.reset(loadTestCrystal("Ni.cif"));
                m_ni = createStructureAdapter(*mcryst_ni);
            }
            if (!m_kbise.get())
            {
                mcryst_kbise.reset(loadTestCrystal("alpha_K2Bi8Se13.cif"));
                m_kbise = createStructureAdapter(*mcryst_kbise);
            }
            if (!m_catio3.get())
            {
                mcryst_catio3.reset(loadTestCrystal("icsd_62149.cif"));
                m_catio3 = createStructureAdapter(*mcryst_catio3);
            }
        }


        void test_typeid()
        {
            TS_ASSERT(typeid(ObjCrystStructureAdapter) == typeid(*m_ni));
        }


        void test_countSites()
        {
            TS_ASSERT_EQUALS(1, m_ni->countSites());
            TS_ASSERT_EQUALS(12, m_kbise->countSites());
            TS_ASSERT_EQUALS(4, m_catio3->countSites());
        }


        void test_totalOccupancy()
        {
            TS_ASSERT_EQUALS(4.0, m_ni->totalOccupancy());
            TS_ASSERT_EQUALS(23.0, m_kbise->totalOccupancy());
            TS_ASSERT_EQUALS(20.0, m_catio3->totalOccupancy());
        }


        void test_numberDensity()
        {
            const double eps = 1.0e-7;
            TS_ASSERT_DELTA(0.0914114, m_ni->numberDensity(), eps);
            TS_ASSERT_DELTA(0.0335565, m_kbise->numberDensity(), eps);
            TS_ASSERT_DELTA(0.0894126, m_catio3->numberDensity(), eps);
        }


        void test_siteCartesianPosition()
        {
            const double eps = 1.0e-5;
            R3::Vector rCa = m_catio3->siteCartesianPosition(0);
            TS_ASSERT_DELTA(5.34323, rCa[0], eps);
            TS_ASSERT_DELTA(0.19603, rCa[1], eps);
            TS_ASSERT_DELTA(1.91003, rCa[2], eps);
        }


        void test_siteAnisotropy()
        {
            for (int i = 0; i < m_ni->countSites(); ++i)
            {
                TS_ASSERT_EQUALS(false, m_ni->siteAnisotropy(i));
            }
            for (int i = 0; i < m_catio3->countSites(); ++i)
            {
                TS_ASSERT_EQUALS(true, m_catio3->siteAnisotropy(i));
            }
        }


        void test_siteCartesianUij()
        {
            // nickel should have all Uij equal zero.
            const double* puij = &(m_ni->siteCartesianUij(0).data()[0]);
            for (int i = 0; i < 9; ++i)
            {
                TS_ASSERT_EQUALS(0.0, puij[i]);
            }
            // check CaTiO3 values
            const R3::Matrix& UTi = m_catio3->siteCartesianUij(1);
            const double eps = 1e-8;
            TS_ASSERT_DELTA(0.0052, UTi(0,0), eps);
            TS_ASSERT_DELTA(0.0049, UTi(1,1), eps);
            TS_ASSERT_DELTA(0.0049, UTi(2,2), eps);
            TS_ASSERT_DELTA(0.00016, UTi(0,1), eps);
            TS_ASSERT_DELTA(0.00016, UTi(1,0), eps);
            TS_ASSERT_DELTA(0.00001, UTi(0,2), eps);
            TS_ASSERT_DELTA(0.00001, UTi(2,0), eps);
            TS_ASSERT_DELTA(0.00021, UTi(1,2), eps);
            TS_ASSERT_DELTA(0.00021, UTi(2,1), eps);
        }


        void test_siteAtomType()
        {
            TS_ASSERT_EQUALS(string("Ni"), m_ni->siteAtomType(0));
            TS_ASSERT_EQUALS(string("K1+"), m_kbise->siteAtomType(0));
            TS_ASSERT_EQUALS(string("Bi3+"), m_kbise->siteAtomType(1));
            TS_ASSERT_EQUALS(string("Se"), m_kbise->siteAtomType(5));
            TS_ASSERT_EQUALS(string("Se"), m_kbise->siteAtomType(11));
        }


        void test_getLattice()
        {
            ObjCrystStructureAdapter* pkbise =
                dynamic_cast<ObjCrystStructureAdapter*>(m_kbise.get());
            TS_ASSERT(pkbise);
            const Lattice& L = pkbise->getLattice();
            const double eps = 1.0e-12;
            TS_ASSERT_DELTA(13.768, L.a(), eps);
            TS_ASSERT_DELTA(12.096, L.b(), eps);
            TS_ASSERT_DELTA(4.1656, L.c(), eps);
            TS_ASSERT_DELTA(89.98, L.alpha(), eps);
            TS_ASSERT_DELTA(98.64, L.beta(), eps);
            TS_ASSERT_DELTA(87.96, L.gamma(), eps);
        }


        void test_serialization()
        {
            stringstream storage(ios::in | ios::out | ios::binary);
            diffpy::serialization::oarchive oa(storage, ios::binary);
            oa << m_kbise;
            diffpy::serialization::iarchive ia(storage, ios::binary);
            StructureAdapterPtr kbise1;
            TS_ASSERT(!kbise1.get());
            ia >> kbise1;
            TS_ASSERT_DIFFERS(m_kbise.get(), kbise1.get());
            const double eps = 1.0e-7;
            TS_ASSERT_EQUALS(12, kbise1->countSites());
            TS_ASSERT_EQUALS(23.0, kbise1->totalOccupancy());
            TS_ASSERT_DELTA(0.0335565, kbise1->numberDensity(), eps);
            TS_ASSERT_EQUALS(string("K1+"), kbise1->siteAtomType(0));
            TS_ASSERT_EQUALS(string("Bi3+"), kbise1->siteAtomType(1));
            TS_ASSERT_EQUALS(string("Se"), kbise1->siteAtomType(5));
            TS_ASSERT_EQUALS(string("Se"), kbise1->siteAtomType(11));
            ObjCrystStructureAdapter* pkbise =
                dynamic_cast<ObjCrystStructureAdapter*>(m_kbise.get());
            ObjCrystStructureAdapter* pkbise1 =
                dynamic_cast<ObjCrystStructureAdapter*>(kbise1.get());
            const Lattice& L = pkbise->getLattice();
            const Lattice& L1 = pkbise1->getLattice();
            TS_ASSERT_EQUALS(L.a(), L1.a());
            TS_ASSERT_EQUALS(L.b(), L1.b());
            TS_ASSERT_EQUALS(L.c(), L1.c());
            TS_ASSERT_EQUALS(L.alpha(), L1.alpha());
            TS_ASSERT_EQUALS(L.beta(), L1.beta());
            TS_ASSERT_EQUALS(L.gamma(), L1.gamma());
        }



};  // class TestObjCrystStructureAdapter

//////////////////////////////////////////////////////////////////////////////
// class TestObjCrystStructureBondGenerator
//////////////////////////////////////////////////////////////////////////////

class TestObjCrystStructureBondGenerator : public CxxTest::TestSuite
{
    private:

        auto_ptr<Crystal> mcryst_ni;
        StructureAdapterPtr m_ni;
        BaseBondGeneratorPtr m_nibnds;

    public:

        void setUp()
        {
            if (!m_ni.get())
            {
                mcryst_ni.reset(loadTestCrystal("Ni.cif"));
                m_ni = createStructureAdapter(*mcryst_ni);
            }
            m_nibnds = m_ni->createBondGenerator();
        }


        void test_typeid()
        {
            TS_ASSERT(typeid(ObjCrystBondGenerator) ==
                    typeid(*m_nibnds));
        }


        void test_bondCountNickel()
        {
            m_nibnds->selectAnchorSite(0);
            m_nibnds->setRmin(0);
            m_nibnds->setRmax(1.0);
            TS_ASSERT_EQUALS(0, countBonds(*m_nibnds));
            m_nibnds->setRmin(-10);
            TS_ASSERT_EQUALS(0, countBonds(*m_nibnds));
            // there are 12 nearest neighbors at 2.49
            m_nibnds->setRmax(3);
            TS_ASSERT_EQUALS(12, countBonds(*m_nibnds));
            // there are no bonds between 2.6 and 3.4
            m_nibnds->setRmin(2.6);
            m_nibnds->setRmax(3.4);
            TS_ASSERT_EQUALS(0, countBonds(*m_nibnds));
            // there are 6 second nearest neighbors at 3.52
            m_nibnds->setRmax(3.6);
            TS_ASSERT_EQUALS(6, countBonds(*m_nibnds));
            // which sums to 18 neigbhors within radius 3.6
            m_nibnds->setRmin(0);
            TS_ASSERT_EQUALS(18, countBonds(*m_nibnds));
        }


        void test_bondCountWurtzite()
        {
            auto_ptr<Crystal> wurtzite(loadTestCrystal("ZnS_wurtzite.cif"));
            StructureAdapterPtr stru = createStructureAdapter(*wurtzite);
            BaseBondGeneratorPtr bnds = stru->createBondGenerator();
            TS_ASSERT_EQUALS(2, stru->countSites());
            bnds->selectAnchorSite(0);
            bnds->selectSiteRange(0, 2);
            bnds->setRmin(0);
            // there should be no bond below the ZnS distance 2.31
            bnds->setRmax(2.2);
            TS_ASSERT_EQUALS(0, countBonds(*bnds));
            // z-neighbor is slightly more distant than 3 in the lower plane
            bnds->setRmax(2.35);
            TS_ASSERT_EQUALS(3, countBonds(*bnds));
            bnds->setRmax(2.5);
            TS_ASSERT_EQUALS(4, countBonds(*bnds));
            // there are 12 second nearest neighbors at 3.81
            bnds->setRmin(3.7);
            bnds->setRmax(3.82);
            TS_ASSERT_EQUALS(12, countBonds(*bnds));
            // and one more at 3.83
            bnds->setRmax(3.85);
            TS_ASSERT_EQUALS(13, countBonds(*bnds));
            // making the total 17
            bnds->setRmin(0);
            TS_ASSERT_EQUALS(17, countBonds(*bnds));
            // and the same happens for all other sites
            bnds->selectAnchorSite(1);
            TS_ASSERT_EQUALS(17, countBonds(*bnds));
        }


        void test_LiTaO3()
        {
            const string lithium = "Li1+";
            const string tantalum = "Ta5+";
            const string oxygen = "O2-";
            const double epsu = 1e-5;
            auto_ptr<Crystal> litao3(loadTestCrystal("LiTaO3.cif"));
            StructureAdapterPtr stru = createStructureAdapter(*litao3);
            TS_ASSERT_EQUALS(3, stru->countSites());
            BaseBondGeneratorPtr bnds = stru->createBondGenerator();
            bnds->selectAnchorSite(0);
            bnds->selectSiteRange(0, 3);
            // there are 3 oxygen neighbors at 2.065
            bnds->setRmax(2.1);
            TS_ASSERT_EQUALS(3, countBonds(*bnds));
            // Li at site 0 is isotropic, oxygens have equal msd-s towards Li
            for (bnds->rewind(); !bnds->finished(); bnds->next())
            {
                TS_ASSERT_EQUALS(lithium, stru->siteAtomType(bnds->site0()));
                TS_ASSERT_EQUALS(oxygen, stru->siteAtomType(bnds->site1()));
                TS_ASSERT_DELTA(0.00265968, testmsd0(stru, bnds), epsu);
                TS_ASSERT_DELTA(0.00710945, testmsd1(stru, bnds), epsu);
            }
            // there are 3 oxygen neighbors at 2.26
            bnds->setRmin(2.2);
            bnds->setRmax(2.3);
            TS_ASSERT_EQUALS(3, countBonds(*bnds));
            for (bnds->rewind(); !bnds->finished(); bnds->next())
            {
                TS_ASSERT_EQUALS(oxygen, stru->siteAtomType(bnds->site1()));
                TS_ASSERT_DELTA(0.00265968, testmsd0(stru, bnds), epsu);
                TS_ASSERT_DELTA(0.00824319, testmsd1(stru, bnds), epsu);
            }
            // finally there are 4 Ta neighbors between 2.8 and 3.1
            bnds->setRmin(2.8);
            bnds->setRmax(3.1);
            TS_ASSERT_EQUALS(4, countBonds(*bnds));
            for (bnds->rewind(); !bnds->finished(); bnds->next())
            {
                TS_ASSERT_DELTA(0.00265968, testmsd0(stru, bnds), epsu);
                TS_ASSERT_EQUALS(tantalum, stru->siteAtomType(bnds->site1()));
                R3::Vector r01xy = bnds->r01();
                r01xy[2] = 0.0;
                // for the tantalum above Li the msd equals U33
                if (R3::norm(r01xy) < 0.1)
                {
                    TS_ASSERT_DELTA(0.00356, testmsd1(stru, bnds), epsu);
                }
                // other 3 tantalums are related by tripple axis and
                // have the same msd towards the central Li
                else
                {
                    TS_ASSERT_DELTA(0.00486942, testmsd1(stru, bnds), epsu);
                }
            }
        }



};  // class TestObjCrystStructureBondGenerator


//////////////////////////////////////////////////////////////////////////////
// class TestObjCrystMoleculeAdapter
//////////////////////////////////////////////////////////////////////////////


class TestObjCrystMoleculeAdapter : public CxxTest::TestSuite
{
    private:

        auto_ptr<Molecule> mmol_c60;
        StructureAdapterPtr m_c60;

    public:

        void setUp()
        {
            if (!mmol_c60.get())
            {
                mmol_c60.reset(makeC60Molecule());
                m_c60 = createStructureAdapter(*mmol_c60);
            }
        }

        void test_typeid()
        {
            TS_ASSERT(typeid(ObjCrystMoleculeAdapter) == typeid(*m_c60));
        }


        void test_countSites()
        {
            TS_ASSERT_EQUALS(60, m_c60->countSites());
        }


        void test_totalOccupancy()
        {
            TS_ASSERT_EQUALS(60, m_c60->totalOccupancy());
        }


        void test_numberDensity()
        {
            const double eps = 1.0e-7;
            TS_ASSERT_DELTA(0.0, m_c60->numberDensity(), eps);
        }


        void test_siteCartesianPosition()
        {
            const double eps = 1.0e-5;
            R3::Vector rC0 = m_c60->siteCartesianPosition(0);
            TS_ASSERT_DELTA(3.45127, rC0[0], eps);
            TS_ASSERT_DELTA(0.68500, rC0[1], eps);
            TS_ASSERT_DELTA(0.00000, rC0[2], eps);
        }


        void test_siteAnisotropy()
        {
            for (int i = 0; i < m_c60->countSites(); ++i)
            {
                TS_ASSERT_EQUALS(false, m_c60->siteAnisotropy(i));
            }
        }


        void test_siteCartesianUij()
        {
            double BtoU = 1.0 / (8 * M_PI * M_PI);
            const R3::Matrix Uij = m_c60->siteCartesianUij(0);
            TS_ASSERT_EQUALS(0.0, Uij(0,1));
            TS_ASSERT_EQUALS(0.0, Uij(0,2));
            TS_ASSERT_EQUALS(0.0, Uij(1,0));
            TS_ASSERT_EQUALS(0.0, Uij(1,2));
            TS_ASSERT_EQUALS(0.0, Uij(2,0));
            TS_ASSERT_EQUALS(0.0, Uij(2,1));
            TS_ASSERT_EQUALS(BtoU, Uij(0,0));
            TS_ASSERT_EQUALS(BtoU, Uij(1,1));
            TS_ASSERT_EQUALS(BtoU, Uij(2,2));
        }


        void test_siteAtomType()
        {
            for (int i = 0; i < m_c60->countSites(); ++i)
            {
                TS_ASSERT_EQUALS(string("C"), m_c60->siteAtomType(i));
            }
        }


        void test_getLattice()
        {
            ObjCrystMoleculeAdapter* padpt =
                dynamic_cast<ObjCrystMoleculeAdapter*>(m_c60.get());
            TS_ASSERT(padpt);
            const Lattice& L = padpt->getLattice();
            const double eps = 1.0e-12;
            TS_ASSERT_DELTA(1, L.a(), eps);
            TS_ASSERT_DELTA(1, L.b(), eps);
            TS_ASSERT_DELTA(1, L.c(), eps);
            TS_ASSERT_DELTA(90, L.alpha(), eps);
            TS_ASSERT_DELTA(90, L.beta(), eps);
            TS_ASSERT_DELTA(90, L.gamma(), eps);
        }


        void test_serialization()
        {
            stringstream storage(ios::in | ios::out | ios::binary);
            diffpy::serialization::oarchive oa(storage, ios::binary);
            oa << m_c60;
            diffpy::serialization::iarchive ia(storage, ios::binary);
            StructureAdapterPtr c60a;
            TS_ASSERT(!c60a.get());
            ia >> c60a;
            TS_ASSERT_DIFFERS(m_c60.get(), c60a.get());
            TS_ASSERT_EQUALS(60, c60a->countSites());
            TS_ASSERT_EQUALS(60.0, c60a->totalOccupancy());
            TS_ASSERT_EQUALS(0.0, c60a->numberDensity());
            R3::Vector dxyz(1, 2, 3);
            dxyz = m_c60->siteCartesianPosition(0) - c60a->siteCartesianPosition(0);
            TS_ASSERT_EQUALS(0.0, R3::norm(dxyz));
            dxyz = m_c60->siteCartesianPosition(37) - c60a->siteCartesianPosition(37);
            TS_ASSERT_EQUALS(0.0, R3::norm(dxyz));
            for (int i = 0; i < c60a->countSites(); ++i)
            {
                TS_ASSERT_EQUALS(string("C"), c60a->siteAtomType(i));
            }
        }

};  // class TestObjCrystMoleculeAdapter


//////////////////////////////////////////////////////////////////////////////
// class TestObjCrystMoleculeBondGenerator
//////////////////////////////////////////////////////////////////////////////

class TestObjCrystMoleculeBondGenerator : public CxxTest::TestSuite
{
    private:

        auto_ptr<Molecule> mmol_c60;
        StructureAdapterPtr m_c60;
        BaseBondGeneratorPtr m_c60bnds;

    public:

        void setUp()
        {
            if (!m_c60.get())
            {
                mmol_c60.reset(makeC60Molecule());
                m_c60 = createStructureAdapter(*mmol_c60);
            }
            m_c60bnds = m_c60->createBondGenerator();
        }

        void test_typeid()
        {
            TS_ASSERT(typeid(ObjCrystMoleculeBondGenerator) ==
                    typeid(*m_c60bnds));
        }


        void test_bondCountC60()
        {
            m_c60bnds->selectAnchorSite(0);
            m_c60bnds->setRmin(0);
            m_c60bnds->setRmax(1.0);
            TS_ASSERT_EQUALS(0, countBonds(*m_c60bnds));
            m_c60bnds->setRmin(-10);
            TS_ASSERT_EQUALS(0, countBonds(*m_c60bnds));
            // there are 3 nearest neighbors at 1.44
            m_c60bnds->setRmax(1.5);
            TS_ASSERT_EQUALS(3, countBonds(*m_c60bnds));
            // There are 59 neighbors total
            m_c60bnds->setRmax(8);
            TS_ASSERT_EQUALS(59, countBonds(*m_c60bnds));
        }


};  // class TestObjCrystMoleculeBondGenerator


namespace {

Molecule* makeC60Molecule()
{

    Crystal* cryst = new Crystal(1, 1, 1, "P1");
    Molecule* mol = new Molecule(*cryst, "c60");
    cryst->AddScatterer(mol);

    // Populate the molecule
    ScatteringPowerAtom* sp = new ScatteringPowerAtom("C", "C");
    cryst->AddScatteringPower(sp);
    mol->AddAtom(3.451266498, 0.685000000, 0.000000000, sp, "C0");
    mol->AddAtom(3.451266498, -0.685000000, 0.000000000, sp, "C1");
    mol->AddAtom(-3.451266498, 0.685000000, 0.000000000, sp, "C2");
    mol->AddAtom(-3.451266498, -0.685000000, 0.000000000, sp, "C3");
    mol->AddAtom(0.685000000, 0.000000000, 3.451266498, sp, "C4");
    mol->AddAtom(-0.685000000, 0.000000000, 3.451266498, sp, "C5");
    mol->AddAtom(0.685000000, 0.000000000, -3.451266498, sp, "C6");
    mol->AddAtom(-0.685000000, 0.000000000, -3.451266498, sp, "C7");
    mol->AddAtom(0.000000000, 3.451266498, 0.685000000, sp, "C8");
    mol->AddAtom(0.000000000, 3.451266498, -0.685000000, sp, "C9");
    mol->AddAtom(0.000000000, -3.451266498, 0.685000000, sp, "C10");
    mol->AddAtom(0.000000000, -3.451266498, -0.685000000, sp, "C11");
    mol->AddAtom(3.003809890, 1.409000000, 1.171456608, sp, "C12");
    mol->AddAtom(3.003809890, 1.409000000, -1.171456608, sp, "C13");
    mol->AddAtom(3.003809890, -1.409000000, 1.171456608, sp, "C14");
    mol->AddAtom(3.003809890, -1.409000000, -1.171456608, sp, "C15");
    mol->AddAtom(-3.003809890, 1.409000000, 1.171456608, sp, "C16");
    mol->AddAtom(-3.003809890, 1.409000000, -1.171456608, sp, "C17");
    mol->AddAtom(-3.003809890, -1.409000000, 1.171456608, sp, "C18");
    mol->AddAtom(-3.003809890, -1.409000000, -1.171456608, sp, "C19");
    mol->AddAtom(1.409000000, 1.171456608, 3.003809890, sp, "C20");
    mol->AddAtom(1.409000000, -1.171456608, 3.003809890, sp, "C21");
    mol->AddAtom(-1.409000000, 1.171456608, 3.003809890, sp, "C22");
    mol->AddAtom(-1.409000000, -1.171456608, 3.003809890, sp, "C23");
    mol->AddAtom(1.409000000, 1.171456608, -3.003809890, sp, "C24");
    mol->AddAtom(1.409000000, -1.171456608, -3.003809890, sp, "C25");
    mol->AddAtom(-1.409000000, 1.171456608, -3.003809890, sp, "C26");
    mol->AddAtom(-1.409000000, -1.171456608, -3.003809890, sp, "C27");
    mol->AddAtom(1.171456608, 3.003809890, 1.409000000, sp, "C28");
    mol->AddAtom(-1.171456608, 3.003809890, 1.409000000, sp, "C29");
    mol->AddAtom(1.171456608, 3.003809890, -1.409000000, sp, "C30");
    mol->AddAtom(-1.171456608, 3.003809890, -1.409000000, sp, "C31");
    mol->AddAtom(1.171456608, -3.003809890, 1.409000000, sp, "C32");
    mol->AddAtom(-1.171456608, -3.003809890, 1.409000000, sp, "C33");
    mol->AddAtom(1.171456608, -3.003809890, -1.409000000, sp, "C34");
    mol->AddAtom(-1.171456608, -3.003809890, -1.409000000, sp, "C35");
    mol->AddAtom(2.580456608, 0.724000000, 2.279809890, sp, "C36");
    mol->AddAtom(2.580456608, 0.724000000, -2.279809890, sp, "C37");
    mol->AddAtom(2.580456608, -0.724000000, 2.279809890, sp, "C38");
    mol->AddAtom(2.580456608, -0.724000000, -2.279809890, sp, "C39");
    mol->AddAtom(-2.580456608, 0.724000000, 2.279809890, sp, "C40");
    mol->AddAtom(-2.580456608, 0.724000000, -2.279809890, sp, "C41");
    mol->AddAtom(-2.580456608, -0.724000000, 2.279809890, sp, "C42");
    mol->AddAtom(-2.580456608, -0.724000000, -2.279809890, sp, "C43");
    mol->AddAtom(0.724000000, 2.279809890, 2.580456608, sp, "C44");
    mol->AddAtom(0.724000000, -2.279809890, 2.580456608, sp, "C45");
    mol->AddAtom(-0.724000000, 2.279809890, 2.580456608, sp, "C46");
    mol->AddAtom(-0.724000000, -2.279809890, 2.580456608, sp, "C47");
    mol->AddAtom(0.724000000, 2.279809890, -2.580456608, sp, "C48");
    mol->AddAtom(0.724000000, -2.279809890, -2.580456608, sp, "C49");
    mol->AddAtom(-0.724000000, 2.279809890, -2.580456608, sp, "C50");
    mol->AddAtom(-0.724000000, -2.279809890, -2.580456608, sp, "C51");
    mol->AddAtom(2.279809890, 2.580456608, 0.724000000, sp, "C52");
    mol->AddAtom(-2.279809890, 2.580456608, 0.724000000, sp, "C53");
    mol->AddAtom(2.279809890, 2.580456608, -0.724000000, sp, "C54");
    mol->AddAtom(-2.279809890, 2.580456608, -0.724000000, sp, "C55");
    mol->AddAtom(2.279809890, -2.580456608, 0.724000000, sp, "C56");
    mol->AddAtom(-2.279809890, -2.580456608, 0.724000000, sp, "C57");
    mol->AddAtom(2.279809890, -2.580456608, -0.724000000, sp, "C58");
    mol->AddAtom(-2.279809890, -2.580456608, -0.724000000, sp, "C59");

    return mol;
}

}   // namespace


// End of file
