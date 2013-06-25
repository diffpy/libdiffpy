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
* class TestDiffPyStructureAdapter -- unit tests for an adapter
*     to Structure class from diffpy.Structure
*
* class TestDiffPyStructureBondGenerator -- unit tests for bond generator
*
*****************************************************************************/

#include <typeinfo>
#include <sstream>
#include <cxxtest/TestSuite.h>

#include <diffpy/serialization.hpp>
#include <diffpy/srreal/PythonStructureAdapter.hpp>
#include <diffpy/srreal/DiffPyStructureAdapter.hpp>
#include <diffpy/srreal/PointsInSphere.hpp>
#include "python_helpers.hpp"

using namespace std;
using namespace boost;
using namespace diffpy::srreal;

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

}   // namespace

//////////////////////////////////////////////////////////////////////////////
// class TestDiffPyStructureAdapter
//////////////////////////////////////////////////////////////////////////////

class TestDiffPyStructureAdapter : public CxxTest::TestSuite
{
    private:

        StructureAdapterPtr m_ni;
        StructureAdapterPtr m_kbise;
        StructureAdapterPtr m_catio3;
        StructureAdapterPtr m_pswt;

    public:

        void setUp()
        {
            CxxTest::setAbortTestOnFail(true);
            python::object stru;
            if (!m_ni.get())
            {
                stru = loadTestStructure("Ni.cif");
                m_ni = createStructureAdapter(stru);
            }
            if (!m_kbise.get())
            {
                stru = loadTestStructure("alpha_K2Bi8Se13.cif");
                m_kbise = createStructureAdapter(stru);
            }
            if (!m_catio3.get())
            {
                stru = loadTestStructure("icsd_62149.cif");
                m_catio3 = createStructureAdapter(stru);
            }
            if (!m_pswt.get())
            {
                stru = loadTestStructure("PbScW25TiO3.stru");
                m_pswt = createStructureAdapter(stru);
            }
            CxxTest::setAbortTestOnFail(false);
        }


        void test_typeid()
        {
            TS_ASSERT(typeid(DiffPyStructureAdapter) == typeid(*m_ni));
        }


        void test_countSites()
        {
            TS_ASSERT_EQUALS(4, m_ni->countSites());
            TS_ASSERT_EQUALS(23, m_kbise->countSites());
            TS_ASSERT_EQUALS(20, m_catio3->countSites());
            TS_ASSERT_EQUALS(56, m_pswt->countSites());
        }


        void test_totalOccupancy()
        {
            TS_ASSERT_EQUALS(4.0, m_ni->totalOccupancy());
            TS_ASSERT_EQUALS(23.0, m_kbise->totalOccupancy());
            TS_ASSERT_EQUALS(20.0, m_catio3->totalOccupancy());
            TS_ASSERT_EQUALS(40.0, m_pswt->totalOccupancy());
        }


        void test_numberDensity()
        {
            const double eps = 1.0e-7;
            TS_ASSERT_DELTA(0.0914114, m_ni->numberDensity(), eps);
            TS_ASSERT_DELTA(0.0335565, m_kbise->numberDensity(), eps);
            TS_ASSERT_DELTA(0.0894126, m_catio3->numberDensity(), eps);
            TS_ASSERT_DELTA(0.0760332, m_pswt->numberDensity(), eps);
        }


        void test_siteCartesianPosition()
        {
            const double eps = 1.0e-5;
            R3::Vector rCa = m_catio3->siteCartesianPosition(0);
            TS_ASSERT_DELTA(2.72617, rCa[0], eps);
            TS_ASSERT_DELTA(2.91718, rCa[1], eps);
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
            const double* puij = m_ni->siteCartesianUij(0).data();
            for (int i = 0; i < 9; ++i)
            {
                TS_ASSERT_EQUALS(0.0, puij[i]);
            }
            // check CaTiO3 values
            const R3::Matrix& UTi = m_catio3->siteCartesianUij(7);
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
            TS_ASSERT_EQUALS(string("Ni"), m_ni->siteAtomType(3));
            TS_ASSERT_EQUALS(string("K1+"), m_kbise->siteAtomType(0));
            TS_ASSERT_EQUALS(string("Bi3+"), m_kbise->siteAtomType(2));
            TS_ASSERT_EQUALS(string("Se"), m_kbise->siteAtomType(10));
            TS_ASSERT_EQUALS(string("Se"), m_kbise->siteAtomType(22));
        }


        void test_getLattice()
        {
            DiffPyStructureAdapter* pkbise =
                dynamic_cast<DiffPyStructureAdapter*>(m_kbise.get());
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
            TS_ASSERT_EQUALS(23, kbise1->countSites());
            TS_ASSERT_EQUALS(23.0, kbise1->totalOccupancy());
            TS_ASSERT_DELTA(0.0335565, kbise1->numberDensity(), eps);
            TS_ASSERT_EQUALS(string("K1+"), kbise1->siteAtomType(0));
            TS_ASSERT_EQUALS(string("Bi3+"), kbise1->siteAtomType(2));
            TS_ASSERT_EQUALS(string("Se"), kbise1->siteAtomType(10));
            TS_ASSERT_EQUALS(string("Se"), kbise1->siteAtomType(22));
            DiffPyStructureAdapter* pkbise =
                dynamic_cast<DiffPyStructureAdapter*>(m_kbise.get());
            DiffPyStructureAdapter* pkbise1 =
                dynamic_cast<DiffPyStructureAdapter*>(kbise1.get());
            const Lattice& L = pkbise->getLattice();
            const Lattice& L1 = pkbise1->getLattice();
            TS_ASSERT_EQUALS(L.a(), L1.a());
            TS_ASSERT_EQUALS(L.b(), L1.b());
            TS_ASSERT_EQUALS(L.c(), L1.c());
            TS_ASSERT_EQUALS(L.alpha(), L1.alpha());
            TS_ASSERT_EQUALS(L.beta(), L1.beta());
            TS_ASSERT_EQUALS(L.gamma(), L1.gamma());
        }


};  // class TestDiffPyStructureAdapter

//////////////////////////////////////////////////////////////////////////////
// class TestDiffPyStructureBondGenerator
//////////////////////////////////////////////////////////////////////////////

class TestDiffPyStructureBondGenerator : public CxxTest::TestSuite
{
    private:

        StructureAdapterPtr m_ni;
        BaseBondGeneratorPtr m_nibnds;

    public:

        void setUp()
        {
            CxxTest::setAbortTestOnFail(true);
            if (!m_ni.get())
            {
                python::object stru;
                stru = loadTestStructure("Ni.cif");
                m_ni.reset(new DiffPyStructureAdapter(stru));
            }
            m_nibnds = m_ni->createBondGenerator();
            CxxTest::setAbortTestOnFail(false);
        }


        void test_typeid()
        {
            TS_ASSERT(typeid(DiffPyStructurePeriodicBondGenerator) ==
                    typeid(*m_nibnds));
        }


        void test_bondCountNickel()
        {
            m_nibnds->selectAnchorSite(0);
            m_nibnds->setRmin(0);
            m_nibnds->setRmax(1.0);
            TS_ASSERT_EQUALS(0, countBonds(*m_nibnds));
            m_nibnds->selectAnchorSite(3);
            TS_ASSERT_EQUALS(0, countBonds(*m_nibnds));
            m_nibnds->setRmin(-10);
            TS_ASSERT_EQUALS(0, countBonds(*m_nibnds));
            // there are 12 nearest neighbors at 2.49
            m_nibnds->setRmax(3);
            TS_ASSERT_EQUALS(12, countBonds(*m_nibnds));
            m_nibnds->selectAnchorSite(0);
            TS_ASSERT_EQUALS(12, countBonds(*m_nibnds));
            // there are no self neighbors below the cell length of 3.52
            m_nibnds->selectAnchorSite(0);
            m_nibnds->selectSiteRange(0, 1);
            TS_ASSERT_EQUALS(0, countBonds(*m_nibnds));
            // and any other unit cell atom would give 4 neighbors
            m_nibnds->selectAnchorSite(0);
            m_nibnds->selectSiteRange(3, 4);
            TS_ASSERT_EQUALS(4, countBonds(*m_nibnds));
            // there are no bonds between 2.6 and 3.4
            m_nibnds->setRmin(2.6);
            m_nibnds->setRmax(3.4);
            m_nibnds->selectSiteRange(0, 4);
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
            python::object wurtzite = loadTestStructure("ZnS_wurtzite.cif");
            StructureAdapterPtr stru = createStructureAdapter(wurtzite);
            BaseBondGeneratorPtr bnds = stru->createBondGenerator();
            TS_ASSERT_EQUALS(4, stru->countSites());
            bnds->selectAnchorSite(0);
            bnds->selectSiteRange(0, 4);
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
            bnds->selectAnchorSite(2);
            TS_ASSERT_EQUALS(17, countBonds(*bnds));
            bnds->selectAnchorSite(3);
            TS_ASSERT_EQUALS(17, countBonds(*bnds));
        }


        void test_LiTaO3()
        {
            const string lithium = "Li1+";
            const string tantalum = "Ta5+";
            const string oxygen = "O2-";
            const double epsu = 1e-5;
            python::object litao3 = loadTestStructure("LiTaO3.cif");
            StructureAdapterPtr stru = createStructureAdapter(litao3);
            TS_ASSERT_EQUALS(30, stru->countSites());
            BaseBondGeneratorPtr bnds = stru->createBondGenerator();
            bnds->selectAnchorSite(0);
            bnds->selectSiteRange(0, 30);
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

};  // class TestDiffPyStructureBondGenerator

// End of file
