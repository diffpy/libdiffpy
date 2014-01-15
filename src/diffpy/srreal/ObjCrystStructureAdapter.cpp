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
* Factory function createStructureAdapter(ObjCryst::Molecule)
*   -- builds AtomicStructureAdapter from the ObjCryst Molecule object
*
*****************************************************************************/

#include <algorithm>
#include <cassert>
#include <cmath>
#include <set>

#include <diffpy/srreal/ObjCrystStructureAdapter.hpp>

using namespace std;

namespace diffpy {
namespace srreal {

namespace {

// Utility functions ---------------------------------------------------------

R3::Matrix
getUij(const ObjCryst::ScatteringPower* sp)
{
    const double BtoU = 1.0 / (8 * M_PI * M_PI);
    R3::Matrix Uij = R3::zeromatrix();
    if (sp->IsIsotropic())
    {
        Uij(0,0) = Uij(1,1) = Uij(2,2) = sp->GetBiso() * BtoU;
    }
    else
    {
        Uij(0,0) = sp->GetBij(1,1) * BtoU;
        Uij(1,1) = sp->GetBij(2,2) * BtoU;
        Uij(2,2) = sp->GetBij(3,3) * BtoU;
        Uij(0,1) = Uij(1,0) = sp->GetBij(1,2) * BtoU;
        Uij(0,2) = Uij(2,0) = sp->GetBij(1,3) * BtoU;
        Uij(1,2) = Uij(2,1) = sp->GetBij(2,3) * BtoU;
    }
    return Uij;
}


CrystalStructureAdapter::SymOpVector
fetchSymmetryOperations(const ObjCryst::SpaceGroup& spacegroup)
{
    CrystalStructureAdapter::SymOpVector rv;
    const int nbsym = spacegroup.GetNbSymmetrics();
    rv.resize(nbsym);
    const vector<ObjCryst::SpaceGroup::SMx>& sgsymops =
        spacegroup.GetSymmetryOperations();
    // copy convert the ObjCryst symmetry operations
    int last = sgsymops.size();
    assert(last <= nbsym);
    for (int i = 0; i != last; ++i)
    {
        CrystalStructureAdapter::SymOpVector::value_type& op = rv[i];
        copy(sgsymops[i].mx, sgsymops[i].mx + 9, op.R.data().begin());
        copy(sgsymops[i].tr, sgsymops[i].tr + 3, op.t.data().begin());
    }
    // apply all origin translations
    const vector<ObjCryst::SpaceGroup::TRx>& sgtrans =
        spacegroup.GetTranslationVectors();
    const int nbtran = spacegroup.GetNbTranslationVectors();
    assert(int(sgtrans.size()) == nbtran);
    assert(nbtran * last <= nbsym);
    for (int nt = 0; nt < nbtran; ++nt)
    {
        const double* pt = sgtrans[nt].tr;
        R3::Vector sgt(pt[0], pt[1], pt[2]);
        for (int i = 0; i < last; ++i)
        {
            int k = i + nt * last;
            assert(k < nbsym);
            rv[k] = rv[i];
            rv[k].t += sgt;
        }
    }
    last = sgsymops.size() * nbtran;
    // apply the inversion center if necessary
    if (last < nbsym)
    {
        assert(spacegroup.IsCentrosymmetric());
        assert(nbsym == 2 * last);
        bool noCenter = false, noTransl = true, noIdentical = true;
        CrystMatrix_REAL xyzinv = spacegroup.GetAllSymmetrics(
                0.0, 0.0, 0.0, noCenter, noTransl, noIdentical);
        // check if symmetry center is away from the origin
        diffpy::mathutils::EpsilonEqual allclose;
        R3::Vector dxyz = R3::zerovector;
        for (int i = 0; i < xyzinv.rows(); ++i)
        {
            using mathutils::eps_eq;
            R3::Vector xyzi(
                    round(xyzinv(i, 0) / 0.5),
                    round(xyzinv(i, 1) / 0.5),
                    round(xyzinv(i, 2) / 0.5));
            if (eps_eq(0.0, R3::norm(xyzi)))  continue;
            // here we found inversion center at non-origin location
            dxyz = xyzi - R3::floor(xyzi);
            break;
        }
        for (int i = 0, j = last; i < last; ++i, ++j)
        {
            rv[j] = rv[i];
            rv[j].R *= -1;
            rv[j].t = dxyz - rv[j].t;
        }
        last *= 2;
    }
    assert(nbsym == last);
    return rv;
}

}   // namespace

//////////////////////////////////////////////////////////////////////////////
// StructureAdapter factory for ObjCryst::Crystal
//////////////////////////////////////////////////////////////////////////////

StructureAdapterPtr
createStructureAdapter(const ObjCryst::Crystal& cryst)
{
    const double radtodeg = 180 / M_PI;
    CrystalStructureAdapterPtr adpt(new CrystalStructureAdapter);
    adpt->setLatPar(
            cryst.GetLatticePar(0),
            cryst.GetLatticePar(1),
            cryst.GetLatticePar(2),
            radtodeg * cryst.GetLatticePar(3),
            radtodeg * cryst.GetLatticePar(4),
            radtodeg * cryst.GetLatticePar(5));
    // find out number of scatterers in the asymmetric unit
    const ObjCryst::ScatteringComponentList& scl =
        cryst.GetScatteringComponentList();
    size_t nbComponent = scl.GetNbComponent();
    adpt->reserve(nbComponent);
    Atom ai;
    const Lattice& L = adpt->getLattice();
    for (size_t i = 0; i < nbComponent; ++i)
    {
        const ObjCryst::ScatteringComponent& sc = scl(i);
        const ObjCryst::ScatteringPower* sp = sc.mpScattPow;
        // Skip over this if it is a dummy atom. A dummy atom has no
        // mpScattPow, and therefore no type. It's just in a structure as a
        // reference position.
        if (sp == NULL) continue;
        ai.occupancy = sc.mOccupancy;
        ai.anisotropy = !(sp->IsIsotropic());
        ai.atomtype = sp->GetSymbol();
        R3::Vector xyz(sc.mX, sc.mY, sc.mZ);
        ai.xyz_cartn = L.cartesian(xyz);
        // Store Uij
        R3::Matrix uijl = getUij(sp);
        ai.uij_cartn = ai.anisotropy ?  L.cartesianMatrix(uijl) : uijl;
        adpt->append(ai);
    }
    const ObjCryst::SpaceGroup& spacegroup = cryst.GetSpaceGroup();
    CrystalStructureAdapter::SymOpVector symops =
        fetchSymmetryOperations(spacegroup);
    CrystalStructureAdapter::SymOpVector::const_iterator op;
    for (op = symops.begin(); op != symops.end(); ++op)  adpt->addSymOp(*op);
    adpt->updateSymmetryPositions();
    return adpt;
}

//////////////////////////////////////////////////////////////////////////////
// StructureAdapter factory for ObjCryst::Molecule
//////////////////////////////////////////////////////////////////////////////

StructureAdapterPtr
createStructureAdapter(const ObjCryst::Molecule& molecule)
{
    AtomicStructureAdapterPtr adpt(new AtomicStructureAdapter);
    size_t nbComponent = molecule.GetNbComponent();
    adpt->reserve(nbComponent);
    Atom ai;
    for (size_t i = 0; i < nbComponent; ++i)
    {
        const ObjCryst::MolAtom& atom = molecule.GetAtom(i);
        if (atom.IsDummy()) continue;
        ai.occupancy = atom.GetOccupancy();
        const ObjCryst::ScatteringPower* sp = &(atom.GetScatteringPower());
        ai.anisotropy = !(sp->IsIsotropic());
        ai.atomtype = sp->GetSymbol();
        ai.xyz_cartn[0] = atom.X();
        ai.xyz_cartn[1] = atom.Y();
        ai.xyz_cartn[2] = atom.Z();
        // Store Uij
        ai.uij_cartn = getUij(sp);
        adpt->append(ai);
    }
    return adpt;
}

}   // namespace srreal
}   // namespace diffpy

// End of file
