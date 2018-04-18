/*****************************************************************************
*
* libdiffpy         by DANSE Diffraction group
*                   Simon J. L. Billinge
*                   (c) 2009 The Trustees of Columbia University
*                   in the City of New York.  All rights reserved.
*
* File coded by:    Chris Farrow
*
* See AUTHORS.txt for a list of people who contributed.
* See LICENSE_DANSE.txt for license information.
*
******************************************************************************
*
* StructureAdapterPtr createStructureAdapter(const ObjCryst::Crystal&)
* StructureAdapterPtr createStructureAdapter(const ObjCryst::Molecule&)
*   -- structure adapter factories for ObjCryst objects.
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


bool isGetInversionCenterDoubled()
{
    static int value_cached = 0;
    const int BIT_CACHED = 1;
    const int BIT_VALUE = 2;
    // short circuit when the flag is already cached
    if (BIT_CACHED & value_cached) {
        const bool rv = BIT_VALUE & value_cached;
        return rv;
    }
    ObjCryst::SpaceGroup sg129("P 4/n m m :1");
    const CrystVector_REAL xyzinv = sg129.GetInversionCenter();
    if (0.5 == xyzinv(0) && 0.5 == xyzinv(1) && 0 == xyzinv(2)) {
        value_cached |= BIT_VALUE;
    }
    else if (0.25 == xyzinv(0) && 0.25 == xyzinv(1) && 0 == xyzinv(2)) {
        value_cached &= ~BIT_VALUE;
    }
    else {
        const char* emsg =
            "Unexpected value of ObjCryst::SpaceGroup::GetInversionCenter()";
        throw logic_error(emsg);
    }
    value_cached |= BIT_CACHED;
    return isGetInversionCenterDoubled();
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
        CrystVector_REAL xyzinv = spacegroup.GetInversionCenter();
        // get tinv the overall translation associated with inversion center.
        // fox-objcryst 2017.2 returns doubled coordinates for the inversion
        // center so the tinv is the same.  Expect fixup which will require
        // doubling the tinv value here.
        R3::Vector tinv(xyzinv(0), xyzinv(1), xyzinv(2));
        tinv *= isGetInversionCenterDoubled() ? 1.0 : 2.0;
        // now apply the inversion center here
        for (int i = 0, j = last; i < last; ++i, ++j)
        {
            rv[j].R = -1 * rv[i].R;
            rv[j].t = tinv - rv[i].t;
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
