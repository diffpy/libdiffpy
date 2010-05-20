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
* class ObjCrystAperiodicBondGenerator 
*   -- Generate bonds from aperiodic ObjCrystStructureAdapter.
* class ObjCrystPeriodicBondGenerator     
*   -- Generate bonds from periodic ObjCrystStructureAdapter.
*
* $Id$
*
*****************************************************************************/

#include <algorithm>
#include <cassert>
#include <cmath>
#include <set>
#include <string>
#include <vector>

#include <diffpy/PythonInterface.hpp>
#include <diffpy/srreal/PythonStructureAdapter.hpp>
#include <diffpy/srreal/Lattice.hpp>
#include <diffpy/srreal/ObjCrystStructureAdapter.hpp>
#include <diffpy/srreal/PDFUtils.hpp>
#include <diffpy/srreal/PointsInSphere.hpp>

using namespace std;
using namespace diffpy::srreal;

namespace {

const double rtod = 180 / M_PI;
const double UtoB = 8 * M_PI * M_PI;
const double BtoU = 1.0 / UtoB;

}


//////////////////////////////////////////////////////////////////////////////
// class ObjCrystStructureAdapter
//////////////////////////////////////////////////////////////////////////////

// static data

const double ObjCrystStructureAdapter::mtoler = 1e-5;

// Constructor ---------------------------------------------------------------

ObjCrystStructureAdapter::
ObjCrystStructureAdapter(const ObjCryst::Crystal& cryst) : mpcryst(&cryst)
{
    using ObjCryst::Crystal;
    mlattice.setLatPar( mpcryst->GetLatticePar(0), 
                        mpcryst->GetLatticePar(1),
                        mpcryst->GetLatticePar(2), 
                        rtod * mpcryst->GetLatticePar(3),
                        rtod * mpcryst->GetLatticePar(4), 
                        rtod * mpcryst->GetLatticePar(5) );

    // The dynamic population correction in ObjCryst is used by its interal
    // calculators, but slows down the calculation of atom positions. This is
    // also required to get the proper ScatteringComponentList for aperiodic
    // structures. Since we're not using any ObjCryst calculators, we turn it
    // off momentarily.
    int usepopcorr = mpcryst->GetUseDynPopCorr();
    const_cast<Crystal*>(mpcryst)->SetUseDynPopCorr(0);
    if (this->isPeriodic())
    {
        this->getPeriodicUnitCell();
    }
    else
    {
        this->getAperiodicUnitCell();
    }
    // Undo change
    const_cast<Crystal*>(mpcryst)->SetUseDynPopCorr(usepopcorr);
}


// Public Methods ------------------------------------------------------------

BaseBondGenerator* 
ObjCrystStructureAdapter::
createBondGenerator() const
{

    BaseBondGenerator* bnds = this->isPeriodic() ?
        new ObjCrystPeriodicBondGenerator(this) :
        new ObjCrystAperiodicBondGenerator(this);
    return bnds;
}


int 
ObjCrystStructureAdapter::
countSites() const
{
    return mvsc.size();
}


double 
ObjCrystStructureAdapter::
numberDensity() const
{
    double rv = this->isPeriodic() ?
        (this->totalOccupancy() / mlattice.volume()) : 0.0;
    return rv;
}


const Lattice& 
ObjCrystStructureAdapter::
getLattice() const
{
    return mlattice;
}


const R3::Vector& 
ObjCrystStructureAdapter::
siteCartesianPosition(int idx) const
{
    assert(0 <= idx && idx < this->countSites());
    return mvsym[idx][0];
}


double 
ObjCrystStructureAdapter::
siteOccupancy(int idx) const
{
    assert(0 <= idx && idx < this->countSites());
    return mvsc[idx].mOccupancy;
}


bool 
ObjCrystStructureAdapter::
siteAnisotropy(int idx) const
{
    assert(0 <= idx && idx < this->countSites());
    return !(mvsc[idx].mpScattPow->IsIsotropic());
}


double
ObjCrystStructureAdapter::
siteMultiplicity(int idx) const
{

    assert(0 <= idx && idx < this->countSites());
    return mvsym[idx].size();
}


const R3::Matrix&
ObjCrystStructureAdapter::
siteCartesianUij(int idx) const
{
    assert(0 <= idx && idx < this->countSites());
    return mvuij[idx][0];
}


const string& 
ObjCrystStructureAdapter::
siteAtomType(int idx) const
{
    assert(0 <= idx && idx < this->countSites());
    return mvsc[idx].mpScattPow->GetSymbol();
}


bool 
ObjCrystStructureAdapter::
isPeriodic() const
{
    const Lattice& L = this->getLattice();
    bool rv = !(R3::EpsEqual(R3::identity(), L.base()));
    return rv;
}


// Private Methods -----------------------------------------------------------

/* Get the conventional unit cell from the crystal. */
void
ObjCrystStructureAdapter::
getPeriodicUnitCell()
{

    // Expand each scattering component in the primitive cell and record the
    // new scatterers.
    const ObjCryst::ScatteringComponentList& scl =
        mpcryst->GetScatteringComponentList();
    size_t nbComponent = scl.GetNbComponent();

    // Various things we need to know about the symmetry operations
    const ObjCryst::SpaceGroup& spacegroup = mpcryst->GetSpaceGroup();
    size_t nbSymmetrics = spacegroup.GetNbSymmetrics();

    double x, y, z, junk;
    CrystMatrix<double> symmetricsCoords;

    mvsc.clear();
    mvsym.clear();
    mvuij.clear();
    mvsc.reserve(nbComponent);
    mvsym.reserve(nbComponent);
    mvuij.reserve(nbComponent);

    typedef std::set<R3::Vector, R3::EpsCompare> SymPosSet;

    // Get the symmetry operations
    const std::vector<ObjCryst::SpaceGroup::SMx>& symops = 
        spacegroup.GetSymmetryOperations();
    const size_t symsize = symops.size();

    // For each scattering component, find its position in the primitive cell
    // and expand that position.
    for (size_t i = 0; i < nbComponent; ++i)
    {

        const ObjCryst::ScatteringPower* sp = scl(i).mpScattPow;

        // Skip over this if it is a dummy atom. A dummy atom has no
        // mpScattPow, and therefore no type. It's just in a structure as a
        // reference position.
        if (sp == NULL) continue;

        mvsc.push_back(scl(i));
        SymPosSet symset = SymPosSet(R3::EpsCompare(mtoler));
        SymPosVec symvec = SymPosVec();
        SymUijVec uijvec = SymUijVec();
        size_t numsym = 0;

        // Store Uij.
        R3::Matrix Uij = objcrystutil::getUij(sp);

        // Get all the symmetric coordinates. Symmetric coordinates are created
        // by translation vector in the outer loop, and rotation matrix in the
        // inner loop. The first translation is [0,0,0]. We need to know this
        // to determine the rotation matrices.
        symmetricsCoords = spacegroup.GetAllSymmetrics(
            scl(i).mX, scl(i).mY, scl(i).mZ);

        // Collect the unique symmetry positions.
        for (size_t j = 0; j < nbSymmetrics; ++j)
        {
            x = modf(symmetricsCoords(j, 0), &junk);
            y = modf(symmetricsCoords(j, 1), &junk);
            z = modf(symmetricsCoords(j, 2), &junk);
            if (fabs(x) < mtoler) x = 0;
            if (fabs(y) < mtoler) y = 0;
            if (fabs(z) < mtoler) z = 0;
            if (x < 0) x += 1.;
            if (y < 0) y += 1.;
            if (z < 0) z += 1.;

            // Record the position
            R3::Vector xyz; 
            xyz = x, y, z;

            // We use this to filter unique positions
            symset.insert(xyz);

            // See if we've got a new position
            if(symset.size() > numsym)
            {
                ++numsym;

                // Get this in Cartesian
                R3::Vector xyzc = mlattice.cartesian(xyz);
                
                // Store this in the symvec so we are assured that the order
                // will not change.
                symvec.push_back(xyzc);

                // Get the symmetry operation used to generate this position.
                R3::Matrix M;
                size_t k = j % symsize;
                M = symops[k].mx[0], symops[k].mx[1], symops[k].mx[2],
                    symops[k].mx[3], symops[k].mx[4], symops[k].mx[5],
                    symops[k].mx[6], symops[k].mx[7], symops[k].mx[8];

                // rotate the Uij matrix
                R3::Matrix Utmp, Urot, UCart; 
                // Get UCart
                if (sp->IsIsotropic())
                {
                    UCart = Uij;
                }
                else
                {
                    Utmp = R3::product(Uij, R3::transpose(M));
                    Urot = R3::product(M, Utmp);
                    UCart = mlattice.cartesianMatrix(Urot);
                }
                uijvec.push_back(UCart);
            }
        }

        assert( uijvec.size() == symvec.size() );

        // Store symmetric positions and Uij matrices
        mvsym.push_back(symvec);
        mvuij.push_back(uijvec);

    }

}


void
ObjCrystStructureAdapter::
getAperiodicUnitCell()
{

    // Aperiodic implies P1 symmetry. We speed things up a bit for this case.

    const ObjCryst::ScatteringComponentList& scl =
        mpcryst->GetScatteringComponentList();
    size_t nbComponent = scl.GetNbComponent();

    // Clear old vectors and make room for the new
    mvsc.clear();
    mvsym.clear();
    mvuij.clear();
    mvsc.reserve(nbComponent);
    mvsym.reserve(nbComponent);
    mvuij.reserve(nbComponent);

    // For each scattering component, find and record its position.
    for (size_t i = 0; i < nbComponent; ++i)
    {

        const ObjCryst::ScatteringPower* sp = scl(i).mpScattPow;

        // Skip over this if it is a dummy atom. A dummy atom has no
        // mpScattPow, and therefore no type. It's just in a structure as a
        // reference position.
        if (sp == NULL) continue;

        // Make the requisite storage vectors
        SymPosVec symvec = SymPosVec();
        SymUijVec uijvec = SymUijVec();

        // Record the ScatteringPower
        mvsc.push_back(scl(i));

        // Record the position
        R3::Vector xyz; 
        xyz = scl(i).mX, scl(i).mY, scl(i).mZ;
        symvec.push_back(xyz);

        // Record UCart
        R3::Matrix Uij = objcrystutil::getUij(sp);
        R3::Matrix UCart;
        if (!sp->IsIsotropic())
        {
            UCart = mlattice.cartesianMatrix(Uij);
        }
        else
        {
            UCart = Uij;
        }
        uijvec.push_back(UCart);

        // Store teh requisite vectors
        mvsym.push_back(symvec);
        mvuij.push_back(uijvec);
    }
}


//////////////////////////////////////////////////////////////////////////////
// class ObjCrystMoleculeAdapter
//////////////////////////////////////////////////////////////////////////////

// static data

// Constructor ---------------------------------------------------------------

ObjCrystMoleculeAdapter::
ObjCrystMoleculeAdapter(const ObjCryst::Molecule& molecule) 
: mpmolecule(&molecule)
{
    using ObjCryst::Molecule;
    using ObjCryst::MolAtom;
    mlattice.setLatPar(1, 1, 1, 90, 90, 90);
    const ObjCryst::ScatteringComponentList& scl = 
        molecule.GetScatteringComponentList();

    size_t nbComponent = scl.GetNbComponent();
    // Initialize positions and uij
    mvatoms.clear();
    mvpos.clear();
    mvuij.clear();
    mvatoms.reserve(nbComponent);
    mvpos.reserve(nbComponent);
    mvuij.reserve(nbComponent);

    for (size_t i = 0; i < nbComponent; ++i)
    {

        const MolAtom& atom = mpmolecule->GetAtom(i);
        if (atom.IsDummy()) continue;

        // Store the atom
        mvatoms.push_back(atom);

        // Store position
        R3::Vector xyz; 
        xyz = atom.X(), atom.Y(), atom.Z();
        mvpos.push_back(xyz);

        // Store Uij
        R3::Matrix Uij = objcrystutil::getUij(&atom.GetScatteringPower());
        mvuij.push_back(Uij);
    }

}

// Public Methods ------------------------------------------------------------

BaseBondGenerator* 
ObjCrystMoleculeAdapter::
createBondGenerator() const
{

    BaseBondGenerator* bnds = new ObjCrystMoleculeBondGenerator(this);
    return bnds;
}


int 
ObjCrystMoleculeAdapter::
countSites() const
{
    return mvatoms.size();
}


double 
ObjCrystMoleculeAdapter::
numberDensity() const
{
    return 0.0;
}


const Lattice& 
ObjCrystMoleculeAdapter::
getLattice() const
{
    return mlattice;
}


const R3::Vector& 
ObjCrystMoleculeAdapter::
siteCartesianPosition(int idx) const
{
    assert(0 <= idx && idx < this->countSites());
    return mvpos[idx];
}


double 
ObjCrystMoleculeAdapter::
siteOccupancy(int idx) const
{
    assert(0 <= idx && idx < this->countSites());
    return mvatoms[idx].GetOccupancy();
}


bool 
ObjCrystMoleculeAdapter::
siteAnisotropy(int idx) const
{
    assert(0 <= idx && idx < this->countSites());
    return !(mvatoms[idx].GetScatteringPower().IsIsotropic());
}


const R3::Matrix&
ObjCrystMoleculeAdapter::
siteCartesianUij(int idx) const
{
    assert(0 <= idx && idx < this->countSites());
    return mvuij[idx];
}


const string& 
ObjCrystMoleculeAdapter::
siteAtomType(int idx) const
{
    assert(0 <= idx && idx < this->countSites());
    return mvatoms[idx].GetScatteringPower().GetSymbol();
}


//////////////////////////////////////////////////////////////////////////////
// class ObjCrystAperiodicBondGenerator
//////////////////////////////////////////////////////////////////////////////

// Constructor ---------------------------------------------------------------

ObjCrystAperiodicBondGenerator::
ObjCrystAperiodicBondGenerator(const ObjCrystStructureAdapter* adpt) 
    : BaseBondGenerator(adpt), mpstructure(adpt), msymidx(0)
{
}

// Public Methods ------------------------------------------------------------

const R3::Vector& 
ObjCrystAperiodicBondGenerator::
r1() const
{
    size_t siteidx = this->site1();
    static R3::Vector rv;
    assert(msymidx < mpstructure->mvsym[siteidx].size());
    rv = mpstructure->mvsym[siteidx][msymidx];
    return rv;
}


double 
ObjCrystAperiodicBondGenerator::
msd0() const
{
    double rv = msd(this->site0(), 0);
    return rv;
}


double 
ObjCrystAperiodicBondGenerator::
msd1() const
{
    double rv = msd(this->site1(), msymidx);
    return rv;
}


bool 
ObjCrystAperiodicBondGenerator::
iterateSymmetry()
{
    // Rewind and iterate the symmetry counter. 
    // If that is finished, then we're done.
    this->uncache();
    if (++msymidx == mpstructure->mvsym[this->site1()].size())
    {
        return false;
    }
    return true;
}


void 
ObjCrystAperiodicBondGenerator::
rewindSymmetry()
{
    this->uncache();
    this->msymidx = 0;
}


double 
ObjCrystAperiodicBondGenerator::
msd(int siteidx, int symidx) const
{
    // Get the proper Uij tensor for the given site, and symmetry indices
    const R3::Matrix& UCart = mpstructure->mvuij[siteidx][symidx];
    const R3::Vector& s = this->r01();
    bool anisotropy = mpstructure->siteAnisotropy(siteidx);
    double rv = meanSquareDisplacement(UCart, s, anisotropy);
    return rv;
}

//////////////////////////////////////////////////////////////////////////////
// class ObjCrystPeriodicBondGenerator
//////////////////////////////////////////////////////////////////////////////

// Constructor ---------------------------------------------------------------

ObjCrystPeriodicBondGenerator::
ObjCrystPeriodicBondGenerator(const ObjCrystStructureAdapter* adpt) 
    : ObjCrystAperiodicBondGenerator(adpt)
{
}


// Public Methods ------------------------------------------------------------

void 
ObjCrystPeriodicBondGenerator::
rewind()
{
    // Delay msphere instantiation to here instead of in constructor,
    // so it is possible to use setRmin, setRmax.
    if (!msphere.get())
    {
        // Make a Lattice instance
        const Lattice& L = mpstructure->getLattice();
        double buffzone = L.ucMaxDiagonalLength();
        double rsphmin = this->getRmin() - buffzone;
        double rsphmax = this->getRmax() + buffzone;
        msphere.reset(new PointsInSphere(rsphmin, rsphmax, L));
    }
    // BaseBondGenerator::rewind calls this->rewindSymmetry,
    // which takes care of msphere and msymidx configuration
    this->ObjCrystAperiodicBondGenerator::rewind();
}


const R3::Vector& 
ObjCrystPeriodicBondGenerator::
r1() const
{
    const Lattice& L = mpstructure->getLattice();
    static R3::Vector rv; 
    rv = L.cartesian(msphere->mno()) + ObjCrystAperiodicBondGenerator::r1();
    return rv;
}


void
ObjCrystPeriodicBondGenerator::
setRmin(double rmin)
{
    // destroy msphere so it will be created on rewind with new rmin
    if (this->getRmin() != rmin)    msphere.reset(NULL);
    this->ObjCrystAperiodicBondGenerator::setRmin(rmin);
}


void 
ObjCrystPeriodicBondGenerator::
setRmax(double rmax)
{
    // destroy msphere so it will be created on rewind with new rmax
    if (this->getRmax() != rmax)    msphere.reset(NULL);
    this->ObjCrystAperiodicBondGenerator::setRmax(rmax);
}


bool 
ObjCrystPeriodicBondGenerator::
iterateSymmetry()
{
    // Iterate the sphere. If it is finished, rewind and iterate the symmetry
    // counter. If that is also finished, then we're done.
    this->uncache();
    msphere->next();
    if (msphere->finished())
    {
        if (++msymidx == mpstructure->mvsym[this->site1()].size())
        {
            return false;
        }
        msphere->rewind();
    }
    return true;
}


void 
ObjCrystPeriodicBondGenerator::
rewindSymmetry()
{
    this->msphere->rewind();
    ObjCrystAperiodicBondGenerator::rewindSymmetry();
}


//////////////////////////////////////////////////////////////////////////////
// class ObjCrystMoleculeBondGenerator
//////////////////////////////////////////////////////////////////////////////

// Constructor ---------------------------------------------------------------

ObjCrystMoleculeBondGenerator::
ObjCrystMoleculeBondGenerator(const ObjCrystMoleculeAdapter* adpt) 
    : BaseBondGenerator(adpt), mpstructure(adpt)
{
}


// Public Methods ------------------------------------------------------------

double 
ObjCrystMoleculeBondGenerator::
msd0() const
{
    double rv = msd(this->site0());
    return rv;
}


double 
ObjCrystMoleculeBondGenerator::
msd1() const
{
    double rv = msd(this->site1());
    return rv;
}


double 
ObjCrystMoleculeBondGenerator::
msd(int siteidx) const
{
    // Get the proper Uij tensor for the given site, and symmetry indices
    const R3::Matrix& UCart = mpstructure->mvuij[siteidx];
    const R3::Vector& s = this->r01();
    bool anisotropy = mpstructure->siteAnisotropy(siteidx);
    double rv = meanSquareDisplacement(UCart, s, anisotropy);
    return rv;
}


// Utility functions --------------------------------------------------------

R3::Matrix 
objcrystutil::
getUij(const ObjCryst::ScatteringPower* sp)
{
    R3::Matrix Uij;
    if (sp->IsIsotropic())
    {
        Uij(0,0) = Uij(1,1) = Uij(2,2) = sp->GetBiso() * BtoU;
        Uij(0,1) = Uij(1,0) = 0;
        Uij(0,2) = Uij(2,0) = 0;
        Uij(2,1) = Uij(1,2) = 0;
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

// Factory Function and its Registration -------------------------------------

StructureAdapterPtr
createPyObjCrystStructureAdapter(const boost::python::object& stru)
{
    using diffpy::importFromPyModule;
    boost::python::object cls_Crystal, None;
    cls_Crystal = importFromPyModule("pyobjcryst.crystal", "Crystal", None);
    StructureAdapterPtr rv;
    if (cls_Crystal.ptr() != Py_None &&
        PyObject_IsInstance(stru.ptr(), cls_Crystal.ptr()) == 1)
    {
        ObjCryst::Crystal* pcryst =
            boost::python::extract<ObjCryst::Crystal*>(stru);
        rv = createStructureAdapter(*pcryst);
    }
    return rv;
}

StructureAdapterPtr
createPyObjCrystMoleculeAdapter(const boost::python::object& stru)
{
    using diffpy::importFromPyModule;
    boost::python::object cls_Molecule, None;
    cls_Molecule = importFromPyModule("pyobjcryst.molecule", "Molecule", None);
    StructureAdapterPtr rv;
    if (cls_Molecule.ptr() != Py_None &&
        PyObject_IsInstance(stru.ptr(), cls_Molecule.ptr()) == 1)
    {
        ObjCryst::Molecule* pmolecule =
            boost::python::extract<ObjCryst::Molecule*>(stru);
        rv = createStructureAdapter(*pmolecule);
    }
    return rv;
}



bool reg_PyObjCrystStructureAdapter =
registerPythonStructureAdapterFactory(createPyObjCrystStructureAdapter);

bool reg_PyObjCrystMoleculeAdapter =
registerPythonStructureAdapterFactory(createPyObjCrystMoleculeAdapter);

// End of file
