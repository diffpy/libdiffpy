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
*****************************************************************************/

#include <diffpy/PythonInterface.hpp>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <set>

#include <diffpy/serialization.hpp>
#include <diffpy/srreal/PythonStructureAdapter.hpp>
#include <diffpy/srreal/ObjCrystStructureAdapter.hpp>
#include <diffpy/srreal/PointsInSphere.hpp>

using namespace std;
using namespace diffpy::srreal;

namespace {

// Local constants -----------------------------------------------------------

const double rtod = 180 / M_PI;
const double UtoB = 8 * M_PI * M_PI;
const double BtoU = 1.0 / UtoB;

// Utility functions ---------------------------------------------------------

R3::Matrix
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

}   // namespace


//////////////////////////////////////////////////////////////////////////////
// class ObjCrystStructureAdapter
//////////////////////////////////////////////////////////////////////////////

// static data

const double ObjCrystStructureAdapter::mtoler = 1e-5;

// Constructor ---------------------------------------------------------------

ObjCrystStructureAdapter::
ObjCrystStructureAdapter(const ObjCryst::Crystal& cryst)
{
    using ObjCryst::Crystal;
    mlattice.setLatPar( cryst.GetLatticePar(0),
                        cryst.GetLatticePar(1),
                        cryst.GetLatticePar(2),
                        rtod * cryst.GetLatticePar(3),
                        rtod * cryst.GetLatticePar(4),
                        rtod * cryst.GetLatticePar(5) );
    // The dynamic population correction in ObjCryst is used by its interal
    // calculators, but slows down the calculation of atom positions. This is
    // also required to get the proper ScatteringComponentList for aperiodic
    // structures. Since we're not using any ObjCryst calculators, we turn it
    // off momentarily.
    //
    // NOTE: 2012-08-20 Pavol - I have not seen this issue and want to use
    // molecule rotations in a crystal structure of an organic compound.
    // Code changed to leave the population correction alone.
    // It can be set as desired using the pyobjcryst interface.
    //
    // int usepopcorr = cryst.GetUseDynPopCorr();
    // const_cast<Crystal&>(cryst).SetUseDynPopCorr(0);
    this->getUnitCell(cryst);
    // Undo change
    // const_cast<Crystal&>(cryst).SetUseDynPopCorr(usepopcorr);
}


// Public Methods ------------------------------------------------------------

BaseBondGeneratorPtr
ObjCrystStructureAdapter::
createBondGenerator() const
{
    BaseBondGeneratorPtr bnds(new ObjCrystBondGenerator(shared_from_this()));
    return bnds;
}


int
ObjCrystStructureAdapter::
countSites() const
{
    return mvsym.size();
}


double
ObjCrystStructureAdapter::
numberDensity() const
{
    double rv = this->totalOccupancy() / mlattice.volume();
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
    return moccupancies[idx];
}


bool
ObjCrystStructureAdapter::
siteAnisotropy(int idx) const
{
    assert(0 <= idx && idx < this->countSites());
    return manisotropies[idx];
}


int
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
    return matomtypes[idx];
}

// Private Methods -----------------------------------------------------------

/* Get the conventional unit cell from the crystal. */
void
ObjCrystStructureAdapter::
getUnitCell(const ObjCryst::Crystal& cryst)
{

    // Expand each scattering component in the primitive cell and record the
    // new scatterers.
    const ObjCryst::ScatteringComponentList& scl =
        cryst.GetScatteringComponentList();
    size_t nbComponent = scl.GetNbComponent();

    // Various things we need to know about the symmetry operations
    const ObjCryst::SpaceGroup& spacegroup = cryst.GetSpaceGroup();
    size_t nbSymmetrics = spacegroup.GetNbSymmetrics();

    double x, y, z, junk;
    CrystMatrix<double> symmetricsCoords;

    moccupancies.clear();
    manisotropies.clear();
    matomtypes.clear();
    mvsym.clear();
    mvuij.clear();
    moccupancies.reserve(nbComponent);
    manisotropies.reserve(nbComponent);
    matomtypes.reserve(nbComponent);
    mvsym.reserve(nbComponent);
    mvuij.reserve(nbComponent);

    using diffpy::mathutils::EpsilonLess;
    typedef std::set<R3::Vector, EpsilonLess> SymPosSet;

    // Get the symmetry operations
    const std::vector<ObjCryst::SpaceGroup::SMx>& symops =
        spacegroup.GetSymmetryOperations();
    const size_t symsize = symops.size();

    // For each scattering component, find its position in the primitive cell
    // and expand that position.
    for (size_t i = 0; i < nbComponent; ++i)
    {

        const ObjCryst::ScatteringComponent& sc = scl(i);
        const ObjCryst::ScatteringPower* sp = sc.mpScattPow;

        // Skip over this if it is a dummy atom. A dummy atom has no
        // mpScattPow, and therefore no type. It's just in a structure as a
        // reference position.
        if (sp == NULL) continue;

        moccupancies.push_back(sc.mOccupancy);
        manisotropies.push_back(!(sp->IsIsotropic()));
        matomtypes.push_back(sp->GetSymbol());
        SymPosSet symset = SymPosSet(EpsilonLess(mtoler));
        SymPosVec symvec = SymPosVec();
        SymUijVec uijvec = SymUijVec();
        size_t numsym = 0;

        // Store Uij.
        R3::Matrix Uij = getUij(sp);

        // Get all the symmetric coordinates. Symmetric coordinates are created
        // by translation vector in the outer loop, and rotation matrix in the
        // inner loop. The first translation is [0,0,0]. We need to know this
        // to determine the rotation matrices.
        symmetricsCoords = spacegroup.GetAllSymmetrics(
            sc.mX, sc.mY, sc.mZ);

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

    // Final sanity checks.
    assert(int(moccupancies.size()) == this->countSites());
    assert(int(manisotropies.size()) == this->countSites());
    assert(int(matomtypes.size()) == this->countSites());
    assert(int(mvsym.size()) == this->countSites());
    assert(int(mvuij.size()) == this->countSites());
}

//////////////////////////////////////////////////////////////////////////////
// class ObjCrystBondGenerator
//////////////////////////////////////////////////////////////////////////////

// Constructor ---------------------------------------------------------------

ObjCrystBondGenerator::
ObjCrystBondGenerator(StructureAdapterConstPtr adpt)
    : BaseBondGenerator(adpt)
{
    mpstructure = dynamic_cast<const ObjCrystStructureAdapter*>(adpt.get());
    assert(mpstructure);
    msymidx = 0;
}

// Public Methods ------------------------------------------------------------

void
ObjCrystBondGenerator::
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
    this->BaseBondGenerator::rewind();
}


const R3::Matrix&
ObjCrystBondGenerator::
Ucartesian1() const
{
    const R3::Matrix& rv = mpstructure->mvuij[this->site1()][msymidx];
    return rv;
}


void
ObjCrystBondGenerator::
setRmin(double rmin)
{
    // destroy msphere so it will be created on rewind with new rmin
    if (this->getRmin() != rmin)    msphere.reset(NULL);
    this->BaseBondGenerator::setRmin(rmin);
}


void
ObjCrystBondGenerator::
setRmax(double rmax)
{
    // destroy msphere so it will be created on rewind with new rmax
    if (this->getRmax() != rmax)    msphere.reset(NULL);
    this->BaseBondGenerator::setRmax(rmax);
}


bool
ObjCrystBondGenerator::
iterateSymmetry()
{
    // Iterate the sphere. If it is finished, rewind and iterate the symmetry
    // counter. If that is also finished, then we're done.
    msphere->next();
    if (msphere->finished())
    {
        if (++msymidx == mpstructure->mvsym[this->site1()].size())
        {
            return false;
        }
        msphere->rewind();
    }
    this->updater1();
    return true;
}


void
ObjCrystBondGenerator::
rewindSymmetry()
{
    this->msphere->rewind();
    msymidx = 0;
    this->updater1();
}

// Private Methods -----------------------------------------------------------

void
ObjCrystBondGenerator::
updater1()
{
    size_t siteidx = this->site1();
    assert(msymidx < mpstructure->mvsym[siteidx].size());
    const Lattice& L = mpstructure->getLattice();
    mr1 = L.cartesian(msphere->mno()) + mpstructure->mvsym[siteidx][msymidx];
    this->updateDistance();
}

//////////////////////////////////////////////////////////////////////////////
// class ObjCrystMoleculeAdapter
//////////////////////////////////////////////////////////////////////////////

// Constructor ---------------------------------------------------------------

ObjCrystMoleculeAdapter::
ObjCrystMoleculeAdapter(const ObjCryst::Molecule& molecule)
{
    using ObjCryst::MolAtom;

    mlattice.setLatPar(1, 1, 1, 90, 90, 90);

    size_t nbComponent = molecule.GetNbComponent();
    // Initialize all cache stores.
    moccupancies.clear();
    manisotropies.clear();
    matomtypes.clear();
    mvpos.clear();
    mvuij.clear();
    moccupancies.reserve(nbComponent);
    manisotropies.reserve(nbComponent);
    matomtypes.reserve(nbComponent);
    mvpos.reserve(nbComponent);
    mvuij.reserve(nbComponent);

    for (size_t i = 0; i < nbComponent; ++i)
    {

        const MolAtom& atom = molecule.GetAtom(i);
        if (atom.IsDummy()) continue;

        // Store the atom
        moccupancies.push_back(atom.GetOccupancy());
        const ObjCryst::ScatteringPower* sp = &(atom.GetScatteringPower());
        manisotropies.push_back(!(sp->IsIsotropic()));
        matomtypes.push_back(sp->GetSymbol());

        // Store position
        R3::Vector xyz;
        xyz = atom.X(), atom.Y(), atom.Z();
        mvpos.push_back(xyz);

        // Store Uij
        R3::Matrix Uij = getUij(&atom.GetScatteringPower());
        mvuij.push_back(Uij);
    }

    // final sanity checks
    assert(int(moccupancies.size()) == this->countSites());
    assert(int(manisotropies.size()) == this->countSites());
    assert(int(matomtypes.size()) == this->countSites());
    assert(int(mvpos.size()) == this->countSites());
    assert(int(mvuij.size()) == this->countSites());

}

// Public Methods ------------------------------------------------------------

BaseBondGeneratorPtr
ObjCrystMoleculeAdapter::
createBondGenerator() const
{

    BaseBondGeneratorPtr bnds(
            new ObjCrystMoleculeBondGenerator(shared_from_this()));
    return bnds;
}


int
ObjCrystMoleculeAdapter::
countSites() const
{
    return mvpos.size();
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
    return moccupancies[idx];
}


bool
ObjCrystMoleculeAdapter::
siteAnisotropy(int idx) const
{
    assert(0 <= idx && idx < this->countSites());
    return manisotropies[idx];
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
    return matomtypes[idx];
}


//////////////////////////////////////////////////////////////////////////////
// class ObjCrystMoleculeBondGenerator
//////////////////////////////////////////////////////////////////////////////

// Constructor ---------------------------------------------------------------

ObjCrystMoleculeBondGenerator::
ObjCrystMoleculeBondGenerator(StructureAdapterConstPtr adpt)
    : BaseBondGenerator(adpt)
{
    mpstructure = dynamic_cast<const ObjCrystMoleculeAdapter*>(adpt.get());
    assert(mpstructure);
}

// Factory Function and its Registration -------------------------------------

StructureAdapterPtr
createPyObjCrystStructureAdapter(boost::python::object stru)
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
createPyObjCrystMoleculeAdapter(boost::python::object stru)
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

// Serialization -------------------------------------------------------------

DIFFPY_INSTANTIATE_SERIALIZE(diffpy::srreal::ObjCrystStructureAdapter)
BOOST_CLASS_EXPORT(diffpy::srreal::ObjCrystStructureAdapter)

DIFFPY_INSTANTIATE_SERIALIZE(diffpy::srreal::ObjCrystMoleculeAdapter)
BOOST_CLASS_EXPORT(diffpy::srreal::ObjCrystMoleculeAdapter)

// End of file
