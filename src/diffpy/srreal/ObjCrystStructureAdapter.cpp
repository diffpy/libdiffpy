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
* ******************************************************************************
*
* class ObjCrystStructureAdapter -- adapter to the Crystal class from
* ObjCryst++.
* class ObjCrystBondGenerator -- Generate bonds from ObjCrystStructureAdapter.
*
* class ObjCrystBondGenerator -- bond generator
*
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
#include <diffpy/srreal/PDFCalculator.hpp>
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

// Constructor ---------------------------------------------------------------

const double ObjCrystStructureAdapter::toler = 1e-5;

ObjCrystStructureAdapter::
ObjCrystStructureAdapter(const ObjCryst::Crystal& cryst) : mpcryst(&cryst)
{
    mlattice.setLatPar( mpcryst->GetLatticePar(0), 
                        mpcryst->GetLatticePar(1),
                        mpcryst->GetLatticePar(2), 
                        rtod * mpcryst->GetLatticePar(3),
                        rtod * mpcryst->GetLatticePar(4), 
                        rtod * mpcryst->GetLatticePar(5) );
    this->getUnitCell();
}

// Public Methods ------------------------------------------------------------

BaseBondGenerator* 
ObjCrystStructureAdapter::
createBondGenerator() const
{
    BaseBondGenerator* bnds = new ObjCrystBondGenerator(this);
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
    double rv = this->totalOccupancy() / mpcryst->GetVolume();
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


// Private Methods -----------------------------------------------------------

/* Get the conventional unit cell from the crystal. */
void
ObjCrystStructureAdapter::
getUnitCell()
{
    // local constants
    R3::Matrix zeros; zeros = 0.0;

    // Expand each scattering component in the primitive cell and record the
    // new scatterers.
    const ObjCryst::ScatteringComponentList& scl =
        mpcryst->GetScatteringComponentList();
    size_t nbComponent = scl.GetNbComponent();

    // Various things we need to know about the symmetry operations
    const ObjCryst::SpaceGroup& spacegroup = mpcryst->GetSpaceGroup();
    size_t nbSymmetrics = spacegroup.GetNbSymmetrics();
    size_t nbTrans = spacegroup.GetNbTranslationVectors();
    size_t cmult = spacegroup.IsCentrosymmetric() ? 2 : 1;
    size_t nbRot = nbSymmetrics / nbTrans / cmult;
    assert( nbSymmetrics = nbTrans * nbRot * cmult );

    double x, y, z, junk;
    CrystMatrix<double> symmetricsCoords;

    mvsc.clear();
    mvsym.clear();
    mvuij.clear();
    mvsc.reserve(nbComponent);
    mvsym.reserve(nbComponent);
    mvuij.reserve(nbComponent);

    typedef std::set<R3::Vector, R3::EpsCompare> SymPosSet;

    // Get the symmetry rotations
    std::vector< R3::Matrix > rotations = getRotations();

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
        SymPosSet symset = SymPosSet(R3::EpsCompare(toler));
        SymPosVec symvec = SymPosVec();
        SymUijVec uijvec = SymUijVec();
        size_t numsym = 0;
        size_t rotidx = 0;

        // Store Uij in cartesian
        R3::Matrix UCart = getUCart(sp);

        // Get all the symmetric coordinates. Symmetric coordinates are created
        // by translation vector in the outer loop, and rotation matrix in the
        // inner loop. The first translation is [0,0,0]. We need to know this
        // to determine the rotation matrices.
        symmetricsCoords = spacegroup.GetAllSymmetrics(
            scl(i).mX, scl(i).mY, scl(i).mZ);

        // Collect the unique symmetry positions.
        for (size_t j = 0; j < nbSymmetrics; ++j, ++rotidx)
        {
            x = modf(symmetricsCoords(j, 0), &junk);
            y = modf(symmetricsCoords(j, 1), &junk);
            z = modf(symmetricsCoords(j, 2), &junk);
            if (fabs(x) < toler) x = 0;
            if (fabs(y) < toler) y = 0;
            if (fabs(z) < toler) z = 0;
            if (x < 0) x += 1.;
            if (y < 0) y += 1.;
            if (z < 0) z += 1.;

            // Get this in cartesian
            mpcryst->FractionalToOrthonormalCoords(x, y, z);

            // Record the position
            R3::Vector xyz; 
            xyz = x, y, z;

            // We use this to filter unique positions
            symset.insert(xyz);

            // Keep track of the rotation index. We use the fact that the
            // translations are in the outer loop when the symmetric
            // coordinates are created. This is a hack.
            if(rotidx >= nbRot) rotidx -= nbRot;

            // See if we've got a new position
            if(symset.size() > numsym)
            {
                ++numsym;
                
                // Store this in the symvec so we are assured that the order
                // will not change.
                symvec.push_back(xyz);

                if(numsym > 0)
                {
                    // Get the symmetry operation used to generate this postion.
                    const R3::Matrix& M = rotations[rotidx];

                    // rotate the UCart matrix
                    R3::Matrix tmp, Urot; 
                    tmp = R3::product(UCart, M);
                    Urot = R3::product(R3::transpose(M), tmp);
                    // Check for centro-symmetry
                    if(j >= nbRot * nbTrans)
                    {
                        Urot = R3::transpose(Urot);
                    }
                    uijvec.push_back(Urot);
                }
                else
                {
                    uijvec.push_back(UCart);
                }
            }
        }

        assert( uijvec.size() == symvec.size() );

        // Store symmetric positions and Uij matrices
        mvsym.push_back(symvec);
        mvuij.push_back(uijvec);

    }

}

R3::Matrix 
ObjCrystStructureAdapter::
getUCart(const ObjCryst::ScatteringPower* sp) const
{
    R3::Matrix UCart;
    if (sp->IsIsotropic())
    {
        UCart(0,0) = UCart(1,1) = UCart(2,2) = sp->GetBiso() * BtoU;
        UCart(0,1) = UCart(1,0) = 0;
        UCart(0,2) = UCart(2,0) = 0;
        UCart(2,1) = UCart(1,2) = 0;
    }
    else
    {
        R3::Matrix Ufrac;
        Ufrac(0,0) = sp->GetBij(1,1) * BtoU;
        Ufrac(1,1) = sp->GetBij(2,2) * BtoU;
        Ufrac(2,2) = sp->GetBij(3,3) * BtoU;
        Ufrac(0,1) = Ufrac(1,0) = sp->GetBij(1,2) * BtoU;
        Ufrac(0,2) = Ufrac(2,0) = sp->GetBij(1,3) * BtoU;
        Ufrac(1,2) = Ufrac(2,1) = sp->GetBij(2,3) * BtoU;
        UCart = mlattice.cartesianMatrix(Ufrac);
    }

    return UCart;

}

std::vector< R3::Matrix >
ObjCrystStructureAdapter::
getRotations() const
{
    const ObjCryst::SpaceGroup& spacegroup = mpcryst->GetSpaceGroup();
    size_t nbSymmetrics = spacegroup.GetNbSymmetrics();
    size_t nbTrans = spacegroup.GetNbTranslationVectors();
    size_t cmult = spacegroup.IsCentrosymmetric() ? 2 : 1;
    size_t nbRot = nbSymmetrics / nbTrans / cmult;

    std::vector< R3::Matrix > rotations(nbRot);
    rotations[0] = R3::identity();

    // Expand a single postion that is not at the origin and use this to
    // determine the symmetry operations. We use the fact that the symmetric
    // postions are expanded using the [0, 0, 0] translation first.
    double x = 0.1, y = 0.11, z = 0.2;
    CrystMatrix<double> coords = spacegroup.GetAllSymmetrics(x, y, z);
    mpcryst->FractionalToOrthonormalCoords(x, y, z);
    R3::Vector p0; 
    p0 = x, y, z;

    for(size_t i = 1; i < nbRot; ++i)
    {
        x = coords(i, 0);
        y = coords(i, 1);
        z = coords(i, 2);
        mpcryst->FractionalToOrthonormalCoords(x, y, z);

        R3::Vector p1;
        p1 = x, y, z;

        // Get _a_ rotation matrix that takes p0 to p1
        R3::Matrix M = R3::rotationfrom(p0, p1);

        rotations[i] = M;
    }

    return rotations;

}

//////////////////////////////////////////////////////////////////////////////
// class ObjCrystBondGenerator
//////////////////////////////////////////////////////////////////////////////

// Constructor ---------------------------------------------------------------

ObjCrystBondGenerator::
ObjCrystBondGenerator(const ObjCrystStructureAdapter* adpt) 
    : BaseBondGenerator(adpt), mpstructure(adpt), msymidx(0)
{
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


const R3::Vector& 
ObjCrystBondGenerator::
r1() const
{
    size_t siteidx = this->site1();
    static R3::Vector rv;
    const Lattice& L = mpstructure->getLattice();
    assert(msymidx < mpstructure->mvsym[siteidx].size());
    rv = mpstructure->mvsym[siteidx][msymidx];
    rv += L.cartesian(msphere->mno());
    return rv;
}


double 
ObjCrystBondGenerator::
msd0() const
{
    double rv = msd(this->site0(), 0);
    return rv;
}


double 
ObjCrystBondGenerator::
msd1() const
{
    double rv = msd(this->site1(), msymidx);
    return rv;
}

bool 
ObjCrystBondGenerator::
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
ObjCrystBondGenerator::
rewindSymmetry()
{
    this->uncache();
    msphere->rewind();
    msymidx = 0;
}


double 
ObjCrystBondGenerator::
msd(int siteidx, int symidx) const
{
    // Get the proper Uij tensor for the given site, and symmetry indices
    const R3::Matrix& UCart = mpstructure->mvuij[siteidx][symidx];
    const R3::Vector& s = this->r01();
    bool anisotropy = mpstructure->siteAnisotropy(siteidx);
    double rv = meanSquareDisplacement(UCart, s, anisotropy);
    return rv;
}


// Factory Function and its Registration -------------------------------------

StructureAdapter*
createPyObjCrystStructureAdapter(const boost::python::object& stru)
{
    using diffpy::importFromPyModule;
    boost::python::object cls_Crystal;
    cls_Crystal = importFromPyModule("pyobjcryst.crystal", "Crystal");
    StructureAdapter* rv = NULL;
    if (cls_Crystal.ptr() != Py_None &&
        PyObject_IsInstance(stru.ptr(), cls_Crystal.ptr()) == 1)
    {
        const ObjCryst::Crystal* pcryst =
            boost::python::extract<ObjCryst::Crystal*>(stru);
        rv = createPQAdapter(*pcryst);
    }
    return rv;
}


bool reg_PyObjCrystStructureAdapter =
registerPythonStructureAdapterFactory(createPyObjCrystStructureAdapter);

// End of file
