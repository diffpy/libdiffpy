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
* class DiffPyStructureAdapter -- adapter to the Structure class from the
*     Python diffpy.Structure package.
*
* class DiffPyStructureBaseBondGenerator -- bond generator for
*     non-periodic structures
*
* class DiffPyStructurePeriodicBondGenerator -- bond generator for
*     periodic structures
*
* $Id$
*
*****************************************************************************/

#include <cassert>
#include <boost/python.hpp>
#include <boost/python/stl_iterator.hpp>

#include <diffpy/PythonInterface.hpp>
#include <diffpy/srreal/DiffPyStructureAdapter.hpp>
#include <diffpy/srreal/PythonStructureAdapter.hpp>
#include <diffpy/srreal/PointsInSphere.hpp>
#include <diffpy/srreal/PDFCalculator.hpp>
#include <diffpy/srreal/ScaleEnvelope.hpp>

using namespace std;
using namespace boost;
using namespace diffpy::srreal;

//////////////////////////////////////////////////////////////////////////////
// class DiffPyStructureAdapter
//////////////////////////////////////////////////////////////////////////////

// Constructor ---------------------------------------------------------------

DiffPyStructureAdapter::DiffPyStructureAdapter(python::object dpstru)
{
    mpystructure = dpstru;
    this->fetchPythonData();
}

// Public Methods ------------------------------------------------------------

BaseBondGenerator* DiffPyStructureAdapter::createBondGenerator() const
{
    // FIXME: hack for handling non-periodic structures
    // should diffpy.Structure get an isPeriodic method?
    BaseBondGenerator* bnds = this->isPeriodic() ?
        new DiffPyStructurePeriodicBondGenerator(this) :
        new DiffPyStructureBaseBondGenerator(this);
    return bnds;
}


int DiffPyStructureAdapter::countSites() const
{
    return mcartesian_positions.size();
}


double DiffPyStructureAdapter::numberDensity() const
{
    double rv = this->isPeriodic() ?
        (this->totalOccupancy() / mlattice.volume()) : 0.0;
    return rv;
}


const Lattice& DiffPyStructureAdapter::getLattice() const
{
    return mlattice;
}


const R3::Vector& DiffPyStructureAdapter::siteCartesianPosition(int idx) const
{
    assert(0 <= idx && idx < this->countSites());
    return mcartesian_positions[idx];
}


double DiffPyStructureAdapter::siteOccupancy(int idx) const
{
    assert(0 <= idx && idx < this->countSites());
    return moccupancies[idx];
}


bool DiffPyStructureAdapter::siteAnisotropy(int idx) const
{
    assert(0 <= idx && idx < this->countSites());
    return manisotropies[idx];
}


const R3::Matrix& DiffPyStructureAdapter::siteCartesianUij(int idx) const
{
    assert(0 <= idx && idx < this->countSites());
    return mcartesian_uijs[idx];
}


const string& DiffPyStructureAdapter::siteAtomType(int idx) const
{
    assert(0 <= idx && idx < this->countSites());
    return matomtypes[idx];
}


void DiffPyStructureAdapter::customPQConfig(PairQuantity* pq) const
{
    PDFCalculator* pdfc = dynamic_cast<PDFCalculator*>(pq);
    if (pdfc)  this->configurePDFCalculator(pdfc);
}

// Protected Methods ---------------------------------------------------------

void DiffPyStructureAdapter::fetchPythonData()
{
    // mlattice
    python::object lattice;
    lattice = mpystructure.attr("lattice");
    double a, b, c, alpha, beta, gamma;
    a = python::extract<double>(lattice.attr("a"));
    b = python::extract<double>(lattice.attr("b"));
    c = python::extract<double>(lattice.attr("c"));
    alpha = python::extract<double>(lattice.attr("alpha"));
    beta = python::extract<double>(lattice.attr("beta"));
    gamma = python::extract<double>(lattice.attr("gamma"));
    mlattice.setLatPar(a, b, c, alpha, beta, gamma);
    // atom properties
    mcartesian_positions.clear();
    moccupancies.clear();
    manisotropies.clear();
    mcartesian_uijs.clear();
    matomtypes.clear();
    int num_atoms = python::len(mpystructure);
    for (int i = 0; i < num_atoms; ++i)
    {
        python::object ai;
        ai = mpystructure[i];
        // mcartesian_positions
        python::object xyz = ai.attr("xyz");
        double x, y, z;
        x = python::extract<double>(xyz[0]);
        y = python::extract<double>(xyz[1]);
        z = python::extract<double>(xyz[2]);
        R3::Vector xyz_frac(x, y, z);
        R3::Vector xyz_cartn;
        xyz_cartn = mlattice.cartesian(xyz_frac);
        mcartesian_positions.push_back(xyz_cartn);
        // moccupancies
        double occupancy = python::extract<double>(ai.attr("occupancy"));
        moccupancies.push_back(occupancy);
        // manisotropies
        bool aniso = ai.attr("anisotropy");
        manisotropies.push_back(aniso);
        // mcartesian_uijs
        R3::Matrix Ufrac;
        python::object uflat = ai.attr("U").attr("flat");
        python::stl_input_iterator<double> ufirst(uflat), ulast;
        std::copy(ufirst, ulast, Ufrac.data());
        R3::Matrix Ucart = mlattice.cartesianMatrix(Ufrac);
        mcartesian_uijs.push_back(Ucart);
        // matomtypes
        string atp = python::extract<string>(ai.attr("element"));
        matomtypes.push_back(atp);
    }
    assert(int(mcartesian_positions.size()) == this->countSites());
    assert(int(mcartesian_positions.size()) == this->countSites());
    assert(int(manisotropies.size()) == this->countSites());
    assert(int(mcartesian_uijs.size()) == this->countSites());
    assert(int(matomtypes.size()) == this->countSites());
    // fetch the pdffit dictionary if present in mpystructure
    mpdffit.clear();
    if (PyObject_HasAttrString(mpystructure.ptr(), "pdffit"))
    {
        // get method of the stru.pdffit dictionary
        python::object stru_pdffit = mpystructure.attr("pdffit");
        python::object pfget = stru_pdffit.attr("get");
        mpdffit["scale"] = python::extract<double>(pfget("scale", 1.0));
        mpdffit["delta1"] = python::extract<double>(pfget("delta1", 0.0));
        mpdffit["delta2"] = python::extract<double>(pfget("delta2", 0.0));
        if (pfget("spdiameter").ptr() != Py_None)
        {
            mpdffit["spdiameter"] =
                python::extract<double>(pfget("spdiameter"));
        }
        if (pfget("stepcut").ptr() != Py_None)
        {
            mpdffit["stepcut"] = python::extract<double>(pfget("stepcut"));
        }
    }
}


void DiffPyStructureAdapter::configurePDFCalculator(PDFCalculator* pdfc) const
{
    // this is only needed if diffpy.Structure instance has pdffit attribute
    // with PDF-related structure properties
    if (mpdffit.empty())  return;
    // scale
    PDFEnvelopePtr envelope;
    envelope = PDFEnvelope::createByType("scale");
    envelope->setDoubleAttr("scale", mpdffit.at("scale"));
    pdfc->addEnvelope(envelope);
    // delta1, delta2 - set these only when using JeongPeakWidth model
    if (pdfc->getPeakWidthModel()->type() == "jeong")
    {
        pdfc->setDoubleAttr("delta1", mpdffit.at("delta1"));
        pdfc->setDoubleAttr("delta2", mpdffit.at("delta2"));
    }
    // spdiameter
    if (mpdffit.count("spdiameter"))
    {
        envelope = PDFEnvelope::createByType("sphericalshape");
        envelope->setDoubleAttr("spdiameter", mpdffit.at("spdiameter"));
        pdfc->addEnvelope(envelope);
    }
    // stepcut
    if (mpdffit.count("stepcut"))
    {
        envelope = PDFEnvelope::createByType("stepcut");
        envelope->setDoubleAttr("stepcut", mpdffit.at("stepcut"));
        pdfc->addEnvelope(envelope);
    }
}


bool DiffPyStructureAdapter::isPeriodic() const
{
    const Lattice& L = this->getLattice();
    bool rv = !(R3::EpsEqual(R3::identity(), L.base()));
    return rv;
}

//////////////////////////////////////////////////////////////////////////////
// class DiffPyStructureBaseBondGenerator
//////////////////////////////////////////////////////////////////////////////

// Constructor ---------------------------------------------------------------

DiffPyStructureBaseBondGenerator::DiffPyStructureBaseBondGenerator(
        const DiffPyStructureAdapter* adpt) : BaseBondGenerator(adpt)
{
    mdpstructure = adpt;
}

//////////////////////////////////////////////////////////////////////////////
// class DiffPyStructurePeriodicBondGenerator
//////////////////////////////////////////////////////////////////////////////

// Constructor ---------------------------------------------------------------

DiffPyStructurePeriodicBondGenerator::DiffPyStructurePeriodicBondGenerator(
        const DiffPyStructureAdapter* adpt
        ) : DiffPyStructureBaseBondGenerator(adpt)
{ }

// Public Methods ------------------------------------------------------------

void DiffPyStructurePeriodicBondGenerator::rewind()
{
    // Delay msphere instantiation to here instead of in constructor,
    // so it is possible to use setRmin, setRmax.
    if (!msphere.get())
    {
        const Lattice& L = mdpstructure->getLattice();
        double buffzone = L.ucMaxDiagonalLength();
        double rsphmin = this->getRmin() - buffzone;
        double rsphmax = this->getRmax() + buffzone;
        msphere.reset(new PointsInSphere(rsphmin, rsphmax, L));
    }
    // BaseBondGenerator::rewind calls this->rewindSymmetry,
    // which takes care of msphere configuration
    this->BaseBondGenerator::rewind();
}


void DiffPyStructurePeriodicBondGenerator::setRmin(double rmin)
{
    // destroy msphere so it will be created on rewind with new rmin
    if (this->getRmin() != rmin)    msphere.reset(NULL);
    this->BaseBondGenerator::setRmin(rmin);
}


void DiffPyStructurePeriodicBondGenerator::setRmax(double rmax)
{
    // destroy msphere so it will be created on rewind with new rmax
    if (this->getRmax() != rmax)    msphere.reset(NULL);
    this->BaseBondGenerator::setRmax(rmax);
}


const R3::Vector& DiffPyStructurePeriodicBondGenerator::r1() const
{
    static R3::Vector rv;
    const Lattice& L = mdpstructure->getLattice();
    rv = this->BaseBondGenerator::r1() + L.cartesian(msphere->mno());
    return rv;
}

// Protected Methods ---------------------------------------------------------

bool DiffPyStructurePeriodicBondGenerator::iterateSymmetry()
{
    this->uncache();
    msphere->next();
    return !msphere->finished();
}


void DiffPyStructurePeriodicBondGenerator::rewindSymmetry()
{
    this->uncache();
    msphere->rewind();
}

// Factory Function and its Registration -------------------------------------

StructureAdapterPtr createDiffPyStructureAdapter(python::object stru)
{
    using diffpy::importFromPyModule;
    python::object cls_Structure, None;
    cls_Structure = importFromPyModule("diffpy.Structure", "Structure", None);
    StructureAdapterPtr rv;
    if (cls_Structure.ptr() != Py_None &&
        PyObject_IsInstance(stru.ptr(), cls_Structure.ptr()) == 1)
    {
        rv.reset(new DiffPyStructureAdapter(stru));
    }
    return rv;
}

bool reg_DiffPyStructureAdapterFactory =
registerPythonStructureAdapterFactory(createDiffPyStructureAdapter);

// End of file
