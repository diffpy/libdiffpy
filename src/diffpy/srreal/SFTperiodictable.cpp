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
* class SFTperiodictableXray
* class SFTElectron
* class SFTperiodictableNeutron
*
* Implementations of x-ray and neutron ScatteringFactorTable using Paul
* Kienzle periodictable library for Python.  The instances can be created
* using the createByType factory method, see the end of this file for
* available type strings.
*
* SFTElectron gives Q-dependent electron scattering factor according to
* the formula in International Tables Volume C, page 224.
* The formula diverges at Q = 0, where SFTElectron returns DOUBLE_MAX.
*
* $Id$
*
*****************************************************************************/

#include <diffpy/PythonInterface.hpp>

#include <stdexcept>

#include <diffpy/srreal/ScatteringFactorTable.hpp>
#include <diffpy/mathutils.hpp>

using namespace std;
namespace python = boost::python;

namespace diffpy {
namespace srreal {

//////////////////////////////////////////////////////////////////////////////
// class SFTperiodictableXray
//////////////////////////////////////////////////////////////////////////////

class SFTperiodictableXray : public ScatteringFactorTable
{
    public:

        // HasClassRegistry methods

        ScatteringFactorTablePtr create() const
        {
            ScatteringFactorTablePtr rv(new SFTperiodictableXray());
            return rv;
        }


        ScatteringFactorTablePtr clone() const
        {
            ScatteringFactorTablePtr rv(new SFTperiodictableXray(*this));
            return rv;
        }

        const string& type() const
        {
            static string rv = "periodictablexray";
            return rv;
        }

        // own methods

        const string& radiationType() const
        {
            static string rv = "X";
            return rv;
        }


        double standardLookup(const string& smbl, double q) const
        {
            diffpy::initializePython();
            static python::object fxrayatq = diffpy::importFromPyModule(
                    "periodictable.cromermann", "fxrayatq");
            const string zp = "0+";
            string::size_type pos = smbl.size() - zp.size();
            string::size_type pe = smbl.rfind(zp, pos);
            string smblnozp = smbl.substr(0, pe);
            double rv;
            try {
                rv = python::extract<double>(fxrayatq(smblnozp, q) + 0.0);
            }
            catch (python::error_already_set e) {
                string emsg = diffpy::getPythonErrorString();
                PyErr_Clear();
                throw invalid_argument(emsg);
            }
            return rv;
        }

};  // class SFTperiodictableXray

//////////////////////////////////////////////////////////////////////////////
// class SFTElectron
//////////////////////////////////////////////////////////////////////////////

class SFTElectron : public SFTperiodictableXray
{
    public:

        // HasClassRegistry methods

        ScatteringFactorTablePtr create() const
        {
            ScatteringFactorTablePtr rv(new SFTElectron());
            return rv;
        }


        ScatteringFactorTablePtr clone() const
        {
            ScatteringFactorTablePtr rv(new SFTElectron(*this));
            return rv;
        }

        const string& type() const
        {
            static string rv = "electron";
            return rv;
        }

        // own methods

        const string& radiationType() const
        {
            static string rv = "E";
            return rv;
        }


        double standardLookup(const string& smbl, double q) const
        {
            using namespace diffpy::mathutils;
            // resolve Z first so that invalid symbol can throw an exception
            double Z = round(
                this->SFTperiodictableXray::standardLookup(smbl, 0.0));
            if (eps_eq(q, 0.0))     return DOUBLE_MAX;
            double stol = q / (4 * M_PI);
            double rv = 0.023934 *
                (Z - this->SFTperiodictableXray::standardLookup(smbl, q)) /
                (stol * stol);
            return rv;
        }

};  // class SFTElectron

//////////////////////////////////////////////////////////////////////////////
// class SFTperiodictableNeutron
//////////////////////////////////////////////////////////////////////////////

class SFTperiodictableNeutron : public ScatteringFactorTable
{
    public:

        // HasClassRegistry methods

        ScatteringFactorTablePtr create() const
        {
            ScatteringFactorTablePtr rv(new SFTperiodictableNeutron());
            return rv;
        }


        ScatteringFactorTablePtr clone() const
        {
            ScatteringFactorTablePtr rv(new SFTperiodictableNeutron(*this));
            return rv;
        }


        const string& type() const
        {
            static string rv = "periodictableneutron";
            return rv;
        }

        // own methods

        const string& radiationType() const
        {
            static string rv = "N";
            return rv;
        }


        double standardLookup(const string& smbl, double q) const
        {
            diffpy::initializePython();
            static python::object isotope = diffpy::importFromPyModule(
                    "periodictable", "elements").attr("isotope");
            double rv;
            string::size_type pe = smbl.find_last_not_of("+-012345678 \t");
            string smblnocharge = smbl.substr(0, pe + 1);
            try {
                python::object el = isotope(smblnocharge);
                python::object b_c = el.attr("neutron").attr("b_c");
                rv = python::extract<double>(b_c);
            }
            catch (python::error_already_set) {
                string emsg = diffpy::getPythonErrorString();
                PyErr_Clear();
                throw invalid_argument(emsg);
            }
            return rv;
        }

};  // class SFTperiodictableNeutron

// Registration --------------------------------------------------------------

bool reg_SFTperiodictableXray = (
        SFTperiodictableXray().registerThisType() &&
        ScatteringFactorTable::aliasType("periodictablexray", "X")
        );

bool reg_SFTElectron = (
        SFTElectron().registerThisType() &&
        ScatteringFactorTable::aliasType("electron", "E")
        );

bool reg_SFTperiodictableNeutron = (
        SFTperiodictableNeutron().registerThisType() &&
        ScatteringFactorTable::aliasType("periodictableneutron", "N")
        );

}   // namespace srreal
}   // namespace diffpy

// End of file
