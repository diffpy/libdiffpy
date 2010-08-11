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
* class SFTperiodictableNeutron
*
* Implementations of x-ray and neutron ScatteringFactorTable using Paul
* Kienzle periodictable library for Python.  The instances can be created
* using the createByType factory method, see the end of this file for
* available type strings.
*
* $Id$
*
*****************************************************************************/

#include <stdexcept>
#include <string>
#include <boost/python.hpp>

#include <diffpy/srreal/ScatteringFactorTable.hpp>
#include <diffpy/PythonInterface.hpp>

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


        double lookupatq(const string& smbl, double q) const
        {
            return this->fetchatq(smbl, q);
        }

    protected:

        // methods

        double fetchatq(const string& smbl, double q) const
        {
            diffpy::initializePython();
            static python::object fxrayatq = diffpy::importFromPyModule(
                    "periodictable.cromermann", "fxrayatq");
            double rv;
            try {
                rv = python::extract<double>(fxrayatq(smbl, q) + 0.0);
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


        double lookupatq(const string& smbl, double q) const
        {
            diffpy::initializePython();
            static python::object isotope = diffpy::importFromPyModule(
                    "periodictable", "elements").attr("isotope");
            double rv;
            try {
                python::object el = isotope(smbl);
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

bool reg_SFTperiodictableNeutron = (
        SFTperiodictableNeutron().registerThisType() &&
        ScatteringFactorTable::aliasType("periodictableneutron", "N")
        );

}   // namespace srreal
}   // namespace diffpy

// End of file
