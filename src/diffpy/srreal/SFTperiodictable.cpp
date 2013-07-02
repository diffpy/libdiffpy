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
* class SFTperiodictableNeutron
*
* Implementations of x-ray and neutron ScatteringFactorTable using Paul
* Kienzle periodictable library for Python.  The instances can be created
* using the createByType factory method, see the end of this file for
* available type strings.
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

bool reg_SFTperiodictableNeutron = (
        SFTperiodictableNeutron().registerThisType() &&
        ScatteringFactorTable::aliasType("periodictableneutron", "N")
        );

}   // namespace srreal
}   // namespace diffpy

// End of file
