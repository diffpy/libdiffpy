/*****************************************************************************
*
* diffpy.srreal     by DANSE Diffraction group
*                   Simon J. L. Billinge
*                   (c) 2010 Trustees of the Columbia University
*                   in the City of New York.  All rights reserved.
*
* File coded by:    Pavol Juhas
*
* See AUTHORS.txt for a list of people who contributed.
* See LICENSE.txt for license information.
*
******************************************************************************
*
* class SFTElectronNumber
*
* ScatteringFactorTable implementation where scattering power equals number
* of valence electrons and is constant with Q.
*
*****************************************************************************/

#include <stdexcept>
#include <string>
#include <cstdlib>
#include <cassert>
#include <fstream>
#include <boost/unordered_map.hpp>
#include <boost/scoped_ptr.hpp>

#include <diffpy/srreal/ScatteringFactorTable.hpp>
#include <diffpy/srreal/StructureAdapter.hpp>
#include <diffpy/runtimepath.hpp>
#include <diffpy/validators.hpp>

using namespace std;

namespace diffpy {
namespace srreal {

//////////////////////////////////////////////////////////////////////////////
// class SFTElectronNumber
//////////////////////////////////////////////////////////////////////////////

class SFTElectronNumber : public ScatteringFactorTable
{
    public:

        // HasClassRegistry methods

        ScatteringFactorTablePtr create() const
        {
            ScatteringFactorTablePtr rv(new SFTElectronNumber());
            return rv;
        }


        ScatteringFactorTablePtr clone() const
        {
            ScatteringFactorTablePtr rv(new SFTElectronNumber(*this));
            return rv;
        }

        const string& type() const
        {
            static string rv = "electronnumber";
            return rv;
        }

        // own methods

        const string& radiationType() const
        {
            static string rv = "EN";
            return rv;
        }


        double standardLookup(const string& smbl, double q) const
        {
            const boost::unordered_map<string,int>&
                entable = this->getElectronNumberTable();
            boost::unordered_map<string,int>::const_iterator en;
            en = entable.find(smbl);
            // try to build standard symbol when not found
            if (en == entable.end())
            {
                ostringstream smbl1;
                smbl1 << atomBareSymbol(smbl);
                int v = atomValence(smbl);
                if (v)  smbl1 << abs(v) << ((v > 0) ? '+' : '-');
                en = entable.find(smbl1.str());
            }
            // throw exception if still not found
            if (en == entable.end())
            {
                ostringstream emsg;
                emsg << "Unknown atom symbol '" << smbl << "'.";
                throw invalid_argument(emsg.str());
            }
            return en->second;
        }

    private:

        // methods

        const boost::unordered_map<string,int>& getElectronNumberTable() const
        {
            using diffpy::runtimepath::datapath;
            using diffpy::validators::ensureFileOK;
            typedef boost::unordered_map<string,int> mapping;
            static boost::scoped_ptr<mapping> entable;
            if (!entable)
            {
                entable.reset(new mapping);
                string ionfile = datapath("ionlist.dat");
                ifstream fp0(ionfile.c_str());
                ensureFileOK(ionfile, fp0);
                string line;
                int lineno = 0;
                while (getline(fp0, line))
                {
                    ++lineno;
                    istringstream fp(line);
                    string element;
                    int z = 0;
                    fp >> element;
                    if (element.empty() || element[0] == '#')  continue;
                    fp >> z;
                    if (!fp)  this->throwBadLine(ionfile, lineno);
                    (*entable)[element] = z;
                    for (int v; fp >> v;)
                    {
                        ostringstream smbl;
                        smbl << element << abs(v) << ((v > 0) ? '+' : '-');
                        (*entable)[smbl.str()] = z - v;
                    }
                }
                assert(436 <= entable->size());
            }
            return *entable;
        }


        void throwBadLine(const string& fname, int lineno) const
        {
            ostringstream emsg;
            emsg << "Invalid data format in " <<
                fname << " line " << lineno << '.';
            throw runtime_error(emsg.str());
        }


};  // class SFTElectronNumber

// Registration --------------------------------------------------------------

bool reg_SFTElectronNumber = SFTElectronNumber().registerThisType() &&
        ScatteringFactorTable::aliasType("electronnumber", "EN");

}   // namespace srreal
}   // namespace diffpy

// End of file
