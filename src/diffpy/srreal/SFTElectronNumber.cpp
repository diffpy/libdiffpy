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
* $Id$
*
*****************************************************************************/

#include <stdexcept>
#include <string>
#include <cstdlib>
#include <cassert>
#include <boost/unordered_map.hpp>

#include <diffpy/srreal/ScatteringFactorTable.hpp>
#include <diffpy/srreal/StructureAdapter.hpp>

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


        double lookupatq(const string& smbl, double q) const
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

        // data
        static const char* mendata_lines[];

        // methods

        const boost::unordered_map<string,int>& getElectronNumberTable() const
        {
            static boost::unordered_map<string,int> entable;
            if (entable.empty())
            {
                const char** line;
                for (line = mendata_lines; *line != NULL; ++line)
                {
                    istringstream fp(*line);
                    string element;
                    int z = 0;
                    fp >> element >> z;
                    assert(fp);
                    entable[element] = z;
                    for (int v = 0; fp >> v;)
                    {
                        ostringstream smbl;
                        smbl << element << abs(v) << ((v > 0) ? '+' : '-');
                        entable[smbl.str()] = z - v;
                    }
                }
                assert(436 == entable.size());
            }
            return entable;
        }

};  // class SFTElectronNumber

// Ions from bvparm2009cif_lines and Python periodictable package ------------

const char* SFTElectronNumber::mendata_lines[] = {
    // element atomic-number list-of-ions
    "H 1 -1 +1",
    "He 2",
    "Li 3 -1 +1",
    "Be 4 -2 +2",
    "B 5 +3",
    "C 6 -4 +2 +4",
    "N 7 -3 -2 +3 +5",
    "O 8 -2 -1 +1 +2",
    "F 9 -1 +1",
    "Ne 10",
    "Na 11 -1 +1",
    "Mg 12 -2 +2",
    "Al 13 -3 +3",
    "Si 14 -4 +4",
    "Siva 14",
    "P 15 -3 +3 +4 +5",
    "S 16 -4 -2 +2 +4 +6",
    "Cl 17 -2 -1 +1 +3 +5 +7",
    "Ar 18",
    "K 19 -1 +1",
    "Ca 20 -2 +2",
    "Sc 21 -3 +3",
    "Ti 22 -4 -3 -2 +2 +3 +4",
    "V 23 -5 -3 -2 +1 +2 +3 +4 +5",
    "Cr 24 -3 -2 +2 +3 +4 +5 +6",
    "Mn 25 -4 -3 -2 +2 +3 +4 +6 +7",
    "Fe 26 -3 -2 +2 +3 +4 +6",
    "Co 27 -3 -2 -1 +1 +2 +3 +4",
    "Ni 28 -3 -2 +2 +3 +4",
    "Cu 29 -2 -1 +1 +2 +3",
    "Zn 30 -2 +2",
    "Ga 31 -3 +1 +3",
    "Ge 32 -4 +4",
    "As 33 -3 +2 +3 +5",
    "Se 34 -2 -1 +2 +4 +6",
    "Br 35 -1 +1 +3 +5 +7",
    "Kr 36 +2",
    "Rb 37 -1 +1",
    "Sr 38 -2 +2",
    "Y 39 +3",
    "Zr 40 -4 +2 +4",
    "Nb 41 -5 -3 +3 +4 +5",
    "Mo 42 -6 -5 -3 +2 +3 +4 +5 +6 +7",
    "Tc 43 +3 +4 +5 +6 +7",
    "Ru 44 -4 -3 +2 +3 +4 +5 +6 +7",
    "Rh 45 -4 -3 +3 +4 +5",
    "Pd 46 -4 -2 +2 +4",
    "Ag 47 -2 -1 +1 +2 +3",
    "Cd 48 -2 +2",
    "In 49 -3 +1 +3",
    "Sn 50 -4 -2 +2 +4",
    "Sb 51 -5 -3 +3 +5",
    "Te 52 -2 +4 +6",
    "I 53 -2 -1 +1 +3 +5 +7",
    "Xe 54 +2 +4 +6 +8",
    "Cs 55 -1 +1",
    "Ba 56 -2 +2",
    "La 57 -3 +3",
    "Ce 58 -4 -3 +3 +4",
    "Pr 59 -4 -3 +3",
    "Nd 60 -3 +2 +3",
    "Pm 61 -3 +3",
    "Sm 62 -3 +2 +3",
    "Eu 63 -3 -2 +2 +3",
    "Gd 64 -3 +2 +3",
    "Tb 65 -3 +3",
    "Dy 66 -3 +2 +3",
    "Ho 67 -3 +3",
    "Er 68 -3 +2 +3",
    "Tm 69 -3 +3",
    "Yb 70 -3 -2 +2 +3",
    "Lu 71 -3 +3",
    "Hf 72 -4 +3 +4",
    "Ta 73 -5 +4 +5",
    "W 74 -6 +4 +5 +6",
    "Re 75 +1 +3 +4 +5 +6 +7",
    "Os 76 -4 +4 +5 +6 +8",
    "Ir 77 -4 -3 +4 +5",
    "Pt 78 -4 -2 +2 +3 +4",
    "Au 79 -3 -1 +1 +3 +5",
    "Hg 80 -2 -1 +1 +2",
    "Tl 81 -3 -1 +1 +3",
    "Pb 82 -4 -2 +2 +4",
    "Bi 83 -5 -3 +3 +5",
    "Po 84 +4",
    "At 85",
    "Rn 86",
    "Fr 87",
    "Ra 88 -2",
    "Ac 89 -3 +3",
    "Th 90 -4 +4",
    "Pa 91 +4 +5",
    "U 92 -6 -4 -3 +2 +3 +4 +5 +6",
    "Np 93 -6 -4 -3 +3 +4 +5 +6 +7",
    "Pu 94 -6 -4 -3 +3 +4 +5 +6 +7",
    "Am 95 +3 +4 +5 +6",
    "Cm 96 +3 +4",
    "Bk 97 +3 +4",
    "Cf 98 +3 +4",
    NULL,
};

// Registration --------------------------------------------------------------

bool reg_SFTElectronNumber = SFTElectronNumber().registerThisType();

}   // namespace srreal
}   // namespace diffpy

// End of file
