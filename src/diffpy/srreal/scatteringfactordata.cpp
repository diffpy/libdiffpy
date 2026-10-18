/*****************************************************************************
*
* libdiffpy         Complex Modeling Initiative
*                   (c) 2013 Brookhaven Science Associates,
*                   Brookhaven National Laboratory.
*                   All rights reserved.
*
* File coded by:    Pavol Juhas
*
* See AUTHORS.txt for a list of people who contributed.
* See LICENSE.txt for license information.
*
******************************************************************************
*
* X-ray scattering factors for ions and neutral atoms obtained from
* f0_WaasKirf.dat and the associated reference
* D. Waasmaier, A. Kirfel, Acta Cryst. (1995). A51, 416-413
* http://dx.doi.org/10.1107/S0108767394013292
*
* Electron scattering factors approximated from the X-rays
# Reference: International Tables Volume C, page 224.
*
*****************************************************************************/

#include <algorithm>
#include <fstream>
#include <sstream>
#include <memory>
#include <unordered_set>
#include <unordered_map>

#include <diffpy/srreal/scatteringfactordata.hpp>
#include <diffpy/srreal/AtomUtils.hpp>
#include <diffpy/runtimepath.hpp>
#include <diffpy/validators.hpp>

using namespace std;

// Helper classes and functions ----------------------------------------------

namespace {

// X-ray scattering factors

const int WKTerms = 5;

class WaasKirfFormula
{
    public:

        // Constructor
        WaasKirfFormula()
        {
            fill(a, a + WKTerms, 0.0);
            fill(b, b + WKTerms, 0.0);
            c = 0.0;
        }

        WaasKirfFormula(const WaasKirfFormula& wk0)
        {
            symbol = wk0.symbol;
            copy(wk0.a, wk0.a + WKTerms, a);
            copy(wk0.b, wk0.b + WKTerms, b);
            c = wk0.c;
        }

        // methods
        double xrayatstol(double stol) const
        {
            double rv = c;
            for (int i = 0; i < WKTerms; ++i)
            {
                rv += a[i] * exp(-b[i] * stol * stol);
            }
            return rv;
        }

        double Z() const
        {
            double rv = c;
            for (int i = 0; i < WKTerms; ++i)  rv += a[i];
            return round(rv);
        }

        double electronatstol(double stol) const
        {
            using diffpy::mathutils::eps_eq;
            using diffpy::mathutils::DOUBLE_MAX;
            if (eps_eq(stol, 0.0))  return DOUBLE_MAX;
            double rv = 0.023934 *
                (this->Z() - this->xrayatstol(stol)) / (stol * stol);
            return rv;
        }

        // data
        string symbol;
        double a[WKTerms];
        double b[WKTerms];
        double c;
};


class wksmbl_equal
{
    public:
        bool operator() (
                const WaasKirfFormula& wk0, const WaasKirfFormula& wk1) const
        {
            return wk0.symbol == wk1.symbol;
        }
};


class wksmbl_hash
{
    public:
        size_t operator()(const WaasKirfFormula& wk) const
        {
            hash<string> hasher;
            return hasher(wk.symbol);
        }
};


typedef unordered_set<WaasKirfFormula,
        wksmbl_hash, wksmbl_equal> SetOfWKFormulas;


const SetOfWKFormulas& getWKFormulasSet()
{
    using namespace diffpy::runtimepath;
    using diffpy::validators::ensureFileOK;
    static unique_ptr<SetOfWKFormulas> the_set;
    if (the_set)  return *the_set;
    the_set.reset(new SetOfWKFormulas);
    string wkfile = datapath("f0_WaasKirf.dat");
    ifstream fp(wkfile.c_str());
    ensureFileOK(wkfile, fp);
    LineReader line;
    line.commentmark = '#';
    for (WaasKirfFormula wk; fp >> line;)
    {
        if (line.iscomment() && (0 == line.words[0].compare(0, 2, "#S")))
        {
            if (line.wcount() < 3)
            {
                throw line.format_error(wkfile,
                        "Expected at least 3 columns of data.");
            }
            wk.symbol = line.words[2];
            continue;
        }
        if (line.iscomment() && (0 == line.words[0].compare(0, 2, "#L")))
        {
            if (wk.symbol.empty())
            {
                throw line.format_error(wkfile,
                        "Missing \"#S\" line with atom symbol definition.");
            }
            fp >> line;
            if (line.wcount() != 11)
            {
                throw line.format_error(wkfile, "Expected 11 values.");
            }
            istringstream fpline(line.line);
            fpline >>
                wk.a[0] >> wk.a[1] >> wk.a[2] >> wk.a[3] >> wk.a[4] >>
                wk.c >> wk.b[0] >> wk.b[1] >> wk.b[2] >> wk.b[3] >> wk.b[4];
            if (!fpline)
            {
                throw line.format_error(wkfile,
                        "Line should contain 11 floating point values.");
            }
            if (the_set->count(wk))
            {
                string emsg = "Duplicate atom symbol \"";
                emsg += wk.symbol + "\".";
                throw line.format_error(wkfile, emsg);
            }
            the_set->insert(wk);
            wk.symbol.clear();
        }
    }
    return getWKFormulasSet();
}


SetOfWKFormulas::const_iterator findWKFormula(const string& smbl)
{
    using diffpy::srreal::atomBareSymbol;
    using diffpy::srreal::atomValence;
    WaasKirfFormula wk;
    wk.symbol = smbl;
    const SetOfWKFormulas& swk = getWKFormulasSet();
    SetOfWKFormulas::const_iterator wkit = swk.find(wk);
    if (wkit != swk.end())  return wkit;
    // try again with an adjusted charge suffix
    wk.symbol = atomBareSymbol(smbl);
    int v = atomValence(smbl);
    if (v == 0)  wkit = swk.find(wk);
    else if (abs(v) == 1)
    {
        wk.symbol += ((v > 0) ? "1+" : "1-");
        wkit = swk.find(wk);
    }
    // raise exception if no match
    if (wkit == swk.end())
    {
        string emsg("Unknown atom or ion symbol '");
        emsg += smbl + "'.";
        throw invalid_argument(emsg);
    }
    return wkit;
}


// Electron numbers for elements and ions

typedef unordered_map<string,int> ElectronNumberStorage;

ElectronNumberStorage& getElectronNumberTable()
{
    using namespace diffpy::runtimepath;
    using diffpy::validators::ensureFileOK;
    static unique_ptr<ElectronNumberStorage> entable;
    typedef ElectronNumberStorage::value_type ENPair;
    if (!entable)
    {
        entable.reset(new ElectronNumberStorage);
        string ionfile = datapath("ionlist.dat");
        ifstream fp0(ionfile.c_str());
        ensureFileOK(ionfile, fp0);
        LineReader line;
        while (fp0 >> line)
        {
            if (line.isignored())  continue;
            istringstream fpline(line.line);
            string element;
            int z = 0;
            fpline >> element >> z;
            if (!fpline)
            {
                throw line.format_error(ionfile,
                        "Expected at least 2 columns for (symbol, Z).");
            }
            entable->insert(ENPair(element, z));
            for (int v; fpline >> v;)
            {
                ostringstream smbl;
                smbl << element << abs(v) << ((v > 0) ? '+' : '-');
                entable->insert(ENPair(smbl.str(), z - v));
            }
        }
        const size_t mintablesize = 447;
        if (entable->size() < mintablesize)
        {
            ostringstream emsg;
            emsg << "Incomplete file.  Expected " << mintablesize <<
                " items loaded " << entable->size() << ".";
            throw line.format_error(ionfile, emsg.str());
        }
    }
    return *entable;
}

// Neutron scattering lengths

typedef unordered_map<string,double> NeutronBCStorage;

const NeutronBCStorage& getNeutronBCTable()
{
    using namespace diffpy::runtimepath;
    using diffpy::validators::ensureFileOK;
    typedef NeutronBCStorage::value_type BCPair;
    static unique_ptr<NeutronBCStorage> bctable;
    if (bctable)  return *bctable;
    bctable.reset(new NeutronBCStorage);
    string nsffile = datapath("nsftable.dat");
    ifstream fp(nsffile.c_str());
    ensureFileOK(nsffile, fp);
    LineReader line;
    line.commentmark = '#';
    line.separator = ',';
    while (fp >> line)
    {
        if (line.isignored())  continue;
        if (line.wcount() != 11)
        {
            throw line.format_error(nsffile,
                    "Expected 11 comma-separated items.");
        }
        if (line.words[3].empty())  continue;
        string smbl = line.words[0];
        size_t p0 = smbl.find_first_not_of("0123456789-");
        if (p0 == string::npos)
        {
            throw line.format_error(nsffile, "Missing or invalid atom symbol.");
        }
        smbl.erase(0, p0);
        size_t p1 = smbl.find_last_of('-');
        if (p1 != string::npos)
        {
            smbl = smbl.substr(p1 + 1) + "-" + smbl.substr(0, p1);
        }
        istringstream fpbc(line.words[3]);
        double bc;
        fpbc >> bc;
        if (!fpbc)
        {
            throw line.format_error(nsffile, "Invalid b_c value.");
        }
        if (bctable->count(smbl))
        {
            string emsg = "Duplicate atom symbol \"";
            emsg += smbl + "\".";
            throw line.format_error(nsffile, emsg);
        }
        bctable->insert(BCPair(smbl, bc));
        // elements are not explicitly included if there is just one isotope
        // or if all isotopes are unstable
        size_t p2 = smbl.find_first_of('-');
        if (p2 != string::npos)
        {
            string el = smbl.substr(p2 + 1);
            // there is just one isotope
            const string& chlf = line.words[1];
            bool addel = (chlf == "100") ||
                (!bctable->count(el) && *chlf.rbegin() == 'Y');
            if (addel)
            {
                if (bctable->count(el))
                {
                    string emsg = "Duplicate element entry for \"";
                    emsg += el + "\".";
                    throw line.format_error(nsffile, emsg);;
                }
                bctable->insert(BCPair(el, bc));
            }
        }
    }
    // define aliases for neutron, deuterium and tritium
    bctable->insert(BCPair("n", bctable->at("1-n")));
    bctable->insert(BCPair("D", bctable->at("2-H")));
    bctable->insert(BCPair("T", bctable->at("3-H")));
    return getNeutronBCTable();
}

}   // namespace

// Implementation ------------------------------------------------------------

namespace diffpy {
namespace srreal {

/// X-ray scattering factor of an element or ion a given Q
double fxrayatq(const string& smbl, double q)
{
    const double stol = q / (4 * M_PI);
    return fxrayatstol(smbl, stol);
}


/// X-ray scattering factor of an element or ion a given sin(theta)/lambda
double fxrayatstol(const string& smbl, double stol)
{
    SetOfWKFormulas::const_iterator wkit = findWKFormula(smbl);
    return wkit->xrayatstol(stol);
}


/// Electron scattering factor of an element or ion a given Q
double felectronatq(const string& smbl, double q)
{
    const double stol = q / (4 * M_PI);
    SetOfWKFormulas::const_iterator wkit = findWKFormula(smbl);
    return wkit->electronatstol(stol);
}


/// Number of electrons for an element or ion
int electronnumber(const string& smbl)
{
    const ElectronNumberStorage& entable = getElectronNumberTable();
    ElectronNumberStorage::const_iterator it;
    it = entable.find(smbl);
    // try to build standard symbol when not found
    if (it == entable.end())
    {
        ostringstream smbl1;
        smbl1 << atomBareSymbol(smbl);
        int v = atomValence(smbl);
        if (v)  smbl1 << abs(v) << ((v > 0) ? '+' : '-');
        it = entable.find(smbl1.str());
    }
    // throw exception if still not found
    if (it == entable.end())
    {
        ostringstream emsg;
        emsg << "Unknown atom symbol '" << smbl << "'.";
        throw invalid_argument(emsg.str());
    }
    return it->second;
}


/// Coherent scattering length of an element or isotope in fm
double bcneutron(const string& smbl)
{
    const NeutronBCStorage& bctable = getNeutronBCTable();
    NeutronBCStorage::const_iterator it = bctable.find(smbl);
    if (it != bctable.end())  return it->second;
    size_t pe = smbl.find_last_not_of("+-012345678 \t");
    string smblnocharge =
        (pe != string::npos) ? smbl.substr(0, pe + 1) : smbl;
    it = bctable.find(smblnocharge);
    if (it == bctable.end())
    {
        string emsg("Unknown atom or isotope symbol '");
        emsg += smbl + "'.";
        throw invalid_argument(emsg);
    }
    return it->second;
}

}   // namespace srreal
}   // namespace diffpy
