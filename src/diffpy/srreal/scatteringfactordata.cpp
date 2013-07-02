/*****************************************************************************
*
* diffpy.srreal     Complex Modeling Initiative
*                   Pavol Juhas
*                   (c) 2013 Brookhaven National Laboratory,
*                   Upton, New York.  All rights reserved.
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
#include <boost/unordered_set.hpp>
#include <boost/scoped_ptr.hpp>
#include <blitz/tinyvec-et.h>

#include <diffpy/srreal/scatteringfactordata.hpp>
#include <diffpy/srreal/StructureAdapter.hpp>
#include <diffpy/runtimepath.hpp>
#include <diffpy/validators.hpp>

using namespace std;

// Helper classes and functions ----------------------------------------------

namespace {

const int WKTerms = 5;

class WaasKirfFormula
{
    public:

        // Constructor
        WaasKirfFormula()
        {
            a = 0.0;
            b = 0.0;
            c = 0.0;
        }

        // methods
        double atstol(double stol) const
        {
            double rv = c + blitz::sum(a * exp(-b * stol * stol));
            return rv;
        }

        // data
        string symbol;
        blitz::TinyVector<double,WKTerms> a;
        blitz::TinyVector<double,WKTerms> b;
        double c;
};


class wksmbl_equal :
    public binary_function<WaasKirfFormula,WaasKirfFormula,bool>
{
    public:
        bool operator() (
                const WaasKirfFormula& wk0, const WaasKirfFormula& wk1) const
        {
            return wk0.symbol == wk1.symbol;
        }
};


size_t hash_value(const WaasKirfFormula& wk)
{
    boost::hash<string> hasher;
    return hasher(wk.symbol);
}


typedef boost::unordered_set<WaasKirfFormula,
        boost::hash<WaasKirfFormula>, wksmbl_equal> SetOfWKFormulas;


const SetOfWKFormulas& getWKFormulasSet()
{
    using namespace diffpy::runtimepath;
    using diffpy::validators::ensureFileOK;
    static boost::scoped_ptr<SetOfWKFormulas> the_set;
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
            assert(line.wcount() >= 3);
            wk.symbol = line.words[2];
            continue;
        }
        if (line.iscomment() && (0 == line.words[0].compare(0, 2, "#L")))
        {
            assert(!wk.symbol.empty());
            assert(fp >> line);
            assert(11 == line.wcount());
            line.linestream.clear();
            line.linestream.seekg(0);
            line.linestream >>
                wk.a[0] >> wk.a[1] >> wk.a[2] >> wk.a[3] >> wk.a[4] >>
                wk.c >> wk.b[0] >> wk.b[1] >> wk.b[2] >> wk.b[3] >> wk.b[4];
            assert(line.linestream);
            assert(!the_set->count(wk));
            the_set->insert(wk);
            wk.symbol.clear();
        }
    }
    return getWKFormulasSet();
}

}   // namespace

// Implementation ------------------------------------------------------------

namespace diffpy {
namespace srreal {

/// X-ray scattering factor of an element or ion a given Q
double fxrayatq(const string& smbl, double q)
{
    double stol = q / (4 * M_PI);
    return fxrayatstol(smbl, stol);
}


/// X-ray scattering factor of an element or ion a given sin(theta)/lambda
double fxrayatstol(const string& smbl, double stol)
{
    WaasKirfFormula wk;
    wk.symbol = smbl;
    const SetOfWKFormulas& swk = getWKFormulasSet();
    SetOfWKFormulas::const_iterator wkit = swk.find(wk);
    if (wkit == swk.end())
    {
        wk.symbol = atomBareSymbol(smbl);
        int v = atomValence(smbl);
        if (v == 0)  wkit = swk.find(wk);
        if (abs(v) == 1)
        {
            wk.symbol += ((v > 0) ? "1+" : "1-");
            wkit = swk.find(wk);
        }
    }
    // raise exception if no match
    if (wkit == swk.end())
    {
        string emsg("Unknown atom or ion symbol '");
        emsg += smbl + "'.";
        throw invalid_argument(emsg);
    }
    double rv = wkit->atstol(stol);
    return rv;
}


/// Electron scattering factor of an element or ion a given Q
double felectronatq(const string& smbl, double q)
{
    using namespace diffpy::mathutils;
    // resolve Z first so that invalid symbol can throw an exception
    double Z = round(fxrayatstol(smbl, 0.0));
    if (eps_eq(q, 0.0))     return DOUBLE_MAX;
    double stol = q / (4 * M_PI);
    double rv = 0.023934 * (Z - fxrayatstol(smbl, stol)) / (stol * stol);
    return rv;
}

}   // namespace srreal
}   // namespace diffpy
