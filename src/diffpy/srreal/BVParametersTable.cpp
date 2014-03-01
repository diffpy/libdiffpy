/*****************************************************************************
*
* libdiffpy         by DANSE Diffraction group
*                   Simon J. L. Billinge
*                   (c) 2010 The Trustees of Columbia University
*                   in the City of New York.  All rights reserved.
*
* File coded by:    Pavol Juhas
*
* See AUTHORS.txt for a list of people who contributed.
* See LICENSE.txt for license information.
*
******************************************************************************
*
* class BVParametersTable -- table of bond valence sum parameters
*
*****************************************************************************/

#include <cassert>
#include <fstream>
#include <boost/serialization/export.hpp>
#include <boost/scoped_ptr.hpp>

#include <diffpy/serialization.ipp>
#include <diffpy/runtimepath.hpp>
#include <diffpy/validators.hpp>
#include <diffpy/srreal/AtomUtils.hpp>
#include <diffpy/srreal/BVParametersTable.hpp>

using namespace std;

namespace diffpy {
namespace srreal {

// Static Methods ------------------------------------------------------------

const BVParam& BVParametersTable::none()
{
    static const BVParam bpnone;
    return bpnone;
}

// Public Methods ------------------------------------------------------------

const BVParam& BVParametersTable::lookup(const BVParam& bpk) const
{
    SetOfBVParam::const_iterator bpit;
    bpit = mcustomtable.find(bpk);
    if (bpit != mcustomtable.end())  return *bpit;
    const SetOfBVParam& stdtable = this->getStandardSetOfBVParam();
    bpit = stdtable.find(bpk);
    if (bpit != stdtable.end())  return *bpit;
    // unspecified cations are marked with valence 9
    // perform this search only if the first atom is indeed a cation
    if (!(bpk.mvalence0 > 0))  return this->none();
    BVParam bpk9 = bpk;
    bpk9.mvalence0 = 9;
    bpit = mcustomtable.find(bpk9);
    if (bpit != mcustomtable.end())  return *bpit;
    bpit = stdtable.find(bpk9);
    if (bpit != stdtable.end())  return *bpit;
    // not found - return blank BVParam
    return this->none();
}


const BVParam&
BVParametersTable::lookup(const string& smbl0, const string& smbl1) const
{
    const BVParam& rv = this->lookup(
            atomBareSymbol(smbl0), atomValence(smbl0),
            atomBareSymbol(smbl1), atomValence(smbl1));
    return rv;
}


const BVParam& BVParametersTable::lookup(const string& atom0, int valence0,
        const string& atom1, int valence1) const
{
    BVParam bpk(atom0, valence0, atom1, valence1);
    const BVParam& rv = this->lookup(bpk);
    return rv;
}


void BVParametersTable::setCustom(const BVParam& bp)
{
    this->resetCustom(bp);
    mcustomtable.insert(bp);
}


void BVParametersTable::setCustom(const string& atom0, int valence0,
        const string& atom1, int valence1,
        double Ro, double b, std::string ref_id)
{
    BVParam bp(atom0, valence0, atom1, valence1, Ro, b, ref_id);
    this->setCustom(bp);
}


void BVParametersTable::resetCustom(const BVParam& bp)
{
    mcustomtable.erase(bp);
}


void BVParametersTable::resetCustom(const string& atom0, int valence0,
        const string& atom1, int valence1)
{
    BVParam bp(atom0, valence0, atom1, valence1);
    this->resetCustom(bp);
}


void BVParametersTable::resetAll()
{
    mcustomtable.clear();
}


BVParametersTable::SetOfBVParam
BVParametersTable::getAll() const
{
    // insert inserts only those items that are not yet in the set, therefore
    // we start with the custom table and then add the standard entries.
    SetOfBVParam rv = mcustomtable;
    const SetOfBVParam& stdtable = this->getStandardSetOfBVParam();
    rv.insert(stdtable.begin(), stdtable.end());
    return rv;
}

// Private Methods -----------------------------------------------------------

const BVParametersTable::SetOfBVParam&
BVParametersTable::getStandardSetOfBVParam() const
{
    using namespace diffpy::runtimepath;
    using diffpy::validators::ensureFileOK;
    static boost::scoped_ptr<SetOfBVParam> the_set;
    if (!the_set)
    {
        the_set.reset(new SetOfBVParam);
        string bvparmfile = datapath("bvparm2011sel.cif");
        ifstream fp(bvparmfile.c_str());
        ensureFileOK(bvparmfile, fp);
        // read the header up to _valence_param_B and then up to an empty line.
        LineReader lnrd;
        lnrd.commentmark = '#';
        while (fp >> lnrd)
        {
            if (lnrd.wcount() && lnrd.words[0] == "_valence_param_B")  break;
        }
        // skip to an empty line
        while (fp >> lnrd && !lnrd.isblank())  { }
        // load data lines skipping the empty or commented entries
        while (fp >> lnrd)
        {
            if (lnrd.isignored())  continue;
            BVParam bp;
            bp.setFromCifLine(lnrd.line);
            assert(!the_set->count(bp));
            the_set->insert(bp);
        }
    }
    return *the_set;
}

}   // namespace srreal
}   // namespace diffpy

// Serialization -------------------------------------------------------------

DIFFPY_INSTANTIATE_SERIALIZATION(diffpy::srreal::BVParametersTable)

// End of file
