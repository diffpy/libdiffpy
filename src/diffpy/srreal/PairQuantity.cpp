/*****************************************************************************
*
* diffpy.srreal     by DANSE Diffraction group
*                   Simon J. L. Billinge
*                   (c) 2008 Trustees of the Columbia University
*                   in the City of New York.  All rights reserved.
*
* File coded by:    Pavol Juhas
*
* See AUTHORS.txt for a list of people who contributed.
* See LICENSE.txt for license information.
*
******************************************************************************
*
* class PairQuantity -- general implementation of pair quantity calculator
*
* $Id$
*
*****************************************************************************/

#include <algorithm>
#include <locale>
#include <diffpy/srreal/PairQuantity.hpp>
#include <diffpy/mathutils.hpp>

using namespace std;

namespace diffpy {
namespace srreal {

// Class Constants -----------------------------------------------------------

const int PairQuantity::ALLATOMSINT = -1;
const string PairQuantity::ALLATOMSSTR = "all";

// Constructor ---------------------------------------------------------------

PairQuantity::PairQuantity()
{
    this->setRmin(0.0);
    this->setRmax(DEFAULT_BONDGENERATOR_RMAX);
    this->setEvaluator(BASIC);
    this->maskAllPairs(true);
    // attributes
    this->registerDoubleAttribute("rmin", this,
            &PairQuantity::getRmin, &PairQuantity::setRmin);
    this->registerDoubleAttribute("rmax", this,
            &PairQuantity::getRmax, &PairQuantity::setRmax);
}

// Public Methods ------------------------------------------------------------

const QuantityType& PairQuantity::eval()
{
    return this->eval(mstructure);
}


void PairQuantity::setStructure(StructureAdapterPtr stru)
{
    mstructure = stru;
    mstructure->customPQConfig(this);
    this->updateMaskData();
    this->resetValue();
}


const StructureAdapterPtr& PairQuantity::getStructure() const
{
    return mstructure;
}


const QuantityType& PairQuantity::value() const
{
    return mvalue;
}


void PairQuantity::mergeParallelValue(const QuantityType& pvalue, int ncpu)
{
    if (mmergedvaluescount >= ncpu)
    {
        const char* emsg = "Number of merged values exceeds NCPU.";
        throw runtime_error(emsg);
    }
    this->executeParallelMerge(pvalue);
    ++mmergedvaluescount;
    if (mmergedvaluescount == ncpu)  this->finishValue();
}


void PairQuantity::setRmin(double rmin)
{
    mrmin = rmin;
}


const double& PairQuantity::getRmin() const
{
    return mrmin;
}


void PairQuantity::setRmax(double rmax)
{
    mrmax = rmax;
}


const double& PairQuantity::getRmax() const
{
    return mrmax;
}


void PairQuantity::setEvaluator(PQEvaluatorType evtp)
{
    if (mevaluator.get() && mevaluator->typeint() == evtp)  return;
    mevaluator = createPQEvaluator(evtp);
    this->resetValue();
}


void PairQuantity::setupParallelRun(int cpuindex, int ncpu)
{
    mevaluator->setupParallelRun(cpuindex, ncpu);
}


void PairQuantity::maskAllPairs(bool mask)
{
    minvertpairmask.clear();
    mtypemask.clear();
    mdefaultpairmask = mask;
}


void PairQuantity::invertMask()
{
    mdefaultpairmask = !mdefaultpairmask;
    TypeMaskStorage::iterator tpmsk;
    for (tpmsk = mtypemask.begin(); tpmsk != mtypemask.end(); ++tpmsk)
    {
        tpmsk->second = !(tpmsk->second);
    }
}


void PairQuantity::setPairMask(int i, int j, bool mask)
{
    mtypemask.clear();
    if (i < 0)  i = ALLATOMSINT;
    if (j < 0)  j = ALLATOMSINT;
    // short circuit for all-all
    if (ALLATOMSINT == i && ALLATOMSINT == j)
    {
        this->maskAllPairs(mask);
        return;
    }
    this->setPairMaskValue(i, j, mask);
}


bool PairQuantity::getPairMask(int i, int j) const
{
    pair<int,int> ij = (i > j) ? make_pair(j, i) : make_pair(i, j);
    bool rv = minvertpairmask.count(ij) ?
        !mdefaultpairmask : mdefaultpairmask;
    return rv;
}


void PairQuantity::
setTypeMask(string smbli, string smblj, bool mask)
{
    static string upcaseall;
    if (upcaseall.empty())
    {
        upcaseall = ALLATOMSSTR;
        transform(upcaseall.begin(), upcaseall.end(),
                upcaseall.begin(), ::toupper);
    }
    if (upcaseall == smbli)  smbli = ALLATOMSSTR;
    if (upcaseall == smblj)  smblj = ALLATOMSSTR;
    pair<string,string> allall(ALLATOMSSTR, ALLATOMSSTR);
    pair<string,string> smblij = (smbli > smblj) ?
        make_pair(smblj, smbli) : make_pair(smbli, smblj);
    // short circuit for all-all
    if (allall == smblij)
    {
        this->maskAllPairs(mask);
        return;
    }
    // when all is used, remove all typemask elements with the other type
    if (ALLATOMSSTR == smbli || ALLATOMSSTR == smblj)
    {
        const string& sk = (ALLATOMSSTR != smbli) ? smbli : smblj;
        TypeMaskStorage::iterator tpmsk;
        for (tpmsk = mtypemask.begin(); tpmsk != mtypemask.end();)
        {
            tpmsk = (sk == tpmsk->first.first || sk == tpmsk->first.second) ?
                mtypemask.erase(tpmsk) : ++tpmsk;
        }
    }
    mtypemask[smblij] = mask;
}


bool PairQuantity::getTypeMask(const string& smbli, const string& smblj) const
{
    pair<string,string> smblij = (smbli > smblj) ?
        make_pair(smblj, smbli) : make_pair(smbli, smblj);
    pair<string,string> alli = (smbli > ALLATOMSSTR) ?
        make_pair(ALLATOMSSTR, smbli) : make_pair(smbli, ALLATOMSSTR);
    pair<string,string> allj = (smblj > ALLATOMSSTR) ?
        make_pair(ALLATOMSSTR, smblj) : make_pair(smblj, ALLATOMSSTR);
    bool rv = mtypemask.count(smblij) ? mtypemask.at(smblij) :
        mtypemask.count(alli) ? mtypemask.at(alli) :
        mtypemask.count(allj) ? mtypemask.at(allj) :
        mdefaultpairmask;
    return rv;
}

// Protected Methods ---------------------------------------------------------

void PairQuantity::resizeValue(size_t sz)
{
    mvalue.resize(sz);
}


void PairQuantity::resetValue()
{
    mmergedvaluescount = 0;
    fill(mvalue.begin(), mvalue.end(), 0.0);
}


void PairQuantity::configureBondGenerator(BaseBondGenerator& bnds) const
{
    bnds.setRmin(this->getRmin());
    bnds.setRmax(this->getRmax());
}


void PairQuantity::executeParallelMerge(const QuantityType& pvalue)
{
    if (pvalue.size() != mvalue.size())
    {
        throw invalid_argument("Merged value array must have the same size.");
    }
    transform(mvalue.begin(), mvalue.end(), pvalue.begin(),
            mvalue.begin(), plus<double>());
}


int PairQuantity::countSites() const
{
    int rv = mstructure.get() ? mstructure->countSites() : 0;
    return rv;
}

// Private Methods -----------------------------------------------------------

void PairQuantity::updateMaskData()
{
    int cntsites = this->countSites();
    // Propagate masks with ALLATOMSINT to all valid indices.
    if (mtypemask.empty())
    {
        for (int i = 0; i < cntsites; ++i)
        {
            bool maskall = minvertpairmask.count(make_pair(ALLATOMSINT, i));
            for (int j = 0; maskall && j < cntsites; ++j)
            {
                this->setPairMaskValue(i, j, !mdefaultpairmask);
            }
        }
    }
    // For type masking propagate atom types to corresponding atom indices.
    else
    {
        // build a list of indices per each unique atom type
        boost::unordered_map< string, list<int> >  siteindices;
        for (int i = 0; i < cntsites; ++i)
        {
            const string& smbl = mstructure->siteAtomType(i);
            siteindices[smbl].push_back(i);
            siteindices[ALLATOMSSTR].push_back(i);
        }
        // rebuild minvertpairmask according to mtypemask
        minvertpairmask.clear();
        TypeMaskStorage::const_iterator tpmsk;
        // build a list of type masks with all-masks at the begining
        list< pair<string,string> > orderedpairs;
        for (tpmsk = mtypemask.begin(); tpmsk != mtypemask.end(); ++tpmsk)
        {
            bool hasall = (ALLATOMSSTR == tpmsk->first.first ||
                    ALLATOMSSTR == tpmsk->first.second);
            if (hasall)  orderedpairs.push_front(tpmsk->first);
            else  orderedpairs.push_back(tpmsk->first);
        }
        list< pair<string,string> >::const_iterator tpp;
        for (tpp = orderedpairs.begin(); tpp != orderedpairs.end(); ++tpp)
        {
            const list<int>& isites = siteindices[tpp->first];
            const list<int>& jsites = siteindices[tpp->second];
            bool msk = mtypemask.at(*tpp);
            list<int>::const_iterator ii, jj;
            for (ii = isites.begin(); ii != isites.end(); ++ii)
            {
                jj = (&isites == &jsites) ? ii : jsites.begin();
                for (; jj != jsites.end(); ++jj)
                {
                    this->setPairMaskValue(*ii, *jj, msk);
                }
            }
        }
    }
}


void PairQuantity::setPairMaskValue(int i, int j, bool mask)
{
    assert(i != ALLATOMSINT || j != ALLATOMSINT);
    pair<int,int> ij = (i > j) ? make_pair(j, i) : make_pair(i, j);
    if (ALLATOMSINT == ij.first && mask == mdefaultpairmask)
    {
        // erase any inverted masks for the second site
        int k = ij.second;
        boost::unordered_set< pair<int,int> >::iterator ii;
        for (ii = minvertpairmask.begin(); ii != minvertpairmask.end();)
        {
            ii = (ii->first == k || ii->second == k) ?
                minvertpairmask.erase(ii) : ++ii;
        }
    }
    if (mask == mdefaultpairmask)  minvertpairmask.erase(ij);
    else    minvertpairmask.insert(ij);
}

}   // namespace srreal
}   // namespace diffpy

// Serialization -------------------------------------------------------------

BOOST_CLASS_EXPORT(diffpy::srreal::PairQuantity)

// End of file
