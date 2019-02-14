/*****************************************************************************
*
* libdiffpy         by DANSE Diffraction group
*                   Simon J. L. Billinge
*                   (c) 2008 The Trustees of Columbia University
*                   in the City of New York.  All rights reserved.
*
* File coded by:    Pavol Juhas
*
* See AUTHORS.txt for a list of people who contributed.
* See LICENSE_DANSE.txt for license information.
*
******************************************************************************
*
* class PairQuantity -- general implementation of pair quantity calculator
*
*****************************************************************************/

#include <algorithm>
#include <locale>
#include <sstream>

#include <diffpy/srreal/PairQuantity.hpp>
#include <diffpy/mathutils.hpp>
#include <diffpy/serialization.ipp>

using namespace std;

namespace diffpy {
namespace srreal {

// Class Constants -----------------------------------------------------------

const int PairQuantity::ALLATOMSINT = -1;
const string PairQuantity::ALLATOMSSTR = "all";

// Constructor ---------------------------------------------------------------

PairQuantity::PairQuantity() :
    mstructure(emptyStructureAdapter()),
    mrmin(0.0),
    mrmax(DEFAULT_BONDGENERATOR_RMAX),
    mdefaultpairmask(true)
{
    this->setEvaluatorType(BASIC);
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


const QuantityType& PairQuantity::eval(StructureAdapterPtr stru)
{
    mevaluator->updateValue(*this, stru);
    this->finishValue();
    return this->value();
}


void PairQuantity::setStructure(StructureAdapterPtr stru)
{
    mstructure = stru.get() ? stru : emptyStructureAdapter();
    mstructure->customPQConfig(this);
    this->updateMaskData();
    this->resetValue();
}


StructureAdapterPtr& PairQuantity::getStructure()
{
    return mstructure;
}


const StructureAdapterPtr& PairQuantity::getStructure() const
{
    return mstructure;
}


const QuantityType& PairQuantity::value() const
{
    return mvalue;
}


void PairQuantity::mergeParallelData(const string& pdata, int ncpu)
{
    if (mmergedvaluescount >= ncpu)
    {
        const char* emsg = "Number of merged values exceeds NCPU.";
        throw runtime_error(emsg);
    }
    this->executeParallelMerge(pdata);
    ++mmergedvaluescount;
    if (mmergedvaluescount == ncpu)  this->finishValue();
}


string PairQuantity::getParallelData() const
{
    ostringstream storage(ios::binary);
    diffpy::serialization::oarchive oa(storage, ios::binary);
    oa << this->value();
    return storage.str();
}


void PairQuantity::setRmin(double rmin)
{
    if (mrmin != rmin)  mticker.click();
    mrmin = rmin;
}


const double& PairQuantity::getRmin() const
{
    return mrmin;
}


void PairQuantity::setRmax(double rmax)
{
    if (mrmax != rmax)  mticker.click();
    mrmax = rmax;
}


const double& PairQuantity::getRmax() const
{
    return mrmax;
}


void PairQuantity::setEvaluatorType(PQEvaluatorType evtp)
{
    if (mevaluator.get() && mevaluator->typeint() == evtp)  return;
    PQEvaluatorPtr pqev = createPQEvaluator(evtp, mevaluator);
    // validate the new evaluator object
    try
    {
        pqev->validate(*this);
    }
    catch (logic_error e)
    {
        string emsg("EvaluatorType not supported.  ");
        emsg += e.what();
        throw invalid_argument(emsg);
    }
    // validator is good here, use it now.
    mevaluator = pqev;
    this->resetValue();
}


PQEvaluatorType PairQuantity::getEvaluatorType() const
{
    return mevaluator->typeint();
}


PQEvaluatorType PairQuantity::getEvaluatorTypeUsed() const
{
    return mevaluator->typeintused();
}


void PairQuantity::setupParallelRun(int cpuindex, int ncpu)
{
    mevaluator->setupParallelRun(cpuindex, ncpu);
}


void PairQuantity::maskAllPairs(bool mask)
{
    bool nochange = minvertpairmask.empty() && msiteallmask.empty() &&
        mtypemask.empty() && (mdefaultpairmask == mask);
    if (!nochange)  mticker.click();
    minvertpairmask.clear();
    msiteallmask.clear();
    mtypemask.clear();
    mdefaultpairmask = mask;
}


void PairQuantity::invertMask()
{
    mticker.click();
    mdefaultpairmask = !mdefaultpairmask;
    unordered_map<int, bool>::iterator mm;
    for (mm = msiteallmask.begin(); mm != msiteallmask.end(); ++mm)
    {
        mm->second = !(mm->second);
    }
    TypeMaskStorage::iterator tpmsk;
    for (tpmsk = mtypemask.begin(); tpmsk != mtypemask.end(); ++tpmsk)
    {
        tpmsk->second = !(tpmsk->second);
    }
}


void PairQuantity::setPairMask(int i, int j, bool mask)
{
    if (i < 0)  i = ALLATOMSINT;
    if (j < 0)  j = ALLATOMSINT;
    // short circuit for all-all
    if (ALLATOMSINT == i && ALLATOMSINT == j)
    {
        this->maskAllPairs(mask);
        return;
    }
    bool modified = false;
    // update ticker if we are switching from type-mask mode
    if (!mtypemask.empty())
    {
        mtypemask.clear();
        modified = true;
    }
    // handle one ALLATOMSINT argument
    if (ALLATOMSINT == i || ALLATOMSINT == j)
    {
        int k = (ALLATOMSINT != i) ? i : j;
        pair<unordered_map<int, bool>::iterator, bool> pmm;
        pmm = msiteallmask.emplace(k, mask);
        if (pmm.second || pmm.first->second != mask)  modified = true;
        pmm.first->second = mask;
        int cntsites = this->countSites();
        for (int l = 0; l < cntsites; ++l)
        {
            this->setPairMaskValue(k, l, mask);
        }
        if (modified)  mticker.click();
        return;
    }
    // here neither i nor j is ALLATOMSINT
    unordered_map<int, bool>::iterator ii, jj;
    ii = msiteallmask.find(i);
    if (ii != msiteallmask.end() && mask != ii->second)
    {
        msiteallmask.erase(ii);
        modified = true;
    }
    jj = msiteallmask.find(j);
    if (jj != msiteallmask.end() && mask != jj->second)
    {
        msiteallmask.erase(jj);
        modified = true;
    }
    modified = this->setPairMaskValue(i, j, mask) ? true : modified;
    if (modified)  mticker.click();
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
    // accept "ALL" (upper ALLATOMSSTR) for smbli and smblj
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
    // msiteallmask is only relevant for index-based masking.
    // clean it up to avoid unexpected effects when switching to index masks.
    msiteallmask.clear();
    bool modified = false;
    // when all is used, remove all typemask elements with the other type
    if (ALLATOMSSTR == smbli || ALLATOMSSTR == smblj)
    {
        const string& sk = (ALLATOMSSTR != smbli) ? smbli : smblj;
        TypeMaskStorage::iterator tpmsk;
        TypeMaskStorage::iterator skip = mtypemask.find(smblij);
        for (tpmsk = mtypemask.begin(); tpmsk != mtypemask.end();)
        {
            const bool eraseitem = (tpmsk != skip) &&
                (sk == tpmsk->first.first || sk == tpmsk->first.second);
            if (eraseitem)
            {
                if (mask != tpmsk->second)  modified = true;
                tpmsk = mtypemask.erase(tpmsk);
            }
            else  ++tpmsk;
        }
    }
    pair<TypeMaskStorage::iterator, bool> pmm;
    pmm = mtypemask.emplace(smblij, mask);
    if (pmm.second || pmm.first->second != mask)  modified = true;
    pmm.first->second = mask;
    if (modified)  mticker.click();
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


void PairQuantity::executeParallelMerge(const string& pdata)
{
    istringstream storage(pdata, ios::binary);
    diffpy::serialization::iarchive ia(storage, ios::binary);
    QuantityType pvalue;
    ia >> pvalue;
    if (pvalue.size() != mvalue.size())
    {
        throw invalid_argument("Merged data array must have the same size.");
    }
    transform(mvalue.begin(), mvalue.end(), pvalue.begin(),
            mvalue.begin(), plus<double>());
}


int PairQuantity::countSites() const
{
    int rv = mstructure.get() ? mstructure->countSites() : 0;
    return rv;
}


bool PairQuantity::hasMask() const
{
    bool rv = !(mdefaultpairmask && minvertpairmask.empty() &&
            msiteallmask.empty() && mtypemask.empty());
    return rv;
}


bool PairQuantity::hasPairMask() const
{
    bool rv = this->hasMask() && !this->hasTypeMask();
    return rv;
}


bool PairQuantity::hasTypeMask() const
{
    return !mtypemask.empty();
}


void PairQuantity::stashPartialValue()
{
    const char* emsg =
        "stashPartialValue() is not defined in the calculator class.";
    throw logic_error(emsg);
}


void PairQuantity::restorePartialValue()
{
    const char* emsg =
        "restorePartialValue() is not defined in the calculator class.";
    throw logic_error(emsg);
}

// Private Methods -----------------------------------------------------------

void PairQuantity::updateMaskData()
{
    int cntsites = this->countSites();
    // Propagate masks with ALLATOMSINT to all valid indices.
    if (mtypemask.empty())
    {
        unordered_map<int, bool>::const_iterator ia;
        for (ia = msiteallmask.begin(); ia != msiteallmask.end(); ++ia)
        {
            for (int j = 0; j < cntsites; ++j)
            {
                this->setPairMaskValue(ia->first, j, ia->second);
            }
        }
    }
    // For type masking propagate atom types to corresponding atom indices.
    else
    {
        // build a list of indices per each unique atom type
        unordered_map< string, list<int> >  siteindices;
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


bool PairQuantity::setPairMaskValue(int i, int j, bool mask)
{
    assert(i >= 0 && j >= 0);
    pair<int,int> ij = (i > j) ? make_pair(j, i) : make_pair(i, j);
    bool rv;
    rv = (mask == mdefaultpairmask) ?
        minvertpairmask.erase(ij) :
        minvertpairmask.insert(ij).second;
    return rv;
}

// Other functions -----------------------------------------------------------

/// The purpose of this function is to support Python pickling of
/// PairQuantity objects that hold Python-derived StructureAdapter classes.
/// Use it only if you absolutely have to and you know what you do.
StructureAdapterPtr
replacePairQuantityStructure(PairQuantity& pq, StructureAdapterPtr stru)
{
    StructureAdapterPtr rv = pq.mstructure;
    pq.mstructure = stru;
    return rv;
}

}   // namespace srreal
}   // namespace diffpy

// Serialization -------------------------------------------------------------

BOOST_CLASS_EXPORT_IMPLEMENT(diffpy::srreal::PairQuantity)
DIFFPY_INSTANTIATE_SERIALIZATION(diffpy::srreal::PairQuantity)
DIFFPY_INSTANTIATE_PTR_SERIALIZATION(diffpy::srreal::PairQuantity)

// End of file
