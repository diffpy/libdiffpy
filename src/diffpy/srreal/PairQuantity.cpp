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

#include <diffpy/srreal/PairQuantity.hpp>
#include <diffpy/mathutils.hpp>

using namespace std;

namespace diffpy {
namespace srreal {

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

void PairQuantity::setStructure(StructureAdapterPtr stru)
{
    mstructure = stru;
    mstructure->customPQConfig(this);
    this->updateMaskData();
    this->resetValue();
}


const QuantityType& PairQuantity::value() const
{
    return mvalue;
}


void PairQuantity::mergeParallelValue(const QuantityType& pvalue)
{
    if (pvalue.size() != mvalue.size())
    {
        throw invalid_argument("Merged value array must have the same size.");
    }
    transform(mvalue.begin(), mvalue.end(), pvalue.begin(),
            mvalue.begin(), plus<double>());
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
setTypeMask(const string& smbli, const string& smblj, bool mask)
{
    pair<string,string> smblij = (smbli > smblj) ?
        make_pair(smblj, smbli) : make_pair(smbli, smblj);
    mtypemask[smblij] = mask;
}


bool PairQuantity::getTypeMask(const string& smbli, const string& smblj) const
{
    pair<string,string> smblij = (smbli > smblj) ?
        make_pair(smblj, smbli) : make_pair(smbli, smblj);
    bool rv = mtypemask.count(smblij) ?
        mtypemask.at(smblij) : mdefaultpairmask;
    return rv;
}

// Protected Methods ---------------------------------------------------------

void PairQuantity::resizeValue(size_t sz)
{
    mvalue.resize(sz);
}


void PairQuantity::resetValue()
{
    fill(mvalue.begin(), mvalue.end(), 0.0);
}


void PairQuantity::configureBondGenerator(BaseBondGenerator& bnds) const
{
    bnds.setRmin(this->getRmin());
    bnds.setRmax(this->getRmax());
}


int PairQuantity::countSites() const
{
    int rv = mstructure.get() ? mstructure->countSites() : 0;
    return rv;
}

// Private Methods -----------------------------------------------------------

void PairQuantity::updateMaskData()
{
    if (mtypemask.empty())  return;
    // build a list of indices per each unique atom type
    boost::unordered_map< string, list<int> >  siteindices;
    int cntsites = this->countSites();
    for (int i = 0; i < cntsites; ++i)
    {
        const string& smbl = mstructure->siteAtomType(i);
        siteindices[smbl].push_back(i);
    }
    // rebuild minvertpairmask according to mtypemask
    minvertpairmask.clear();
    TypeMaskStorage::const_iterator tpmsk;
    for (tpmsk = mtypemask.begin(); tpmsk != mtypemask.end(); ++tpmsk)
    {
        const list<int>& isites = siteindices[tpmsk->first.first];
        const list<int>& jsites = siteindices[tpmsk->first.second];
        list<int>::const_iterator ii, jj;
        for (ii = isites.begin(); ii != isites.end(); ++ii)
        {
            jj = (&isites == &jsites) ? ii : jsites.begin();
            for (; jj != jsites.end(); ++jj)
            {
                this->setPairMaskValue(*ii, *jj, tpmsk->second);
            }
        }
    }
}


void PairQuantity::setPairMaskValue(int i, int j, bool mask)
{
    pair<int,int> ij = (i > j) ? make_pair(j, i) : make_pair(i, j);
    if (mask == mdefaultpairmask)  minvertpairmask.erase(ij);
    else    minvertpairmask.insert(ij);
}

}   // namespace srreal
}   // namespace diffpy

// Serialization -------------------------------------------------------------

BOOST_CLASS_EXPORT(diffpy::srreal::PairQuantity)

// End of file
