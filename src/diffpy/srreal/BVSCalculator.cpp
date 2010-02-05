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
* class BVSCalculator -- concrete counter of pairs in a structure.
*
* $Id$
*
*****************************************************************************/

#include <cmath>
#include <cassert>

#include <diffpy/srreal/BVSCalculator.hpp>
#include <diffpy/srreal/StructureAdapter.hpp>
#include <diffpy/srreal/BaseBondGenerator.hpp>
#include <diffpy/srreal/BVParametersTable.hpp>

using namespace std;

namespace diffpy {
namespace srreal {

// Constructor ---------------------------------------------------------------

BVSCalculator::BVSCalculator()
{
    // default configuration
    const double valence_precision = 1e-5;
    BVParametersTable bvtb;
    this->setBVParamTable(bvtb);
    this->setValencePrecision(valence_precision);
}

// Public Methods ------------------------------------------------------------

// results

double BVSCalculator::bvmsdiff() const
{
    // FIXME
    return 0.0;
}


double BVSCalculator::bvrmsdiff() const
{
    double rvsq = this->bvmsdiff();
    return sqrt(rvsq);
}


void BVSCalculator::setBVParamTable(const BVParametersTable& bvpt)
{
    if (mbvptable.get() == &bvpt)   return;
    mbvptable.reset(bvpt.copy());
}


const BVParametersTable& BVSCalculator::getBVParamTable() const
{
    assert(mbvptable.get());
    return *mbvptable;
}


const BVParam& BVSCalculator::bvpar(int idx0, int idx1) const
{
    // FIXME
    static BVParam rv;
    return rv;
}


void BVSCalculator::setValencePrecision(double eps)
{
    BVParametersTable::SetOfBVParam allbp = this->getBVParamTable().getAll();
    double dmx = 0.0;
    BVParametersTable::SetOfBVParam::const_iterator bpit;
    for (bpit = allbp.begin(); bpit != allbp.end(); ++bpit)
    {
        dmx = max(dmx, bpit->bonddistance(eps));
    }
    this->setRmax(dmx);
}


double BVSCalculator::getValencePrecision() const
{
    BVParametersTable::SetOfBVParam allbp = this->getBVParamTable().getAll();
    double veps = 0.0;
    BVParametersTable::SetOfBVParam::const_iterator bpit;
    for (bpit = allbp.begin(); bpit != allbp.end(); ++bpit)
    {
        veps = max(veps, bpit->bondvalence(this->getRmax()));
    }
    return veps;
}

// Protected Methods ---------------------------------------------------------

void BVSCalculator::resetValue()
{ }


void BVSCalculator::configureBondGenerator(BaseBondGenerator& bnds)
{ }


void BVSCalculator::addPairContribution(const BaseBondGenerator& bnds)
{ }


}   // namespace srreal
}   // namespace diffpy

// End of file
