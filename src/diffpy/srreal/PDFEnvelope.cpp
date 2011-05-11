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
* class PDFEnvelope -- abstract base class for PDF envelope functions
*
* $Id$
*
*****************************************************************************/

#include <sstream>
#include <stdexcept>
#include <boost/serialization/export.hpp>

#include <diffpy/srreal/PDFEnvelope.hpp>
#include <diffpy/HasClassRegistry.ipp>

using namespace std;
using diffpy::srreal::PDFEnvelope;

// Unique instantiation of the template registry base class.
template class HasClassRegistry<PDFEnvelope>;

namespace diffpy {
namespace srreal {

// class PDFEnvelopeOwner ----------------------------------------------------

// public methods

QuantityType PDFEnvelopeOwner::applyEnvelopes(
        const QuantityType& x, const QuantityType& y) const
{
    assert(x.size() == y.size());
    QuantityType z = y;
    EnvelopeStorage::const_iterator evit;
    for (evit = menvelope.begin(); evit != menvelope.end(); ++evit)
    {
        PDFEnvelope& fenvelope = *(evit->second);
        QuantityType::const_iterator xi = x.begin();
        QuantityType::iterator zi = z.begin();
        for (; xi != x.end(); ++xi, ++zi)
        {
            *zi *= fenvelope(*xi);
        }
    }
    return z;
}


void PDFEnvelopeOwner::addEnvelope(PDFEnvelopePtr envlp)
{
    assert(envlp.get());
    menvelope[envlp->type()] = envlp;
}


void PDFEnvelopeOwner::addEnvelopeByType(const string& tp)
{
    // this throws invalid_argument for invalid type
    PDFEnvelopePtr envlp = PDFEnvelope::createByType(tp);
    // we get here only when createByType was successful
    menvelope[envlp->type()] = envlp;
}


void PDFEnvelopeOwner::popEnvelope(PDFEnvelopePtr envlp)
{
    EnvelopeStorage::iterator evit = menvelope.find(envlp->type());
    if (evit != menvelope.end() && evit->second.get() == envlp.get())
    {
        menvelope.erase(evit);
    }
}


void PDFEnvelopeOwner::popEnvelopeByType(const string& tp)
{
    menvelope.erase(tp);
}


const PDFEnvelopePtr& PDFEnvelopeOwner::getEnvelopeByType(const string& tp) const
{
    // call non-constant method
    const PDFEnvelopePtr& rv =
        const_cast<PDFEnvelopeOwner*>(this)->getEnvelopeByType(tp);
    return rv;
}


PDFEnvelopePtr& PDFEnvelopeOwner::getEnvelopeByType(const string& tp)
{
    if (!menvelope.count(tp))
    {
        ostringstream emsg;
        emsg << "Invalid or missing PDFEnvelope type '" << tp << "'.";
        throw invalid_argument(emsg.str());
    }
    PDFEnvelopePtr& rv = menvelope[tp];
    return rv;
}


set<string> PDFEnvelopeOwner::usedEnvelopeTypes() const
{
    set<string> rv;
    EnvelopeStorage::const_iterator evit;
    for (evit = menvelope.begin(); evit != menvelope.end(); ++evit)
    {
        rv.insert(rv.end(), evit->first);
    }
    return rv;
}


void PDFEnvelopeOwner::clearEnvelopes()
{
    menvelope.clear();
}

}   // namespace srreal
}   // namespace diffpy

// End of file
