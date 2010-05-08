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
*     A concrete instance of PDFEnvelope is a functor, that calculates
*     PDF scaling coefficients at a given pair distance r.  Several functors
*     can be defined and applied in PDFCalculator.
*
* $Id$
*
*****************************************************************************/

#ifndef PDFENVELOPE_HPP_INCLUDED
#define PDFENVELOPE_HPP_INCLUDED

#include <string>
#include <set>

#include <diffpy/Attributes.hpp>
#include <diffpy/HasClassRegistry.hpp>
#include <diffpy/srreal/PairQuantity.hpp>

namespace diffpy {
namespace srreal {


/// @class PDFEnvelope
/// @brief abstract base class for PDF envelope scaling function

class PDFEnvelope :
    public diffpy::Attributes,
    public diffpy::HasClassRegistry<PDFEnvelope>
{
    public:
        virtual double operator()(const double& r) const = 0;
};


typedef PDFEnvelope::SharedPtr PDFEnvelopePtr;


/// @class PDFEnvelopeOwner
/// @brief storage of one or more PDFEnvelope functions and their
/// application to unscaled x, y arrays

class PDFEnvelopeOwner
{
    public:

        // application on (x, y) data
        QuantityType applyEnvelopes(const QuantityType& x, const QuantityType& y) const;

        // access and configuration of PDF envelope functions
        // configuration of envelopes
        void addEnvelope(PDFEnvelopePtr);
        void addEnvelopeByType(const std::string& tp);
        void popEnvelope(PDFEnvelopePtr);
        void popEnvelopeByType(const std::string& tp);
        const PDFEnvelopePtr getEnvelopeByType(const std::string& tp) const;
        PDFEnvelopePtr getEnvelopeByType(const std::string& tp);
        std::set<std::string> usedEnvelopeTypes() const;
        void clearEnvelopes();

    private:

        // types
        typedef std::map<std::string, PDFEnvelopePtr> EnvelopeStorage;

        // data
        EnvelopeStorage menvelope;
};

}   // namespace srreal
}   // namespace diffpy

#endif  // PDFENVELOPE_HPP_INCLUDED
