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
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/map.hpp>
#include <boost/serialization/assume_abstract.hpp>
#include <boost/serialization/export.hpp>

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
        const PDFEnvelopePtr& getEnvelopeByType(const std::string& tp) const;
        PDFEnvelopePtr& getEnvelopeByType(const std::string& tp);
        std::set<std::string> usedEnvelopeTypes() const;
        void clearEnvelopes();

    private:

        // types
        typedef std::map<std::string, PDFEnvelopePtr> EnvelopeStorage;

        // data
        EnvelopeStorage menvelope;

        // serialization
        friend class boost::serialization::access;
        BOOST_SERIALIZATION_SPLIT_MEMBER()

        template<class Archive>
        void save(Archive & ar, const unsigned int version) const
        {
            using namespace std;
            using namespace diffpy::attributes;
            map<string, AttributesDataMap> thedata;
            EnvelopeStorage::const_iterator evit;
            for (evit = menvelope.begin(); evit != menvelope.end(); ++evit)
            {
                thedata[evit->first] = saveAttributesData(*(evit->second));
            }
            ar & thedata;
        }

        template<class Archive>
        void load(Archive & ar, const unsigned int version)
        {
            using namespace std;
            using namespace diffpy::attributes;
            map<string, AttributesDataMap> thedata;
            ar & thedata;
            map<string, AttributesDataMap>::const_iterator dt;
            this->clearEnvelopes();
            for (dt = thedata.begin(); dt != thedata.end(); ++dt)
            {
                this->addEnvelopeByType(dt->first);
                PDFEnvelopePtr e = this->getEnvelopeByType(dt->first);
                loadAttributesData(*e, dt->second);
            }
        }

};

}   // namespace srreal
}   // namespace diffpy

// Serialization -------------------------------------------------------------

BOOST_SERIALIZATION_ASSUME_ABSTRACT(diffpy::srreal::PDFEnvelopeOwner)

#endif  // PDFENVELOPE_HPP_INCLUDED
