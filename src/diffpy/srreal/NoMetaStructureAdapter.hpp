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
* class NoMetaStructureAdapter -- StructureAdapter proxy class that disables
*     customPQConfig of another StructureAdapter instance.  This may be used
*     for preventing the DiffPyStructureAdapter from applying scale,
*     spdiameter and such metadata to the PDFCalculator.
*
* nometa -- factory function that returns a NoMetaStructureAdapter
*     instance wrapped in StructureAdapterPtr
*
* $Id$
*
*****************************************************************************/

#ifndef NOMETASTRUCTUREADAPTER_HPP_INCLUDED
#define NOMETASTRUCTUREADAPTER_HPP_INCLUDED

#include <diffpy/srreal/StructureAdapter.hpp>

namespace diffpy {
namespace srreal {

class NoMetaStructureAdapter : public StructureAdapter
{
    public:

        // constructors
        NoMetaStructureAdapter(StructureAdapterPtr);

        // methods - overloaded
        virtual BaseBondGenerator* createBondGenerator() const;
        virtual int countSites() const;
        virtual double numberDensity() const;
        virtual const R3::Vector& siteCartesianPosition(int idx) const;
        virtual double siteOccupancy(int idx) const;
        virtual bool siteAnisotropy(int idx) const;
        virtual const R3::Matrix& siteCartesianUij(int idx) const;
        virtual const std::string& siteAtomType(int idx) const;
        virtual void customPQConfig(PairQuantity* pq) const;

    private:

        // data
        StructureAdapterPtr msrcstructure;
};

// Routines ------------------------------------------------------------------

/// create NoMetaStructureAdapter from an existing StructureAdapter
StructureAdapterPtr nometa(StructureAdapterPtr stru);


/// create NoMetaStructureAdapter from an adaptable structure object
template <class T>
StructureAdapterPtr nometa(const T& stru)
{
    StructureAdapterPtr bstru = createStructureAdapter(stru);
    return nometa(bstru);
}

}   // namespace srreal
}   // namespace diffpy

#endif  // NOMETASTRUCTUREADAPTER_HPP_INCLUDED
