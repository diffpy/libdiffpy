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
* Shared types and small functions used by PairQuantity calculators.
*
* QuantityType -- type that stores PairQuantity results, an array of doubles
*
* $Id$
*
*****************************************************************************/

#ifndef PAIRQUANTITYUTILS_HPP_INCLUDED
#define PAIRQUANTITYUTILS_HPP_INCLUDED

#include <vector>

namespace diffpy {
namespace srreal {

typedef std::vector<double> QuantityType;

}   // namespace srreal
}   // namespace diffpy

#endif  // PAIRQUANTITYUTILS_HPP_INCLUDED
