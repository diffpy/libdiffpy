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
* class BVParam -- bond valence parameters for a cation-anion pair
*
* $Id$
*
*****************************************************************************/

#ifndef BVPARAM_HPP_INCLUDED
#define BVPARAM_HPP_INCLUDED

#include <string>
#include <map>
#include <set>

namespace diffpy {
namespace srreal {

class BVParam
{
    public:

        // constructors
        BVParam();

        // methods
        /// obtain data from a cif record in bvparm.cif
        void setFromCifLine(const std::string&);
        /// Return bond valence at a specified distance
        double atd(double distance);

        // data
        std::string matom0;
        std::string matom1;
        int mvalence0;
        int mvalence1;
        double mRo;
        double mB;
        std::string mref_id;
};  // class BVParam


}   // namespace srreal
}   // namespace diffpy

#endif  // BVPARAM_HPP_INCLUDED
