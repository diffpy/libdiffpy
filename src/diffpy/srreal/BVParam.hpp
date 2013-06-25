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
*****************************************************************************/

#ifndef BVPARAM_HPP_INCLUDED
#define BVPARAM_HPP_INCLUDED

#include <string>
#include <boost/serialization/base_object.hpp>

namespace diffpy {
namespace srreal {

class BVParam
{
    public:

        // constructors
        BVParam();
        BVParam(const std::string& atom0, int valence0,
                const std::string& atom1, int valence1,
                double Ro=0.0, double B=0.0, std::string ref_id="");

        // methods
        /// comparison binary_function
        bool operator==(const BVParam& other) const;
        /// Return bond valence at a specified distance
        double bondvalence(double distance) const;
        /// Return distance corresponding  to a specified bond valence
        double bondvalenceToDistance(double bvalence) const;
        /// obtain data from a cif record in bvparm.cif
        void setFromCifLine(const std::string&);

        // data
        std::string matom0;
        int mvalence0;
        std::string matom1;
        int mvalence1;
        double mRo;
        double mB;
        std::string mref_id;

        // comparison binary_function for unordered_set hasher
        class HashEqual : public std::binary_function<BVParam,BVParam,bool>
        {
            public:
                bool operator()(const BVParam& bp0, const BVParam& bp1) const;
        };

    private:

        // serialization
        friend class boost::serialization::access;
        template<class Archive>
            void serialize(Archive& ar, const unsigned int version)
        {
            ar & matom0;
            ar & mvalence0;
            ar & matom1;
            ar & mvalence1;
            ar & mRo;
            ar & mB;
            ar & mref_id;
        }

};  // class BVParam

// Functions -----------------------------------------------------------------

size_t hash_value(const BVParam& bp);

}   // namespace srreal
}   // namespace diffpy

#endif  // BVPARAM_HPP_INCLUDED
