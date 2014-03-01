/*****************************************************************************
*
* libdiffpy         by DANSE Diffraction group
*                   Simon J. L. Billinge
*                   (c) 2010 The Trustees of Columbia University
*                   in the City of New York.  All rights reserved.
*
* File coded by:    Pavol Juhas
*
* See AUTHORS.txt for a list of people who contributed.
* See LICENSE_DANSE.txt for license information.
*
******************************************************************************
*
* class BVParametersTable -- table of bond valence sum parameters
*
*****************************************************************************/

#ifndef BVPARAMETERSTABLE_HPP_INCLUDED
#define BVPARAMETERSTABLE_HPP_INCLUDED

#include <string>
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/shared_ptr.hpp>
#include <diffpy/boostextensions/serialize_unordered_set.hpp>
#include <diffpy/srreal/BVParam.hpp>

namespace diffpy {
namespace srreal {

typedef boost::shared_ptr<class BVParametersTable> BVParametersTablePtr;

class BVParametersTable
{

    public:

        // types
        typedef boost::unordered_set<BVParam,
                boost::hash<BVParam>, BVParam::HashEqual> SetOfBVParam;

        // static methods
        static const BVParam& none();

        // methods
        const BVParam& lookup(const BVParam&) const;
        const BVParam& lookup(
                const std::string& smbl0, const std::string& smbl1) const;
        const BVParam& lookup(const std::string& atom0, int valence0,
                const std::string& atom1, int valence1) const;
        void setCustom(const BVParam&);
        void setCustom(const std::string& atom0, int valence0,
                const std::string& atom1, int valence1,
                double Ro, double b, std::string ref_id="");
        void resetCustom(const BVParam&);
        void resetCustom(const std::string& atom0, int valence0,
                const std::string& atom1, int valence1);
        void resetAll();
        SetOfBVParam getAll() const;

    private:

        // data
        SetOfBVParam mcustomtable;

        // methods
        const SetOfBVParam& getStandardSetOfBVParam() const;

        // serialization
        friend class boost::serialization::access;
        template<class Archive>
            void serialize(Archive& ar, const unsigned int version)
        {
            ar & mcustomtable;
        }

};  // class BVParametersTable

}   // namespace srreal
}   // namespace diffpy

#endif  // BVPARAMETERSTABLE_HPP_INCLUDED
