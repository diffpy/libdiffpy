/*****************************************************************************
*
* libdiffpy         by DANSE Diffraction group
*                   Simon J. L. Billinge
*                   (c) 2009 The Trustees of Columbia University
*                   in the City of New York.  All rights reserved.
*
* File coded by:    Pavol Juhas
*
* See AUTHORS.txt for a list of people who contributed.
* See LICENSE_DANSE.txt for license information.
*
******************************************************************************
*
* class BaseBondGenerator -- semi-abstract class for a generation
*     of all atom pairs containing specified anchor atom.
*
*****************************************************************************/

#ifndef BASEBONDGENERATOR_HPP_INCLUDED
#define BASEBONDGENERATOR_HPP_INCLUDED

#include <diffpy/srreal/R3linalg.hpp>
#include <diffpy/srreal/forwardtypes.hpp>

namespace diffpy {
namespace srreal {

/// Default rmax is very large, but still within the integer type limits.
const double DEFAULT_BONDGENERATOR_RMAX = 1.0e6;

class BaseBondGenerator
{
    public:

        // constructor
        BaseBondGenerator(StructureAdapterConstPtr);
        virtual ~BaseBondGenerator()  { }

        // methods
        // loop control
        virtual void rewind();
        bool finished() const;
        void next();

        // configuration
        virtual void selectAnchorSite(int);
        void selectSiteRange(int first, int last);
        void selectSites(const SiteIndices&);
        void selectSites(
                SiteIndices::const_iterator first,
                SiteIndices::const_iterator last);
        virtual void setRmin(double);
        virtual void setRmax(double);

        // get data
        const double& getRmin() const;
        const double& getRmax() const;
        int site0() const;
        int site1() const;
        int multiplicity() const;
        const R3::Vector& r0() const;
        const R3::Vector& r1() const;
        const double& distance() const;
        const R3::Vector& r01() const;
        virtual const R3::Matrix& Ucartesian0() const;
        virtual const R3::Matrix& Ucartesian1() const;
        double msd() const;

    protected:

        // data
        int msite_anchor;
        SiteIndices::const_iterator msite_first;
        SiteIndices::const_iterator msite_last;
        SiteIndices::const_iterator msite_current;
        double mrmin;
        double mrmax;
        StructureAdapterConstPtr mstructure;
        R3::Vector mr0;
        R3::Vector mr1;
        R3::Vector mr01;
        double mdistance;
        SiteIndices msite_all;
        SiteIndices msite_selection;

        // methods
        virtual bool iterateSymmetry();
        virtual void rewindSymmetry();
        virtual void getNextBond();
        void updateDistance();

    private:

        // methods
        void advanceWhileInvalid();
        bool bondOutOfRange() const;
        bool atSelfPair() const;
        void setFinishedFlag();

};


}   // namespace srreal
}   // namespace diffpy

#endif  // BASEBONDGENERATOR_HPP_INCLUDED
