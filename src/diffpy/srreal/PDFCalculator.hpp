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
* class PDFCalculator -- real space PDF calculator
*
* $Id$
*
*****************************************************************************/

#ifndef PDFCALCULATOR_HPP_INCLUDED
#define PDFCALCULATOR_HPP_INCLUDED

#include <set>
#include <boost/shared_ptr.hpp>

#include <diffpy/srreal/PairQuantity.hpp>
#include <diffpy/srreal/PeakProfile.hpp>
#include <diffpy/srreal/PeakWidthModel.hpp>
#include <diffpy/srreal/PDFBaseline.hpp>
#include <diffpy/srreal/PDFEnvelope.hpp>
#include <diffpy/srreal/ScatteringFactorTable.hpp>

namespace diffpy {
namespace srreal {

class PDFCalculator :
    public PairQuantity,
    public PeakWidthModelOwner,
    public ScatteringFactorTableOwner,
    public PDFEnvelopeOwner
{
    public:

        // constructor
        PDFCalculator();

        // results
        QuantityType getPDF() const;
        QuantityType getRDF() const;
        QuantityType getRDFperR() const;
        QuantityType getF() const;

        /// PDF on an r-range extended for termination ripples
        QuantityType getExtendedPDF() const;
        /// RDF on an r-range extended for termination ripples
        QuantityType getExtendedRDF() const;
        /// RDF divided by r on an r-range extended for termination ripples
        QuantityType getExtendedRDFperR() const;
        /// F(Q) on a zero-padded grid that reaches r-sampling Qmax = PI/dr
        QuantityType getExtendedF() const;
        /// r-grid extended for termination ripples
        QuantityType getExtendedRgrid() const;

        // Q-range methods
        QuantityType getQgrid() const;
        // Q-range configuration
        void setQmin(double);
        const double& getQmin() const;
        void setQmax(double);
        const double& getQmax() const;
        const double& getQstep() const;

        // R-range methods
        QuantityType getRgrid() const;
        // R-range configuration
        virtual void setRmin(double);
        virtual void setRmax(double);
        void setRstep(double);
        const double& getRstep() const;
        /// maximum total extension of the r-range accounting for both
        /// termination ripples and peak tails
        void setMaxExtension(double);
        /// maximum total extension of the r-range accounting for both
        /// termination ripples and peak tails
        const double& getMaxExtension() const;
        /// lower bound for the r-range extended for termination ripples
        double getExtendedRmin() const;
        /// upper bound of the r-range extended for termination ripples
        double getExtendedRmax() const;

        // PDF profile configuration
        void setPeakProfile(PeakProfilePtr);
        void setPeakProfileByType(const std::string& tp);
        PeakProfilePtr& getPeakProfile();
        const PeakProfilePtr& getPeakProfile() const;

        // PDF baseline configuration
        // application on an array
        QuantityType applyBaseline(const QuantityType& x, const QuantityType& y) const;
        void setBaseline(PDFBaselinePtr);
        void setBaselineByType(const std::string& tp);
        PDFBaselinePtr& getBaseline();
        const PDFBaselinePtr& getBaseline() const;

    protected:

        // Attributes overload to direct visitors around data structures
        virtual void accept(diffpy::BaseAttributesVisitor& v);
        virtual void accept(diffpy::BaseAttributesVisitor& v) const;

        // PairQuantity overloads
        virtual void resetValue();
        virtual void configureBondGenerator(BaseBondGenerator&) const;
        virtual void addPairContribution(const BaseBondGenerator&, int);

    private:

        // methods - calculation specific
        /// complete lower bound extension of the calculated grid
        double rcalclo() const;
        /// complete upper bound extension of the calculated grid
        double rcalchi() const;
        /// r-range extension to allow propagation of termination ripples
        double extFromTerminationRipples() const;
        /// r-range extension to account for tails from out-of-range peaks
        double extFromPeakTails() const;
        /// number of dr steps at rcalclo from 0
        int rcalcloSteps() const;
        /// number of dr steps at rhalchi from 0
        int rcalchiSteps() const;
        /// number of dr steps at lower extension due to termination ripples
        int extendedRminSteps() const;
        /// number of dr steps at upper extension due to termination ripples
        int extendedRmaxSteps() const;
        /// number of points in the r-grid extended with termination ripples
        int countExtendedPoints() const;
        /// number of points in the complete calculated r-grid
        int countCalcPoints() const;
        /// index of a nearby point in the complete calculated r-grid
        int calcIndex(double r) const;
        /// reduce extended grid to user-requested results grid
        /// by cutting away the points for termination ripples
        void cutRipplePoints(QuantityType& y) const;

        // structure factors - fast lookup by site index
        /// effective scattering factor at a given site scaled by occupancy
        const double& sfSite(int) const;
        /// average scattering factor
        double sfAverage() const;
        void cacheStructureData();
        void cacheRlimitsData();

        // data
        // configuration
        double mqmin;
        double mqmax;
        double mrstep;
        double mmaxextension;
        PeakProfilePtr mpeakprofile;
        PDFBaselinePtr mbaseline;
        struct {
            std::vector<double> sfsite;
            double sfaverage;
            double totaloccupancy;
        } mstructure_cache;
        struct {
            int extendedrminsteps;
            int extendedrmaxsteps;
            int rcalclosteps;
            int rcalchisteps;
        } mrlimits_cache;

};  // class PDFCalculator

}   // namespace srreal
}   // namespace diffpy

#endif  // PDFCALCULATOR_HPP_INCLUDED
