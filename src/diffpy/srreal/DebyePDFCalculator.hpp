/*****************************************************************************
*
* diffpy.srreal     by DANSE Diffraction group
*                   Simon J. L. Billinge
*                   (c) 2010 Trustees of the Columbia University
*                   in the City of New York.  All rights reserved.
*
* File coded by:    Pavol Juhas, Chris Farrow
*
* See AUTHORS.txt for a list of people who contributed.
* See LICENSE.txt for license information.
*
******************************************************************************
*
* class DebyePDFCalculator -- calculate PDF from the Debye equation.
*
* $Id$
*
*****************************************************************************/

// Comments from Chris Farrow's sources for PDFNanoProfile:
//
// This calculates PDF from the Deybe scattering equation.  
// See chapter 10 of "x-ray diffraction" by B. E. Warren.
//
// I(Q) = 1/N sum_i sum_j f_i f_j sin(Q r_ij)/(Q r_ij)
// From "underneath the bragg peaks"
// S(Q) = (I(Q) - N*(<f^2> + <f>^2))/(N*<f>^2)
//      = ([1/N sum_i sum_j f_i f_j sin(Q r_ij)/(Q r_ij)] - <f^2>)/<f>^2 + 1
// Equation 2.9 from "underneath the bragg peaks"
// F(Q) = Q(S(Q)-1)
// F(Q) = Q ([1/N sum_i sum_j f_i f_j sin(Q r_ij)/(Q r_ij)] - <f^2>)/<f>^2
// The sums where i==j result in N<f^2>. We can skip these terms and avoid
// subtracting N<f^2>.
// F(Q) = 1/N (sum_i sum_(j!=i) f_i f_j sin(Q r_ij)) / (r_ij <f>^2)
//      = 1/N sum_i sum_(j!=i) (f_i/<f>) (f_j/<f>) sin(Q r_ij)/r_ij
//      = 2/N sum_i sum_(j<i) (f_i/<f>) (f_j/<f>) sin(Q r_ij)/r_ij
//
// We then add a thermal motion by multiplying a each term by a Debye-Waller
// factor: exp( -0.5 rij^2 MSDij )
//
// According to Warren, this gives us g(r) = 4 pi r rho(r).  When we
// calculate G(r) from reciprocal space data into real space, we see:  
// G(r) = 4 pi r (rho(r) - rho0).
// In the theory (Warren), the rho0 is added ad hoc. We account for this by
// cutting off the scattering at a low-Q value. The exact value is not
// precise. See the setQmin function. This is described in
// Acta Cryst. A65, 232-239 (2009)

#ifndef DEBYEPDFCALCULATOR_HPP_INCLUDED
#define DEBYEPDFCALCULATOR_HPP_INCLUDED

#include <diffpy/srreal/BaseDebyeSum.hpp>
#include <diffpy/srreal/ScatteringFactorTable.hpp>

namespace diffpy {
namespace srreal {

class DebyePDFCalculator : public BaseDebyeSum,
    public ScatteringFactorTableOwner
{
    public:

        // constructor
        DebyePDFCalculator();

        // results
        /// PDF on the specified r-grid
        QuantityType getPDF() const;
        QuantityType getRgrid() const;

        // Q-range configuration
        virtual void setQstep(double);
        void setOptimumQstep();

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

    protected:

        // Attributes overload to direct visitors around data structures
        virtual void accept(diffpy::BaseAttributesVisitor& v);
        virtual void accept(diffpy::BaseAttributesVisitor& v) const;

        // BaseDebyeSum overloads
        virtual void resetValue();
        virtual void configureBondGenerator(BaseBondGenerator&);
        virtual double sfSiteAtQ(int, const double& Q) const;

    private:

        // methods
        void updateQstep();
        /// complete lower bound extension of the calculated grid
        double rcalclo() const;
        /// complete upper bound extension of the calculated grid
        double rcalchi() const;
        /// r-range extension to allow propagation of termination ripples
        double extFromTerminationRipples() const;
        /// r-range extension to account for tails from out-of-range peaks
        double extFromPeakTails() const;
        void cacheRlimitsData();

        // data
        bool moptimumqstep;
        double mrstep;
        double mmaxextension;
        struct RLimitsCache {
            double totalextension;
            RLimitsCache() :
                totalextension(0.0)
            { }
        };
        RLimitsCache mrlimits_cache;

};  // class DebyePDFCalculator

}   // namespace srreal
}   // namespace diffpy

#endif  // DEBYEPDFCALCULATOR_HPP_INCLUDED