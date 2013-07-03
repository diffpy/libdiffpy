/*****************************************************************************
*
* diffpy.srreal     Complex Modeling Initiative
*                   Pavol Juhas
*                   (c) 2013 Brookhaven National Laboratory,
*                   Upton, New York.  All rights reserved.
*
* File coded by:    Pavol Juhas
*
* See AUTHORS.txt for a list of people who contributed.
* See LICENSE.txt for license information.
*
******************************************************************************
*
* X-ray scattering factors for ions and neutral atoms obtained from
* f0_WaasKirf.dat and the associated reference
* D. Waasmaier, A. Kirfel, Acta Cryst. (1995). A51, 416-413
* http://dx.doi.org/10.1107/S0108767394013292
*
* Electron scattering factors approximated from the X-rays
# Reference: International Tables Volume C, page 224.
*
*****************************************************************************/

#ifndef SCATTERINGFACTORS_HPP_INCLUDED
#define SCATTERINGFACTORS_HPP_INCLUDED

#include <string>

namespace diffpy {
namespace srreal {

/// X-ray scattering factor of an element or ion a given Q
double fxrayatq(const std::string& smbl, double q);


/// X-ray scattering factor of an element or ion a given sin(theta)/lambda
double fxrayatstol(const std::string& smbl, double stol);

/// Electron scattering factor of an element or ion a given Q
double felectronatq(const std::string& smbl, double q);

/// Number of electrons for an element or ion
int electronnumber(const std::string& smbl);

/// Coherent scattering length of an element or isotope in fm
double bcneutron(const std::string& smbl);

}   // namespace srreal
}   // namespace diffpy

#endif  // SCATTERINGFACTORS_HPP_INCLUDED
