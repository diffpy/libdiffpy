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
* class ZeroBaseline -- linear PDF baseline
*
*****************************************************************************/

#include <diffpy/srreal/ZeroBaseline.hpp>
#include <diffpy/serialization.ipp>

using namespace std;

namespace diffpy {
namespace srreal {

// Constructors --------------------------------------------------------------

PDFBaselinePtr ZeroBaseline::create() const
{
    PDFBaselinePtr rv(new ZeroBaseline());
    return rv;
}


PDFBaselinePtr ZeroBaseline::clone() const
{
    PDFBaselinePtr rv(new ZeroBaseline(*this));
    return rv;
}

// Public Methods ------------------------------------------------------------

const string& ZeroBaseline::type() const
{
    static string rv = "zero";
    return rv;
}


double ZeroBaseline::operator()(const double& r) const
{
    return 0.0;
}

// Registration --------------------------------------------------------------

bool reg_ZeroBaseline = ZeroBaseline().registerThisType();

}   // namespace srreal
}   // namespace diffpy

// Serialization -------------------------------------------------------------

DIFFPY_INSTANTIATE_SERIALIZATION(diffpy::srreal::ZeroBaseline)
BOOST_CLASS_EXPORT_IMPLEMENT(diffpy::srreal::ZeroBaseline)

// End of file
