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
* class PDFBaseline -- abstract base class for PDF baseline functions
*
*****************************************************************************/

#include <diffpy/srreal/PDFBaseline.hpp>
#include <diffpy/HasClassRegistry.ipp>
#include <diffpy/serialization.ipp>

namespace diffpy {

// Unique instantiation of the template registry base class.
template class HasClassRegistry<srreal::PDFBaseline>;

}   // namespace diffpy

// Serialization -------------------------------------------------------------

DIFFPY_INSTANTIATE_SERIALIZATION(diffpy::srreal::PDFBaseline)
DIFFPY_INSTANTIATE_PTR_SERIALIZATION(diffpy::srreal::PDFBaseline)

// End of file
