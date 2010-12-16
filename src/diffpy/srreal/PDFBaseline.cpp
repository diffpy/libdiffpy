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
* class PDFBaseline -- abstract base class for PDF baseline functions
*
* $Id$
*
*****************************************************************************/

#include <boost/serialization/export.hpp>

#include <diffpy/srreal/PDFBaseline.hpp>
#include <diffpy/HasClassRegistry.ipp>

// Unique instantiation of the template registry base class.
template class HasClassRegistry<diffpy::srreal::PDFBaseline>;

// Serialization -------------------------------------------------------------

BOOST_CLASS_EXPORT(diffpy::srreal::PDFBaselinePtr)

// End of file
