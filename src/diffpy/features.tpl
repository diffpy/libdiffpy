/*****************************************************************************
*
* libdiffpy         Complex Modeling Initiative
*                   (c) 2014 Brookhaven Science Associates,
*                   Brookhaven National Laboratory.
*                   All rights reserved.
*
* File coded by:    Pavol Juhas
*
* See AUTHORS.txt for a list of people who contributed.
* See LICENSE.txt for license information.
*
******************************************************************************
*
* Macro definitions for optional features in libdiffpy.
*
*****************************************************************************/

#ifndef FEATURES_HPP_INCLUDED
#define FEATURES_HPP_INCLUDED

// Optional Features ---------------------------------------------------------

#if ${DIFFPY_HAS_OBJCRYST}
# define DIFFPY_HAS_OBJCRYST
#endif

// FIXME -- temporary features, remove when released
#define DIFFPY_DEV_PEAKWIDTHMODEL_SERIALIZATION
#define DIFFPY_DEV_DIRECT_SERIALIZATION
#define DIFFPY_DEV_CONSTANTPEAKWIDTH_UISOWIDTH

#endif  // FEATURES_HPP_INCLUDED

// vim:ft=cpp:
