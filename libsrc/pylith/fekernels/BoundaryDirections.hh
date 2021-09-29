/* -*- C++ -*-
 *
 * ----------------------------------------------------------------------
 *
 * Brad T. Aagaard, U.S. Geological Survey
 * Charles A. Williams, GNS Science
 * Matthew G. Knepley, University at Buffalo
 *
 * This code was developed as part of the Computational Infrastructure
 * for Geodynamics (http:*geodynamics.org).
 *
 * Copyright (c) 2010-2021 University of California, Davis
 *
 * See LICENSE.md for license information.
 *
 * ----------------------------------------------------------------------
 */

/** @file libsrc/fekernels/BoundaryDirections.hh
 *
 * Kernels for computing tangential directions for boundaries.
 */

#if !defined(pylith_fekernels_boundarydirections_hh)
#define pylith_fekernels_boundarydirections_hh

#include "fekernelsfwd.hh" // forward declarations

#include "pylith/utils/types.hh"

class pylith::fekernels::BoundaryDirections {
    // PUBLIC MEMBERS /////////////////////////////////////////////////////////////////////////////
public:

    /* Compute tangential directions for boundary in 3-D.
     *
     * @param[out] tanDir1 First tangential direction.
     * @param[out] tanDIr2 Second tangential direction.
     * @param[in] refDir1 First choice for reference direction.
     * @param[in] refDir2 Second choice for reference direction if first fails.
     * @param[in] normDir Normal direction.
     */
    static
    void tangential_directions(PylithScalar tanDir1[],
                               PylithScalar tanDir2[],
                               const PylithScalar refDir1[],
                               const PylithScalar refDir2[],
                               const PylithScalar normDir[]);

}; // BoundaryDirections

#endif // pylith_fekernels_boundarydirections_hh

/* End of file */
