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
 * Copyright (c) 2010-2022 University of California, Davis
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

#include <cassert> // USES assert()

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
    static inline
    void tangential_directions(PylithScalar tanDir1[],
                               PylithScalar tanDir2[],
                               const PylithScalar refDir1[],
                               const PylithScalar refDir2[],
                               const PylithScalar normalDir[]) {
        assert(tanDir1);
        assert(tanDir2);
        assert(refDir1);
        assert(refDir2);
        assert(normalDir);

        const int dim = 3;

        // Choose reference direction 1 unless it nearly coincides with normal direction.
        PylithScalar refDir[dim] = { refDir1[0], refDir1[1], refDir1[2] };
        if (fabs(refDir[0]*normalDir[0] + refDir[1]*normalDir[1] + refDir[2]*normalDir[2]) > 0.999) {
            for (int i = 0; i < dim; ++i) {
                refDir[i] = refDir2[i];
            } // for
        } // if

        // tanDir1 = refDir x normalDir
        _unitCrossProduct(tanDir1, refDir, normalDir);

        // tranDir2 = normalDir x tanDir1
        _unitCrossProduct(tanDir2, normalDir, tanDir1);
    }

private:

    static inline
    void _unitCrossProduct(PylithScalar result[],
                           const PylithScalar vec1[],
                           const PylithScalar vec2[]) {
        assert(result);
        assert(vec1);
        assert(vec2);

        result[0] = +vec1[1]*vec2[2] - vec1[2]*vec2[1];
        result[1] = +vec1[2]*vec2[0] - vec1[0]*vec2[2];
        result[2] = +vec1[0]*vec2[1] - vec1[1]*vec2[0];

        const int dim = 3;
        PylithScalar resultMag = 0.0;
        for (int i = 0; i < dim; ++i) {
            resultMag += result[i]*result[i];
        } // for
        resultMag = sqrt(resultMag);
        assert(resultMag >= 1.0e-4);
        for (int i = 0; i < dim; ++i) {
            result[i] /= resultMag;
        } // for
    } // _unitCrossProduct

}; // BoundaryDirections

#endif // pylith_fekernels_boundarydirections_hh

/* End of file */
