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

/** @file libsrc/fekernels/FaultCohesive.hh
 *
 * General kernels for faults.
 *
 * Kernel interface.
 *
 * @param[in] dim Spatial dimension.
 * @param[in] numS Number of registered subfields in solution field.
 * @param[in] numA Number of registered subfields in auxiliary field.
 * @param[in] sOff Offset of registered subfields in solution field [numS].
 * @param[in] sOff_x Offset of registered subfields in gradient of the solution field [numS].
 * @param[in] s Solution field with all subfields.
 * @param[in] s_t Time derivative of solution field.
 * @param[in] s_x Gradient of solution field.
 * @param[in] aOff Offset of registered subfields in auxiliary field [numA]
 * @param[in] aOff_x Offset of registered subfields in gradient of auxiliary field [numA]
 * @param[in] a Auxiliary field with all subfields.
 * @param[in] a_t Time derivative of auxiliary field.
 * @param[in] a_x Gradient of auxiliary field.
 * @param[in] t Time for residual evaluation.
 * @param[in] x Coordinates of point evaluation.
 * @param[in] n Outward normal vector for boundary.
 * @param[in] numConstants Number of registered constants.
 * @param[in] constants Array of registered constants.
 * @param[out] f0 [dim].
 */

#if !defined(pylith_fekernels_faultcohesive_hh)
#define pylith_fekernels_faultcohesive_hh

// Include directives ---------------------------------------------------
#include "fekernelsfwd.hh" // forward declarations

#include "pylith/fekernels/Tensor.hh" // USES Tensor
#include "pylith/fekernels/BoundaryDirections.hh" // USES BoundaryDirections

#include "pylith/utils/types.hh"
#include "pylith/utils/journals.hh"

class pylith::fekernels::FaultCohesive {
    // PUBLIC METHODS /////////////////////////////////////////////////////////////////////////////
public:

    // Kernels for diagnostic subfields

    // ------------------------------------------------------------------------------------------------
    /** Fault normal direction.
     */
    static inline
    void normalDir(const PylithInt dim,
                   const PylithInt numS,
                   const PylithInt numA,
                   const PylithInt sOff[],
                   const PylithInt sOff_x[],
                   const PylithScalar s[],
                   const PylithScalar s_t[],
                   const PylithScalar s_x[],
                   const PylithInt aOff[],
                   const PylithInt aOff_x[],
                   const PylithScalar a[],
                   const PylithScalar a_t[],
                   const PylithScalar a_x[],
                   const PylithReal t,
                   const PylithScalar x[],
                   const PylithScalar n[],
                   const PylithInt numConstants,
                   const PylithScalar constants[],
                   PylithScalar normalDir[]) {
        const PylithInt spaceDim = dim; // dim passed in is embedded dimension (spaceDim).

        for (PylithInt i = 0; i < spaceDim; ++i) {
            normalDir[i] = n[i];
        } // for
    } // normalDir

    // ------------------------------------------------------------------------------------------------
    /** Direction along fault strike.
     */
    static inline
    void strikeDir(const PylithInt dim,
                   const PylithInt numS,
                   const PylithInt numA,
                   const PylithInt sOff[],
                   const PylithInt sOff_x[],
                   const PylithScalar s[],
                   const PylithScalar s_t[],
                   const PylithScalar s_x[],
                   const PylithInt aOff[],
                   const PylithInt aOff_x[],
                   const PylithScalar a[],
                   const PylithScalar a_t[],
                   const PylithScalar a_x[],
                   const PylithReal t,
                   const PylithScalar x[],
                   const PylithScalar n[],
                   const PylithInt numConstants,
                   const PylithScalar constants[],
                   PylithScalar strikeDir[]) {
        const PylithInt spaceDim = dim; // dim passed in is embedded dimension (spaceDim).

        switch (spaceDim) {
        case 2: {
            const PylithInt _dim = 2;
            const PylithScalar tanDir[2] = {-n[1], n[0] };
            for (PylithInt i = 0; i < _dim; ++i) {
                strikeDir[i] = tanDir[i];
            } // for
            break;
        } // case 2
        case 3: {
            const PylithInt _dim = 3;
            const PylithScalar* refDir1 = &constants[0];
            const PylithScalar* refDir2 = &constants[3];
            PylithScalar tanDir1[3], tanDir2[3];
            pylith::fekernels::BoundaryDirections::tangential_directions(tanDir1, tanDir2, refDir1, refDir2, n);

            for (PylithInt i = 0; i < _dim; ++i) {
                strikeDir[i] = tanDir1[i];
            } // for
            break;
        } // case 3
        default:
            assert(0);
        } // switch

    } // strikeDir

    // ------------------------------------------------------------------------------------------------
    /** Direction up dip (only valid in 3D).
     */
    static inline
    void dipDir(const PylithInt dim,
                const PylithInt numS,
                const PylithInt numA,
                const PylithInt sOff[],
                const PylithInt sOff_x[],
                const PylithScalar s[],
                const PylithScalar s_t[],
                const PylithScalar s_x[],
                const PylithInt aOff[],
                const PylithInt aOff_x[],
                const PylithScalar a[],
                const PylithScalar a_t[],
                const PylithScalar a_x[],
                const PylithReal t,
                const PylithScalar x[],
                const PylithScalar n[],
                const PylithInt numConstants,
                const PylithScalar constants[],
                PylithScalar dipDir[]) {
        const PylithInt spaceDim = dim; // dim passed in is embedded dimension (spaceDim).

        switch (spaceDim) {
        case 2: {
            PYLITH_JOURNAL_LOGICERROR("Dip direction is not defined in 2D.");
            break;
        } // case 2
        case 3: {
            const PylithInt _dim = 3;
            const PylithScalar* refDir1 = &constants[0];
            const PylithScalar* refDir2 = &constants[3];
            PylithScalar tanDir1[3], tanDir2[3];
            pylith::fekernels::BoundaryDirections::tangential_directions(tanDir1, tanDir2, refDir1, refDir2, n);

            for (PylithInt i = 0; i < _dim; ++i) {
                dipDir[i] = tanDir2[i];
            } // for
            break;
        } // case 3
        default:
            assert(0);
        } // switch

    } // strikeDir

}; // FaultCohesive

#endif // pylith_fekernels_faultcohesive_hh

// End of file
