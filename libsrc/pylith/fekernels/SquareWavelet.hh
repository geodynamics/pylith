/* -*- C++ -*-
 *
 * ----------------------------------------------------------------------
 *
 * Brad T. Aagaard, U.S. Geological Survey
 * Charles A. Williams, GNS Science
 * Matthew G. Knepley, University of Chicago
 *
 * This code was developed as part of the Computational Infrastructure
 * for Geodynamics (http:*geodynamics.org).
 *
 * Copyright (c) 2010-2015 University of California, Davis
 *
 * See COPYING for license information.
 *
 * ----------------------------------------------------------------------
 */

/** @file libsrc/fekernels/SquareWavelet.hh
 *
 */

#if !defined(pylith_fekernels_squarewavelet_hh)
#define pylith_fekernels_squarewavelet_hh

// Include directives ---------------------------------------------------
#include "fekernelsfwd.hh" // forward declarations

#include "pylith/utils/types.hh"

#include <cassert> // USES assert()

// =====================================================================================================================
// Kernels for the Square Source Time Function in 2D.
// =====================================================================================================================

class pylith::fekernels::SquareWaveletPlaneStrain {
    // PUBLIC MEMBERS ///////////////////////////////////////////////////////
public:

    /** Kernel interface.
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
     * @param[in] numConstants Number of registered constants.
     * @param[in] constants Array of registered constants.
     * @param[out] f0 [dim].
     */

    /** g1 function for velocity equation
     *
     */
    static
    void g1v(const PylithInt dim,
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
             const PylithInt numConstants,
             const PylithScalar constants[],
             PylithScalar g1[]) {
        assert(sOff);
        assert(s);
        assert(g1);

        PylithInt _dim = 2;

        // Incoming re-packed solution field.

        // Incoming re-packed auxiliary field.
        const PylithInt i_momentTensor = 0;
        const PylithInt i_timeDelay = 1;
        const PylithInt i_squarewaveletCenterFrequency = numA - 1;

        const PylithScalar* momentTensor = &a[aOff[i_momentTensor]];
        const PylithScalar timeDelay = a[aOff[i_timeDelay]];
        const PylithScalar squarewaveletCenterFrequency = a[aOff[i_squarewaveletCenterFrequency]];

        // SquareWavelet source time function (time domain)

        PylithScalar rt = t - timeDelay;
        PylithScalar squarewavelet = (rt >= 0.0) ? 1.0 : 0.0;

        // PetscPrintf(PETSC_COMM_WORLD, "timeDelay %f\n", (double)timeDelay);
        // PetscPrintf(PETSC_COMM_WORLD, "t %f\n", (double)t);
        // PetscPrintf(PETSC_COMM_WORLD, "Center Freq %f\n", (double)squarewaveletCenterFrequency);
        // PetscPrintf(PETSC_COMM_WORLD, "squareWavelet %f\n", (double)squarewavelet);

        for (PylithInt i = 0; i < dim*dim; ++i) {
            g1[i] -= momentTensor[i] * squarewavelet;
            // PetscPrintf(PETSC_COMM_WORLD, "g1[%i]: %f - square\n", (int)i, (double)g1[i]);
        } // for
    } // g1v
}; // SquareWaveletPlaneStrain

// =====================================================================================================================
// Kernels for the Square Source Time Function in 3D.
// =====================================================================================================================

class pylith::fekernels::SquareWavelet3D {
    // PUBLIC MEMBERS ///////////////////////////////////////////////////////
public:

    /** Kernel interface.
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
     * @param[in] numConstants Number of registered constants.
     * @param[in] constants Array of registered constants.
     * @param[out] f0 [dim].
     */

    /** g1 function for velocity equation
     *
     */
    static
    void g1v(const PylithInt dim,
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
             const PylithInt numConstants,
             const PylithScalar constants[],
             PylithScalar g1[]) {
        assert(sOff);
        assert(s);
        assert(g1);

        PylithInt _dim = 3;

        // Incoming re-packed solution field.

        // Incoming re-packed auxiliary field.
        const PylithInt i_momentTensor = 0;
        const PylithInt i_timeDelay = 1;
        const PylithInt i_squarewaveletCenterFrequency = numA - 1;

        const PylithScalar* momentTensor = &a[aOff[i_momentTensor]];
        const PylithScalar timeDelay = a[aOff[i_timeDelay]];
        const PylithScalar squarewaveletCenterFrequency = a[aOff[i_squarewaveletCenterFrequency]];

        // SquareWavelet source time function (time domain)

        PylithScalar rt = t - timeDelay;
        PylithScalar squarewavelet = (rt >= 0.0) ? 1.0 : 0.0;

        for (PylithInt i = 0; i < dim*dim; ++i) {
            g1[i] -= momentTensor[i] * squarewavelet;
        } // for
    } // g1v
}; // SquareWaveletPlaneStrain

#endif /* pylith_fekernels_squarewavelet_hh */

/* End of file */
