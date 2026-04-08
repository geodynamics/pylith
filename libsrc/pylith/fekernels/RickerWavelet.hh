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

/** @file libsrc/fekernels/RickerWavelet.hh
 *
 */

#if !defined(pylith_fekernels_rickerwavelet_hh)
#define pylith_fekernels_rickerwavelet_hh

// Include directives ---------------------------------------------------
#include "fekernelsfwd.hh" // forward declarations

#include "pylith/utils/types.hh"

#include <cassert> // USES assert()

// =====================================================================================================================
// Kernels for the Ricker Source Time Function in 2D.
// =====================================================================================================================

class pylith::fekernels::RickerWaveletPlaneStrain {
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
        const PylithInt i_rickerwaveletCenterFrequency = numA - 1;

        const PylithScalar* momentTensor = &a[aOff[i_momentTensor]];
        const PylithScalar timeDelay = a[aOff[i_timeDelay]];
        const PylithScalar rickerwaveletCenterFrequency = a[aOff[i_rickerwaveletCenterFrequency]];

        // RickerWavelet source time function (time domain)

        PylithScalar rt = t - timeDelay;
        PylithScalar rickerwavelet = (1.0 - 2.0*PETSC_PI*PETSC_PI*rickerwaveletCenterFrequency*rickerwaveletCenterFrequency*rt*rt) *
                                     PetscExpReal(-PETSC_PI*PETSC_PI*rickerwaveletCenterFrequency*rickerwaveletCenterFrequency*rt*rt);

        // PetscPrintf(PETSC_COMM_WORLD, "timeDelay %f\n", (double)timeDelay);
        // PetscPrintf(PETSC_COMM_WORLD, "t %f\n", (double)t);
        // PetscPrintf(PETSC_COMM_WORLD, "Center Freq %f\n", (double)rickerwaveletCenterFrequency);
        // PetscPrintf(PETSC_COMM_WORLD, "rickerWavelet %f\n", (double)rickerwavelet);

        for (PylithInt i = 0; i < dim*dim; ++i) {
            g1[i] -= momentTensor[i] * rickerwavelet;
            // PetscPrintf(PETSC_COMM_WORLD, "g1[%i]: %f - ricker\n", (int)i, (double)g1[i]);
        } // for
    } // g1v

}; // RickerWaveletPlaneStrain

// =====================================================================================================================
// Kernels for the Ricker Source Time Function in 3D.
// =====================================================================================================================

class pylith::fekernels::RickerWavelet3D {
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
    static inline
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
        const PylithInt i_rickerwaveletCenterFrequency = numA - 1;

        const PylithScalar* momentTensor = &a[aOff[i_momentTensor]];
        const PylithScalar timeDelay = a[aOff[i_timeDelay]];
        const PylithScalar rickerwaveletCenterFrequency = a[aOff[i_rickerwaveletCenterFrequency]];

        // RickerWavelet source time function (time domain)

        PylithScalar rt = t - timeDelay;
        PylithScalar rickerwavelet = (1.0 - 2.0*PETSC_PI*PETSC_PI*rickerwaveletCenterFrequency*rickerwaveletCenterFrequency*rt*rt) *
                                     PetscExpReal(-PETSC_PI*PETSC_PI*rickerwaveletCenterFrequency*rickerwaveletCenterFrequency*rt*rt);

        for (PylithInt i = 0; i < dim*dim; ++i) {
            g1[i] -= momentTensor[i] * rickerwavelet;
        } // for
    } // g1v

}; // RickerWaveletPlaneStrain

#endif /* pylith_fekernels_rickerwavelet_hh */

/* End of file */
