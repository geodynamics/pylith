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

/** @file libsrc/fekernels/TimeHistoryWavelet.hh
 *
 */

#if !defined(pylith_fekernels_timehistorywavelet_hh)
#define pylith_fekernels_timehistorywavelet_hh

// Include directives ---------------------------------------------------
#include "fekernelsfwd.hh" // forward declarations

#include "pylith/utils/types.hh"

#include <cassert> // USES assert()

// =====================================================================================================================
// Kernels for the TimeHistory Source Time Function in 2D.
// =====================================================================================================================

class pylith::fekernels::TimeHistoryWaveletPlaneStrain {
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
        const PylithInt _numA = 3;
        assert(_numA == numA);
        assert(aOff);
        assert(a);

        // const PylithInt i_amplitude = aOff[numA - 3];
        const PylithInt i_moment_tensor = aOff[0];
        const PylithInt i_start_time = aOff[numA - 2];
        const PylithInt i_value = aOff[numA - 1];

        PylithInt _dim = 2;

        // PetscPrintf(PETSC_COMM_WORLD, "timeDelay %f\n", (double)timeDelay);
        // PetscPrintf(PETSC_COMM_WORLD, "t %f\n", (double)t);
        // PetscPrintf(PETSC_COMM_WORLD, "Center Freq %f\n", (double)timehistorywaveletCenterFrequency);
        // PetscPrintf(PETSC_COMM_WORLD, "timehistoryWavelet %f\n", (double)timehistorywavelet);

        const PylithScalar* momentTensor = &a[aOff[i_moment_tensor]];

        for (PylithInt i = 0; i < dim*dim; ++i) {
            g1[i] -= momentTensor[i] * a[i_value];
            // PetscPrintf(PETSC_COMM_WORLD, "g1[%i]: %f - timehistory\n", (int)i, (double)g1[i]);
        } // for
    } // g1v

}; // TimeHistoryWaveletPlaneStrain

// =====================================================================================================================
// Kernels for the TimeHistory Source Time Function in 3D.
// =====================================================================================================================

class pylith::fekernels::TimeHistoryWavelet3D {
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
        const PylithInt _numA = 3;
        assert(_numA == numA);
        assert(aOff);
        assert(a);

        // const PylithInt i_amplitude = aOff[numA - 3];
        const PylithInt i_moment_tensor = aOff[0];
        // const PylithInt i_start_time = aOff[numA - 2];
        const PylithInt i_value = aOff[numA - 1];

        PylithInt _dim = 3;

        // PetscPrintf(PETSC_COMM_WORLD, "timeDelay %f\n", (double)timeDelay);
        // PetscPrintf(PETSC_COMM_WORLD, "t %f\n", (double)t);
        // PetscPrintf(PETSC_COMM_WORLD, "Center Freq %f\n", (double)timehistorywaveletCenterFrequency);
        // PetscPrintf(PETSC_COMM_WORLD, "timehistoryWavelet %f\n", (double)timehistorywavelet);

        const PylithScalar* momentTensor = &a[aOff[i_moment_tensor]];

        for (PylithInt i = 0; i < dim*dim; ++i) {
            g1[i] -= momentTensor[i] * a[i_value];
            // PetscPrintf(PETSC_COMM_WORLD, "g1[%i]: %f - timehistory\n", (int)i, (double)g1[i]);
        } // for
    } // g1v

}; // g1v

#endif /* pylith_fekernels_timehistorywavelet_hh */

/* End of file */
