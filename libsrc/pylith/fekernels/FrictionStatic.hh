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

/** @file libsrc/fekernels/FrictionStatic`.hh
 *
 * Kernels for static friction on a fault.
 *
 * Solution fields: [disp(dim), vel(dim, optional), lagrange(dim)]
 *
 * Auxiliary fields: [...] (not used)
 *
 * LHS Residual
 *
 *
 * LHS Jacobian
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

#if !defined(pylith_fekernels_frictionstatic_hh)
#define pylith_fekernels_frictionstatic_hh

// Include directives ---------------------------------------------------
#include "fekernelsfwd.hh" // forward declarations

#include "pylith/fekernels/FaultRheology.hh" // USES FaultRheology

#include "pylith/utils/types.hh"

class pylith::fekernels::FrictionStatic {
    // PUBLIC STRUCTS /////////////////////////////////////////////////////////////////////////////
public:

    struct Context {
        PylithReal staticCoefficient;
    };

    // PUBLIC METHODS /////////////////////////////////////////////////////////////////////////////
public:

    // --------------------------------------------------------------------------------------------
    static inline
    void setContext(Context* context,
                    const PylithInt dim,
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
                    const pylith::fekernels::TensorOps& tensorOps) {
        assert(context);

        const PylithInt i_static = numA-1;

        assert(a);
        assert(aOff);
        assert(aOff[i_static] >= 0);

        context->staticCoefficient = a[aOff[i_static]];assert(context->staticCoefficient >= 0.0);
    } // setContext

    // --------------------------------------------------------------------------------------------
    /** Compute coefficient of friction.
     *
     */
    static inline
    void frictionCoefficient(const pylith::fekernels::FaultRheology::SlipContext& slipContext,
                             void* rheologyContext,
                             PylithReal* coefficient) {
        Context* context = (Context*)(rheologyContext);
        assert(context);
        assert(coefficient);

        *coefficient = context->staticCoefficient;
    } // friction_coefficient

}; // FrictionStatic`

#endif // pylith_fekernels_frictionstatic_hh

/* End of file */
