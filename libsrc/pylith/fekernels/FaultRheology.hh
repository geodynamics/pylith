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

/** @file libsrc/fekernels/FaultRheology.hh
 *
 * Kernels for fault rheologies.
 *
 * Solution fields: [disp(dim), vel(dim, optional), lagrange(dim)]
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

#if !defined(pylith_fekernels_faultrheology_hh)
#define pylith_fekernels_faultrheology_hh

// Include directives ---------------------------------------------------
#include "fekernelsfwd.hh" // forward declarations

#include "pylith/utils/types.hh"

class pylith::fekernels::FaultRheology {
    // PUBLIC TYPEDEFS ////////////////////////////////////////////////////////////////////////////
public:

    struct SlipContext {
        PylithReal t;
        PylithReal slip;
        PylithReal slipRate;
    };

    struct FrictionContext {
        PylithInt dim;
        PylithReal cohesion;
        PylithReal tractionNormal;
        PylithReal slipRate[2];
        PylithReal opening;
        bool openFreeSurface;
    };

    typedef void (*frictioncoeffn_type)(const SlipContext&,
                                        void*, // rheologyContext
                                        PylithReal*); // tractionFriction

    typedef void (*frictiondirfn_type)(const FrictionContext&,
                                       const PylithReal, // frictionMag
                                       PylithReal*); // tractionFriction

    // PUBLIC METHODS /////////////////////////////////////////////////////////////////////////////
public:

    // --------------------------------------------------------------------------------------------
    /** Compute friction traction.
     */
    static inline
    void friction(const SlipContext& slipContext,
                  const FrictionContext& frictionContext,
                  void* rheologyContext,
                  frictioncoeffn_type frictionCoefFn,
                  frictiondirfn_type frictionDirFn,
                  PylithReal* tractionFriction) {
        assert(rheologyContext);
        assert(frictionCoefFn);
        assert(frictionDirFn);
        assert(tractionFriction);

        PylithReal frictionCoef = 0.0;
        frictionCoefFn(slipContext, rheologyContext, &frictionCoef);

        PylithReal frictionMag;
        if (0.0 == frictionContext.opening) {
            frictionMag = frictionContext.cohesion + frictionCoef * frictionContext.tractionNormal;
        } else if (!frictionContext.openFreeSurface) {
            frictionMag = frictionContext.cohesion;
        } else {
            frictionMag = 0.0;
        } // if/else

        frictionDirFn(frictionContext, frictionMag, tractionFriction);
    } // friction

    // --------------------------------------------------------------------------------------------
    /** Compute friction direction in 2D.
     */
    static inline
    void frictionDir2D(const FrictionContext& frictionContext,
                       const PylithReal frictionMag,
                       PylithReal* tractionFriction) {
        assert(tractionFriction);
        assert(frictionContext.slipRate);

        const PylithReal slipRate = frictionContext.slipRate[0];
        tractionFriction[0] = -slipRate/fabs(slipRate) * frictionMag;
        tractionFriction[1] = 0.0;
    } // frictionDir2D

    // --------------------------------------------------------------------------------------------
    /** Compute friction direction in 3D.
     */
    static inline
    void frictionDir3D(const FrictionContext& frictionContext,
                       const PylithReal frictionMag,
                       PylithReal* tractionFriction) {
        assert(tractionFriction);
        assert(frictionContext.slipRate);

        const PylithReal slipRate1 = frictionContext.slipRate[0];
        const PylithReal slipRate2 = frictionContext.slipRate[2];
        const PylithReal slipRateMag = sqrt(slipRate1*slipRate1 + slipRate2*slipRate2);
        tractionFriction[0] = -slipRate1/slipRateMag * frictionMag;
        tractionFriction[1] = -slipRate2/slipRateMag * frictionMag;
        tractionFriction[2] = 0.0;
    } // frictionDir3D

}; // FaultRheology

#endif // pylith_fekernels_faultrheology_hh

/* End of file */
