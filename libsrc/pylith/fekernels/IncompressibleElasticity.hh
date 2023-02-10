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

/** @file libsrc/fekernels/IncompressibleElasticity.hh
 *
 * Kernels for incompressible elasticity independent of rheology.
 *
 * Solution fields: [disp(dim), pressure(1)]
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
 * @param[in] numConstants Number of registered constants.
 * @param[in] constants Array of registered constants.
 * @param[out] f0 [dim].
 */

#if !defined(pylith_fekernels_incompressibleelasticity_hh)
#define pylith_fekernels_incompressibleelasticity_hh

#include "fekernelsfwd.hh" // forward declarations

#include "pylith/fekernels/Elasticity.hh" // USES Elasticity kernels

#include "pylith/utils/types.hh"

#include <cassert> // USES assert()

// ------------------------------------------------------------------------------------------------
class pylith::fekernels::IncompressibleElasticity {
public:

    // Function interface for computing incompressible term.
    typedef void (*incompressiblefn_type) (const PylithInt,
                                           const PylithInt,
                                           const PylithInt,
                                           const PylithInt[],
                                           const PylithInt[],
                                           const PylithScalar[],
                                           const PylithScalar[],
                                           const PylithScalar[],
                                           const PylithInt[],
                                           const PylithInt[],
                                           const PylithScalar[],
                                           const PylithScalar[],
                                           const PylithScalar[],
                                           const PylithReal t,
                                           const PylithScalar[],
                                           const PylithInt,
                                           const PylithScalar[],
                                           const pylith::fekernels::Tensor&,
                                           const pylith::fekernels::TensorOps&,
                                           PylithScalar*);

    // PUBLIC MEMBERS //////////////////////////////////////////////////////////////////////////////////////////////////
public:

    // --------------------------------------------------------------------------------------------
    // f0p helper function.
    static inline
    void f0p(const PylithInt dim,
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
             pylith::fekernels::Elasticity::strainfn_type strainFn,
             pylith::fekernels::IncompressibleElasticity::incompressiblefn_type incompressibleFn,
             const pylith::fekernels::TensorOps& tensorOps,
             PylithScalar f0[]) {
        assert(f0);

        pylith::fekernels::Tensor strain;
        strainFn(dim, numS, sOff, sOff, s, s_t, s_x, x, &strain);

        PylithScalar value = 0.0;
        incompressibleFn(dim, numS, numA, sOff, sOff, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x, t, x,
                         numConstants, constants, strain, tensorOps, &value);

        f0[0] += value;
    } // f0p

    // --------------------------------------------------------------------------------------------
    /** Jf1pu entry function for pressure equation for incompressible elasticity.
     *
     * Solution fields: [disp(dim), pressure(1)]
     */
    static inline
    void Jf1pu(const PylithInt dim,
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
               const PylithReal s_tshift,
               const PylithScalar x[],
               const PylithInt numConstants,
               const PylithScalar constants[],
               PylithScalar Jf1[]) {
        assert(Jf1);

        /* j(f,g,dg), f=0, g=0..dim, dg=0..dim
         *
         * j == 1 if g==dg, otherwise 0.
         *
         * 3-D
         * 0: j000 = 1
         * 1: j001 = 0
         * 2: j002 = 0
         * 3: j010 = 0
         * 4: j011 = 1
         * 5: j012 = 0
         * 6: j020 = 0
         * 7: j021 = 0
         * 8: j022 = 1
         */

        for (PylithInt i = 0; i < dim; ++i) {
            Jf1[i*dim+i] += 1.0;
        } // for
    } // Jf1pu

    // --------------------------------------------------------------------------------------------
    /** Jf2up entry function for elasticity equation for incompressible elasticity.
     *
     * Solution fields: [disp(dim), pressure(1)]
     */
    static inline
    void Jf2up(const PylithInt dim,
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
               const PylithReal s_tshift,
               const PylithScalar x[],
               const PylithInt numConstants,
               const PylithScalar constants[],
               PylithScalar Jf2[]) {
        assert(Jf2);

        /* j(f,g,df), f=0..dim, df=0..dim, g=0
         *
         * j == 1 if f==df, otherwise 0.
         *
         * 3-D
         * 0: j000 = 1
         * 1: j001 = 0
         * 2: j002 = 0
         * 3: j100 = 0
         * 4: j101 = 1
         * 5: j102 = 0
         * 6: j200 = 0
         * 7: j201 = 0
         * 8: j202 = 1
         */

        for (PylithInt i = 0; i < dim; ++i) {
            Jf2[i*dim+i] += 1.0;
        } // for
    } // Jf2up

    // --------------------------------------------------------------------------------------------
    /** Calculate mean stress for isotropic linear incompressible elasticity WITHOUT reference stress
     * and strain.
     */
    static inline
    void meanStress(const PylithReal pressure,
                    pylith::fekernels::Tensor* stress) {
        assert(stress);

        stress->xx -= pressure;
        stress->yy -= pressure;
        stress->zz -= pressure;
    } // meanStress

    // --------------------------------------------------------------------------------------------
    /** Calculate mean stress for isotropic linear incompressible elasticity WITH reference stress
     * and strain.
     */
    static inline
    void meanStress_refState(const PylithReal pressure,
                             const pylith::fekernels::Tensor& refStress,
                             pylith::fekernels::Tensor* stress) {
        assert(stress);

        const PylithReal meanRefStress = (refStress.xx + refStress.yy + refStress.zz) / 3.0;
        const PylithScalar meanStress = meanRefStress - pressure;

        stress->xx += meanStress;
        stress->yy += meanStress;
        stress->zz += meanStress;
    } // meanStress_refState

}; // IncompressibleElasticity

#endif // pylith_fekernels_incompressibleelasticity_hh

// End of file
