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

/** @file libsrc/fekernels/IsotropicLinearIncompElasticityPlaneStrain.hh
 *
 * Kernels for linear incompressible elasticity plane strain WITHOUT inertia.
 *
 * Solution fields: [disp(dim), pressure(1)]
 *
 * Auxiliary fields:
 * - 0: density(1)
 * - 1: shear_modulus(1)
 * - 2: bulk_modulus(1)
 * - 3: gravity_field (2, optional)
 * - 4: body_force(2,optional)
 * - 5: reference_stress(4,optional) (stress_xx, stress_yy, stress_xy, stress_zz)
 * - 6: reference_strain(4,optional) (strain_xx, strain_yy, strain_xy, strain_zz)
 *
 * :TODO: @charles Add equation here
 *
 * ======================================================================
 */

#if !defined(pylith_fekernels_isotropiclinearincompelasticityplanestrain_hh)
#define pylith_fekernels_isotropiclinearincompelasticityplanestrain_hh

// Include directives ---------------------------------------------------
#include "fekernelsfwd.hh" // forward declarations

#include "pylith/utils/types.hh"

class pylith::fekernels::IsotropicLinearIncompElasticityPlaneStrain {

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


    /** g1 function for isotropic linear incompressible elasticity plane strain
     *  WITHOUT reference stress and reference strain.
     *
     * Solution fields: [disp(dim), pressure(1)]
     * Auxiliary fields: [density(1), shear_modulus(1), bulk_modulus(1), ...]
     */
    static
    void g1u(const PylithInt dim,
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
             PylithScalar g1[]);


    /** g1 function for isotropic linear incompressible elasticity plane strain
     *  WITH reference stress and reference strain.
     *
     * Solution fields: [disp(dim), pressure(1)]
     * Auxiliary fields: [density(1), shear_modulus(1), bulk_modulus(1), ..., refstress(4), refstrain(4)]
     */
    static
    void g1u_refstate(const PylithInt dim,
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
                      PylithScalar g1[]);


    /** Jg1_pu entry function for 2-D plane strain incompressible isotropic
     * linear elasticity.
     *
     */
    static
    void Jg1pu(const PylithInt dim,
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
               PylithScalar Jg1[]);


    /** Calculate mean stress for 2-D plane strain isotropic linear
     * elasticity WITHOUT reference stress and strain.
     *
     * Solution fields: [disp(dim)]
     * Auxiliary fields: [bulk_modulus(1)]
     */
    static
    void meanStress(const PylithInt dim,
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
                    PylithScalar stress[]);


    /** Calculate mean stress for 2-D plane strain isotropic linear
     * elasticity WITH reference stress and reference strain.
     *
     * Solution fields: [disp(dim)]
     * Auxiliary fields: [bulk_modulus(1), refstress(4), refstrain(4)]
     */
    static
    void meanStress_refstate(const PylithInt dim,
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
                             PylithScalar stress[]);


}; // IsotropicLinearIncompElasticityPlaneStrain

#endif // pylith_fekernels_isotropiclinearincompelasticityplanestrain_hh


// End of file
