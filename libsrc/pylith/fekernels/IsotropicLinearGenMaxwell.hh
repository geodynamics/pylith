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

/** @file libsrc/fekernels/IsotropicLinearGenMaxwellPlaneStrain.hh
 *
 * Kernels for linear Generalized Maxwell viscoelastic plane strain
 * with 3 Maxwell elements.
 *
 * Solution fields: [disp(dim), vel(dim, optional)]
 *
 * Isotropic, linear Generalized Maxwell viscoelastic plane strain without reference stress/strain.
 * Auxiliary fields:
 * - 0: density(1)
 * - 1: shear_modulus(1)
 * - 2: bulk_modulus(1)
 * - 3: maxwell_time(3) (maxwell_time_1, maxwell_time_2, maxwell_time_3)
 * - 4: shear_modulus_ratio(3) (shear_modulus_ratio_1, shear_modulus_ratio_2, shear_modulus_ratio_3)
 * - 5: total_strain
 *     2D: 4 components (strain_xx, strain_yy, strain_zz, strain_xy)
 *     3D: 6 components (strain_xx, strain_yy, strain_zz, strain_xy, strain_yz, strain_xz)
 * - 6: viscous_strain
 *     2D: 3*4 components (strain1_xx, strain1_yy, strain1_zz, strain1_xy, ...)
 *     3D: 3*6 components (strain1_xx, strain1_yy, strain1_zz, strain1_xy, strain1_yz, strain1_xz, ...)
 * - 7: gravity_field (2, optional)
 * - 8: body_force(2,optional)
 * - 9: reference_stress(optional) (stress_xx, stress_yy, stress_zz, stress_xy)
 *     2D: 4 components (stress_xx, stress_yy, stress_zz, stress_xy)
 *     3D: 6 components (stress_xx, stress_yy, stress_zz, stress_xy, stress_yz, stress_xz)
 * -10: reference_strain(optional) (strain_xx, strain_yy, strain_zz, strain_xy)
 *     2D: 4 components (strain_xx, strain_yy, strain_zz, strain_xy)
 *     3D: 6 components (strain_xx, strain_yy, strain_zz, strain_xy, strain_yz, strain_xz)
 *
 * \int_V \vec{\phi}_u \cdot \left( \rho \frac{\partial \vec{v}(t)}{\partial t} \right) \, dV =
 *   \int_V \vec{\phi}_u \cdot \vec{f}(t) - \nabla \vec{\phi}_u : \tensor{\sigma}(\vec{u}) \, dV +
 *   \int_{S_\tau} \vec{\phi}_u \cdot \vec{\tau}(t) \, dS.
 *
 * =====================================================================================================================
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

#if !defined(pylith_fekernels_isotropiclineargenmaxwell_hh)
#define pylith_fekernels_isotropiclineargenmaxwell_hh

#include "fekernelsfwd.hh" // forward declarations

#include "pylith/utils/types.hh"

// ---------------------------------------------------------------------------------------------------------------------
/// Kernels for isotropic, linearly generalized Maxwell viscoelastic that are independent of spatial dimension.
class pylith::fekernels::IsotropicLinearGenMaxwell {
    // PUBLIC MEMBERS //////////////////////////////////////////////////////////////////////////////////////////////////
public:

    /** f0 function for isotropic linear generalized Maxwell.
     *
     * Solution fields: [disp(dim), vel(dim)]
     * Auxiliary fields: [density(1), ...]
     */
    static
    void f0v(const PylithInt dim,
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
             PylithScalar f0[]);

    /** Jf0 function for isotropoc linear generalized Maxwell.
     *
     * Solution fields: [...]
     * Auxiliary fields: [density(1), ...]
     */
    static
    void Jf0vv(const PylithInt dim,
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
               PylithScalar Jf0[]);

    /** g0 function for isotropic linear Generalized Maxwell with both gravity and body forces.
     *
     * Solution fields: [...]
     * Auxiliary fields: [0: density(1),
     *                    1: shear_modulus(1),
     *                    2: bulk_modulus(1),
     *                    3: maxwell_time(3),
     *                    4: shear_modulus_ratio(3),
     *                    5: total_strain(4),
     *                    6: viscous_strain(12),
     *                    7: gravity_field (2, optional),
     *                    8: body_force(2,optional),
     *                    9: reference_stress(4,optional),
     *                   10: reference_strain(4,optional)]
     */
    static
    void g0v_gravbodyforce(const PylithInt dim,
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
                           PylithScalar g0[]);

    /** g0 function for isotropic linear generalized Maxwell with gravity.
     *
     * Solution fields: [...]
     * Auxiliary fields: [0: density(1),
     *                    1: shear_modulus(1),
     *                    2: bulk_modulus(1),
     *                    3: maxwell_time(3),
     *                    4: shear_modulus_ratio(3),
     *                    5: total_strain(4),
     *                    6: viscous_strain(12),
     *                    7: gravity_field (2, optional),
     *                    8: body_force(2,optional),
     *                    9: reference_stress(4,optional),
     *                   10: reference_strain(4,optional)]
     */
    static
    void g0v_grav(const PylithInt dim,
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
                  PylithScalar g0[]);

    /** g0 function for isotropic linear generalized Maxwell with body forces.
     *
     * Solution fields: [...]
     * Auxiliary fields: [0: density(1),
     *                    1: shear_modulus(1),
     *                    2: bulk_modulus(1),
     *                    3: maxwell_time(3),
     *                    4: shear_modulus_ratio(3),
     *                    5: total_strain(4),
     *                    6: viscous_strain(12),
     *                    7: gravity_field (2, optional),
     *                    8: body_force(2,optional),
     *                    9: reference_stress(4,optional),
     *                   10: reference_strain(4,optional)]
     */
    static
    void g0v_bodyforce(const PylithInt dim,
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
                       PylithScalar g0[]);

}; // IsotropicLinearGenMaxwell

// ---------------------------------------------------------------------------------------------------------------------
/// Kernels for isotropic, linear generalized Maxwell plane strain.
class pylith::fekernels::IsotropicLinearGenMaxwellPlaneStrain {
    // PUBLIC MEMBERS //////////////////////////////////////////////////////////////////////////////////////////////////
public:

    /** g1 function for isotropic linear generalized Maxwell plane strain WITHOUT reference stress and reference strain.
     *
     * Solution fields: [disp(dim), ...]
     * Auxiliary fields: [0: density(1),
     *                    1: shear_modulus(1),
     *                    2: bulk_modulus(1),
     *                    3: maxwell_time(3),
     *                    4: shear_modulus_ratio(3),
     *                    5: total_strain(4),
     *                    6: viscous_strain(12),
     *                    7: gravity_field (2, optional),
     *                    8: body_force(2,optional),
     *                    9: reference_stress(4,optional),
     *                   10: reference_strain(4,optional)]
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
             PylithScalar g1[]);

    /** g1 function for isotropic linear generalized Maxwell plane strain WITH reference stress and reference strain.
     *
     * Solution fields: [disp(dim), ...]
     * Auxiliary fields: [0: density(1),
     *                    1: shear_modulus(1),
     *                    2: bulk_modulus(1),
     *                    3: maxwell_time(3),
     *                    4: shear_modulus_ratio(3),
     *                    5: total_strain(4),
     *                    6: viscous_strain(12),
     *                    7: gravity_field (2, optional),
     *                    8: body_force(2,optional),
     *                    9: reference_stress(4,optional),
     *                   10: reference_strain(4,optional)]
     */
    static
    void g1v_refstate(const PylithInt dim,
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

    /** Jg3_vu entry function for 2-D plane strain isotropic linear generalized Maxwell viscoelasticity.
     *
     * Solution fields: [...]
     * Auxiliary fields: [0: density(1),
     *                    1: shear_modulus(1),
     *                    2: bulk_modulus(1),
     *                    3: maxwell_time(3),
     *                    4: shear_modulus_ratio(3),
     *                    5: total_strain(4),
     *                    6: viscous_strain(12),
     *                    7: gravity_field (2, optional),
     *                    8: body_force(2,optional),
     *                    9: reference_stress(4,optional),
     *                   10: reference_strain(4,optional)]
     */
    static
    void Jg3vu(const PylithInt dim,
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
               PylithScalar Jg3[]);

    /** Calculate deviatoric stress for 2-D plane strain isotropic linear
     * generalized Maxwell viscoelasticity WITHOUT reference stress and strain.
     *
     * Solution fields: [disp(dim)]
     * Auxiliary fields: [0: shear_modulus(1),
     *                    1: maxwell_time(3),
     *                    2: shear_modulus_ratio(3),
     *                    3: total_strain(4),
     *                    4: viscous_strain(12),
     */
    static
    void deviatoricStress(const PylithInt dim,
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

    /** Calculate deviatoric stress for 2-D plane strain isotropic linear
     * generalized Maxwell viscoelasticity WITH reference stress and strain.
     *
     * Solution fields: [disp(dim)]
     * Auxiliary fields: [0: shear_modulus(1),
     *                    1: maxwell_time(3),
     *                    2: shear_modulus_ratio(3),
     *                    3: total_strain(4),
     *                    4: viscous_strain(12),
     *                    5: reference_stress(4),
     *                    6: reference_strain(4)]
     */
    static
    void deviatoricStress_refstate(const PylithInt dim,
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

    /** Calculate viscous strain at t+dt for 2-D plane strain isotropic linear
     * generalized Maxwell viscoelasticity.
     *
     * Solution fields: [disp(dim)]
     * Auxiliary fields: [maxwell_time(3), total_strain(4), viscous_strain(12)]
     */
    static
    void computeViscousStrain(const PylithInt dim,
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
                              PylithScalar visStrainTpdt[]);

    /** Update total strain for 2-D plane strain isotropic linear
     * generalized Maxwell viscoelasticity.
     *
     * Solution fields: [disp(dim)]
     * Auxiliary fields: [...]
     */
    static
    void updateTotalStrain(const PylithInt dim,
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
                           PylithScalar totalStrainTpdt[]);

    /** Update viscous strain for generalized Maxwell.
     *
     * Solution fields: [disp(dim)]
     * Auxiliary fields: [0: density(1),
     *                    1: shear_modulus(1),
     *                    2: bulk_modulus(1),
     *                    3: maxwell_time(3),
     *                    4: shear_modulus_ratio(3),
     *                    5: total_strain(4),
     *                    6: viscous_strain(12),
     *                   ...]
     */
    static
    void updateViscousStrain(const PylithInt dim,
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
                             PylithScalar visStrainTpdt[]);

}; // IsotropicLinearGenMaxwellPlaneStrain

// ---------------------------------------------------------------------------------------------------------------------
/// Kernels for isotropic, linear generalized Maxwell viscoelastic in 3D.
class pylith::fekernels::IsotropicLinearGenMaxwell3D {
    // PUBLIC MEMBERS //////////////////////////////////////////////////////////////////////////////////////////////////
public:

    /** g1 function for isotropic linear generalized Maxwell 3D WITHOUT reference stress and reference strain.
     *
     * Solution fields: [disp(dim), ...]
     * Auxiliary fields: [0: density(1),
     *                    1: shear_modulus(1),
     *                    2: bulk_modulus(1),
     *                    3: maxwell_time(3),
     *                    4: shear_modulus_ratio(3),
     *                    5: viscous_strain(18),
     *                    6: total_strain(6),
     *                    7: gravity_field (3, optional),
     *                    8: body_force(3,optional),
     *                    9: reference_stress(6,optional),
     *                   10: reference_strain(6,optional)]
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
             PylithScalar g1[]);

    /** g1 function for isotropic linear generalized Maxwell 3D WITH reference stress and reference strain.
     *
     * Solution fields: [disp(dim), ...]
     * Auxiliary fields: [0: density(1),
     *                    1: shear_modulus(1),
     *                    2: bulk_modulus(1),
     *                    3: maxwell_time(3),
     *                    4: shear_modulus_ratio(3),
     *                    5: viscous_strain(18),
     *                    6: total_strain(6),
     *                    7: gravity_field (3, optional),
     *                    8: body_force(3,optional),
     *                    9: reference_stress(6,optional),
     *                   10: reference_strain(6,optional)]
     */
    static
    void g1v_refstate(const PylithInt dim,
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

    /** Jg3_vu entry function for 3-D isotropic linear generalized Maxwell viscoelasticity.
     *
     * Solution fields: [...]
     * Auxiliary fields: [0: density(1),
     *                    1: shear_modulus(1),
     *                    2: bulk_modulus(1),
     *                    3: maxwell_time(3),
     *                    4: shear_modulus_ratio(3),
     *                    5: viscous_strain(18),
     *                    6: total_strain(6),
     *                    7: gravity_field (3, optional),
     *                    8: body_force(3,optional),
     *                    9: reference_stress(6,optional),
     *                   10: reference_strain(6,optional)]
     */
    static
    void Jg3vu(const PylithInt dim,
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
               const PylithReal utshift,
               const PylithScalar x[],
               const PylithInt numConstants,
               const PylithScalar constants[],
               PylithScalar Jg3[]);

    /** Calculate deviatoric stress for 3-D isotropic linear
     * generalized Maxwell viscoelasticity WITHOUT reference stress and strain.
     *
     * Solution fields: [disp(dim)]
     * Auxiliary fields: [0: shear_modulus(1),
     *                    1: maxwell_time(3),
     *                    2: shear_modulus_ratio(3),
     *                    3: viscous_strain(18),
     *                    4: total_strain(6),
     */
    static
    void deviatoricStress(const PylithInt dim,
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

    /** Calculate deviatoric stress for 3-D isotropic linear
     * generalized Maxwell viscoelasticity WITH reference stress and strain.
     *
     * Solution fields: [disp(dim)]
     * Auxiliary fields: [0: shear_modulus(1),
     *                    1: maxwell_time(3),
     *                    2: shear_modulus_ratio(3),
     *                    3: viscous_strain(18),
     *                    4: total_strain(6),
     *                    5: reference_stress(6),
     *                    6: reference_strain(6)]
     */
    static
    void deviatoricStress_refstate(const PylithInt dim,
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

    /** Calculate viscous strain at t+dt for 3-D isotropic linear
     * generalized Maxwell viscoelasticity.
     *
     * Solution fields: [disp(dim)]
     * Auxiliary fields: [maxwell_time(3), viscous_strain(18), total_strain(6)]
     */
    static
    void computeViscousStrain(const PylithInt dim,
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
                              PylithScalar visStrainTpdt[]);

    /** Update total strain for 3-D isotropic linear
     * generalized Maxwell viscoelasticity.
     *
     * Solution fields: [disp(dim)]
     * Auxiliary fields: [...]
     */
    static
    void updateTotalStrain(const PylithInt dim,
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
                           PylithScalar totalStrainTpdt[]);

    /** Update viscous strain for generalized Maxwell.
     *
     * Solution fields: [disp(dim)]
     * Auxiliary fields: [0: density(1),
     *                    1: shear_modulus(1),
     *                    2: bulk_modulus(1),
     *                    3: maxwell_time(3),
     *                    4: shear_modulus_ratio(3),
     *                    4: viscous_strain(18),
     *                    6: total_strain(6),
     *                   ...]
     */
    static
    void updateViscousStrain(const PylithInt dim,
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
                             PylithScalar visStrainTpdt[]);

}; // IsotropicLinearGenMaxwell3D

#endif // pylith_fekernels_isotropiclineargenmaxwell_hh

// End of file
