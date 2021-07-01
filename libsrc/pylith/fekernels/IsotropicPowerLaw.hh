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
 * Copyright (c) 2010-2021 University of California, Davis
 *
 * See LICENSE.md for license information.
 *
 * ----------------------------------------------------------------------
 */

/** @file libsrc/fekernels/IsotropicPowerLaw.hh
 *
 * Kernels for power-law viscoelastic material.
 *
 * Solution fields: [disp(dim), ...]
 *
 * Isotropic, power-law viscoelastic material.
 * Auxiliary fields:
 * - 0: density(1)
 * - 1: body_force(2,optional)
 * - 2: gravity_field (2, optional)
 * - 3: reference_stress(optional)
 *     2D: 4 components (stress_xx, stress_yy, stress_zz, stress_xy)
 *     3D: 6 components (stress_xx, stress_yy, stress_zz, stress_xy, stress_yz, stress_xz)
 * - 4: reference_strain(optional)
 *     2D: 4 components (strain_xx, strain_yy, strain_zz, strain_xy)
 *     3D: 6 components (strain_xx, strain_yy, strain_zz, strain_xy, strain_yz, strain_xz)
 * - 5: shear_modulus(1)
 * - 6: bulk_modulus(1)
 * - 7: power_law_reference_strain_rate(1)
 * - 8: power_law_reference_stress(1)
 * - 9: power_law_exponent(1)
 * -10: viscous_strain
 *     2D: 4 components (strain_xx, strain_yy, strain_zz, strain_xy)
 *     3D: 6 components (strain_xx, strain_yy, strain_zz, strain_xy, strain_yz, strain_xz)
 * -11: stress
 *     2D: 4 components (stress_xx, stress_yy, stress_zz, stress_xy)
 *     3D: 6 components (stress_xx, stress_yy, stress_zz, stress_xy, stress_yz, stress_xz)
 *
 * The elasticity subfields come first (with required ones before optional ones) followed by the rheology subfields
 * (optional ones before required ones). The rheology fields have required fields last because we index from the back.
 *
 * Viscous strain must be before total strain, because viscous strain at t+dt depends on total strain at t.
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

#if !defined(pylith_fekernels_isotropicpowerlaw_hh)
#define pylith_fekernels_isotropicpowerlaw_hh

#include "fekernelsfwd.hh" // forward declarations

#include "pylith/utils/types.hh"

// ---------------------------------------------------------------------------------------------------------------------
/// Kernels for isotropic power-law plane strain.
class pylith::fekernels::IsotropicPowerLawPlaneStrain {
    // PUBLIC MEMBERS //////////////////////////////////////////////////////////////////////////////////////////////////
public:

    /** f1 function for isotropic power-law plane strain WITHOUT reference stress and reference strain.
     *
     * Solution fields: [disp(dim), ...]
     * Auxiliary fields: [..., shear_modulus(1), bulk_modulus(1), power_law_reference_strain_rate(1),
     * power_law_reference_stress(1), power_law_exponent(1), viscous_strain(4), stress(4)]
     */
    static
    void f1v(const PylithInt dim,
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
             PylithScalar f1[]);

    /** f1 function for isotropic power-law plane strain WITH reference stress and reference strain.
     *
     * Solution fields: [disp(dim), ...]
     * Auxiliary fields: [..., reference_stress(4), reference_strain(4), shear_modulus(1), bulk_modulus(1),
     *                    power_law_reference_strain_rate(1), power_law_reference_stress(1), power_law_exponent(1),
     *                    viscous_strain(4), stress(4)]
     */
    static
    void f1v_refstate(const PylithInt dim,
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
                      PylithScalar f1[]);

    /** Jf3_vu entry function for plane strain isotropic power-law viscoelasticity WITHOUT reference stress/strain.
     *
     * Solution fields: [...]
     * Auxiliary fields: [..., shear_modulus(1), bulk_modulus(1), power_law_reference_strain_rate(1),
     * power_law_reference_stress(1), power_law_exponent(1), viscous_strain(4), stress(4)]
     */
    static
    void Jf3vu(const PylithInt dim,
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
               PylithScalar Jf3[]);

    /** Jf3_vu entry function for plane strain isotropic power-law viscoelasticity WITH reference stress/strain.
     *
     * Solution fields: [...]
     * Auxiliary fields: [..., reference_stress(4), reference_strain(4), shear_modulus(1), bulk_modulus(1),
     *                    power_law_reference_strain_rate(1), power_law_reference_stress(1), power_law_exponent(1),
     *                    viscous_strain(4), stress(4)]
     */
    static
    void Jf3vu_refstate(const PylithInt dim,
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
                        PylithScalar Jf3[]);

    /** Calculate deviatoric stress for 2-D plane strain isotropic power-law
     * viscoelasticity WITHOUT reference stress and strain.
     *
     * Solution fields: [disp(dim)]
     * Auxiliary fields: [..., shear_modulus(1), bulk_modulus(1), power_law_reference_strain_rate(1),
     * power_law_reference_stress(1), power_law_exponent(1), viscous_strain(4), stress(4)]
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

    /** Calculate deviatoric stress for 2-D plane strain isotropic power-law
     * viscoelasticity WITH reference stress and strain.
     *
     * Solution fields: [disp(dim)]
     * Auxiliary fields: [..., reference_stress(4), reference_strain(4), shear_modulus(1), bulk_modulus(1),
     *                    power_law_reference_strain_rate(1), power_law_reference_stress(1), power_law_exponent(1),
     *                    viscous_strain(4), stress(4)]
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

    /** Calculate deviatoric stress including stress_zz for 2-D plane strain isotropic power-law
     * viscoelasticity WITHOUT reference stress and strain.
     *
     * Solution fields: [disp(dim)]
     * Auxiliary fields: [..., shear_modulus(1), bulk_modulus(1), power_law_reference_strain_rate(1),
     * power_law_reference_stress(1), power_law_exponent(1), viscous_strain(4), stress(4)]
     */
    static
    void deviatoricStress4(const PylithInt dim,
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
                           PylithScalar devStressTpdt[]);

    /** Calculate deviatoric stress including stress_zz for 2-D plane strain isotropic power-law
     * viscoelasticity WITH reference stress and strain.
     *
     * Solution fields: [disp(dim)]
     * Auxiliary fields: [..., reference_stress(4), reference_strain(4), shear_modulus(1), bulk_modulus(1),
     *                    power_law_reference_strain_rate(1), power_law_reference_stress(1), power_law_exponent(1),
     *                    viscous_strain(4), stress(4)]
     */
    static
    void deviatoricStress4_refstate(const PylithInt dim,
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

    /** Update stress for 2-D plane strain isotropic power-law
     * viscoelasticity WITHOUT reference stress/strain.
     *
     * Solution fields: [disp(dim)]
     * Auxiliary fields: [...]
     */
    static
    void updateStress(const PylithInt dim,
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

    /** Update stress for 2-D plane strain isotropic power-law
     * viscoelasticity WITH reference stress/strain.
     *
     * Solution fields: [disp(dim)]
     * Auxiliary fields: [...]
     */
    static
    void updateStress_refstate(const PylithInt dim,
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

    /** Update viscous strain for plane strain isotropic power-law viscoelasticity WITHOUT reference stress/strain.
     *
     * Solution fields: [disp(dim)]
     * Auxiliary fields: [..., shear_modulus(1), bulk_modulus(1), power_law_reference_strain_rate(1),
     * power_law_reference_stress(1), power_law_exponent(1), viscous_strain(4), stress(4)]
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
                             PylithScalar visStrain[]);

    /** Update viscous strain for plane strain isotropic power-law viscoelasticity WITH reference stress/strain.
     *
     * Solution fields: [disp(dim)]
     * Auxiliary fields: [..., reference_stress(4), reference_strain(4), shear_modulus(1), bulk_modulus(1),
     *                    power_law_reference_strain_rate(1), power_law_reference_stress(1), power_law_exponent(1),
     *                    viscous_strain(4), stress(4)]
     */
    static
    void updateViscousStrain_refstate(const PylithInt dim,
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
                                      PylithScalar visStrain[]);

    /** Calculate stress for 2-D plane strain isotropic power-law
     * WITHOUT a reference stress and strain.
     *
     * Used in outputing the stress field.
     *
     * Solution fields: [disp(dim)]
     * Auxiliary fields: [..., shear_modulus(1), bulk_modulus(1), power_law_reference_strain_rate(1),
     *                    power_law_reference_stress(1), power_law_exponent(1), viscous_strain(4), stress(4)]
     */
    static
    void cauchyStress(const PylithInt dim,
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
                      PylithScalar stressVector[]);

    /** Calculate stress for 2-D plane strain isotropic power-law
     * WITH a reference stress/strain.
     *
     * Used in outputing the stress field.
     *
     * Solution fields: [disp(dim)]
     * Auxiliary fields: [..., reference_stress(4), reference_strain(4), shear_modulus(1), bulk_modulus(1),
     *                    power_law_reference_strain_rate(1), power_law_reference_stress(1), power_law_exponent(1),
     *                    viscous_strain(4), stress(4)]
     */
    static
    void cauchyStress_refstate(const PylithInt dim,
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
                               PylithScalar stressVector[]);

}; // IsotropicPowerLawPlaneStrain

// ---------------------------------------------------------------------------------------------------------------------
/// Kernels for isotropic power-law viscoelastic material in 3D.
class pylith::fekernels::IsotropicPowerLaw3D {
    // PUBLIC MEMBERS //////////////////////////////////////////////////////////////////////////////////////////////////
public:

    /** f1 function for isotropic power-law viscoelastic material in 3D WITHOUT reference stress and reference strain.
     *
     * Solution fields: [disp(dim), ...]
     * Auxiliary fields: [..., shear_modulus(1), bulk_modulus(1), power_law_reference_strain_rate(1),
     *                    power_law_reference_stress(1), power_law_exponent(1), viscous_strain(4), stress(4)]
     */
    static
    void f1v(const PylithInt dim,
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
             PylithScalar f1[]);

    /** f1 function for isotropic power-law viscoelastic material in 3D WITH reference stress and reference strain.
     *
     * Solution fields: [disp(dim), ...]
     * Auxiliary fields: [..., reference_stress(4), reference_strain(4), shear_modulus(1), bulk_modulus(1),
     *                    power_law_reference_strain_rate(1), power_law_reference_stress(1), power_law_exponent(1),
     *                    viscous_strain(4), stress(4)]
     */
    static
    void f1v_refstate(const PylithInt dim,
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
                      PylithScalar f1[]);

    /** Jf3_vu entry function for 3-D isotropic power-law viscoelasticity WITHOUT reference stress and
     * reference strain.
     *
     * Solution fields: [...]
     * Auxiliary fields: [..., shear_modulus(1), bulk_modulus(1), power_law_reference_strain_rate(1),
     *                    power_law_reference_stress(1), power_law_exponent(1), viscous_strain(4), stress(4)]
     */
    static
    void Jf3vu(const PylithInt dim,
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
               PylithScalar Jf3[]);

    /** Jf3_vu entry function for 3-D isotropic power-law viscoelasticity WITH reference stress and
     * reference strain.
     *
     * Solution fields: [...]
     * Auxiliary fields: [..., reference_stress(4), reference_strain(4), shear_modulus(1), bulk_modulus(1),
     *                    power_law_reference_strain_rate(1), power_law_reference_stress(1), power_law_exponent(1),
     *                    viscous_strain(4), stress(4)]
     */
    static
    void Jf3vu_refstate(const PylithInt dim,
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
                        PylithScalar Jf3[]);

    /** Calculate deviatoric stress for 3-D isotropic power-law
     * viscoelasticity WITHOUT reference stress and strain.
     *
     * Solution fields: [disp(dim)]
     * Auxiliary fields: [..., shear_modulus(1), bulk_modulus(1), power_law_reference_strain_rate(1),
     *                    power_law_reference_stress(1), power_law_exponent(1), viscous_strain(4), stress(4)]
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

    /** Calculate deviatoric stress for 3-D isotropic power-law
     * viscoelasticity WITH reference stress and strain.
     *
     * Solution fields: [disp(dim)]
     * Auxiliary fields: [..., reference_stress(4), reference_strain(4), shear_modulus(1), bulk_modulus(1),
     *                    power_law_reference_strain_rate(1), power_law_reference_stress(1), power_law_exponent(1),
     *                    viscous_strain(4), stress(4)]
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

    /** Update stress for 3-D isotropic power-law viscoelasticity WITHOUT reference stress and strain.
     *
     * Solution fields: [disp(dim)]
     * Auxiliary fields: [...]
     */
    static
    void updateStress(const PylithInt dim,
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

    /** Update stress for 3-D isotropic power-law viscoelasticity WITH reference stress and strain.
     *
     * Solution fields: [disp(dim)]
     * Auxiliary fields: [...]
     */
    static
    void updateStress_refstate(const PylithInt dim,
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

    /** Update viscous strain for 3-D isotropic power-law viscoelasticity WITHOUT reference stress and strain.
     *
     * Solution fields: [disp(dim)]
     * Auxiliary fields: [...]
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
                             PylithScalar visStrain[]);

    /** Update viscous strain for 3-D isotropic power-law viscoelasticity WITH reference stress and strain.
     *
     * Solution fields: [disp(dim)]
     * Auxiliary fields: [...]
     */
    static
    void updateViscousStrain_refstate(const PylithInt dim,
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
                                      PylithScalar visStrain[]);

    /** Calculate stress for 3-D isotropic power-law viscoelasticity WITHOUT a reference stress and strain.
     *
     * Used in outputing the stress field.
     *
     * Solution fields: [disp(dim)]
     * Auxiliary fields: [..., shear_modulus(1), bulk_modulus(1), power_law_reference_strain_rate(1),
     *                    power_law_reference_stress(1), power_law_exponent(1), viscous_strain(4), stress(4)]
     */
    static
    void cauchyStress(const PylithInt dim,
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
                      PylithScalar stressVector[]);

    /** Calculate stress for 3-D isotropic power-law viscoelasticity WITH a reference stress/strain.
     *
     * Used in outputing the stress field.
     *
     * Solution fields: [disp(dim)]
     * Auxiliary fields: [..., reference_stress(4), reference_strain(4), shear_modulus(1), bulk_modulus(1),
     *                    power_law_reference_strain_rate(1), power_law_reference_stress(1), power_law_exponent(1),
     *                    viscous_strain(4), stress(4)]
     */
    static
    void cauchyStress_refstate(const PylithInt dim,
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
                               PylithScalar stressVector[]);

}; // IsotropicPowerLaw3D
// ---------------------------------------------------------------------------------------------------------------------
/// Kernels for isotropic power-law viscoelastic material.
class pylith::fekernels::IsotropicPowerLawEffectiveStress {
    // PUBLIC MEMBERS //////////////////////////////////////////////////////////////////////////////////////////////////
public:

    /** Compute effective stress for power-law material, given an initial guess and the current parameters.
     *
     * Used to compute stress and viscous strain.
     *
     */

    static
    PylithScalar computeEffectiveStress(const PylithScalar j2InitialGuess,
                                        const PylithScalar stressScale,
                                        const PylithScalar ae,
                                        const PylithScalar b,
                                        const PylithScalar c,
                                        const PylithScalar d,
                                        const PylithScalar powerLawAlpha,
                                        const PylithScalar dt,
                                        const PylithScalar j2T,
                                        const PylithScalar powerLawExponent,
                                        const PylithScalar powerLawReferenceStrainRate,
                                        const PylithScalar powerLawReferenceStress);

private:

    /** Calculate effective stress function for a power-law viscoelastic material.
     *
     * Used for bracketing.
     *
     */

    static
    PylithScalar _effStressFunc(const PylithScalar j2Tpdt,
                                const PylithScalar ae,
                                const PylithScalar b,
                                const PylithScalar c,
                                const PylithScalar d,
                                const PylithScalar powerLawAlpha,
                                const PylithScalar dt,
                                const PylithScalar j2T,
                                const PylithScalar powerLawExponent,
                                const PylithScalar powerLawReferenceStrainRate,
                                const PylithScalar powerLawReferenceStress);

    /** Calculate effective stress function and its derivative for a power-law viscoelastic material.
     *
     * Used for root-finding.
     *
     */

    static
    void _effStressFuncDerivFunc(PylithScalar* func,
                                 PylithScalar* dfunc,
                                 const PylithScalar j2Tpdt,
                                 const PylithScalar ae,
                                 const PylithScalar b,
                                 const PylithScalar c,
                                 const PylithScalar d,
                                 const PylithScalar powerLawAlpha,
                                 const PylithScalar dt,
                                 const PylithScalar j2T,
                                 const PylithScalar powerLawExponent,
                                 const PylithScalar powerLawReferenceStrainRate,
                                 const PylithScalar powerLawReferenceStress);

    /** Bracket effective stress root.
     *
     * Used to place bounds on effective stress.
     *
     */

    static
    void _bracket(PylithScalar* px1,
                  PylithScalar* px2,
                  const PylithScalar ae,
                  const PylithScalar b,
                  const PylithScalar c,
                  const PylithScalar d,
                  const PylithScalar powerLawAlpha,
                  const PylithScalar dt,
                  const PylithScalar j2T,
                  const PylithScalar powerLawExponent,
                  const PylithScalar powerLawReferenceStrainRate,
                  const PylithScalar powerLawReferenceStress);

    /** Find zero of effective stress function using Newton's method with bisection.
     *
     * Used to find the effective stress.
     *
     */

    static
    PylithScalar _search(PylithScalar x1,
                         PylithScalar x2,
                         const PylithScalar ae,
                         const PylithScalar b,
                         const PylithScalar c,
                         const PylithScalar d,
                         const PylithScalar powerLawAlpha,
                         const PylithScalar dt,
                         const PylithScalar j2T,
                         const PylithScalar powerLawExponent,
                         const PylithScalar powerLawReferenceStrainRate,
                         const PylithScalar powerLawReferenceStress);

}; // IsotropicPowerLaw

#endif // pylith_fekernels_isotropicpowerlaw_hh

// End of file
