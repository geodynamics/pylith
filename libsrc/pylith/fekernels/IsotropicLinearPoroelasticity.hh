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

/** @file libsrc/fekernels/IsotropicLinearPoroelasticity.hh
 *
 * Kernels for linear poroelasticity plane strain.
 *
 * Solution fields: [disp(dim), pressure(1),trace_strain(1) ] (QS)
 * OR
 * Solution fields: [disp(dim), pressure(1),velocity(dim) ] (QS)
 *
 * Auxiliary fields:
 * -- numA : number of auxiliary fields
 ***** Required fields(govening equations) + option fields + required fields (rheology)
 * - 0: solid_density(1)
 * - 1: fluid_density(1)
 * - 2: fluid_viscosity(1)
 * - 3: porosity(1)
 *
 ** Optional fields
 * - +1: body_force(dim,optional)
 * - +1: gravity_field (dim, optional)
 * - +1: source_density(1,optional)
 * - +1: constant_pressure_source(1,optional)*
 * - +1: reference_stress(4,optional) (stress_xx, stress_yy, stress_xy, stress_zz)
 * - +1: reference_strain(4,optional) (strain_xx, strain_yy, strain_xy, strain_zz)
 *
 ** Rheological fields
 * - numA - 5: addShearModulus(1)
 * - numA - 4: addBulkModulus(1)
 * - numA - 3: addBiotCoefficient(1)
 * - numA - 2: addIsotropicPermeability(1)
 * - numA - 1: addFluidBulkModulus(1)
 *
 *
 *
 *
 * \int_V \vec{\phi}_u \cdot \left( \rho \frac{\partial \vec{v}(t)}{\partial t} \right) \, dV =
 *   \int_V \vec{\phi}_u \cdot \vec{f}(t) - \nabla \vec{\phi}_u : \tensor{\sigma}(\vec{u}) \, dV +
 *   \int_{S_\tau} \vec{\phi}_u \cdot \vec{\tau}(t) \, dS.
 *
 * ======================================================================
 */

#if !defined(pylith_fekernels_isotropiclinearporoelasticity_hh)
#define pylith_fekernels_isotropiclinearporoelasticity_hh

// Include directives ---------------------------------------------------
#include "fekernelsfwd.hh" // forward declarations

#include "pylith/utils/types.hh"
// ---------------------------------------------------------------------------------------------------------------------
/// Kernels specific to isotropic, linearl poroelasticity in Plane Strain.
class pylith::fekernels::IsotropicLinearPoroelasticityPlaneStrain {
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

    // ================================= MMS =======================================

    // ----------------------------------------------------------------------
    // f0u function for quadratic space and linear time MMS.
    static void f0_mms_ql_u(const PylithInt dim,
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

    // ----------------------------------------------------------------------
    // f0p function for quadratic space and linear time MMS.
    static void f0_mms_ql_p(const PylithInt dim,
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

    // ----------------------------------------------------------------------
    // f0u function for quadratic space and trigonometric time MMS.
    static void f0_mms_qt_u(const PylithInt dim,
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

    // ----------------------------------------------------------------------
    // f0p function for quadratic space and trigonometric time MMS.
    static void f0_mms_qt_p(const PylithInt dim,
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

    // ----------------------------------------------------------------------
    // f0u function for trigonometric space and linear time MMS.
    static void f0_mms_tl_u(const PylithInt dim,
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

    // ----------------------------------------------------------------------
    // f0p function for trigonometric space and linear time MMS.
    static void f0_mms_tl_p(const PylithInt dim,
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

    // ================================= STD =======================================

    // ================================= LHS =======================================

    // ----------------------------------------------------------------------
    // f0p function for generic poroelasticity terms (source density).
    static void f0p_explicit(const PylithInt dim,
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

    // ----------------------------------------------------------------------
    // f0p function for generic poroelasticity terms (source density).
    static void f0p_implicit(const PylithInt dim,
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

    // ----------------------------------------------------------------------
    // f0p function for generic poroelasticity terms (source density).
    static void f0p_implicit_source(const PylithInt dim,
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

    // ----------------------------------------------------------------------
    // f0p function for generic poroelasticity terms (source density).
    static void f0p_implicit_source_body(const PylithInt dim,
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

    // ----------------------------------------------------------------------
    // f0p function for generic poroelasticity terms (source density).
    static void f0p_implicit_source_grav(const PylithInt dim,
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

    // ----------------------------------------------------------------------
    // f0p function for generic poroelasticity terms (source density).
    static void f0p_implicit_source_grav_body(const PylithInt dim,
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

    // -----------------------------------------------------------------------------
    /** f1u function for isotropic linear poroelasticity plane strain WITHOUT reference stress and reference strain.
     * Quasi - Static Case
     * Solution fields: [disp(dim), pore_pres(dim), vel(dim, optional)]
     * Auxiliary fields: [density(1), shear_modulus(1), bulk_modulus(1), ...]
     */
    static void f1u(const PylithInt dim,
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

    // -----------------------------------------------------------------------------
    /** f1u function for isotropic linear poroelasticity plane strain WITH reference stress and reference strain.
     *
     * Solution fields: [disp(dim), pres(dim), vel(dim, optional)]
     * Auxiliary fields: [density(1), shear_modulus(1), bulk_modulus(1), ..., refstress(4), refstrain(4)]
     */
    static void f1u_refstate(const PylithInt dim,
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

    // -----------------------------------------------------------------------------
    /** f1p / darcy flow / without gravity
     *
     * Solution fields: [...]
     * Auxiliary fields: [density(1), shear_modulus(1), bulk_modulus(1), other poroelastic related param ...]
     */
    static void f1p(const PylithInt dim,
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

    // -----------------------------------------------------------------------------
    /** f1p / darcy flow / without gravity, tensor permeability
     *
     * Solution fields: [...]
     * Auxiliary fields: [density(1), shear_modulus(1), bulk_modulus(1), other poroelastic related param ...]
     */
    static void f1p_tensor_permeability(const PylithInt dim,
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

    // -----------------------------------------------------------------------------
    /** f1p / darcy flow / including body forces
     *
     * Solution fields: [...]
     * Auxiliary fields: [density(1), shear_modulus(1), bulk_modulus(1), other poroelastic related param ...]
     */
    static void f1p_body(const PylithInt dim,
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

    // -----------------------------------------------------------------------------
    /** f1p / darcy flow / including body forces, tensor permeability
     *
     * Solution fields: [...]
     * Auxiliary fields: [density(1), shear_modulus(1), bulk_modulus(1), other poroelastic related param ...]
     */
    static void f1p_body_tensor_permeability(const PylithInt dim,
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

    // -----------------------------------------------------------------------------
    /** f1p / darcy flow / including gravity
     *
     * Solution fields: [...]
     * Auxiliary fields: [density(1), shear_modulus(1), bulk_modulus(1), other poroelastic related param ...]
     */
    static void f1p_gravity(const PylithInt dim,
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

    // -----------------------------------------------------------------------------
    /** f1p / darcy flow / including gravity, tensor permeability
     *
     * Solution fields: [...]
     * Auxiliary fields: [density(1), shear_modulus(1), bulk_modulus(1), other poroelastic related param ...]
     */
    static void f1p_gravity_tensor_permeability(const PylithInt dim,
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

    // -----------------------------------------------------------------------------
    /** f1p / darcy flow / including body forces and gravity
     *
     * Solution fields: [...]
     * Auxiliary fields: [density(1), shear_modulus(1), bulk_modulus(1), other poroelastic related param ...]
     */
    static void f1p_body_gravity(const PylithInt dim,
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

    // -----------------------------------------------------------------------------
    /** f1p / darcy flow / including body forces and gravity, tensor permeability
     *
     * Solution fields: [...]
     * Auxiliary fields: [density(1), shear_modulus(1), bulk_modulus(1), other poroelastic related param ...]
     */
    static void f1p_body_gravity_tensor_permeability(const PylithInt dim,
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

    // -----------------------------------------------------------------------------
    /** f1p / darcy flow / including gravity
     *
     * Solution fields: [...]
     * Auxiliary fields: [density(1), shear_modulus(1), bulk_modulus(1), other poroelastic related param ...]
     */
    static void f1p_cpsource(const PylithInt dim,
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

    // -----------------------------------------------------------------------------
    /** f1p / darcy flow / including gravity, tensor permeability
     *
     * Solution fields: [...]
     * Auxiliary fields: [density(1), shear_modulus(1), bulk_modulus(1), other poroelastic related param ...]
     */
    static void f1p_cpsource_tensor_permeability(const PylithInt dim,
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

    // -----------------------------------------------------------------------------
    /** f1p / darcy flow / including gravity
     *
     * Solution fields: [...]
     * Auxiliary fields: [density(1), shear_modulus(1), bulk_modulus(1), other poroelastic related param ...]
     */
    static void f1p_body_cpsource(const PylithInt dim,
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

    // -----------------------------------------------------------------------------
    /** f1p / darcy flow / including gravity, tensor permeability
     *
     * Solution fields: [...]
     * Auxiliary fields: [density(1), shear_modulus(1), bulk_modulus(1), other poroelastic related param ...]
     */
    static void f1p_body_cpsource_tensor_permeability(const PylithInt dim,
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

    // -----------------------------------------------------------------------------
    /** f1p / darcy flow / including gravity
     *
     * Solution fields: [...]
     * Auxiliary fields: [density(1), shear_modulus(1), bulk_modulus(1), other poroelastic related param ...]
     */
    static void f1p_gravity_cpsource(const PylithInt dim,
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

    // -----------------------------------------------------------------------------
    /** f1p / darcy flow / including gravity, tensor permeability
     *
     * Solution fields: [...]
     * Auxiliary fields: [density(1), shear_modulus(1), bulk_modulus(1), other poroelastic related param ...]
     */
    static void f1p_gravity_cpsource_tensor_permeability(const PylithInt dim,
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

    // -----------------------------------------------------------------------------
    /** f1p / darcy flow / including gravity
     *
     * Solution fields: [...]
     * Auxiliary fields: [density(1), shear_modulus(1), bulk_modulus(1), other poroelastic related param ...]
     */
    static void f1p_body_gravity_cpsource(const PylithInt dim,
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

    // -----------------------------------------------------------------------------
    /** f1p / darcy flow / including gravity, tensor permeability
     *
     * Solution fields: [...]
     * Auxiliary fields: [density(1), shear_modulus(1), bulk_modulus(1), other poroelastic related param ...]
     */
    static void f1p_body_gravity_cpsource_tensor_permeability(const PylithInt dim,
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

    // =========================== LHS Jacobian ============================

    // ----------------------------------------------------------------------
    /* Jf3_uu entry function for isotropic linear poroelasticity.
     */
    static void Jf3uu(const PylithInt dim,
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
                      PylithScalar Jf3[]);

    // ----------------------------------------------------------------------
    /** Jf2_up entry function for isotropic linear poroelasticity.
     * Solution fields: [...]
     * Auxiliary fields: [density(1), shear_modulus(1), bulk_modulus(1), other poroelastic related param ...]
     */
    static void Jf2up(const PylithInt dim,
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
                      PylithScalar Jf2[]);

    // -----------------------------------------------------------------------------
    // Jf2ue function for isotropic linear poroelasticity.
    static void Jf2ue(const PylithInt dim,
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
                      PylithScalar Jf2[]);

    // ----------------------------------------------------------------------
    /** Jf3pp entry function for isotropic linear poroelasticity.
     *
     * Solution fields: [...]
     * Auxiliary fields: [density(1), shear_modulus(1), bulk_modulus(1), other poroelastic related param ...]
     */
    static void Jf3pp(const PylithInt dim,
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
                      PylithScalar Jf3[]);

    // ----------------------------------------------------------------------
    /** Jf3pp entry function for isotropic linear poroelasticity.
     *
     * Solution fields: [...]
     * Auxiliary fields: [density(1), shear_modulus(1), bulk_modulus(1), other poroelastic related param ...]
     */
    static void Jf3pp_tensor_permeability(const PylithInt dim,
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
                                          PylithScalar Jf3[]);

    // ----------------------------------------------------------------------
    /** Jf0_pp entry function for isotropic linear poroelasticity.
     *
     * Solution fields: [...]
     * Auxiliary fields: [density(1), shear_modulus(1), bulk_modulus(1), other poroelastic related param ...]
     */
    static void Jf0pp(const PylithInt dim,
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
                      PylithScalar Jf0[]);

    // ----------------------------------------------------------------------
    /** Jf0_pe entry function for isotropic linear poroelasticity.
     *
     * Solution fields: [...]
     * Auxiliary fields: [density(1), shear_modulus(1), bulk_modulus(1), other poroelastic related param ...]
     */
    static void Jf0pe(const PylithInt dim,
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
                      PylithScalar Jf0[]);

    // ----------------------------------------------------------------------
    /** Jf0_ppdot entry function for isotropic linear poroelasticity.
     *
     * Solution fields: [...]
     * Auxiliary fields: [density(1), shear_modulus(1), bulk_modulus(1), other poroelastic related param ...]
     */
    static
    void Jf0ppdot(const PylithInt dim,
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
                  PylithScalar Jf0[]);

    // ----------------------------------------------------------------------
    /** Jf0_pedot entry function for isotropic linear poroelasticity.
     *
     * Solution fields: [...]
     * Auxiliary fields: [density(1), shear_modulus(1), bulk_modulus(1), other poroelastic related param ...]
     */
    static
    void Jf0pedot(const PylithInt dim,
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
                  PylithScalar Jf0[]);

    // ============================== RHS Residual =================================

    // ----------------------------------------------------------------------
    // g0p function for generic poroelasticity terms.
    // \left( \alpha \frac{\partial \epsilon_{v}}{\partial t} + \frac{1}{M} \frac{\partial p_{f}}{\partial t} \right)
    // \frac{\partial \epsilon)_{v}}{\partial t} = \frac{\nabla \cdot \vec{u}}{dt}
    static void g0p(const PylithInt dim,
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
                    PylithScalar g0p[]);

    // ----------------------------------------------------------------------
    // g0p function for generic poroelasticity terms.
    // \left( \alpha \frac{\partial \epsilon_{v}}{\partial t} + \frac{1}{M} \frac{\partial p_{f}}{\partial t} \right)
    // \frac{\partial \epsilon)_{v}}{\partial t} = \frac{\nabla \cdot \vec{u}}{dt}
    static void g0p_source(const PylithInt dim,
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

    // ----------------------------------------------------------------------
    // g0p function for generic poroelasticity terms.
    // \left( \alpha \frac{\partial \epsilon_{v}}{\partial t} + \frac{1}{M} \frac{\partial p_{f}}{\partial t} \right)
    // \frac{\partial \epsilon)_{v}}{\partial t} = \frac{\nabla \cdot \vec{u}}{dt}
    static void g0p_source_body(const PylithInt dim,
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

    // ----------------------------------------------------------------------
    // g0p function for generic poroelasticity terms.
    // \left( \alpha \frac{\partial \epsilon_{v}}{\partial t} + \frac{1}{M} \frac{\partial p_{f}}{\partial t} \right)
    // \frac{\partial \epsilon)_{v}}{\partial t} = \frac{\nabla \cdot \vec{u}}{dt}
    static void g0p_source_grav(const PylithInt dim,
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

    // ----------------------------------------------------------------------
    // g0p function for generic poroelasticity terms.
    // \left( \alpha \frac{\partial \epsilon_{v}}{\partial t} + \frac{1}{M} \frac{\partial p_{f}}{\partial t} \right)
    // \frac{\partial \epsilon)_{v}}{\partial t} = \frac{\nabla \cdot \vec{u}}{dt}
    static void g0p_source_grav_body(const PylithInt dim,
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

    // -----------------------------------------------------------------------------
    /** g1p / darcy flow / including gravity
     *
     * Solution fields: [...]
     * Auxiliary fields: [density(1), shear_modulus(1), bulk_modulus(1), other poroelastic related param ...]
     */
    static void g1p_gravity(const PylithInt dim,
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

    // -----------------------------------------------------------------------------
    /** g1p / darcy flow / including gravity, tensor permeability
     *
     * Solution fields: [...]
     * Auxiliary fields: [density(1), shear_modulus(1), bulk_modulus(1), other poroelastic related param ...]
     */
    static void g1p_gravity_tensor_permeability(const PylithInt dim,
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

    // -----------------------------------------------------------------------------
    /** g1p / darcy flow / without gravity
     *
     * Solution fields: [...]
     * Auxiliary fields: [density(1), shear_modulus(1), bulk_modulus(1), other poroelastic related param ...]
     */
    static void g1p(const PylithInt dim,
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

    // -----------------------------------------------------------------------------
    /** g1p / darcy flow / without gravity, tensor permeability
     *
     * Solution fields: [...]
     * Auxiliary fields: [density(1), shear_modulus(1), bulk_modulus(1), other poroelastic related param ...]
     */
    static void g1p_tensor_permeability(const PylithInt dim,
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

    // -----------------------------------------------------------------------------
    /** g1v function for isotropic linear poroelasticity plane strain WITHOUT reference stress and reference strain.
     * Dynamic Case
     * Solution fields: [disp(dim), pore_pres(dim), vel(dim, optional)]
     * Auxiliary fields: [density(1), shear_modulus(1), bulk_modulus(1), ...]
     */
    static void g1v(const PylithInt dim,
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

    // -----------------------------------------------------------------------------
    /** g1v function for isotropic linear poroelasticity plane strain WITH reference stress and reference strain.
     *
     * Solution fields: [disp(dim), pres(dim), vel(dim, optional)]
     * Auxiliary fields: [density(1), shear_modulus(1), bulk_modulus(1), ..., refstress(4), refstrain(4)]
     */
    static void g1v_refstate(const PylithInt dim,
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

    // -----------------------------------------------------------------------------
    /** Calculate stress for isotropic linear elasticity WITHOUT a reference stress and strain.
     *
     * Used to output the stress field.
     *
     * Solution fields: [disp(dim)]
     * Auxiliary fields: [..., shear_modulus(1), bulk_modulus(1)]
     */
    static void cauchyStress(const PylithInt dim,
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

    /** Calculate stress for isotropic linear elasticity WITH a reference stress/strain.
     *
     * Used to output the stress field.
     *
     * Solution fields: [disp(dim)]
     * Auxiliary fields: [..., refstress(6), refstrain(6), shear_modulus(1), bulk_modulus(1)]
     */
    static void cauchyStress_refstate(const PylithInt dim,
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

    // ========================== Update Kernels ===================================

    // ---------------------------------------------------------------------------------------------------------------------
    /* Update porosity for a linear poroelastic material.
     */
    static
    void updatePorosity(const PylithInt dim,
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
                        PylithScalar porosity[]);

}; // IsotropicLinearPoroelasticityPlaneStrain

// ---------------------------------------------------------------------------------------------------------------------
/// Kernels specific to isotropic, linearl poroelasticity in 3D.
class pylith::fekernels::IsotropicLinearPoroelasticity3D {
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

    // ================================= MMS =======================================

    // ----------------------------------------------------------------------
    // f0u function for quadratic space and linear time MMS.
    static void f0_mms_ql_u(const PylithInt dim,
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

    // ----------------------------------------------------------------------
    // f0p function for quadratic space and linear time MMS.
    static void f0_mms_ql_p(const PylithInt dim,
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

    // ----------------------------------------------------------------------
    // f0u function for quadratic space and trigonometric time MMS.
    static void f0_mms_qt_u(const PylithInt dim,
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

    // ----------------------------------------------------------------------
    // f0p function for quadratic space and trigonometric time MMS.
    static void f0_mms_qt_p(const PylithInt dim,
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

    // ----------------------------------------------------------------------
    // f0u function for trigonometric space and linear time MMS.
    static void f0_mms_tl_u(const PylithInt dim,
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

    // ----------------------------------------------------------------------
    // f0p function for trigonometric space and linear time MMS.
    static void f0_mms_tl_p(const PylithInt dim,
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

    // ================================= STD =======================================

    // ================================= LHS =======================================

    // ----------------------------------------------------------------------
    // f0p function for generic poroelasticity terms (source density).
    static void f0p_explicit(const PylithInt dim,
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

    // ----------------------------------------------------------------------
    // f0p function for generic poroelasticity terms (source density).
    static void f0p_implicit(const PylithInt dim,
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

    // ----------------------------------------------------------------------
    // f0p function for generic poroelasticity terms (source density).
    static void f0p_implicit_source(const PylithInt dim,
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

    // ----------------------------------------------------------------------
    // f0p function for generic poroelasticity terms (source density).
    static void f0p_implicit_source_body(const PylithInt dim,
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

    // ----------------------------------------------------------------------
    // f0p function for generic poroelasticity terms (source density).
    static void f0p_implicit_source_grav(const PylithInt dim,
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

    // ----------------------------------------------------------------------
    // f0p function for generic poroelasticity terms (source density).
    static void f0p_implicit_source_grav_body(const PylithInt dim,
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

    // -----------------------------------------------------------------------------
    /** f1u function for isotropic linear poroelasticity plane strain WITHOUT reference stress and reference strain.
     * Quasi - Static Case
     * Solution fields: [disp(dim), pore_pres(dim), vel(dim, optional)]
     * Auxiliary fields: [density(1), shear_modulus(1), bulk_modulus(1), ...]
     */
    static void f1u(const PylithInt dim,
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

    // -----------------------------------------------------------------------------
    /** f1u function for isotropic linear poroelasticity plane strain WITH reference stress and reference strain.
     *
     * Solution fields: [disp(dim), pres(dim), vel(dim, optional)]
     * Auxiliary fields: [density(1), shear_modulus(1), bulk_modulus(1), ..., refstress(4), refstrain(4)]
     */
    static void f1u_refstate(const PylithInt dim,
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

    // -----------------------------------------------------------------------------
    /** f1p / darcy flow / without gravity
     *
     * Solution fields: [...]
     * Auxiliary fields: [density(1), shear_modulus(1), bulk_modulus(1), other poroelastic related param ...]
     */
    static void f1p(const PylithInt dim,
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

    // -----------------------------------------------------------------------------
    /** f1p / darcy flow / without gravity, tensor permeability
     *
     * Solution fields: [...]
     * Auxiliary fields: [density(1), shear_modulus(1), bulk_modulus(1), other poroelastic related param ...]
     */
    static void f1p_tensor_permeability(const PylithInt dim,
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

    // -----------------------------------------------------------------------------
    /** f1p / darcy flow / including body forces
     *
     * Solution fields: [...]
     * Auxiliary fields: [density(1), shear_modulus(1), bulk_modulus(1), other poroelastic related param ...]
     */
    static void f1p_body(const PylithInt dim,
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

    // -----------------------------------------------------------------------------
    /** f1p / darcy flow / including body forces, tensor permeability
     *
     * Solution fields: [...]
     * Auxiliary fields: [density(1), shear_modulus(1), bulk_modulus(1), other poroelastic related param ...]
     */
    static void f1p_body_tensor_permeability(const PylithInt dim,
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

    // -----------------------------------------------------------------------------
    /** f1p / darcy flow / including gravity
     *
     * Solution fields: [...]
     * Auxiliary fields: [density(1), shear_modulus(1), bulk_modulus(1), other poroelastic related param ...]
     */
    static void f1p_gravity(const PylithInt dim,
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

    // -----------------------------------------------------------------------------
    /** f1p / darcy flow / including gravity, tensor permeability
     *
     * Solution fields: [...]
     * Auxiliary fields: [density(1), shear_modulus(1), bulk_modulus(1), other poroelastic related param ...]
     */
    static void f1p_gravity_tensor_permeability(const PylithInt dim,
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

    // -----------------------------------------------------------------------------
    /** f1p / darcy flow / including body forces and gravity
     *
     * Solution fields: [...]
     * Auxiliary fields: [density(1), shear_modulus(1), bulk_modulus(1), other poroelastic related param ...]
     */
    static void f1p_body_gravity(const PylithInt dim,
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

    // -----------------------------------------------------------------------------
    /** f1p / darcy flow / including body forces and gravity, tensor permeability
     *
     * Solution fields: [...]
     * Auxiliary fields: [density(1), shear_modulus(1), bulk_modulus(1), other poroelastic related param ...]
     */
    static void f1p_body_gravity_tensor_permeability(const PylithInt dim,
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

    // -----------------------------------------------------------------------------
    /** f1p / darcy flow / including gravity
     *
     * Solution fields: [...]
     * Auxiliary fields: [density(1), shear_modulus(1), bulk_modulus(1), other poroelastic related param ...]
     */
    static void f1p_cpsource(const PylithInt dim,
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

    // -----------------------------------------------------------------------------
    /** f1p / darcy flow / including gravity, tensor permeability
     *
     * Solution fields: [...]
     * Auxiliary fields: [density(1), shear_modulus(1), bulk_modulus(1), other poroelastic related param ...]
     */
    static void f1p_cpsource_tensor_permeability(const PylithInt dim,
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

    // -----------------------------------------------------------------------------
    /** f1p / darcy flow / including gravity
     *
     * Solution fields: [...]
     * Auxiliary fields: [density(1), shear_modulus(1), bulk_modulus(1), other poroelastic related param ...]
     */
    static void f1p_body_cpsource(const PylithInt dim,
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

    // -----------------------------------------------------------------------------
    /** f1p / darcy flow / including gravity, tensor permeability
     *
     * Solution fields: [...]
     * Auxiliary fields: [density(1), shear_modulus(1), bulk_modulus(1), other poroelastic related param ...]
     */
    static void f1p_body_cpsource_tensor_permeability(const PylithInt dim,
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

    // -----------------------------------------------------------------------------
    /** f1p / darcy flow / including gravity
     *
     * Solution fields: [...]
     * Auxiliary fields: [density(1), shear_modulus(1), bulk_modulus(1), other poroelastic related param ...]
     */
    static void f1p_gravity_cpsource(const PylithInt dim,
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

    // -----------------------------------------------------------------------------
    /** f1p / darcy flow / including gravity, tensor permeability
     *
     * Solution fields: [...]
     * Auxiliary fields: [density(1), shear_modulus(1), bulk_modulus(1), other poroelastic related param ...]
     */
    static void f1p_gravity_cpsource_tensor_permeability(const PylithInt dim,
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

    // -----------------------------------------------------------------------------
    /** f1p / darcy flow / including gravity
     *
     * Solution fields: [...]
     * Auxiliary fields: [density(1), shear_modulus(1), bulk_modulus(1), other poroelastic related param ...]
     */
    static void f1p_body_gravity_cpsource(const PylithInt dim,
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

    // -----------------------------------------------------------------------------
    /** f1p / darcy flow / including gravity, tensor permeability
     *
     * Solution fields: [...]
     * Auxiliary fields: [density(1), shear_modulus(1), bulk_modulus(1), other poroelastic related param ...]
     */
    static void f1p_body_gravity_cpsource_tensor_permeability(const PylithInt dim,
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

    // =========================== LHS Jacobian ============================

    // ----------------------------------------------------------------------
    /* Jf3_uu entry function for isotropic linear poroelasticity.
     */
    static void Jf3uu(const PylithInt dim,
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
                      PylithScalar Jf3[]);

    // ----------------------------------------------------------------------
    /** Jf2_up entry function for isotropic linear poroelasticity.
     * Solution fields: [...]
     * Auxiliary fields: [density(1), shear_modulus(1), bulk_modulus(1), other poroelastic related param ...]
     */
    static void Jf2up(const PylithInt dim,
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
                      PylithScalar Jf2[]);

    // -----------------------------------------------------------------------------
    // Jf2ue function for isotropic linear poroelasticity.
    static void Jf2ue(const PylithInt dim,
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
                      PylithScalar Jf2[]);

    // ----------------------------------------------------------------------
    /** Jf3pp entry function for isotropic linear poroelasticity.
     *
     * Solution fields: [...]
     * Auxiliary fields: [density(1), shear_modulus(1), bulk_modulus(1), other poroelastic related param ...]
     */
    static void Jf3pp(const PylithInt dim,
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
                      PylithScalar Jf3[]);

    // ----------------------------------------------------------------------
    /** Jf3pp entry function for isotropic linear poroelasticity.
     *
     * Solution fields: [...]
     * Auxiliary fields: [density(1), shear_modulus(1), bulk_modulus(1), other poroelastic related param ...]
     */
    static void Jf3pp_tensor_permeability(const PylithInt dim,
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
                                          PylithScalar Jf3[]);

    // ----------------------------------------------------------------------
    /** Jf0_pp entry function for isotropic linear poroelasticity.
     *
     * Solution fields: [...]
     * Auxiliary fields: [density(1), shear_modulus(1), bulk_modulus(1), other poroelastic related param ...]
     */
    static void Jf0pp(const PylithInt dim,
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
                      PylithScalar Jf0[]);

    // ----------------------------------------------------------------------
    /** Jf0_pe entry function for isotropic linear poroelasticity.
     *
     * Solution fields: [...]
     * Auxiliary fields: [density(1), shear_modulus(1), bulk_modulus(1), other poroelastic related param ...]
     */
    static void Jf0pe(const PylithInt dim,
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
                      PylithScalar Jf0[]);

    // ----------------------------------------------------------------------
    /** Jf0_ppdot entry function for isotropic linear poroelasticity.
     *
     * Solution fields: [...]
     * Auxiliary fields: [density(1), shear_modulus(1), bulk_modulus(1), other poroelastic related param ...]
     */
    static
    void Jf0ppdot(const PylithInt dim,
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
                  PylithScalar Jf0[]);

    // ----------------------------------------------------------------------
    /** Jf0_pedot entry function for isotropic linear poroelasticity.
     *
     * Solution fields: [...]
     * Auxiliary fields: [density(1), shear_modulus(1), bulk_modulus(1), other poroelastic related param ...]
     */
    static
    void Jf0pedot(const PylithInt dim,
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
                  PylithScalar Jf0[]);

    // ============================== RHS Residual =================================

    // ----------------------------------------------------------------------
    // g0p function for generic poroelasticity terms.
    // \left( \alpha \frac{\partial \epsilon_{v}}{\partial t} + \frac{1}{M} \frac{\partial p_{f}}{\partial t} \right)
    // \frac{\partial \epsilon)_{v}}{\partial t} = \frac{\nabla \cdot \vec{u}}{dt}
    static void g0p(const PylithInt dim,
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

    // ----------------------------------------------------------------------
    // g0p function for generic poroelasticity terms.
    // \left( \alpha \frac{\partial \epsilon_{v}}{\partial t} + \frac{1}{M} \frac{\partial p_{f}}{\partial t} \right)
    // \frac{\partial \epsilon)_{v}}{\partial t} = \frac{\nabla \cdot \vec{u}}{dt}
    static void g0p_source(const PylithInt dim,
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

    // ----------------------------------------------------------------------
    // g0p function for generic poroelasticity terms.
    // \left( \alpha \frac{\partial \epsilon_{v}}{\partial t} + \frac{1}{M} \frac{\partial p_{f}}{\partial t} \right)
    // \frac{\partial \epsilon)_{v}}{\partial t} = \frac{\nabla \cdot \vec{u}}{dt}
    static void g0p_source_body(const PylithInt dim,
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

    // ----------------------------------------------------------------------
    // g0p function for generic poroelasticity terms.
    // \left( \alpha \frac{\partial \epsilon_{v}}{\partial t} + \frac{1}{M} \frac{\partial p_{f}}{\partial t} \right)
    // \frac{\partial \epsilon)_{v}}{\partial t} = \frac{\nabla \cdot \vec{u}}{dt}
    static void g0p_source_grav(const PylithInt dim,
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

    // ----------------------------------------------------------------------
    // g0p function for generic poroelasticity terms.
    // \left( \alpha \frac{\partial \epsilon_{v}}{\partial t} + \frac{1}{M} \frac{\partial p_{f}}{\partial t} \right)
    // \frac{\partial \epsilon)_{v}}{\partial t} = \frac{\nabla \cdot \vec{u}}{dt}
    static void g0p_source_grav_body(const PylithInt dim,
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

    // -----------------------------------------------------------------------------
    /** g1p / darcy flow / including gravity
     *
     * Solution fields: [...]
     * Auxiliary fields: [density(1), shear_modulus(1), bulk_modulus(1), other poroelastic related param ...]
     */
    static void g1p_gravity(const PylithInt dim,
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

    // -----------------------------------------------------------------------------
    /** g1p / darcy flow / including gravity, tensor permeability
     *
     * Solution fields: [...]
     * Auxiliary fields: [density(1), shear_modulus(1), bulk_modulus(1), other poroelastic related param ...]
     */
    static void g1p_gravity_tensor_permeability(const PylithInt dim,
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

    // -----------------------------------------------------------------------------
    /** g1p / darcy flow / without gravity
     *
     * Solution fields: [...]
     * Auxiliary fields: [density(1), shear_modulus(1), bulk_modulus(1), other poroelastic related param ...]
     */
    static void g1p(const PylithInt dim,
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

    // -----------------------------------------------------------------------------
    /** g1p / darcy flow / without gravity, tensor permeability
     *
     * Solution fields: [...]
     * Auxiliary fields: [density(1), shear_modulus(1), bulk_modulus(1), other poroelastic related param ...]
     */
    static void g1p_tensor_permeability(const PylithInt dim,
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

    // -----------------------------------------------------------------------------
    /** g1v function for isotropic linear poroelasticity plane strain WITHOUT reference stress and reference strain.
     * Dynamic Case
     * Solution fields: [disp(dim), pore_pres(dim), vel(dim, optional)]
     * Auxiliary fields: [density(1), shear_modulus(1), bulk_modulus(1), ...]
     */
    static void g1v(const PylithInt dim,
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

    // -----------------------------------------------------------------------------
    /** g1v function for isotropic linear poroelasticity plane strain WITH reference stress and reference strain.
     *
     * Solution fields: [disp(dim), pres(dim), vel(dim, optional)]
     * Auxiliary fields: [density(1), shear_modulus(1), bulk_modulus(1), ..., refstress(4), refstrain(4)]
     */
    static void g1v_refstate(const PylithInt dim,
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

    // -----------------------------------------------------------------------------
    /** Calculate stress for isotropic linear elasticity WITHOUT a reference stress and strain.
     *
     * Used to output the stress field.
     *
     * Solution fields: [disp(dim)]
     * Auxiliary fields: [..., shear_modulus(1), bulk_modulus(1)]
     */
    static void cauchyStress(const PylithInt dim,
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

    /** Calculate stress for isotropic linear elasticity WITH a reference stress/strain.
     *
     * Used to output the stress field.
     *
     * Solution fields: [disp(dim)]
     * Auxiliary fields: [..., refstress(6), refstrain(6), shear_modulus(1), bulk_modulus(1)]
     */
    static void cauchyStress_refstate(const PylithInt dim,
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

    // ========================== Update Kernels ===================================

    // ---------------------------------------------------------------------------------------------------------------------
    /* Update porosity for a linear poroelastic material.
     */
    static
    void updatePorosity(const PylithInt dim,
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
                        PylithScalar porosity[]);

}; // IsotropicLinearPoroelasticity3D

#endif // pylith_fekernels_isotropiclinearporoelasticity_hh

// End of file
