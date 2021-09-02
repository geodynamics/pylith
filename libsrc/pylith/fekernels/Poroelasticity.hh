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

/** @file libsrc/fekernels/Poroelasticity.hh
 *
 * Generic poroelasticity pointwise functions
 *
 */

#if !defined(pylith_fekernels_poroelasticity_hh)
#define pylith_fekernels_poroelasticity_hh

// Include directives ---------------------------------------------------
#include "fekernelsfwd.hh" // forward declarations

#include "pylith/utils/types.hh"

class pylith::fekernels::Poroelasticity
{
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
     * @param[out] Storage coefficient at constant strain.
     */
    // =============================================================================
    // Displacement
    // =============================================================================
    // ----------------------------------------------------------------------
    // f0u
    static void f0u(const PylithInt dim,
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

    // =============================================================================
    // Velocity
    // =============================================================================
    static void f0v_implicit(const PylithInt dim,
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

    static void f0v_explicit(const PylithInt dim,
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

    // =============================================================================
    // Volumetric Strain
    // =============================================================================
    // ----------------------------------------------------------------------
    // f0e function for isotropic linear Poroelasticity.
    static void f0e(const PylithInt dim,
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

    // =============================================================================
    // Time Derivative of Pressure
    // =============================================================================
    // ----------------------------------------------------------------------
    static void f0pdot(const PylithInt dim,
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

    // =============================================================================
    // Time Derivative of Volumetric Strain
    // =============================================================================
    // ----------------------------------------------------------------------
    static void f0edot(const PylithInt dim,
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

    /* -------------------------------------------------------------------------- */
    /*                           RHS Residuals                                    */
    /* -------------------------------------------------------------------------- */
    // Quasi-Static

    // =============================================================================
    // Displacement
    // =============================================================================

    // ----------------------------------------------------------------------
    // g0 function for displacement equation: g0u = v.
    static void g0u(const PylithInt dim,
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

    // ---------------------------------------------------------------------------------------------------------------------
    // g0v_grav - g0 function for generic poroelasticity terms ( + grav body forces).
    static void g0v_grav(const PylithInt dim,
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

    // ---------------------------------------------------------------------------------------------------------------------
    // g0v_bodyforce - g0 function for generic poroelasticity terms ( + body forces).
    static void g0v_bodyforce(const PylithInt dim,
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
    // g0v_gravbodyforce - g0 function for isotropic linear Poroelasticity with both gravity and body forces.
    static void g0v_grav_bodyforce(const PylithInt dim,
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

    // =============================================================================
    // Pressure
    // =============================================================================

    // ----------------------------------------------------------------------
    // g0p_sourceDensity - g0p function for generic poroelasticity terms (source density).
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
    // g0p_sourceDensity - g0p function for generic poroelasticity terms (source density).
    static void g0p_sourceDensity(const PylithInt dim,
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

    // ------------------------------------------------------------------------------
    // g0p function for isotropic linear Poroelasticity plane strain with source density, gravity
    static void g0p_sourceDensity_grav(const PylithInt dim,
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

    // ------------------------------------------------------------------------------
    // g0p function for isotropic linear Poroelasticity plane strain with source density, and body force.
    static void g0p_sourceDensity_body(const PylithInt dim,
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

    // ------------------------------------------------------------------------------
    // g0p function for isotropic linear Poroelasticity plane strain with source density, gravity, and body force.
    static void g0p_sourceDensity_grav_body(const PylithInt dim,
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

    /* -------------------------------------------------------------------------- */
    /*                           LHS Jacobian                                     */
    /* -------------------------------------------------------------------------- */

    // -----------------------------------------------------------------------------
    // Jg0ee - Jf0 function for isotropic linear poroelasticity plane strain.
    static void Jf0ee(const PylithInt dim,
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

    // -----------------------------------------------------------------------------
    // Jf1eu - Jf1 function for isotropic linear poroelasticity plane strain.
    static void Jf1eu(const PylithInt dim,
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
                      PylithScalar Jf1[]);

    // ---------------------------------------------------------------------------------------------------------------------
    // Jf0vu function for poroelasticity equation, quasistatic.
    static void Jf0vu_implicit(const PylithInt dim,
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

    // ---------------------------------------------------------------------------------------------------------------------
    // Jf0vv function for poroelasticity equation, quasistatic.
    static void Jf0vv_implicit(const PylithInt dim,
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

    // ---------------------------------------------------------------------------------------------------------------------
    // Jf0vv function for poroelasticity equation, dynamic
    static void Jf0vv_explicit(const PylithInt dim,
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

    // ---------------------------------------------------------------------------------------------------------------------
    // Jf0pdotp function for poroelasticity equation, quasistatic.
    static void Jf0pdotp(const PylithInt dim,
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

    // ---------------------------------------------------------------------------------------------------------------------
    // Jf0pdotpdot function for poroelasticity equation, quasistatic.
    static void Jf0pdotpdot(const PylithInt dim,
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

    // ---------------------------------------------------------------------------------------------------------------------
    // Jf0edote function for poroelasticity equation, quasistatic.
    static void Jf0edote(const PylithInt dim,
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

    // ---------------------------------------------------------------------------------------------------------------------
    // Jf0edotedot function for poroelasticity equation, quasistatic.
    static void Jf0edotedot(const PylithInt dim,
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

}; // Poroelasticity

// ---------------------------------------------------------------------------------------------------------------------

/// Kernels specific to poroelasticity plane strain.
class pylith::fekernels::PoroelasticityPlaneStrain
{
    // PUBLIC MEMBERS //////////////////////////////////////////////////////////////////////////////////////////////////
public:
    /** Calculate Cauchy strain for 2-D plane strain poroelasticity.
     *
     * Order of output components is xx, yy, zz, xy.
     *
     * Solution fields: [disp(dim)]
     */
    static void cauchyStrain(const PylithInt dim,
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
                             PylithScalar strain[]);

}; // PoroelasticityPlaneStrain

// ---------------------------------------------------------------------------------------------------------------------
/// Kernels specific to poroelasticity in 3D.
class pylith::fekernels::Poroelasticity3D
{
    // PUBLIC MEMBERS //////////////////////////////////////////////////////////////////////////////////////////////////
public:
    /** Calculate Cauchy strain for 3-D poroelasticity.
     *
     * Order of output components is xx, yy, zz, xy, yz, xz.
     *
     * Solution fields: [disp(dim)]
     */
    static void cauchyStrain(const PylithInt dim,
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
                             PylithScalar strain[]);

}; // Poroelasticity3D

#endif // pylith_fekernels_poroelasticity_hh

// End of file
