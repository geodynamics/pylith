
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

#include <portinfo>

#include "pylith/fekernels/Poroelasticity.hh"

#include <cassert> // USES assert()

// =====================================================================================================================
// Generic poroelasticity kernels for inertia and body forces.
// =====================================================================================================================

/* -------------------------------------------------------------------------- */
/*                           LHS Residuals                                    */
/* -------------------------------------------------------------------------- */

// ---------------------------------------------------------------------------------------------------------------------
// f0v function for poroelasticity equation, explicit time stepping, dynamic.
void
pylith::fekernels::Poroelasticity::f0v(const PylithInt dim,
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
                                       PylithScalar f0[]) {
    // Incoming solution fields.
    const PylithInt i_velocity = 2;

    // Incoming auxiliary fields.
    const PylithInt i_solid_density = 0;
    const PylithInt i_fluid_density = 1;
    const PylithInt i_porosity = 3;

    const PylithScalar bulkDensity = (1 - a[aOff[i_porosity]]) * a[aOff[i_solid_density]] + a[aOff[i_porosity]] * a[aOff[i_fluid_density]];
    const PylithScalar* velocity_t = &s_t[sOff[i_velocity]]; // acceleration

    for (PylithInt i = 0; i < dim; ++i) {
        f0[i] += velocity_t[i] * bulkDensity;
    } // for
} // f0v


// =============================================================================
// Volumetric Strain
// =============================================================================
// ----------------------------------------------------------------------
// f0e function for isotropic linear Poroelasticity plane strain.
void
pylith::fekernels::Poroelasticity::f0e(const PylithInt dim,
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
                                       PylithScalar f0[]) {
    // Incoming solution fields.
    const PylithInt i_displacement = 0;
    const PylithInt i_trace_strain = 2;

    // Incoming auxiliary fields.

    const PylithScalar* displacement = &s[sOff[i_displacement]];
    const PylithScalar* displacement_x = &s_x[sOff_x[i_displacement]];
    const PylithScalar trace_strain = s[sOff[i_trace_strain]];

    for (PylithInt d = 0; d < dim; ++d) {
        f0[0] += displacement_x[d*dim+d];
    }
    f0[0] -= trace_strain;
} // f0e


/* -------------------------------------------------------------------------- */
/*                           RHS Residuals                                    */
/* -------------------------------------------------------------------------- */
// Quasi-Static

// =============================================================================
// Displacement
// =============================================================================

// ----------------------------------------------------------------------
// g0 function for displacement equation: g0u = v.
void
pylith::fekernels::Poroelasticity::g0u(const PylithInt dim,
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
                                       PylithScalar g0[]) {
    const PylithInt i_velocity = 2;
    const PylithScalar* velocity = &s[sOff[i_velocity]];

    for (PylithInt i = 0; i < dim; ++i) {
        g0[i] += velocity[i];
    } // for
} // g0u


// ---------------------------------------------------------------------------------------------------------------------
// g0v_grav - g0 function for generic poroelasticity terms ( + grav body forces).
void
pylith::fekernels::Poroelasticity::g0v_grav(const PylithInt dim,
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
                                            PylithScalar g0[]) {
    // Incoming auxililary fields.

    // Poroelasticity
    const PylithInt i_solid_density = 0;
    const PylithInt i_fluid_density = 1;
    const PylithInt i_porosity = 3;

    // 3 + n
    const PylithInt i_gravityField = 4;

    const PylithScalar bulkDensity = (1 - a[aOff[i_porosity]]) * a[aOff[i_solid_density]] + a[aOff[i_porosity]] * a[aOff[i_fluid_density]];
    const PylithScalar* gravityField = &a[aOff[i_gravityField]];

    for (PylithInt i = 0; i < dim; ++i) {
        g0[i] += bulkDensity * gravityField[i];
    } // for

} // g0v_grav


// ---------------------------------------------------------------------------------------------------------------------
// g0v_bodyforce - g0 function for generic poroelasticity terms ( + body forces).
void
pylith::fekernels::Poroelasticity::g0v_bodyforce(const PylithInt dim,
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
                                                 PylithScalar g0[]) {
    // Incoming auxiliary fields

    // Poroelasticity

    // 3 + n
    const PylithInt i_bodyForce = 4;

    const PylithScalar* bodyForce = &a[aOff[i_bodyForce]];

    for (PylithInt i = 0; i < dim; ++i) {
        g0[i] += bodyForce[i];
    } // for
} // g0v_bodyforce


// ----------------------------------------------------------------------
// g0v_gravbodyforce - g0 function for isotropic linear Poroelasticity with both gravity and body forces.
void
pylith::fekernels::Poroelasticity::g0v_grav_bodyforce(const PylithInt dim,
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
                                                      PylithScalar g0[]) {
    // Incoming auxiliary fields.

    // Poroelasticity
    const PylithInt i_solid_density = 0;
    const PylithInt i_fluid_density = 1;
    const PylithInt i_porosity = 3;

    // 3 + n
    const PylithInt i_bodyForce = 4;
    const PylithInt i_gravityField = 5;

    const PylithScalar bulkDensity = (1 - a[aOff[i_porosity]]) * a[aOff[i_solid_density]] + a[aOff[i_porosity]] * a[aOff[i_fluid_density]];
    const PylithScalar* gravityField = &a[aOff[i_gravityField]];
    const PylithScalar* bodyForce = &a[aOff[i_bodyForce]];

    // gravity field
    for (PylithInt i = 0; i < dim; ++i) {
        g0[i] += bulkDensity * gravityField[i];
    } // for

    // body force
    for (PylithInt i = 0; i < dim; ++i) {
        g0[i] += bodyForce[i];
    } // for
} // g0v_gravbodyforce


// =============================================================================
// Pressure
// =============================================================================

// ----------------------------------------------------------------------
// g0p_sourceDensity - g0p function for generic poroelasticity terms (source density).
void
pylith::fekernels::Poroelasticity::g0p_source(const PylithInt dim,
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
                                              PylithScalar g0[]) {
    // Incoming auxiliary fields.

    // Poroelasticity

    const PylithInt i_sourceDensity = 0;
    const PylithScalar* sourceDensity = &a[aOff[i_sourceDensity]];

    for (PylithInt i = 0; i < dim; ++i) {
        g0[i] += sourceDensity[i];
    } // for
} // g0p_source


// ------------------------------------------------------------------------------
// g0p function for isotropic linear Poroelasticity plane strain with source density, gravity, and body force.
void
pylith::fekernels::Poroelasticity::g0p_sourceDensity(const PylithInt dim,
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
                                                     PylithScalar g0[]) {
    // Incoming auxiliary fields.

    // Poroelasticity

    // 3 + n
    const PylithInt i_sourceDensity = 3;

    const PylithInt _numS = 1; // Number passed on to g0p_source.

    const PylithInt numASource = 1; // Number passed on to g0p_source.
    const PylithInt aOffSource[1] = { aOff[i_sourceDensity] };
    const PylithInt aOffSource_x[1] = { aOff_x[i_sourceDensity] };

    pylith::fekernels::Poroelasticity::g0p_source(dim, _numS, numASource,
                                                  NULL, NULL, NULL, NULL, NULL,
                                                  aOffSource, aOffSource_x, a, a_t, a_x,
                                                  t, x, numConstants, constants, g0);
} // g0p_sourceDensity


// ------------------------------------------------------------------------------
// g0p function for isotropic linear Poroelasticity plane strain with source density, gravity, and body force.
void
pylith::fekernels::Poroelasticity::g0p_sourceDensity_grav(const PylithInt dim,
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
                                                          PylithScalar g0[]) {
    // Incoming auxiliary fields.

    // Poroelasticity

    // 2 + n
    const PylithInt i_sourceDensity = 4;

    const PylithInt _numS = 1; // Number passed on to g0p_source.

    const PylithInt numASource = 1; // Number passed on to g0p_source.
    const PylithInt aOffSource[1] = { aOff[i_sourceDensity] };
    const PylithInt aOffSource_x[1] = { aOff_x[i_sourceDensity] };

    pylith::fekernels::Poroelasticity::g0p_source(dim, _numS, numASource,
                                                  NULL, NULL, NULL, NULL, NULL,
                                                  aOffSource, aOffSource_x, a, a_t, a_x,
                                                  t, x, numConstants, constants, g0);
} // g0p_sourceDensity_grav


// ------------------------------------------------------------------------------
// g0p function for Poroelasticity with source density, gravity, and body force.
void
pylith::fekernels::Poroelasticity::g0p_sourceDensity_body(const PylithInt dim,
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
                                                          PylithScalar g0[]) {
    // Incoming auxiliary fields.

    // Poroelasticity

    // 3 + n
    const PylithInt i_sourceDensity = 4;

    const PylithInt _numS = 1; // Number passed on to g0p_source.

    const PylithInt numASource = 1; // Number passed on to g0p_source.
    const PylithInt aOffSource[1] = { aOff[i_sourceDensity] };
    const PylithInt aOffSource_x[1] = { aOff_x[i_sourceDensity] };

    pylith::fekernels::Poroelasticity::g0p_source(dim, _numS, numASource,
                                                  NULL, NULL, NULL, NULL, NULL,
                                                  aOffSource, aOffSource_x, a, a_t, a_x,
                                                  t, x, numConstants, constants, g0);
} // g0p_sourceDensity_body


// ------------------------------------------------------------------------------
// g0p function for Poroelasticity with source density, gravity, and body force.
void
pylith::fekernels::Poroelasticity::g0p_sourceDensity_grav_body(const PylithInt dim,
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
                                                               PylithScalar g0[]) {
    // Incoming auxiliary fields.

    // Poroelasticity

    // 3 + n
    const PylithInt i_sourceDensity = 5;

    const PylithInt _numS = 1; // Number passed on to g0p_source.

    const PylithInt numASource = 1; // Number passed on to g0p_source.
    const PylithInt aOffSource[1] = { aOff[i_sourceDensity] };
    const PylithInt aOffSource_x[1] = { aOff_x[i_sourceDensity] };

    pylith::fekernels::Poroelasticity::g0p_source(dim, _numS, numASource,
                                                  NULL, NULL, NULL, NULL, NULL,
                                                  aOffSource, aOffSource_x, a, a_t, a_x,
                                                  t, x, numConstants, constants, g0);
} // g0p_sourceDensity_grav_body


/* -------------------------------------------------------------------------- */
/*                           LHS Jacobian                                     */
/* -------------------------------------------------------------------------- */

// -----------------------------------------------------------------------------
// Jf0ee - Jf0 function for isotropic linear poroelasticity plane strain.
void
pylith::fekernels::Poroelasticity::Jf0ee(const PylithInt dim,
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
                                         PylithScalar Jf0[]) {
    assert(aOff);
    assert(a);

    Jf0[0] = -1.0;
} // Jg0ee


// -----------------------------------------------------------------------------
// Jf1eu - Jf1 function for isotropic linear poroelasticity plane strain.
void
pylith::fekernels::Poroelasticity::Jf1eu(const PylithInt dim,
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
                                         PylithScalar Jf1[]) {
    for (PylithInt d = 0; d < dim; ++d) {
        Jf1[d*dim+d] = 1.0;
    } // for
} // Jf1eu


// ---------------------------------------------------------------------------------------------------------------------
// Jf0 function for poroelasticity equation.
void
pylith::fekernels::Poroelasticity::Jf0vv(const PylithInt dim,
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
                                         PylithScalar Jf0[]) {
    const PylithInt _numA = 1;

    // Incoming auxiliary fields.
    const PylithInt i_solid_density = 0;
    const PylithInt i_fluid_density = 1;
    const PylithInt i_porosity = 3;

    const PylithScalar bulkDensity = (1 - a[aOff[i_porosity]]) * a[aOff[i_solid_density]] + a[aOff[i_porosity]] * a[aOff[i_fluid_density]];

    for (PetscInt i = 0; i < dim; ++i) {
        Jf0[i*dim+i] += s_tshift * bulkDensity;
    } // for
} // Jf0vv


// =====================================================================================================================
// Kernels for poroelasticity plane strain.
// =====================================================================================================================

void
pylith::fekernels::PoroelasticityPlaneStrain::cauchyStrain(const PylithInt dim,
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
                                                           PylithScalar strain[]) {
    const PylithInt _dim = 2;

    assert(_dim == dim);
    assert(numS >= 1);
    assert(sOff_x);
    assert(s_x);
    assert(strain);

    // Incoming solution field.
    const PylithInt i_displacement = 0;
    const PylithScalar* displacement_x = &s_x[sOff_x[i_displacement]];

    const PylithScalar strain_xx = displacement_x[0*_dim+0];
    const PylithScalar strain_yy = displacement_x[1*_dim+1];
    const PylithScalar strain_zz = 0.0;
    const PylithScalar strain_xy = 0.5*(displacement_x[0*_dim+1] + displacement_x[1*_dim+0]);

    strain[0] = strain_xx;
    strain[1] = strain_yy;
    strain[2] = strain_zz;
    strain[3] = strain_xy;
} // cauchyStrain


// =====================================================================================================================
// Kernels for poroelasticity in 3D
// =====================================================================================================================

void
pylith::fekernels::Poroelasticity3D::cauchyStrain(const PylithInt dim,
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
                                                  PylithScalar strain[]) {
    const PylithInt _dim = 3;

    assert(_dim == dim);
    assert(numS >= 1);
    assert(sOff_x);
    assert(s_x);
    assert(strain);

    // Incoming solution field.
    const PylithInt i_displacement = 0;
    const PylithScalar* displacement_x = &s_x[sOff_x[i_displacement]];

    const PylithScalar strain_xx = displacement_x[0*_dim+0];
    const PylithScalar strain_yy = displacement_x[1*_dim+1];
    const PylithScalar strain_zz = displacement_x[2*_dim+2];
    const PylithScalar strain_xy = 0.5*(displacement_x[0*_dim+1] + displacement_x[1*_dim+0]);
    const PylithScalar strain_yz = 0.5*(displacement_x[1*_dim+2] + displacement_x[2*_dim+1]);
    const PylithScalar strain_xz = 0.5*(displacement_x[0*_dim+2] + displacement_x[2*_dim+0]);

    strain[0] = strain_xx;
    strain[1] = strain_yy;
    strain[2] = strain_zz;
    strain[3] = strain_xy;
    strain[4] = strain_yz;
    strain[5] = strain_xz;
} // cauchyStrain
