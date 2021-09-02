
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
// f0u placeholder function for poroelasticity equation
void pylith::fekernels::Poroelasticity::f0u(const PylithInt dim,
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
                                            PylithScalar f0[])
{
    // Incoming solution fields.

    // Incoming auxiliary fields.

    for (PylithInt i = 0; i < dim; ++i)
    {
        f0[i] += 0.0;
        // PetscPrintf(PETSC_COMM_WORLD, "f0u[%i]: %f\n",i, f0[i]);
    } // for
} // f0u

// ---------------------------------------------------------------------------------------------------------------------
// f0v function for poroelasticity equation, implicit time stepping, quasistatic.
void pylith::fekernels::Poroelasticity::f0v_implicit(const PylithInt dim,
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
                                                     PylithScalar f0[])
{
    // Incoming solution fields.
    const PylithInt i_displacement = 0;
    const PylithInt i_velocity = 3;

    assert(sOff);
    assert(sOff[i_displacement] >= 0);
    assert(sOff[i_velocity] >= 0);
    assert(s_t);
    assert(aOff);
    assert(s);

    const PylithScalar *displacement_t = &s_t[sOff[i_displacement]]; // disp_t
    const PylithScalar *velocity = &s[sOff[i_velocity]];             // vel

    for (PylithInt i = 0; i < dim; ++i)
    {
        f0[i] += displacement_t[i];
        f0[i] -= velocity[i];
    } // for
} // f0v_implicit

// ---------------------------------------------------------------------------------------------------------------------
// f0v function for poroelasticity equation, explicit time stepping, dynamic.
void pylith::fekernels::Poroelasticity::f0v_explicit(const PylithInt dim,
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
                                                     PylithScalar f0[])
{
    // Incoming solution fields.
    const PylithInt i_velocity = 2;

    // Incoming auxiliary fields.
    const PylithInt i_solid_density = 0;
    const PylithInt i_fluid_density = 1;
    const PylithInt i_porosity = 3;

    assert(sOff);
    assert(sOff[i_velocity] >= 0);
    assert(s_t);
    assert(aOff);
    assert(aOff[i_solid_density] >= 0);
    assert(aOff[i_fluid_density] >= 0);
    assert(aOff[i_porosity] >= 0);
    assert(a);

    const PylithScalar bulkDensity = (1 - a[aOff[i_porosity]]) * a[aOff[i_solid_density]] + a[aOff[i_porosity]] * a[aOff[i_fluid_density]];
    const PylithScalar *velocity_t = &s_t[sOff[i_velocity]]; // acceleration

    for (PylithInt i = 0; i < dim; ++i)
    {
        f0[i] += velocity_t[i] * bulkDensity;
    } // for
} // f0v_explicit

// =============================================================================
// Volumetric Strain
// =============================================================================
// ----------------------------------------------------------------------
// f0e function for isotropic linear Poroelasticity plane strain.
void pylith::fekernels::Poroelasticity::f0e(const PylithInt dim,
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
                                            PylithScalar f0[])
{
    // Incoming solution fields.
    const PylithInt i_displacement = 0;
    const PylithInt i_trace_strain = 2;

    assert(sOff);
    assert(sOff[i_displacement] >= 0);
    assert(sOff[i_trace_strain] >= 0);
    assert(sOff_x);
    assert(sOff_x[i_displacement] >= 0);
    assert(s);

    const PylithScalar *displacement = &s[sOff[i_displacement]];
    const PylithScalar *displacement_x = &s_x[sOff_x[i_displacement]];
    const PylithScalar trace_strain = s[sOff[i_trace_strain]];

    for (PylithInt d = 0; d < dim; ++d)
    {
        f0[0] += displacement_x[d * dim + d];
    }
    f0[0] -= trace_strain;
} // f0e

// ---------------------------------------------------------------------------------------------------------------------
// f0pdot function for poroelasticity equation, implicit time stepping, quasistatic.
void pylith::fekernels::Poroelasticity::f0pdot(const PylithInt dim,
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
                                               PylithScalar f0[])
{
    // Incoming solution fields.
    const PylithInt i_pressure = 1;
    const PylithInt i_pdot = 4;

    assert(sOff);
    assert(sOff[i_pressure] >= 0);
    assert(sOff[i_pdot] >= 0);
    assert(s);
    assert(s_t);

    const PylithScalar pressure_t = s_t[sOff[i_pressure]]; // disp_t
    const PylithScalar pdot = s[sOff[i_pdot]];             // vel

    f0[0] += pressure_t;
    f0[0] -= pdot;
} // f0pdot

// ---------------------------------------------------------------------------------------------------------------------
// f0edot function for poroelasticity equation, implicit time stepping, quasistatic.
void pylith::fekernels::Poroelasticity::f0edot(const PylithInt dim,
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
                                               PylithScalar f0[])
{
    // Incoming solution fields.
    const PylithInt i_trace_strain = 2;
    const PylithInt i_edot = 5;

    assert(sOff);
    assert(sOff[i_trace_strain] >= 0);
    assert(sOff[i_edot] >= 0);

    const PylithScalar trace_strain_t = s_t[sOff[i_trace_strain]]; // disp_t
    const PylithScalar edot = s[sOff[i_edot]];                     // vel

    f0[0] += trace_strain_t - edot;
} // f0pdot

// =============================================================================
// Displacement
// =============================================================================

// ----------------------------------------------------------------------
// g0 function for displacement equation: g0u = v.
void pylith::fekernels::Poroelasticity::g0u(const PylithInt dim,
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
                                            PylithScalar g0[])
{
    const PylithInt i_velocity = 2;

    assert(sOff);
    assert(sOff[i_velocity] >= 0);

    const PylithScalar *velocity = &s[sOff[i_velocity]];

    for (PylithInt i = 0; i < dim; ++i)
    {
        g0[i] += velocity[i];
    } // for
} // g0u

// ---------------------------------------------------------------------------------------------------------------------
// g0v_grav - g0 function for generic poroelasticity terms ( + grav body forces).
void pylith::fekernels::Poroelasticity::g0v_grav(const PylithInt dim,
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
                                                 PylithScalar g0[])
{
    const PylithInt i_velocity = 2;

    // Incoming auxililary fields.

    // Poroelasticity
    const PylithInt i_solid_density = 0;
    const PylithInt i_fluid_density = 1;
    const PylithInt i_porosity = 3;

    // 3 + n
    const PylithInt i_gravity_field = 4;

    assert(aOff);
    assert(aOff[i_solid_density] >= 0);
    assert(aOff[i_fluid_density] >= 0);
    assert(aOff[i_porosity] >= 0);
    assert(aOff[i_gravity_field] >= 0);
    assert(a);

    const PylithScalar bulkDensity = (1 - a[aOff[i_porosity]]) * a[aOff[i_solid_density]] + a[aOff[i_porosity]] * a[aOff[i_fluid_density]];
    const PylithScalar *gravity_field = &a[aOff[i_gravity_field]];

    for (PylithInt i = 0; i < dim; ++i)
    {
        g0[i] += bulkDensity * gravity_field[i];
    } // for
} // g0u

// ---------------------------------------------------------------------------------------------------------------------
// g0v_bodyforce - g0 function for generic poroelasticity terms ( + body forces).
void pylith::fekernels::Poroelasticity::g0v_bodyforce(const PylithInt dim,
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
                                                      PylithScalar g0[])
{
    // Incoming auxiliary fields

    // Poroelasticity

    // 3 + n
    const PylithInt i_body_force = 4;

    assert(aOff);
    assert(aOff[i_body_force] >= 0);
    assert(a);

    const PylithScalar *body_force = &a[aOff[i_body_force]];

    for (PylithInt i = 0; i < dim; ++i)
    {
        g0[i] += body_force[i];
    } // for
} // g0v_bodyforce

// ----------------------------------------------------------------------
// g0v_gravbodyforce - g0 function for isotropic linear Poroelasticity with both gravity and body forces.
void pylith::fekernels::Poroelasticity::g0v_grav_bodyforce(const PylithInt dim,
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
                                                           PylithScalar g0[])
{
    // Incoming auxiliary fields.

    // Poroelasticity
    const PylithInt i_solid_density = 0;
    const PylithInt i_fluid_density = 1;
    const PylithInt i_porosity = 3;

    // 3 + n
    const PylithInt i_body_force = 4;
    const PylithInt i_gravity_field = 5;

    assert(aOff);
    assert(aOff[i_solid_density] >= 0);
    assert(aOff[i_fluid_density] >= 0);
    assert(aOff[i_porosity] >= 0);
    assert(aOff[i_gravity_field] >= 0);
    assert(aOff[i_body_force] >= 0);
    assert(a);

    const PylithScalar bulkDensity = (1 - a[aOff[i_porosity]]) * a[aOff[i_solid_density]] + a[aOff[i_porosity]] * a[aOff[i_fluid_density]];
    const PylithScalar *gravity_field = &a[aOff[i_gravity_field]];
    const PylithScalar *body_force = &a[aOff[i_body_force]];

    // gravity field
    for (PylithInt i = 0; i < dim; ++i)
    {
        g0[i] += bulkDensity * gravity_field[i];
    } // for

    // body force
    for (PylithInt i = 0; i < dim; ++i)
    {
        g0[i] += body_force[i];
    } // for
} // g0v_gravbodyforce

// =============================================================================
// Pressure
// =============================================================================

// ----------------------------------------------------------------------
// g0p_source - g0p function for generic poroelasticity terms (source density).
void pylith::fekernels::Poroelasticity::g0p_source(const PylithInt dim,
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
                                                   PylithScalar g0[])
{
    // Incoming auxiliary fields.

    // Poroelasticity

    const PylithInt i_source_density = 0;

    assert(aOff);
    assert(aOff[i_source_density] >= 0);
    assert(a);

    const PylithScalar *source_density = &a[aOff[i_source_density]];

    for (PylithInt i = 0; i < dim; ++i)
    {
        g0[i] += source_density[i];
    } // for
} // g0p_source

// ------------------------------------------------------------------------------
// g0p function for isotropic linear Poroelasticity plane strain with source density, gravity, and body force.
void pylith::fekernels::Poroelasticity::g0p_sourceDensity(const PylithInt dim,
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
                                                          PylithScalar g0[])
{
    // Incoming auxiliary fields.

    // Poroelasticity

    // 3 + n
    const PylithInt i_source_density = 3;

    assert(aOff);
    assert(aOff[i_source_density] >= 0);
    assert(a);

    const PylithInt _numS = 1; // Number passed on to g0p_source.

    const PylithInt numASource = 1; // Number passed on to g0p_source.
    const PylithInt aOffSource[1] = {aOff[i_source_density]};
    const PylithInt aOffSource_x[1] = {aOff_x[i_source_density]};

    pylith::fekernels::Poroelasticity::g0p_source(dim, _numS, numASource,
                                                  NULL, NULL, NULL, NULL, NULL,
                                                  aOffSource, aOffSource_x, a, a_t, a_x,
                                                  t, x, numConstants, constants, g0);
} // g0p_sourceDensity

// ------------------------------------------------------------------------------
// g0p function for isotropic linear Poroelasticity plane strain with source density, gravity, and body force.
void pylith::fekernels::Poroelasticity::g0p_sourceDensity_grav(const PylithInt dim,
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
                                                               PylithScalar g0[])
{
    // Incoming auxiliary fields.

    // Poroelasticity

    // 2 + n
    const PylithInt i_source_density = 4;

    assert(aOff);
    assert(aOff[i_source_density] >= 0);
    assert(a);

    const PylithInt _numS = 1; // Number passed on to g0p_source.

    const PylithInt numASource = 1; // Number passed on to g0p_source.
    const PylithInt aOffSource[1] = {aOff[i_source_density]};
    const PylithInt aOffSource_x[1] = {aOff_x[i_source_density]};

    pylith::fekernels::Poroelasticity::g0p_source(dim, _numS, numASource,
                                                  NULL, NULL, NULL, NULL, NULL,
                                                  aOffSource, aOffSource_x, a, a_t, a_x,
                                                  t, x, numConstants, constants, g0);
} // g0p_sourceDensity_grav

// ------------------------------------------------------------------------------
// g0p function for Poroelasticity with source density, gravity, and body force.
void pylith::fekernels::Poroelasticity::g0p_sourceDensity_body(const PylithInt dim,
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
                                                               PylithScalar g0[])
{
    // Incoming auxiliary fields.

    // Poroelasticity

    // 3 + n
    const PylithInt i_source_density = 4;

    assert(aOff);
    assert(aOff[i_source_density] >= 0);
    assert(a);

    const PylithInt _numS = 1; // Number passed on to g0p_source.

    const PylithInt numASource = 1; // Number passed on to g0p_source.
    const PylithInt aOffSource[1] = {aOff[i_source_density]};
    const PylithInt aOffSource_x[1] = {aOff_x[i_source_density]};

    pylith::fekernels::Poroelasticity::g0p_source(dim, _numS, numASource,
                                                  NULL, NULL, NULL, NULL, NULL,
                                                  aOffSource, aOffSource_x, a, a_t, a_x,
                                                  t, x, numConstants, constants, g0);
} // g0p_sourceDensity_body

// ------------------------------------------------------------------------------
// g0p function for Poroelasticity with source density, gravity, and body force.
void pylith::fekernels::Poroelasticity::g0p_sourceDensity_grav_body(const PylithInt dim,
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
                                                                    PylithScalar g0[])
{
    // Incoming auxiliary fields.

    // Poroelasticity

    // 3 + n
    const PylithInt i_source_density = 5;

    assert(aOff);
    assert(aOff[i_source_density] >= 0);
    assert(a);

    const PylithInt _numS = 1; // Number passed on to g0p_source.

    const PylithInt numASource = 1; // Number passed on to g0p_source.
    const PylithInt aOffSource[1] = {aOff[i_source_density]};
    const PylithInt aOffSource_x[1] = {aOff_x[i_source_density]};

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
void pylith::fekernels::Poroelasticity::Jf0ee(const PylithInt dim,
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
                                              PylithScalar Jf0[])
{
    assert(aOff);
    assert(a);

    Jf0[0] = -1.0;
} // Jg0ee

// -----------------------------------------------------------------------------
// Jf1eu - Jf1 function for isotropic linear poroelasticity plane strain.
void pylith::fekernels::Poroelasticity::Jf1eu(const PylithInt dim,
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
                                              PylithScalar Jf1[])
{
    assert(aOff);
    assert(a);

    for (PylithInt d = 0; d < dim; ++d)
    {
        Jf1[d * dim + d] = 1.0;
    } // for
} // Jf1eu

// ---------------------------------------------------------------------------------------------------------------------
// Jf0vu function for poroelasticity equation, quasistatic.
void pylith::fekernels::Poroelasticity::Jf0vu_implicit(const PylithInt dim,
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
                                                       PylithScalar Jf0[])
{
    assert(aOff);
    assert(a);

    // Incoming auxiliary fields.

    for (PylithInt d = 0; d < dim; ++d)
    {
        Jf0[d * dim + d] += s_tshift;
    } // for
} // Jf0vu_implicit

// ---------------------------------------------------------------------------------------------------------------------
// Jf0vv function for poroelasticity equation, quasistatic.
void pylith::fekernels::Poroelasticity::Jf0vv_implicit(const PylithInt dim,
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
                                                       PylithScalar Jf0[])
{
    assert(aOff);
    assert(a);

    // Incoming auxiliary fields.

    for (PylithInt d = 0; d < dim; ++d)
    {
        Jf0[d * dim + d] -= 1.0;
    } // for
} // Jf0vv_implicit

// ---------------------------------------------------------------------------------------------------------------------
// Jf0vv function for poroelasticity equation, dynamic.
void pylith::fekernels::Poroelasticity::Jf0vv_explicit(const PylithInt dim,
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
                                                       PylithScalar Jf0[])
{
    // Incoming auxiliary fields.
    const PylithInt i_solid_density = 0;
    const PylithInt i_fluid_density = 1;
    const PylithInt i_porosity = 3;

    assert(aOff);
    assert(aOff[i_solid_density] >= 0);
    assert(aOff[i_fluid_density] >= 0);
    assert(aOff[i_porosity] >= 0);
    assert(a);

    const PylithScalar bulkDensity = (1 - a[aOff[i_porosity]]) * a[aOff[i_solid_density]] + a[aOff[i_porosity]] * a[aOff[i_fluid_density]];

    for (PetscInt i = 0; i < dim; ++i)
    {
        Jf0[i * dim + i] += s_tshift * bulkDensity;
    } // for
} // Jf0vv_explicit

// -----------------------------------------------------------------------------
// Jf0pdotp - Jf0 function for isotropic linear poroelasticity plane strain.
void pylith::fekernels::Poroelasticity::Jf0pdotp(const PylithInt dim,
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
                                                 PylithScalar Jf0[])
{
    assert(aOff);
    assert(a);

    Jf0[0] += s_tshift;
} // Jg0pdotp

// -----------------------------------------------------------------------------
// Jf0pdotpdot - Jf0 function for isotropic linear poroelasticity plane strain.
void pylith::fekernels::Poroelasticity::Jf0pdotpdot(const PylithInt dim,
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
                                                    PylithScalar Jf0[])
{
    assert(aOff);
    assert(a);

    Jf0[0] -= 1.0;
} // Jg0pdotpdot

// -----------------------------------------------------------------------------
// Jf0edote - Jf0 function for isotropic linear poroelasticity plane strain.
void pylith::fekernels::Poroelasticity::Jf0edote(const PylithInt dim,
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
                                                 PylithScalar Jf0[])
{
    assert(aOff);
    assert(a);

    Jf0[0] += s_tshift;
} // Jg0edote

// -----------------------------------------------------------------------------
// Jf0edotedot - Jf0 function for isotropic linear poroelasticity plane strain.
void pylith::fekernels::Poroelasticity::Jf0edotedot(const PylithInt dim,
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
                                                    PylithScalar Jf0[])
{
    assert(aOff);
    assert(a);

    Jf0[0] -= 1.0;
} // Jg0edotedot

// =====================================================================================================================
// Kernels for poroelasticity plane strain.
// =====================================================================================================================

void pylith::fekernels::PoroelasticityPlaneStrain::cauchyStrain(const PylithInt dim,
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
                                                                PylithScalar strain[])
{
    const PylithInt _dim = 2;

    assert(_dim == dim);
    assert(numS >= 1);
    assert(sOff_x);
    assert(s_x);
    assert(strain);

    // Incoming solution field.
    const PylithInt i_displacement = 0;
    const PylithScalar *displacement_x = &s_x[sOff_x[i_displacement]];

    const PylithScalar strain_xx = displacement_x[0 * _dim + 0];
    const PylithScalar strain_yy = displacement_x[1 * _dim + 1];
    const PylithScalar strain_zz = 0.0;
    const PylithScalar strain_xy = 0.5 * (displacement_x[0 * _dim + 1] + displacement_x[1 * _dim + 0]);

    strain[0] = strain_xx;
    strain[1] = strain_yy;
    strain[2] = strain_zz;
    strain[3] = strain_xy;
} // cauchyStrain

// =====================================================================================================================
// Kernels for poroelasticity in 3D
// =====================================================================================================================

void pylith::fekernels::Poroelasticity3D::cauchyStrain(const PylithInt dim,
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
                                                       PylithScalar strain[])
{
    const PylithInt _dim = 3;

    assert(_dim == dim);
    assert(numS >= 1);
    assert(sOff_x);
    assert(s_x);
    assert(strain);

    // Incoming solution field.
    const PylithInt i_displacement = 0;
    const PylithScalar *displacement_x = &s_x[sOff_x[i_displacement]];

    const PylithScalar strain_xx = displacement_x[0 * _dim + 0];
    const PylithScalar strain_yy = displacement_x[1 * _dim + 1];
    const PylithScalar strain_zz = displacement_x[2 * _dim + 2];
    const PylithScalar strain_xy = 0.5 * (displacement_x[0 * _dim + 1] + displacement_x[1 * _dim + 0]);
    const PylithScalar strain_yz = 0.5 * (displacement_x[1 * _dim + 2] + displacement_x[2 * _dim + 1]);
    const PylithScalar strain_xz = 0.5 * (displacement_x[0 * _dim + 2] + displacement_x[2 * _dim + 0]);

    strain[0] = strain_xx;
    strain[1] = strain_yy;
    strain[2] = strain_zz;
    strain[3] = strain_xy;
    strain[4] = strain_yz;
    strain[5] = strain_xz;
} // cauchyStrain
