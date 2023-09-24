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

/** @file libsrc/fekernels/Poroelasticity.hh
 *
 * Solution fields: [disp(dim), pres, vol_strain]
 *
 * Displacement
 * \int_V \vec{\phi}_u \cdot \left( \rho \frac{\partial \vec{v}(t)}{\partial t} \right) \, dV =
 *   \int_V \vec{\phi}_u \cdot \vec{f}(t) - \nabla \vec{\phi}_u : \tensor{\sigma}(\vec{u}) \, dV +
 *   \int_{S_\tau} \vec{\phi}_u \cdot \vec{\tau}(t) \, dS.
 *
 *
 *** Kernel interface.
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

#if !defined(pylith_fekernels_poroelasticity_hh)
#define pylith_fekernels_poroelasticity_hh

// Include directives ---------------------------------------------------
#include "fekernelsfwd.hh" // forward declarations

#include "pylith/fekernels/Elasticity.hh" // USES Elasticity kernels

#include "pylith/fekernels/Tensor.hh"

#include "pylith/utils/types.hh"

#include <cassert> // USES assert()

class pylith::fekernels::Poroelasticity {
public:

    struct Context {
        PylithInt dim;
        const PylithReal* displacement;
        const PylithReal* displacement_t;
        const PylithReal* displacement_x;
        PylithReal pressure;
        PylithReal pressure_t;
        const PylithReal* pressure_x;
        PylithReal trace_strain;
        PylithReal trace_strain_t;
        const PylithReal* trace_strain_x;
        const PylithReal* velocity;
        const PylithReal* velocity_t;
        const PylithReal* velocity_x;
        PylithReal pressure_dot;
        PylithReal pressure_dot_t;
        const PylithReal* pressure_dot_x;
        PylithReal trace_strain_dot;
        PylithReal trace_strain_dot_t;
        const PylithReal* trace_strain_dot_x;
        const PylithReal* x;
        PylithReal solidDensity;
        PylithReal fluidDensity;
        PylithReal fluidViscosity;
        PylithReal porosity;
        PylithReal bulkDensity;
        PylithReal gravityField[3];
        PylithReal bodyForce[3];
        PylithReal sourceDensity;

        Context(void) :
            dim(0),
            displacement(NULL),
            displacement_t(NULL),
            displacement_x(NULL),
            pressure_x(NULL),
            trace_strain_x(NULL),
            velocity(NULL),
            velocity_t(NULL),
            velocity_x(NULL),
            pressure_dot_x(NULL),
            trace_strain_dot_x(NULL),
            x(NULL),
            sourceDensity(0.0) {
            const size_t _dim = 3;
            for (size_t i = 0; i < _dim; ++i) {
                gravityField[i] = 0.0;
                bodyForce[i] = 0.0;
            } // for
        } // default constructor

    };

    // Interface for functions computing time variation in fluid content term.
    typedef void (*timevarfluidcontentfn_type) (const Context& context,
                                                pylith::fekernels::Tensor*);

    // Interface for functions computing flux term (Darcy).
    typedef void (*fluxratefn_type) (const Context& Context,
                                     void*,
                                     const pylith::fekernels::TensorOps&,
                                     pylith::fekernels::Tensor*);

    // PUBLIC MEMBERS ///////////////////////////////////////////////////////
public:

    // ----------------------------------------------------------------------
    /*
     * Poroelastic context for standard QS
     */
    static inline
    void setContextQuasistatic(Context* context,
                               const PylithInt dim,
                               const PylithInt numS,
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
                               const PylithScalar x[]) {
        assert(context);
        assert(numS >= 3);

        // Incoming solution fields.
        const PylithInt i_displacement = 0;
        const PylithInt i_pressure = 1;
        const PylithInt i_trace_strain = 2;

        // Incoming auxiliary fields.
        const PylithInt i_solidDensity = 0;
        const PylithInt i_fluidDensity = 1;
        const PylithInt i_fluidViscosity = 2;
        const PylithInt i_porosity = 3;

        assert(sOff[i_displacement] >= 0);
        assert(sOff[i_pressure] >= 0);
        assert(sOff[i_trace_strain] >= 0);

        assert(aOff[i_solidDensity] >= 0);
        assert(aOff[i_fluidDensity] >= 0);
        assert(aOff[i_porosity] >= 0);
        assert(aOff[i_fluidViscosity] >= 0);

        context->dim = dim;
        context->displacement = &s[sOff[i_displacement]];
        context->displacement_t = &s_t[sOff[i_displacement]];
        context->displacement_x = &s_x[sOff_x[i_displacement]];
        context->pressure = s[sOff[i_pressure]];
        context->pressure_t = s_t ? s_t[sOff[i_pressure]] : 0.0;
        context->pressure_x = &s_x[sOff_x[i_pressure]];
        context->trace_strain = s[sOff[i_trace_strain]];
        context->trace_strain_t = s_t ? s_t[sOff[i_trace_strain]] : 0.0;
        context->trace_strain_x = &s_x[sOff_x[i_trace_strain]];
        context->x = x;

        // Poroelastic auxiliaries
        context->solidDensity = a[aOff[i_solidDensity]];
        context->fluidDensity = a[aOff[i_fluidDensity]];
        context->porosity = a[aOff[i_porosity]];
        context->fluidViscosity = a[aOff[i_fluidViscosity]];
        context->bulkDensity = (1 - a[aOff[i_porosity]]) * a[aOff[i_solidDensity]] + a[aOff[i_porosity]] * a[aOff[i_fluidDensity]];

    } // setContextQuasistatic

    // ----------------------------------------------------------------------
    /*
     * Poroelastic context for six field QS
     */
    static inline
    void setContextQuasistatic_statevars(Context* context,
                                         const PylithInt dim,
                                         const PylithInt numS,
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
                                         const PylithScalar x[]) {
        assert(context);
        assert(numS >= 6);

        // Incoming solution fields.
        const PylithInt i_velocity = 3;
        const PylithInt i_pressure_dot = 4;
        const PylithInt i_trace_strain_dot = 5;

        assert(sOff[i_velocity] >= 0);
        assert(sOff[i_pressure_dot] >= 0);
        assert(sOff[i_trace_strain_dot] >= 0);

        // Variables
        context->velocity = &s[sOff[i_velocity]];
        context->velocity_t = &s_t[sOff[i_velocity]];
        context->velocity_x = &s_x[sOff_x[i_velocity]];
        context->pressure_dot = s ? s[sOff[i_pressure_dot]] : 0.0;
        context->pressure_dot_t = s_t[sOff[i_pressure_dot]];
        context->pressure_dot_x = &s_x[sOff_x[i_pressure_dot]];
        context->trace_strain_dot = s ? s[sOff[i_trace_strain_dot]] : 0.0;
        context->trace_strain_dot_t = s_t[sOff[i_trace_strain_dot]];
        context->trace_strain_dot_x = &s_x[sOff_x[i_trace_strain_dot]];

    } // setContextQuasistatic_statevars

    // ----------------------------------------------------------------------
    /*
     * Poroelastic context for standard dynamic
     */
    static inline
    void setContextDynamic(Context* context,
                           const PylithInt dim,
                           const PylithInt numS,
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
                           const PylithScalar x[]) {
        assert(context);
        assert(numS >= 3);

        // Incoming solution fields.
        const PylithInt i_displacement = 0;
        const PylithInt i_pressure = 1;
        const PylithInt i_velocity = 2;

        // Incoming auxiliary fields.
        const PylithInt i_solidDensity = 0;
        const PylithInt i_fluidDensity = 1;
        const PylithInt i_fluidViscosity = 2;
        const PylithInt i_porosity = 3;

        assert(sOff[i_displacement] >= 0);
        assert(sOff[i_pressure] >= 0);
        assert(sOff[i_velocity] >= 0);

        assert(aOff[i_solidDensity] >= 0);
        assert(aOff[i_fluidDensity] >= 0);
        assert(aOff[i_porosity] >= 0);
        assert(aOff[i_fluidViscosity] >= 0);

        // Variables &c
        context->dim = dim;
        context->displacement = &s[sOff[i_displacement]];
        context->displacement_t = &s_t[sOff[i_displacement]];
        context->displacement_x = &s_x[sOff_x[i_displacement]];
        context->pressure = s[sOff[i_pressure]];
        context->pressure_t = s_t[sOff[i_pressure]];
        context->pressure_x = &s_x[sOff_x[i_pressure]];
        context->velocity = &s[sOff[i_velocity]];
        context->velocity_t = &s_t[sOff[i_velocity]];
        context->velocity_x = &s_x[sOff_x[i_velocity]];
        context->x = x;

        // Poroelastic auxiliaries
        context->solidDensity = a[aOff[i_solidDensity]];
        context->fluidDensity = a[aOff[i_fluidDensity]];
        context->porosity = a[aOff[i_porosity]];
        context->fluidViscosity = a[aOff[i_fluidViscosity]];
        context->bulkDensity = (1 - a[aOff[i_porosity]]) * a[aOff[i_solidDensity]] + a[aOff[i_porosity]] * a[aOff[i_fluidDensity]];

        // Improvised values
        const PylithScalar *displacement_x = &s_x[sOff_x[i_displacement]];
        const PylithScalar *velocity_x = &s_x[sOff_x[i_velocity]];

        PylithScalar trace_strain = 0.0;
        for (PylithInt d = 0; d < dim; ++d) {
            trace_strain += displacement_x[d * dim + d];
        }

        PylithScalar trace_strain_t = 0.0;
        for (PylithInt d = 0; d < dim; ++d) {
            trace_strain_t += velocity_x[d * dim + d];
        }

        context->trace_strain = trace_strain;
        context->trace_strain_t = trace_strain_t;

    } // setContextDynamic

    // ----------------------------------------------------------------------
    /*
     * Add gravity field to context
     */
    static inline
    void setContextGravity(Context* context,
                           const PylithInt dim,
                           const PylithInt numS,
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
                           const PylithScalar x[]) {
        assert(context);
        assert(numS >= 3);

        const PylithInt i_gravityField = 4;

        assert(aOff[i_gravityField] >= 0);

        // Gravity Field
        const PylithScalar* gravityField = &a[aOff[i_gravityField]];
        for (PylithInt i = 0; i < dim; ++i) {
            context->gravityField[i] = gravityField[i];
        } // for

    } // setContextGravity

    // ----------------------------------------------------------------------
    /*
     * Add source density to context
     */
    static inline
    void setContextBodyForce(Context* context,
                             const PylithInt dim,
                             const PylithInt numS,
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
                             const PylithScalar x[]) {
        assert(context);
        assert(numS >= 3);

        const PylithInt i_bodyForce = 4;

        assert(aOff[i_bodyForce] >= 0);

        const PylithScalar* bodyForce = &a[aOff[i_bodyForce]];
        for (PylithInt i = 0; i < dim; ++i) {
            context->bodyForce[i] = bodyForce[i];
        } // for

    } // setContextBodyForce

    // ----------------------------------------------------------------------
    /*
     * Add source density to context
     */
    static inline
    void setContextSourceDensity(Context* context,
                                 const PylithInt dim,
                                 const PylithInt numS,
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
                                 const PylithScalar x[]) {
        assert(context);
        assert(numS >= 3);

        const PylithInt i_sourceDensity = 4;

        assert(aOff[i_sourceDensity] >= 0);
        context->sourceDensity = a[aOff[i_sourceDensity]];

    } // setContextSourceDensity

    // ----------------------------------------------------------------------
    /*
     * Add gravity field and body force to context
     */
    static inline
    void setContextGravityBodyForce(Context* context,
                                    const PylithInt dim,
                                    const PylithInt numS,
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
                                    const PylithScalar x[]) {
        assert(context);
        assert(numS >= 3);

        const PylithInt i_gravityField = 4;
        const PylithInt i_bodyForce = 5;

        assert(aOff[i_gravityField] >= 0);
        assert(aOff[i_bodyForce] >= 0);

        const PylithScalar* bodyForce = &a[aOff[i_bodyForce]];
        const PylithScalar* gravityField = &a[aOff[i_gravityField]];

        for (PylithInt i = 0; i < dim; ++i) {
            context->gravityField[i] = gravityField[i];
            context->bodyForce[i] = bodyForce[i];
        } // for

    } // setContextGravityBodyForce

    // ----------------------------------------------------------------------
    /*
     * Add gravity field, and source density to context
     */
    static inline
    void setContextGravitySourceDensity(Context* context,
                                        const PylithInt dim,
                                        const PylithInt numS,
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
                                        const PylithScalar x[]) {
        assert(context);
        assert(numS >= 3);

        const PylithInt i_gravityField = 4;
        const PylithInt i_sourceDensity = 5;

        assert(aOff[i_gravityField] >= 0);
        assert(aOff[i_sourceDensity] >= 0);

        context->sourceDensity = a[aOff[i_sourceDensity]];

        const PylithScalar* gravityField = &a[aOff[i_gravityField]];

        for (PylithInt i = 0; i < dim; ++i) {
            context->gravityField[i] = gravityField[i];
        } // for

    } // setContextGravityBodyForceSourceDensity

    // ----------------------------------------------------------------------
    /*
     * Add gravity field, body force, and source density to context
     */
    static inline
    void setContextGravityBodyForceSourceDensity(Context* context,
                                                 const PylithInt dim,
                                                 const PylithInt numS,
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
                                                 const PylithScalar x[]) {
        assert(context);
        assert(numS >= 3);

        const PylithInt i_gravityField = 4;
        const PylithInt i_bodyForce = 5;
        const PylithInt i_sourceDensity = 6;

        assert(aOff[i_gravityField] >= 0);
        assert(aOff[i_bodyForce] >= 0);
        assert(aOff[i_sourceDensity] >= 0);

        context->sourceDensity = a[aOff[i_sourceDensity]];

        const PylithScalar* bodyForce = &a[aOff[i_bodyForce]];
        const PylithScalar* gravityField = &a[aOff[i_gravityField]];

        for (PylithInt i = 0; i < dim; ++i) {
            context->gravityField[i] = gravityField[i];
            context->bodyForce[i] = bodyForce[i];
        } // for

    } // setContextGravityBodyForceSourceDensity

    // ----------------------------------------------------------------------
    /*
     * Add body force, and source density to context
     */
    static inline
    void setContextBodyForceSourceDensity(Context* context,
                                          const PylithInt dim,
                                          const PylithInt numS,
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
                                          const PylithScalar x[]) {
        assert(context);
        assert(numS >= 3);

        const PylithInt i_bodyForce = 4;
        const PylithInt i_sourceDensity = 5;

        assert(aOff[i_bodyForce] >= 0);
        assert(aOff[i_sourceDensity] >= 0);

        const PylithScalar* bodyForce = &a[aOff[i_bodyForce]];

        context->sourceDensity = a[aOff[i_sourceDensity]];

        // Set others to zero
        for (PylithInt i = 0; i < dim; ++i) {
            context->bodyForce[i] = bodyForce[i];
        } // for

    } // setContextBodyForceSourceDensity

#if 0
    // =============================================================================
    // Displacement
    // =============================================================================
    // ----------------------------------------------------------------------
    /** f0 function for poroelasticity equation for displacement field
     *
     * Placeholder function
     *
     */
    static inline
    void f0u(const PylithInt dim,
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
        for (PylithInt i = 0; i < dim; ++i) {
            f0[i] += 0.0;
        } // for
    } // f0u

#endif

    // Calculate bulk density
    static inline
    void bulkDensity_asScalar(const PylithInt dim,
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
                              PylithReal* bulkDensity) {
        // Poroelastic Context
        pylith::fekernels::Poroelasticity::Context Context;
        pylith::fekernels::Poroelasticity::setContextQuasistatic(
            &Context, dim, numS, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x, t, x);

        *bulkDensity = Context.bulkDensity;

    } // bulkDensity_asScalar

    // =============================================================================
    // Velocity
    // =============================================================================
    // ----------------------------------------------------------------------
    /** f0 function for implicit time stepping poroelasticity equation for velocity field
     *
     */

    static inline
    void f0v_implicit(const PylithInt dim,
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
        // Poroelastic Context
        pylith::fekernels::Poroelasticity::Context Context;
        pylith::fekernels::Poroelasticity::setContextQuasistatic(
            &Context, dim, numS, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x, t, x);
        pylith::fekernels::Poroelasticity::setContextQuasistatic_statevars(
            &Context, dim, numS, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x, t, x);

        const PylithScalar* displacement_t = Context.displacement_t; // disp_t
        const PylithScalar* velocity = Context.velocity; // vel
        assert(displacement_t);
        assert(velocity);

        for (PylithInt i = 0; i < dim; ++i) {
            f0[i] += displacement_t[i];
            f0[i] -= velocity[i];
        } // for
    } // f0v_implicit

    // ---------------------------------------------------------------------------------------------------------------------
    /*
     * f0v function for poroelasticity equation, explicit time stepping, dynamic.
     *
     */

    static inline
    void f0v_explicit(const PylithInt dim,
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
        pylith::fekernels::Poroelasticity::Context Context;
        pylith::fekernels::Poroelasticity::setContextDynamic(
            &Context, dim, numS, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x, t, x);

        const PylithScalar bulkDensity = Context.bulkDensity; // Bulk Density
        const PylithScalar* velocity_t = Context.velocity_t; // acceleration

        for (PylithInt i = 0; i < dim; ++i) {
            f0[i] += velocity_t[i] * bulkDensity;
        } // for
    } // f0v_explicit

    // =============================================================================
    // Volumetric Strain
    // =============================================================================
    // ----------------------------------------------------------------------
    // f0e function for isotropic linear Poroelasticity.
    static inline
    void f0e(const PylithInt dim,
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
        pylith::fekernels::Poroelasticity::Context Context;
        pylith::fekernels::Poroelasticity::setContextQuasistatic(
            &Context, dim, numS, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x, t, x);

        const PylithScalar* displacement_x = Context.displacement_x;
        const PylithScalar trace_strain = Context.trace_strain;

        for (PylithInt d = 0; d < dim; ++d) {
            f0[0] += displacement_x[d*dim+d];
        }
        f0[0] -= trace_strain;
    } // f0e

    // =============================================================================
    // Time Derivative of Pressure
    // =============================================================================
    // ----------------------------------------------------------------------
    /*
     * f0 function for poroelasticity equation, implicit time stepping, quasistatic, for
     * time derivative of pressure.
     */
    static inline
    void f0pdot(const PylithInt dim,
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
        pylith::fekernels::Poroelasticity::Context Context;
        pylith::fekernels::Poroelasticity::setContextQuasistatic(
            &Context, dim, numS, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x, t, x);
        pylith::fekernels::Poroelasticity::setContextQuasistatic_statevars(
            &Context, dim, numS, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x, t, x);

        // PylithScalar pressure_t = Context.pressure_t;
        // PylithScalar pressure_dot = Context.pressure_dot;

        PylithScalar pressure_t = (Context.pressure_t ? Context.pressure_t : 0.0);
        PylithScalar pressure_dot = (Context.pressure_dot ? Context.pressure_dot : 0.0);

        f0[0] += pressure_t;
        f0[0] -= pressure_dot;

    } // f0pdot

    // =============================================================================
    // Time Derivative of Volumetric Strain
    // =============================================================================
    // ----------------------------------------------------------------------
    /*
     * f0 function for poroelasticity equation, implicit time stepping, quasistatic, for
     * time derivative of volumetric strain.
     */
    static inline
    void f0edot(const PylithInt dim,
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
        // Poroelastic Context
        pylith::fekernels::Poroelasticity::Context Context;
        pylith::fekernels::Poroelasticity::setContextQuasistatic(
            &Context, dim, numS, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x, t, x);
        pylith::fekernels::Poroelasticity::setContextQuasistatic_statevars(
            &Context, dim, numS, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x, t, x);

        PylithScalar trace_strain_t = (Context.trace_strain_t ? Context.trace_strain_t : 0.0);
        PylithScalar trace_strain_dot = (Context.trace_strain_dot ? Context.trace_strain_dot : 0.0);

        f0[0] += trace_strain_t - trace_strain_dot;
    } // f0edot

    // --------------------------------------------------------------------------------------------
    // f1 function for darcy term.
    static inline
    void f1p(const Context& Context,
             void* rheologyContext,
             Poroelasticity::fluxratefn_type fluxRateFn,
             const TensorOps& tensorOps,
             PylithScalar f1[]) {
        Tensor fluxRate;
        fluxRateFn(Context, rheologyContext, tensorOps, &fluxRate);
        PylithScalar fluxRateTensor[9] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
        tensorOps.toTensor(fluxRate, fluxRateTensor);
        PylithInt dim = Context.dim;

        for (PylithInt i = 0; i < dim; ++i) {
            f1[i] += fluxRateTensor[i*dim+i];
        } // for
    } // f1p

    /* -------------------------------------------------------------------------- */
    /*                           RHS Residuals                                    */
    /* -------------------------------------------------------------------------- */
    // Quasi-Static

    // =============================================================================
    // Displacement
    // =============================================================================

    // ----------------------------------------------------------------------
    /*
     * g0 function for displacement equation: g0u = v.
     */
    static inline
    void g0u(const PylithInt dim,
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
        pylith::fekernels::Poroelasticity::Context Context;
        pylith::fekernels::Poroelasticity::setContextDynamic(
            &Context, dim, numS, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x, t, x);

        const PylithScalar *velocity = Context.velocity;

        for (PylithInt i = 0; i < dim; ++i) {
            g0[i] += velocity[i];
        } // for
    } // g0u

    // ---------------------------------------------------------------------------------------------------------------------
    /*
     * g0v_grav - g0 function for generic poroelasticity terms ( + grav body forces).
     */
    static inline
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
                  PylithScalar g0[]) {
        // Poroelastic Context
        pylith::fekernels::Poroelasticity::Context Context;
        pylith::fekernels::Poroelasticity::setContextDynamic(
            &Context, dim, numS, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x, t, x);
        pylith::fekernels::Poroelasticity::setContextGravity(
            &Context, dim, numS, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x, t, x);

        // Poroelastic Auxiliariies
        const PylithScalar bulkDensity = Context.bulkDensity;
        const PylithScalar* gravityField = Context.gravityField;

        for (PylithInt i = 0; i < dim; ++i) {
            g0[i] += bulkDensity * gravityField[i];
        } // for
    } // g0v_grav

    // ---------------------------------------------------------------------------------------------------------------------
    /*
     * g0v_bodyforce - g0 function for generic poroelasticity terms ( + body forces).
     */
    static inline
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
                       PylithScalar g0[]) {
        // Poroelastic Context
        pylith::fekernels::Poroelasticity::Context Context;
        pylith::fekernels::Poroelasticity::setContextDynamic(
            &Context, dim, numS, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x, t, x);
        pylith::fekernels::Poroelasticity::setContextBodyForce(
            &Context, dim, numS, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x, t, x);

        // Poroelastic Auxiliariies
        const PylithScalar* bodyForce = Context.bodyForce;

        for (PylithInt i = 0; i < dim; ++i) {
            g0[i] += bodyForce[i];
        } // for
    } // g0v_bodyforce

    // ----------------------------------------------------------------------
    /*
     * g0v_gravbodyforce - g0 function for isotropic linear Poroelasticity with both gravity and body forces.
     */
    static inline
    void g0v_grav_bodyforce(const PylithInt dim,
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
        // Poroelastic Context
        pylith::fekernels::Poroelasticity::Context Context;
        pylith::fekernels::Poroelasticity::setContextDynamic(
            &Context, dim, numS, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x, t, x);
        pylith::fekernels::Poroelasticity::setContextGravityBodyForce(
            &Context, dim, numS, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x, t, x);

        // Poroelastic Auxiliariies
        const PylithScalar bulkDensity = Context.bulkDensity;
        const PylithScalar *gravityField = Context.gravityField;
        const PylithScalar *bodyForce = Context.bodyForce;

        for (PylithInt i = 0; i < dim; ++i) {
            g0[i] += bulkDensity * gravityField[i] + bodyForce[i];
        } // for
    } // g0v_grav_bodyforce

    // =============================================================================
    // Pressure
    // =============================================================================

    /* -------------------------------------------------------------------------- */
    /*                           LHS Jacobian                                     */
    /* -------------------------------------------------------------------------- */

    // -----------------------------------------------------------------------------
    /*
     * Jg0ee - Jf0 function for isotropic linear poroelasticity.
     */
    static inline
    void Jf0ee(const PylithInt dim,
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
        Jf0[0] = -1.0;
    } // Jg0ee

    // -----------------------------------------------------------------------------
    /*
     * Jf1eu - Jf1 function for isotropic linear poroelasticity.
     */
    static inline
    void Jf1eu(const PylithInt dim,
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
    /*
     * Jf0vu function for poroelasticity equation, quasistatic.
     */
    static inline
    void Jf0vu_implicit(const PylithInt dim,
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
        for (PylithInt d = 0; d < dim; ++d) {
            Jf0[d * dim + d] += s_tshift;
        } // for
    } // Jf0vu_implicit

    // ---------------------------------------------------------------------------------------------------------------------
    /*
     * Jf0vv function for poroelasticity equation, quasistatic.
     */
    static inline
    void Jf0vv_implicit(const PylithInt dim,
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
        for (PylithInt d = 0; d < dim; ++d) {
            Jf0[d * dim + d] -= 1.0;
        } // for
    } // Jf0vv_implicit

    // ---------------------------------------------------------------------------------------------------------------------
    /*
     * Jf0vv function for poroelasticity equation, dynamic
     */
    static inline
    void Jf0vv_explicit(const PylithInt dim,
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
        pylith::fekernels::Poroelasticity::Context Context;
        pylith::fekernels::Poroelasticity::setContextDynamic(
            &Context, dim, numS, sOff, sOff_x, s, s_t, s_x, aOff, aOff_x, a, a_t, a_x, t, x);

        const PylithScalar bulkDensity = Context.bulkDensity;

        for (PylithInt i = 0; i < dim; ++i) {
            Jf0[i*dim+i] += s_tshift * bulkDensity;
        } // for
    } // Jf0vv_explicit

    // ---------------------------------------------------------------------------------------------------------------------
    /*
     * Jf0pdotp function for poroelasticity equation, quasistatic.
     */
    static inline
    void Jf0pdotp(const PylithInt dim,
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
        Jf0[0] += s_tshift;
    } // Jf0pdotp

    // ---------------------------------------------------------------------------------------------------------------------
    /*
     * Jf0pdotpdot function for poroelasticity equation, quasistatic.
     */

    static inline
    void Jf0pdotpdot(const PylithInt dim,
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
        Jf0[0] -= 1.0;
    } // Jg0pdotpdot

    // ---------------------------------------------------------------------------------------------------------------------
    /*
     * Jf0edote function for poroelasticity equation, quasistatic.
     */
    static inline
    void Jf0edote(const PylithInt dim,
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
        Jf0[0] += s_tshift;
    } // Jf0edote

    // ---------------------------------------------------------------------------------------------------------------------
    // Jf0edotedot function for poroelasticity equation, quasistatic.
    static inline
    void Jf0edotedot(const PylithInt dim,
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
        Jf0[0] -= 1.0;
    } // Jg0edotedot

}; // Poroelasticity

#endif // pylith_fekernels_poroelasticity_hh

// End of file
