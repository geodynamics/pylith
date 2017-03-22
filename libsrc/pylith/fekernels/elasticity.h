/* -*- C -*-
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

/** @file libsrc/fekernels/elasticity.h
 *
 * Generic elasticity kernels for inertia and body forces.
 *
 * Solution fields: [disp(dim), vel(dim)]
 *
 * \int_V \vec{\phi}_u \cdot \left( \rho \frac{\partial \vec{v}(t)}{\partial t} \right) \, dV =
 *   \int_V \vec{\phi}_u \cdot \vec{f}(t) - \nabla \vec{\phi}_u : \tensor{\sigma}(\vec{u}) \, dV +
 *   \int_{S_\tau} \vec{\phi}_u \cdot \vec{\tau}(t) \, dS.
 *
 * Note: Ignore stress term because it depends on the constitutive model.
 */

#if !defined(pylith_fekernels_elasticity_h)
#define pylith_fekernels_elasticity_h

/* Include directives ---------------------------------------------------
 */
#include <portinfo>

#include "pylith/utils/types.hh"
#include "pylith/utils/error.h"

/** f0 function for generic elasticity terms (inertia).
 *
 * \vec{f0} = \rho \frac{\partial \vec{v}(t)}{\partial t}
 *
 * Auxiliary fields: [density(1)]
 *
 * @param[in] dim Spatial dimension.
 * @param[in] numS Number of registered subfields in solution field [2].
 * @param[in] numA Number of registered subfields in auxiliary field [2].
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
 * @param[out] f0 [dim].
 */
void
pylith_fekernels_Elasticity_f0v_inertia(const PylithInt dim,
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
                                        PylithScalar f0[]);


/** g0 function for generic elasticity terms (gravitational body force).
 *
 * g0_Elasticity: \vec{g0} = \vec{f}(t)
 *
 * Solution fields: []
 * Auxiliary fields: [density, gravity_field(dim)]
 *
 * @param[in] dim Spatial dimension.
 * @param[in] numS Number of registered subfields in solution field [2].
 * @param[in] numA Number of registered subfields in auxiliary field [2].
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
 * @param[out] g0 [dim].
 */
void
pylith_fekernels_Elasticity_g0v_grav(const PylithInt dim,
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
                                     PylithScalar g0[]);


/** g0 function for generic elasticity terms (body force).
 *
 * g0_Elasticity: \vec{g0} = \vec{f}(t)
 *
 * Auxiliary fields: [body_force(dim)]
 *
 * @param[in] dim Spatial dimension.
 * @param[in] numS Number of registered subfields in solution field [2].
 * @param[in] numA Number of registered subfields in auxiliary field [2].
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
 * @param[out] g0 [dim].
 */
void
pylith_fekernels_Elasticity_g0v_bodyforce(const PylithInt dim,
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
                                          PylithScalar g0[]);


/** Jf0 function for generic elasticity terms (inertia) with implicit time stepping.
 *
 * Auxiliary fields: [density(1)]
 *
 * @param[in] dim Spatial dimension.
 * @param[in] numS Number of registered subfields in solution field [2].
 * @param[in] numA Number of registered subfields in auxiliary field [0].
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
 * @param[in] utshift Coefficient for dF/ds_t term in Jacobian.
 * @param[in] x Coordinates of point evaluation.
 * @param[out] Jf0 [dim*dim].
 */
void
pylith_fekernels_Elasticity_Jf0vv_inertiaimplicit(const PylithInt dim,
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
                                                  PylithScalar Jf0[]);




/** Jf0 function for generic elasticity terms (inertia) with explicit time stepping.
 *
 * Auxiliary fields: [density(1)]
 *
 * @param[in] dim Spatial dimension.
 * @param[in] numS Number of registered subfields in solution field [2].
 * @param[in] numA Number of registered subfields in auxiliary field [0].
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
 * @param[in] utshift Coefficient for dF/ds_t term in Jacobian.
 * @param[in] x Coordinates of point evaluation.
 * @param[out] Jf0 [dim*dim].
 */
void
pylith_fekernels_Elasticity_Jf0vv_inertiaexplicit(const PylithInt dim,
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
                                                  PylithScalar Jf0[]);

#endif /* pylith_fekernels_elasticity_h */

/* End of file */
