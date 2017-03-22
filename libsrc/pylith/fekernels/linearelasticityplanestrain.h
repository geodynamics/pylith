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

/** @file libsrc/fekernels/linearelasticityplanestrain.h
 *
 * Kernels for linear elasticity plane strain.
 *
 * Solution fields: [disp(dim), vel(dim, optional)]
 *
 * Isotropic, linear elasticity plane strain without reference stress/strain.
 * Auxiliary fields:
 * - 0: density(1)
 * - 1: shear_odulus(1)
 * - 2: bulk_odulus(1)
 * - 3: gravity_field (2, optional)
 * - 4: body_force(2,optional)
 * - 5: reference_stress(4,optional) (stress_xx, stress_yy, stress_xy, stress_zz)
 * - 6: reference_strain(4,optional) (strain_xx, strain_yy, strain_xy, strain_zz)
 *
 * \int_V \vec{\phi}_u \cdot \left( \rho \frac{\partial \vec{v}(t)}{\partial t} \right) \, dV =
 *   \int_V \vec{\phi}_u \cdot \vec{f}(t) - \nabla \vec{\phi}_u : \tensor{\sigma}(\vec{u}) \, dV +
 *   \int_{S_\tau} \vec{\phi}_u \cdot \vec{\tau}(t) \, dS.
 *
 * ======================================================================
 */

#if !defined(pylith_fekernels_linearelasticityplanestrain_h)
#define pylith_fekernels_linearelasticityplanestrain_h

/* Include directives ---------------------------------------------------
 */
#include <portinfo>

#include "pylith/utils/types.hh"
#include "pylith/utils/error.h"

/** f0 function for isotropic linear elasticity plane strain.
 *
 * Solution fields: [disp(dim), vel(dim)]
 * Auxiliary fields: [density(1), ...]
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
pylith_fekernels_IsotropicLinearElasticityPlaneStrain_f0v(const PylithInt dim,
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


/** g0 function for isotropic linear elasticity plane strain with both gravity and body forces.
 *
 * Solution fields: [...]
 * Auxiliary fields: [density(1), shear_modulus(1), bulk_modulus(1), gravity_field(dim), body_force(dim), ...]
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
pylith_fekernels_IsotropicLinearElasticityPlaneStrain_g0v_gravbodyforce(const PylithInt dim,
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


/** g0 function for isotropic linear elasticity plane strain with gravity.
 *
 * Solution fields: [...]
 * Auxiliary fields: [density(1), shear_modulus(1), bulk_modulus(1), gravity_field(dim), ...]
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
pylith_fekernels_IsotropicLinearElasticityPlaneStrain_g0v_grav(const PylithInt dim,
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


/** g0 function for isotropic linear elasticity plane strain with body forces.
 *
 * Solution fields: [...]
 * Auxiliary fields: [density(1), shear_modulus(1), bulk_modulus(1), body_force(dim), ...]
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
pylith_fekernels_IsotropicLinearElasticityPlaneStrain_g0v_bodyforce(const PylithInt dim,
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


/** g1 function for isotropic linear elasticity plane strain WITHOUT reference stress and reference strain.
 *
 * Solution fields: [disp(dim), ...]
 * Auxiliary fields: [density(1), shear_modulus(1), bulk_modulus(1), ...]
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
 * @param[out] g1 [dim*dim].
 */
void
pylith_fekernels_IsotropicLinearElasticityPlaneStrain_g1v(const PylithInt dim,
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
                                                          PylithScalar g1[]);


/** g1 function for isotropic linear elasticity plane strain WITH reference stress and reference strain.
 *
 * Solution fields: [disp(dim), ...]
 * Auxiliary fields: [density(1), shear_modulus(1), bulk_modulus(1), ..., refstress(4), refstrain(4)]
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
 * @param[out] g1 [dim*dim].
 */
void
pylith_fekernels_IsotropicLinearElasticityPlaneStrain_g1v_refstate(const PylithInt dim,
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
                                                                   PylithScalar g1[]);


/** Jf0 function for isotropoc linear elasticity plane strain with implicit time stepping.
 *
 * Solution fields: [...]
 * Auxiliary fields: [density(1), ...]
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
 * @param[in] Jf0 Result [dim*dim].
 */
void
pylith_fekernels_IsotropicLinearElasticityPlaneStrain_Jf0vv_implicit(const PylithInt dim,
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


/** Jf0 function for isotropoc linear elasticity plane strain with explicit time stepping.
 *
 * Solution fields: [...]
 * Auxiliary fields: [density(1), shear_modulus(1), bulk_modulus(1), ..., refstress(4), refstrain(4)]
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
 * @param[in] Jf0 Result [dim*dim].
 */
void
pylith_fekernels_IsotropicLinearElasticityPlaneStrain_Jf0vv_explicit(const PylithInt dim,
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


/** Jg3_vu entry function for 2-D plane strain isotropic linear elasticity.
 *
 * Solution fields: [...]
 * Auxiliary fields: [shear_modulus(1), bulk_modulus(1), ...]
 *
 * @param[in] dim Spatial dimension [3].
 * @param[in] numS Number of registered subfields in solution field [1].
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
 * @param[in] utshift Coefficient for dF/ds_t term in Jacobian.
 * @param[in] x Coordinates of point evaluation.
 * @param[in] Jg3 Result [dim*dim*dim*dim].
 */
void
pylith_fekernels_IsotropicLinearElasticityPlaneStrain_Jg3vu(const PylithInt dim,
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
                                                            PylithScalar Jg3[]);


/** Calculate mean stress for 2-D plane strain isotropic linear
 * elasticity WITHOUT reference stress and reference strain.
 *
 * Solution fields: [disp(dim)]
 * Auxiliary fields: [bulk_modulus(1)]
 *
 * @param[in] dim Spatial dimension.
 * @param[in] numS Number of registered subfields in solution field [1].
 * @param[in] numA Number of registered subfields in auxiliary field [1].
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
 * @param[out] stress Result [dim*dim].
 */
void
pylith_fekernels_IsotropicLinearElasticityPlaneStrain_meanStress(const PylithInt dim,
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
                                                                 PylithScalar stress[]);


/** Calculate mean stress for 2-D plane strain isotropic linear
 * elasticity WITH reference stress and reference strain.
 *
 * Solution fields: [disp(dim)]
 * Auxiliary fields: [bulk_modulus(1), reference_stress(4), reference_strain(4)]
 *
 * @param[in] dim Spatial dimension.
 * @param[in] numS Number of registered subfields in solution field [1].
 * @param[in] numA Number of registered subfields in auxiliary field [1].
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
 * @param[out] stress Result [dim*dim].
 */
void
pylith_fekernels_IsotropicLinearElasticityPlaneStrain_meanStress_refstate(const PylithInt dim,
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
                                                                          PylithScalar stress[]);


/** Calculate deviatoric stress for 2-D plane strain isotropic linear
 * elasticity WITHOUT reference stress and strain.
 *
 * Solution fields: [disp(dim)]
 * Auxiliary fields: [shear_modulus(1)]
 *
 * @param[in] dim Spatial dimension.
 * @param[in] numS Number of registered subfields in solution field [1].
 * @param[in] numA Number of registered subfields in auxiliary field [1].
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
 * @param[out] stress Result [dim*dim].
 */
void
pylith_fekernels_IsotropicLinearElasticityPlaneStrain_deviatoricStress(const PylithInt dim,
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
                                                                       PylithScalar stress[]);


/** Calculate deviatoric stress for 2-D plane strain isotropic linear
 * elasticity WITH reference stress and strain.
 *
 * Solution fields: [disp(dim)]
 * Auxiliary fields: [shear_modulus(1), reference_stress(4), reference_strain(4)]
 *
 * @param[in] dim Spatial dimension.
 * @param[in] numS Number of registered subfields in solution field [1].
 * @param[in] numA Number of registered subfields in auxiliary field [1].
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
 * @param[out] stress Result [dim*dim].
 *
 */
void
pylith_fekernels_IsotropicLinearElasticityPlaneStrain_deviatoricStress_refstate(const PylithInt dim,
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
                                                                                PylithScalar stress[]);


#endif /* pylith_fekernels_linearelasticityplanestrain_h */


/* End of file */
