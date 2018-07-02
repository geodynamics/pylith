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
 * Copyright (c) 2010-2017 University of California, Davis
 *
 * See COPYING for license information.
 *
 * ----------------------------------------------------------------------
 */

/** @file libsrc/fekernels/LinearElasticity3D.hh
 *
 * Kernels for linear elasticity 3D.
 *
 * Solution fields: [disp(dim), vel(dim, optional)]
 *
 * Auxiliary fields:
 * - 0: density(1)
 * - 1: shear_modulus(1)
 * - 2: bulk_modulus(1)
 * - 3: gravity_field (3, optional)
 * - 4: body_force(3,optional)
 * - 5: reference_stress(6,optional) (stress_xx, stress_yy, stress_zz, stress_xy, stress_yz, stress_xz)
 * - 6: reference_strain(6,optional) (strain_xx, strain_yy, strain_zz, strain_xy, strain_yz, strain_xz)
 *
 * \int_V \vec{\phi}_u \cdot \left( \rho \frac{\partial \vec{v}(t)}{\partial t} \right) \, dV =
 *   \int_V \vec{\phi}_u \cdot \vec{f}(t) - \nabla \vec{\phi}_u : \tensor{\sigma}(\vec{u}) \, dV +
 *   \int_{S_\tau} \vec{\phi}_u \cdot \vec{\tau}(t) \, dS.
 *
 * ======================================================================
 */

#if !defined(pylith_fekernels_isotropiclinearelasticity3d_hh)
#define pylith_fekernels_isotropiclinearelasticity3d_hh

// Include directives ---------------------------------------------------
#include "fekernelsfwd.hh" // forward declarations

#include "pylith/utils/types.hh"

class pylith::fekernels::IsotropicLinearElasticity3D {

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


    /** g0 function for isotropic linear elasticity 3D with both gravity and body forces.
     *
     * Solution fields: [...]
     * Auxiliary fields: [density(1), shear_modulus(1), bulk_modulus(1), gravity_field(dim), body_force(dim), ...]
     */
    // MOVE TO materials::IsotropicLinearElasticity3D
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


    /** g0 function for isotropic linear elasticity 3D with gravity.
     *
     * Solution fields: [...]
     * Auxiliary fields: [density(1), shear_modulus(1), bulk_modulus(1), gravity_field(dim), ...]
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


    /** g0 function for isotropic linear elasticity 3D with body forces.
     *
     * Solution fields: [...]
     * Auxiliary fields: [density(1), shear_modulus(1), bulk_modulus(1), body_force(dim), ...]
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


    /** g1 function for isotropic linear elasticity 3D WITHOUT reference stress and reference strain.
     *
     * Solution fields: [disp(dim), ...]
     * Auxiliary fields: [density(1), shear_modulus(1), bulk_modulus(1), ...]
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


    /** g1 function for isotropic linear elasticity 3D WITH reference stress and reference strain.
     *
     * Solution fields: [disp(dim), ...]
     * Auxiliary fields: [density(1), shear_modulus(1), bulk_modulus(1), ..., refstress(6), refstrain(6)]
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


    /** Jg3_vu entry function for 3-D isotropic linear elasticity.
     *
     * Solution fields: [...]
     * Auxiliary fields: [density(1), shear_modulus(1), bulk_modulus(1), ...]
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

    /** Calculate stress for 3-D isotropic linear
     * elasticity WITHOUT a reference stress and strain.
     *
     * Used to output the stress field.
     *
     * Solution fields: [disp(dim)]
     * Auxiliary fields: [density(1), shear_modulus(1), bulk_modulus(1), ...]
     */
    static
    void stress(const PylithInt dim,
                const PylithInt numS,
                const PylithInt numA,
                const PylithInt sOff[],
                const PylithInt sOff_x[],
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

    /** Calculate stress for 3-D isotropic linear
     * elasticity WITH a reference stress/strain.
     *
     * Used to output the stress field.
     *
     * Solution fields: [disp(dim)]
     * Auxiliary fields: [density(1), shear_modulus(1), bulk_modulus(1), ..., refstress(6), refstrain(6)]
     */
    static
    void stress_refstate(const PylithInt dim,
                         const PylithInt numS,
                         const PylithInt numA,
                         const PylithInt sOff[],
                         const PylithInt sOff_x[],
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

}; // IsotropicLinearElasticity3D

#endif // pylith_fekernels_isotropiclinearelasticity3d_hh


// End of file
