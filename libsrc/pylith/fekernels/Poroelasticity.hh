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

class pylith::fekernels::Poroelasticity {

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
    /** g0 function for generic poroelasticity terms (gravitational body force).
     *
     * g0_Poroelasticity: \vec{g0} = \vec{f}(t)
     *
     * Solution fields: []
     * Auxiliary fields: [density, gravity_field(dim)]
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


    /** g0 function for generic poroelasticity terms (body force).
     *
     * g0_Poroelasticity: \vec{g0} = \vec{f}(t)
     *
     * Auxiliary fields: [body_force(dim)]
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

     // ----------------------------------------------------------------------
     //g0v_gravbodyforce - g0 function for isotropic linear Poroelasticity plane strain with both gravity and body forces.
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


     /** ----------------------------------------------------------------------
     * g0p function for isotropic linear poroelasticity plane strain without gravity
     * or body force.
     */
     static
     void g0p_source(const PylithInt dim,
                 const PylithInt numS,
                 const PylithInt numA,
                 const PylithInt sOff[],
                 const PylithInt sOff_x[],
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

     // ------------------------------------------------------------------------------
     // g0p function for isotropic linear Poroelasticity plane strain with source density
     static
     void g0p_sourceDensity(const PylithInt dim,
                                      const PylithInt numS,
                                      const PylithInt numA,
                                      const PylithInt sOff[],
                                      const PylithInt sOff_x[],
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

     // ------------------------------------------------------------------------------
     // g0p function for isotropic linear Poroelasticity plane strain with source density, gravity
     static
     void g0p_sourceDensity_grav(const PylithInt dim,
                                      const PylithInt numS,
                                      const PylithInt numA,
                                      const PylithInt sOff[],
                                      const PylithInt sOff_x[],
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




     // ------------------------------------------------------------------------------
     // g0p function for isotropic linear Poroelasticity plane strain with source density, and body force.
     static
     void g0p_sourceDensity_body(const PylithInt dim,
                                      const PylithInt numS,
                                      const PylithInt numA,
                                      const PylithInt sOff[],
                                      const PylithInt sOff_x[],
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

      // ------------------------------------------------------------------------------
      // g0p function for isotropic linear Poroelasticity plane strain with source density, gravity, and body force.
      static
      void g0p_sourceDensity_grav_body(const PylithInt dim,
                                       const PylithInt numS,
                                       const PylithInt numA,
                                       const PylithInt sOff[],
                                       const PylithInt sOff_x[],
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
     // g0e function for isotropic linear Poroelasticity.
     static
     void g0e(const PylithInt dim,
              const PylithInt numS,
              const PylithInt numA,
              const PylithInt sOff[],
              const PylithInt sOff_x[],
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
              PylithScalar g0e[]);

      // -----------------------------------------------------------------------------
      //Jg0ee - Jg0 function for isotropic linear poroelasticity plane strain.
      static
        void Jg0ee(const PylithInt dim,
                  const PylithInt numS,
                  const PylithInt numA,
                  const PylithInt sOff[],
                  const PylithInt sOff_x[],
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
                  PylithScalar Jg0[]);

      // -----------------------------------------------------------------------------
      // Jg1eu - Jg1 function for isotropic linear poroelasticity plane strain.
      static
        void Jg1eu(const PylithInt dim,
                  const PylithInt numS,
                  const PylithInt numA,
                  const PylithInt sOff[],
                  const PylithInt sOff_x[],
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
                  PylithScalar Jg1[]);
};

// Poroelasticity

// ---------------------------------------------------------------------------------------------------------------------

/// Kernels specific to poroelasticity plane strain.
class pylith::fekernels::PoroelasticityPlaneStrain {
  // PUBLIC MEMBERS //////////////////////////////////////////////////////////////////////////////////////////////////
public:

  /** Calculate Cauchy strain for 2-D plane strain poroelasticity.
   *
   * Order of output components is xx, yy, zz, xy.
   *
   * Solution fields: [disp(dim)]
   */
static
  void cauchyStrain(const PylithInt dim,
                    const PylithInt numS,
                    const PylithInt numA,
                    const PylithInt sOff[],
                    const PylithInt sOff_x[],
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
class pylith::fekernels::Poroelasticity3D {
  // PUBLIC MEMBERS //////////////////////////////////////////////////////////////////////////////////////////////////
public:

  /** Calculate Cauchy strain for 3-D poroelasticity.
   *
   * Order of output components is xx, yy, zz, xy, yz, xz.
   *
   * Solution fields: [disp(dim)]
   */
static
  void cauchyStrain(const PylithInt dim,
                    const PylithInt numS,
                    const PylithInt numA,
                    const PylithInt sOff[],
                    const PylithInt sOff_x[],
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
