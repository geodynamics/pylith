// -*- C++ -*-
//
// ----------------------------------------------------------------------
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University of Chicago
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2015 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

/** @file libsrc/materials/IsotropicLinearIncompElasticityPlaneStrain.hh
 *
 * @brief C++ class for isotropic linear incompressible elastic plane strain material.
 */

#if !defined(pylith_materials_isotropiclinearincompelasticityplanestrain_hh)
#define pylith_materials_isotropiclinearincompelasticityplanestrain_hh

// Include directives ---------------------------------------------------
#include "materialsfwd.hh" // forward declarations

#include "pylith/materials/MaterialNew.hh" // ISA Material

// Material -------------------------------------------------------------
/** @brief C++ class for isotropic linear incompressible elastic plane strain material.
 */

class pylith::materials::IsotropicLinearIncompElasticityPlaneStrain : public pylith::materials::MaterialNew
{ // class IsotropicLinearIncompElasticityPlaneStrain
  friend class TestIsotropicLinearElasticityPlaneStrain; // unit testing

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Default constructor.
  IsotropicLinearIncompElasticityPlaneStrain(void);

  /// Destructor.
  ~IsotropicLinearIncompElasticityPlaneStrain(void);

  /** Include inertia?
   *
   * @param value Flag indicating to include inertial term.
   */
  void useInertia(const bool value);

  /** Include body force?
   *
   * @param value Flag indicating to include body force term.
   */
  void useBodyForce(const bool value);

  /** Include initial stress/strain?
   *
   * @param value Flag indicating to include initial stress/strain.
   */
  void useInitialState(const bool value);

  // PROTECTED METHODS //////////////////////////////////////////////////
protected :

  /// Setup auxiliary subfields (discretization and query fns).
  void _auxFieldsSetup(void);
 
  /** Set kernels for RHS residual G(t,u).
   *
   * @param[in] solution Solution field.
   */
  void _setFEKernelsRHSResidual(const topology::Field& solution) const;

  /** Set kernels for RHS Jacobian G(t,u).
   *
   * @param[in] solution Solution field.
   */
  void _setFEKernelsRHSJacobian(const topology::Field& solution) const ;

  /** Set kernels for LHS residual F(t,u,\dot{u}).
   *
   * @param[in] solution Solution field.
   */
  void _setFEKernelsLHSResidual(const topology::Field& solution) const;


  /** Set kernels for LHS Jacobian F(t,u,\dot{u}) when implicit time-stepping.
   *
   * @param[in] solution Solution field.
   */
  void _setFEKernelsLHSJacobianImplicit(const topology::Field& solution) const;


  /** Set kernels for LHS Jacobian F(t,u,\dot{u}) when explicit time-stepping.
   *
   * @param[in] solution Solution field.
   */
  void _setFEKernelsLHSJacobianExplicit(const topology::Field& solution) const;


  // PRIVATE MEMBERS ////////////////////////////////////////////////////
private :

  bool _useInertia; ///< Flag to include inertial term.
  bool _useBodyForce; ///< Flag to include body force term.
  bool _useInitialState; ///< Flag to include initial stress/strain.

  // NOT IMPLEMENTED ////////////////////////////////////////////////////
private :

  IsotropicLinearIncompElasticityPlaneStrain(const IsotropicLinearIncompElasticityPlaneStrain&); ///< Not implemented.
  const IsotropicLinearIncompElasticityPlaneStrain& operator=(const IsotropicLinearIncompElasticityPlaneStrain&); ///< Not implemented

}; // class IsotropicLinearIncompElasticityPlaneStrain

#endif // pylith_materials_isotropiclinearincompelasticityplanestrain_hh


// End of file 
