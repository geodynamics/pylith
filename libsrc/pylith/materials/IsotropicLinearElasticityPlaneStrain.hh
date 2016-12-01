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

/** @file libsrc/materials/IsotropicLinearElasticityPlaneStrain.hh
 *
 * @brief C++ class for isotropic linear elastic plane strain material.
 */

#if !defined(pylith_materials_isotropiclinearelasticityplanestrain_hh)
#define pylith_materials_isotropiclinearelasticityplanestrain_hh

// Include directives ---------------------------------------------------
#include "materialsfwd.hh" // forward declarations

#include "pylith/materials/MaterialNew.hh" // ISA Material

// Material -------------------------------------------------------------
/** @brief C++ class for isotropic linear elastic plane strain material.
 */

class pylith::materials::IsotropicLinearElasticityPlaneStrain : public pylith::materials::MaterialNew
{ // class IsotropicLinearElasticityPlaneStrain
friend class TestIsotropicLinearElasticityPlaneStrain;   // unit testing

// PUBLIC METHODS /////////////////////////////////////////////////////
public:

/// Default constructor.
IsotropicLinearElasticityPlaneStrain(void);

/// Destructor.
~IsotropicLinearElasticityPlaneStrain(void);

/** Include inertia?
 *
 * @param[in] value Flag indicating to include inertial term.
 */
void useInertia(const bool value);

/** Include inertia?
 *
 * @returns True if including inertial term, false otherwise.
 */
bool useInertia(void) const;

/** Include body force?
 *
 * @param[in] value Flag indicating to include body force term.
 */
void useBodyForce(const bool value);

/** Include body force?
 *
 * @returns True if including body force term, false otherwise.
 */
bool useBodyForce(void) const;

/** Use reference stress and strain in computation of stress and
 * strain?
 *
 * @param[in] value Flag indicating to include reference stress and strain.
 */
void useReferenceState(const bool value);

/** Use reference stress and strain in computation of stress and
 * strain?
 *
 * @returns True if using reference stress and strain, false otherwise.
 */
bool useReferenceState(void) const;

// PROTECTED METHODS //////////////////////////////////////////////////
protected:

/// Setup auxiliary subfields (discretization and query fns).
void _auxFieldsSetup(void);

/** Set kernels for RHS residual G(t,u).
 *
 * @param[in] solution Solution field.
 */
void _setFEKernelsRHSResidual(const topology::Field& solution);

/** Set kernels for RHS Jacobian G(t,u).
 *
 * @param[in] solution Solution field.
 */
void _setFEKernelsRHSJacobian(const topology::Field& solution);

/** Set kernels for LHS residual F(t,u,\dot{u}).
 *
 * @param[in] solution Solution field.
 */
void _setFEKernelsLHSResidual(const topology::Field& solution);


/** Set kernels for LHS Jacobian F(t,u,\dot{u}) when implicit time-stepping.
 *
 * @param[in] solution Solution field.
 */
void _setFEKernelsLHSJacobianImplicit(const topology::Field& solution);


/** Set kernels for LHS Jacobian F(t,u,\dot{u}) when explicit time-stepping.
 *
 * @param[in] solution Solution field.
 */
void _setFEKernelsLHSJacobianExplicit(const topology::Field& solution);


// PRIVATE MEMBERS ////////////////////////////////////////////////////
private:

bool _useInertia;   ///< Flag to include inertial term.
bool _useBodyForce;   ///< Flag to include body force term.
bool _useReferenceState;   ///< Flag to use reference stress and strain.

static const char* _journal; ///< Name for journals.

// NOT IMPLEMENTED ////////////////////////////////////////////////////
private:

IsotropicLinearElasticityPlaneStrain(const IsotropicLinearElasticityPlaneStrain&);   ///< Not implemented.
const IsotropicLinearElasticityPlaneStrain& operator=(const IsotropicLinearElasticityPlaneStrain&);   ///< Not implemented

}; // class IsotropicLinearElasticityPlaneStrain

#endif // pylith_materials_isotropiclinearelasticityplanestrain_hh


// End of file
