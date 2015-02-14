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

/** @file libsrc/friction/StaticFriction.hh
 *
 * @brief C++ static friction fault constitutive model.
 */

#if !defined(pylith_friction_staticfriction_hh)
#define pylith_friction_staticfriction_hh

// Include directives ---------------------------------------------------
#include "FrictionModel.hh" // ISA FrictionModel

// StaticFriction -------------------------------------------------------
/** @brief C++ static friction fault constitutive model.
 *
 * Friction is equal to the product of a coefficient of friction and
 * the normal traction.
 */

class pylith::friction::StaticFriction : public FrictionModel
{ // class StaticFriction
  friend class TestStaticFriction; // unit testing

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Default constructor.
  StaticFriction(void);

  /// Destructor.
  ~StaticFriction(void);

  // PROTECTED METHODS //////////////////////////////////////////////////
protected :

  /// These methods should be implemented by every constitutive model.

  /** Compute properties from values in spatial database.
   *
   * @param propValues Array of property values.
   * @param dbValues Array of database values.
   */
  void _dbToProperties(PylithScalar* const propValues,
		       const scalar_array& dbValues) const;

  /** Nondimensionalize properties.
   *
   * @param values Array of property values.
   * @param nvalues Number of values.
   */
  void _nondimProperties(PylithScalar* const values,
			 const int nvalues) const;

  /** Dimensionalize properties.
   *
   * @param values Array of property values.
   * @param nvalues Number of values.
   */
  void _dimProperties(PylithScalar* const values,
		      const int nvalues) const;

  /** Compute friction from properties and state variables.
   *
   * @param t Time in simulation.
   * @param slip Current slip at location.
   * @param slipRate Current slip rate at location.
   * @param normalTraction Normal traction at location.
   * @param properties Properties at location.
   * @param numProperties Number of properties.
   * @param stateVars State variables at location.
   * @param numStateVars Number of state variables.
   */
  PylithScalar _calcFriction(const PylithScalar t,
			     const PylithScalar slip,
			     const PylithScalar slipRate,
			     const PylithScalar normalTraction,
			     const PylithScalar* properties,
			     const int numProperties,
			     const PylithScalar* stateVars,
			     const int numStateVars);

  /** Compute derivative of friction with slip from properties and
   * state variables.
   *
   * @param t Time in simulation.
   * @param slip Current slip at location.
   * @param slipRate Current slip rate at location.
   * @param normalTraction Normal traction at location.
   * @param properties Properties at location.
   * @param numProperties Number of properties.
   * @param stateVars State variables at location.
   * @param numStateVars Number of state variables.
   *
   * @returns Derivative of friction (magnitude of shear traction) at vertex.
   */
  PylithScalar _calcFrictionDeriv(const PylithScalar t,
				  const PylithScalar slip,
				  const PylithScalar slipRate,
				  const PylithScalar normalTraction,
				  const PylithScalar* properties,
				  const int numProperties,
				  const PylithScalar* stateVars,
				  const int numStateVars);

  // PRIVATE MEMBERS ////////////////////////////////////////////////////
private :

  static const int p_coef;
  static const int p_cohesion;

  static const int db_coef;
  static const int db_cohesion;

  // NOT IMPLEMENTED ////////////////////////////////////////////////////
private :

  StaticFriction(const StaticFriction&); ///< Not implemented.
  const StaticFriction& operator=(const StaticFriction&); ///< Not implemented

}; // class StaticFriction

#endif // pylith_friction_staticfriction_hh


// End of file 
