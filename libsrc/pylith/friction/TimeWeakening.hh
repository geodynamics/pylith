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

/** @file libsrc/friction/TimeWeakening.hh
 *
 * @brief C++ time weakening fault constitutive model.
 */

#if !defined(pylith_friction_timeweakening_hh)
#define pylith_friction_timeweakening_hh

// Include directives ---------------------------------------------------
#include "FrictionModel.hh" // ISA FrictionModel

// TimeWeakening -------------------------------------------------------
/** @brief C++ time weakening fault constitutive model.
 *
 * Friction is equal to the product of a coefficient of friction
 * (function of time since slipping started) and the normal traction.
 */

class pylith::friction::TimeWeakening : public FrictionModel
{ // class TimeWeakening
  friend class TestTimeWeakening; // unit testing

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Default constructor.
  TimeWeakening(void);

  /// Destructor.
  ~TimeWeakening(void);

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
   * @param slip Current slip at location.
   * @param slipRate Current slip rate at location.
   * @param normalTraction Normal traction at location.
   * @param properties Properties at location.
   * @param numProperties Number of properties.
   * @param stateVars State variables at location.
   * @param numStateVars Number of state variables.
   */
  void _dbToStateVars(PylithScalar* const stateValues,
		      const scalar_array& dbValues) const;

  /** Nondimensionalize state variables.
   *
   * @param values Array of initial state values.
   * @param nvalues Number of values.
   */
  void _nondimStateVars(PylithScalar* const values,
			   const int nvalues) const;
  
  /** Dimensionalize state variables.
   *
   * @param values Array of initial state values.
   * @param nvalues Number of values.
   */
  void _dimStateVars(PylithScalar* const values,
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

  /** Update state variables (for next time step).
   *
   * @param slip Current slip at location.
   * @param slipRate Current slip rate at location.
   * @param normalTraction Normal traction at location.
   * @param stateVars State variables at location.
   * @param numStateVars Number of state variables.
   * @param properties Properties at location.
   * @param numProperties Number of properties.
   */
  void _updateStateVars(const PylithScalar t,
			const PylithScalar slip,
			const PylithScalar slipRate,
			const PylithScalar normalTraction,
			PylithScalar* const stateVars,
			const int numStateVars,
			const PylithScalar* properties,
			const int numProperties);

  // PRIVATE MEMBERS ////////////////////////////////////////////////////
private :

  /// Indices for properties in section and spatial database.
  static const int p_coefS;
  static const int p_coefD;
  static const int p_Tc;
  static const int p_cohesion;

  static const int db_coefS;
  static const int db_coefD;
  static const int db_Tc;
  static const int db_cohesion;

  /// Indices for state variables in section and spatial database.
  static const int s_time;

  static const int db_time;

  // NOT IMPLEMENTED ////////////////////////////////////////////////////
private :

  TimeWeakening(const TimeWeakening&); ///< Not implemented.
  const TimeWeakening& operator=(const TimeWeakening&); ///< Not implemented

}; // class TimeWeakening

#endif // pylith_friction_timeweakening_hh


// End of file 
