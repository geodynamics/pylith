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

/** @file libsrc/friction/SlipWeakeningTimeStable.hh
 *
 * @brief C++ slip + time weakening fault constitutive model.
 *
 * This constitutive model is used in the SCEC spontaneous rupture
 * benchmarks suite.
 *
 * This version of slip + time weakening contains a time constant to
 * make it more stable then SlipWeakeningTime.
 */

#if !defined(pylith_friction_slipweakeningtimestable_hh)
#define pylith_friction_slipweakeningtimestable_hh

// Include directives ---------------------------------------------------
#include "FrictionModel.hh" // ISA FrictionModel

// SlipWeakeningTimeStable -------------------------------------------------------
/** @brief C++ slip + time weakening fault constitutive model.
 *
 * Friction is equal to the product of a coefficient of friction (function
 * of slip path length) and the normal traction.
 *
 * This version of slip + time weakening contains a time constant to
 * make it more stable then SlipWeakeningTime.
 */

class pylith::friction::SlipWeakeningTimeStable : public FrictionModel
{ // class SlipWeakeningTimeStable
  friend class TestSlipWeakeningTimeStable; // unit testing

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Default constructor.
  SlipWeakeningTimeStable(void);

  /// Destructor.
  ~SlipWeakeningTimeStable(void);

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
   * @param t Time in simulation.
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
  static const int p_d0;
  static const int p_cohesion;
  static const int p_weaktime;
  static const int p_t0;

  static const int db_coefS;
  static const int db_coefD;
  static const int db_d0;
  static const int db_cohesion;
  static const int db_weaktime;
  static const int db_t0;

  /// Indices for state variables in section and spatial database.
  static const int s_slipCum;
  static const int s_slipPrev;

  static const int db_slipCum;
  static const int db_slipPrev;

  // NOT IMPLEMENTED ////////////////////////////////////////////////////
private :

  SlipWeakeningTimeStable(const SlipWeakeningTimeStable&); ///< Not implemented.
  const SlipWeakeningTimeStable& operator=(const SlipWeakeningTimeStable&); ///< Not implemented

}; // class SlipWeakeningTimeStable

#endif // pylith_friction_slipweakeningtimestable_hh


// End of file 
