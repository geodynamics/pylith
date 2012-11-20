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
// Copyright (c) 2010-2011 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

/** @file libsrc/friction/SlipWeakeningStress.hh
 *
 * @brief C++ slip weakening fault constitutive model with
 * stress-based friction and additional time parameter for nucleation.
 */

#if !defined(pylith_friction_slipweakeningstress_hh)
#define pylith_friction_slipweakeningstress_hh

// Include directives ---------------------------------------------------
#include "FrictionModel.hh" // ISA FrictionModel

// SlipWeakeningStress -------------------------------------------------------
/** @brief C++ slip weakening fault constitutive model.
 *
 * Friction is equal stress as a function of slip.
 */

class pylith::friction::SlipWeakeningStress : public FrictionModel
{ // class SlipWeakeningStress
  friend class TestSlipWeakeningStress; // unit testing

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Default constructor.
  SlipWeakeningStress(void);

  /// Destructor.
  ~SlipWeakeningStress(void);

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
  static const int p_stressStatic;
  static const int p_stressDyn;
  static const int p_d0;
  static const int p_weaktime;

  static const int db_stressStatic;
  static const int db_stressDyn;
  static const int db_d0;
  static const int db_weaktime;

  /// Indices for state variables in section and spatial database.
  static const int s_slipCum;
  static const int s_slipPrev;

  static const int db_slipCum;
  static const int db_slipPrev;

  // NOT IMPLEMENTED ////////////////////////////////////////////////////
private :

  SlipWeakeningStress(const SlipWeakeningStress&); ///< Not implemented.
  const SlipWeakeningStress& operator=(const SlipWeakeningStress&); ///< Not implemented

}; // class SlipWeakeningStress

#endif // pylith_friction_slipweakeningstress_hh


// End of file 
