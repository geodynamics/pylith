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
// Copyright (c) 2010-2017 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

/* @brief C++ ViscousFriction object that implements friction
 * correlated with slip rate.
 *
 * This objects demonstrates how to extend PyLith by adding fault
 * constitutive models.
 *
 * Friction model in which the perturbation from the static
 * coefficient of friction is proportional to slip rate. The physical
 * properties are specified using cohesion, static coefficient of
 * friction, and a reference slip rate.
 *
 * $\mu_f = \mu_s (1 + \dot{D} / v_0)
 */

#if !defined(pylith_friction_viscousfriction_hh)
#define pylith_friction_viscousfriction_hh

// Include directives ---------------------------------------------------
#include "pylith/friction/FrictionModel.hh" // ISA FrictionModel

// Forward declarations
namespace contrib {
  namespace friction {
    class ViscousFriction;
  } // friction
} // pylith

// ViscousFriction -------------------------------------------------------
class contrib::friction::ViscousFriction : public pylith::friction::FrictionModel
{ // class ViscousFriction
  friend class TestViscousFriction; // unit testing

  // --------------------------------------------------------------------
  // All of these functions are required to satisfy the PyLith
  // interface for a fault constitutive model.
  // --------------------------------------------------------------------

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Default constructor.
  ViscousFriction(void);

  /// Destructor.
  ~ViscousFriction(void);

  // PROTECTED METHODS //////////////////////////////////////////////////
protected :

  /// These methods should be implemented by every constitutive model.

  /** Compute properties from values in spatial database.
   *
   * @param propValues Array of property values.
   * @param dbValues Array of database values.
   */
  void _dbToProperties(PylithScalar* const propValues,
		       const pylith::scalar_array& dbValues) const;

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
		      const pylith::scalar_array& dbValues) const;

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

  /** Compute derivative friction with slip from properties and state
   * variables.
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
  
  // --------------------------------------------------------------------
  // Optional function in the PyLith interface for a fault
  // constitutive model. Even though this function is optional, for it
  // to be used it the interface must exactly matched the one
  // specified in FrictionModel.
  // --------------------------------------------------------------------

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

  // --------------------------------------------------------------------
  // We use these constants for consistent access into the arrays of
  // physical properties and state variables.
  // --------------------------------------------------------------------

  static const int p_coefS;
  static const int p_v0;
  static const int p_cohesion;
  static const int db_coefS;
  static const int db_v0;
  static const int db_cohesion;

  static const int s_slipRate;
  static const int db_slipRate;

  // NOT IMPLEMENTED ////////////////////////////////////////////////////
private :

  ViscousFriction(const ViscousFriction&); ///< Not implemented.
  const ViscousFriction& operator=(const ViscousFriction&); ///< Not implemented

}; // class ViscousFriction

#endif // pylith_friction_viscousfriction_hh


// End of file 
