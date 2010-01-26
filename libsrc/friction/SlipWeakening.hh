// -*- C++ -*-
//
// ----------------------------------------------------------------------
//
//                           Brad T. Aagaard
//                        U.S. Geological Survey
//
// {LicenseText}
//
// ----------------------------------------------------------------------
//

/** @file libsrc/friction/SlipWeakening.hh
 *
 * @brief C++ slip weakening fault constitutive model.
 */

#if !defined(pylith_friction_slipweakening_hh)
#define pylith_friction_slipweakening_hh

// Include directives ---------------------------------------------------
#include "FrictionModel.hh" // ISA FrictionModel

// SlipWeakening -------------------------------------------------------
/** @brief C++ slip weakening fault constitutive model.
 *
 * Friction is equal to the product of a coefficient of friction (function
 * of slip path length) and the normal traction.
 */

class pylith::friction::SlipWeakening : public FrictionModel
{ // class SlipWeakening
  friend class TestSlipWeakening; // unit testing

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Default constructor.
  SlipWeakening(void);

  /// Destructor.
  ~SlipWeakening(void);

  // PROTECTED METHODS //////////////////////////////////////////////////
protected :

  /// These methods should be implemented by every constitutive model.

  /** Compute properties from values in spatial database.
   *
   * @param propValues Array of property values.
   * @param dbValues Array of database values.
   */
  void _dbToProperties(double* const propValues,
		       const double_array& dbValues) const;

  /** Nondimensionalize properties.
   *
   * @param values Array of property values.
   * @param nvalues Number of values.
   */
  void _nondimProperties(double* const values,
			 const int nvalues) const;

  /** Dimensionalize properties.
   *
   * @param values Array of property values.
   * @param nvalues Number of values.
   */
  void _dimProperties(double* const values,
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
  virtual
  void _dbToStateVars(double* const stateValues,
		      const double_array& dbValues) const;

  /** Nondimensionalize state variables.
   *
   * @param values Array of initial state values.
   * @param nvalues Number of values.
   */
  virtual
  void _nondimStateVars(double* const values,
			   const int nvalues) const;
  
  /** Dimensionalize state variables.
   *
   * @param values Array of initial state values.
   * @param nvalues Number of values.
   */
  virtual
  void _dimStateVars(double* const values,
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
  double _calcFriction(const double slip,
		       const double slipRate,
		       const double normalTraction,
		       const double* properties,
		       const int numProperties,
		       const double* stateVars,
		       const int numStateVars);
  virtual
  void _updateStateVars(const double slip,
			const double slipRate,
			double* const stateVars,
			const int numStateVars,
			const double* properties,
			const int numProperties);

  // PRIVATE MEMBERS ////////////////////////////////////////////////////
private :

  /// Indices for properties in section and spatial database.
  static const int p_coefS;
  static const int p_coefD;
  static const int p_d0;

  static const int db_coefS;
  static const int db_coefD;
  static const int db_d0;

  /// Indices for state variables in section and spatial database.
  static const int s_slipCum;
  static const int s_slipPrev;

  static const int db_slipCum;
  static const int db_slipPrev;

  // NOT IMPLEMENTED ////////////////////////////////////////////////////
private :

  SlipWeakening(const SlipWeakening&); ///< Not implemented.
  const SlipWeakening& operator=(const SlipWeakening&); ///< Not implemented

}; // class SlipWeakening

#endif // pylith_friction_slipweakening_hh


// End of file 
