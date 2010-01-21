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
  double _calcFriction(const double slip,
		       const double slipRate,
		       const double normalTraction,
		       const double* properties,
		       const int numProperties,
		       const double* stateVars,
		       const int numStateVars);

  // PRIVATE MEMBERS ////////////////////////////////////////////////////
private :

  static const int p_coef;
  static const int db_coef;

  // NOT IMPLEMENTED ////////////////////////////////////////////////////
private :

  StaticFriction(const StaticFriction&); ///< Not implemented.
  const StaticFriction& operator=(const StaticFriction&); ///< Not implemented

}; // class StaticFriction

#endif // pylith_friction_staticfriction_hh


// End of file 
