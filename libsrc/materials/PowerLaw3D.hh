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

/** @file libsrc/materials/PowerLaw3D.h
 *
 * @brief C++ PowerLaw3D object
 *
 * 3-D, isotropic, power-law Maxwell viscoelastic material. The
 * physical properties are specified using density, shear-wave speed,
 * viscosity coefficient, power-law exponent, and compressional-wave speed.
 * The physical properties are stored internally using density, lambda, mu,
 * which are directly related to the elasticity constants used in the
 * finite-element integration. The viscosity information is retained as
 * specified.
 */

#if !defined(pylith_materials_powerlaw3d_hh)
#define pylith_materials_powerlaw3d_hh

#include "ElasticMaterial.hh"

/// Namespace for pylith package
namespace pylith {
  namespace materials {
    class PowerLaw3D;
    class TestPowerLaw3D; // unit testing
  } // materials
} // pylith

/// 3-D, isotropic, linear Maxwell viscoelastic material.
class pylith::materials::PowerLaw3D : public ElasticMaterial
{ // class PowerLaw3D
  friend class TestPowerLaw3D; // unit testing

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Default constructor
  PowerLaw3D(void);

  /// Destructor
  ~PowerLaw3D(void);

  /** Set current time step.
   *
   * @param dt Current time step.
   */
  void timeStep(const double dt);

  /** Set whether elastic or inelastic constitutive relations are used.
   *
   * @param flag True to use elastic, false to use inelastic.
   */
  void useElasticBehavior(const bool flag);

  /** Get flag indicating whether material implements an empty
   * _updateProperties() method.
   *
   * @returns False if _updateProperties() is empty, true otherwise.
   */
  bool usesUpdateProperties(void) const;

  // PROTECTED METHODS //////////////////////////////////////////////////
protected :

  /** Compute properties from values in spatial database.
   *
   * Order of values in arrays matches order used in dbValues() and
   * parameterNames().
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

  /** Nondimensionalize initial state.
   *
   * @param values Array of initial state values.
   * @param nvalues Number of values.
   */
  void _nondimInitState(double* const values,
			const int nvalues) const;

  /** Dimensionalize initial state.
   *
   * @param values Array of initial state values.
   * @param nvalues Number of values.
   */
  void _dimInitState(double* const values,
		     const int nvalues) const;

  /** Compute density from properties.
   *
   * @param density Array for density.
   * @param properties Properties at location.
   * @param numProperties Number of properties.
   */
  void _calcDensity(double* const density,
		    const double* properties,
		    const int numProperties);

  /** Compute stress tensor from properties. If the state variables
   * are from the previous time step, then the computeStateVars flag
   * should be set to true so that the state variables are updated
   * (but not stored) when computing the stresses.
   *
   * @param stress Array for stress tensor.
   * @param stressSize Size of stress tensor.
   * @param properties Properties at location.
   * @param numProperties Number of properties.
   * @param totalStrain Total strain at location.
   * @param strainSize Size of strain tensor.
   * @param initialState Initial state values.
   * @param initialStateSize Size of initial state array.
   */
  void _calcStress(double* const stress,
		   const int stressSize,
		   const double* properties,
		   const int numProperties,
		   const double* totalStrain,
		   const int strainSize,
		   const double* initialState,
		   const int initialStateSize,
		   const bool computeStateVars);

  /** Compute derivatives of elasticity matrix from properties.
   *
   * @param elasticConsts Array for elastic constants.
   * @param numElasticConsts Number of elastic constants.
   * @param properties Properties at location.
   * @param numProperties Number of properties.
   * @param totalStrain Total strain at location.
   * @param strainSize Size of strain tensor.
   * @param initialState Initial state values.
   * @param initialStateSize Size of initial state array.
   */
  void _calcElasticConsts(double* const elasticConsts,
			  const int numElasticConsts,
			  const double* properties,
			  const int numProperties,
			  const double* totalStrain,
			  const int strainSize,
		          const double* initialState,
		          const int initialStateSize);

  /** Get stable time step for implicit time integration.
   *
   * @returns Time step
   */
  double _stableTimeStepImplicit(const double* properties,
				 const int numProperties) const;

  /** Update properties (for next time step).
   *
   * @param properties Properties at location.
   * @param numProperties Number of properties.
   * @param totalStrain Total strain at location.
   * @param strainSize Size of strain tensor.
   * @param initialState Initial state values.
   * @param initialStateSize Size of initial state array.
   */
  void _updateProperties(double* const properties,
		    const int numProperties,
		    const double* totalStrain,
		    const int strainSize,
		    const double* initialState,
		    const int initialStateSize);

  // PRIVATE TYPEDEFS ///////////////////////////////////////////////////
private :

  /// Member prototype for _calcStress()
  typedef void (pylith::materials::PowerLaw3D::*calcStress_fn_type)
    (double* const,
     const int,
     const double*,
     const int,
     const double*,
     const int,
     const double*,
     const int,
     const bool);

  /// Member prototype for _calcElasticConsts()
  typedef void (pylith::materials::PowerLaw3D::*calcElasticConsts_fn_type)
    (double* const,
     const int,
     const double*,
     const int,
     const double*,
     const int,
     const double*,
     const int);

  /// Member prototype for _updateProperties()
  typedef void (pylith::materials::PowerLaw3D::*updateProperties_fn_type)
    (double* const,
     const int,
     const double*,
     const int,
     const double*,
     const int);

  // PRIVATE METHODS ////////////////////////////////////////////////////
private :

/** Compute viscous strains (state variables) for the current time step.
   *
   * @param properties Properties at location.
   * @param numProperties Number of properties.
   * @param totalStrain Total strain at location.
   * @param strainSize Size of strain tensor.
   * @param initialState Initial state values.
   * @param initialStateSize Size of initial state array.
   */
  void _computeStateVars(const double* properties,
			 const int numProperties,
			 const double* totalStrain,
			 const int strainSize,
			 const double* initialState,
			 const int initialStateSize);

  /** Compute stress tensor from properties as an elastic material.
   *
   * @param stress Array for stress tensor.
   * @param stressSize Size of stress tensor.
   * @param properties Properties at locations.
   * @param numProperties Number of properties.
   * @param totalStrain Total strain at locations.
   * @param strainSize Size of strain tensor.
   * @param initialState Initial state values.
   * @param initialStateSize Size of initial state array.
   * @param computeStateVars Flag indicating to compute updated state vars.
   */
  void _calcStressElastic(double* const stress,
			  const int stressSize,
			  const double* properties,
			  const int numProperties,
			  const double* totalStrain,
			  const int strainSize,
			  const double* initialState,
			  const int initialStateSize,
			  const bool computeStateVars);

  /** Compute stress tensor from properties as an viscoelastic material.
   *
   * @param stress Array for stress tensor.
   * @param stressSize Size of stress tensor.
   * @param properties Properties at locations.
   * @param numProperties Number of properties.
   * @param totalStrain Total strain at locations.
   * @param strainSize Size of strain tensor.
   * @param initialState Initial state values.
   * @param initialStateSize Size of initial state array.
   * @param computeStateVars Flag indicating to compute updated state vars.
   */
  void _calcStressViscoelastic(double* const stress,
			       const int stressSize,
			       const double* properties,
			       const int numProperties,
			       const double* totalStrain,
			       const int strainSize,
			       const double* initialState,
			       const int initialStateSize,
			       const bool computeStateVars);

  /** Compute derivatives of elasticity matrix from properties as an
   * elastic material.
   *
   * @param elasticConsts Array for elastic constants.
   * @param numElasticConsts Number of elastic constants.
   * @param properties Properties at location.
   * @param numProperties Number of properties.
   * @param totalStrain Total strain at location.
   * @param strainSize Size of strain tensor.
   * @param initialState Initial state values.
   * @param initialStateSize Size of initial state array.
   */
  void _calcElasticConstsElastic(double* const elasticConsts,
				 const int numElasticConsts,
				 const double* properties,
				 const int numProperties,
				 const double* totalStrain,
				 const int strainSize,
			         const double* initialState,
			         const int initialStateSize);

  /** Compute derivatives of elasticity matrix from properties as a
   * viscoelastic material.
   *
   * @param elasticConsts Array for elastic constants.
   * @param numElasticConsts Number of elastic constants.
   * @param properties Properties at location.
   * @param numProperties Number of properties.
   * @param totalStrain Total strain at location.
   * @param strainSize Size of strain tensor.
   * @param initialState Initial state values.
   * @param initialStateSize Size of initial state array.
   */
  void _calcElasticConstsViscoelastic(double* const elasticConsts,
				      const int numElasticConsts,
				      const double* properties,
				      const int numProperties,
				      const double* totalStrain,
				      const int strainSize,
			              const double* initialState,
			              const int initialStateSize);

  /** Update state variables after solve as an elastic material.
   *
   * @param properties Properties at location.
   * @param numProperties Number of properties.
   * @param totalStrain Total strain at location.
   * @param strainSize Size of strain tensor.
   * @param initialState Initial state values.
   * @param initialStateSize Size of initial state array.
   */
  void _updatePropertiesElastic(double* const properties,
			   const int numProperties,
			   const double* totalStrain,
			   const int strainSize,
			   const double* initialState,
			   const int initialStateSize);

  /** Update state variables after solve as a viscoelastic material.
   *
   * @param properties Properties at location.
   * @param numProperties Number of properties.
   * @param totalStrain Total strain at location.
   * @param strainSize Size of strain tensor.
   * @param initialState Initial state values.
   * @param initialStateSize Size of initial state array.
   */
  void _updatePropertiesViscoelastic(double* const properties,
				const int numProperties,
				const double* totalStrain,
				const int strainSize,
			        const double* initialState,
			        const int initialStateSize);

  // NOT IMPLEMENTED ////////////////////////////////////////////////////
private :

  /// Not implemented
  PowerLaw3D(const PowerLaw3D& m);

  /// Not implemented
  const PowerLaw3D& operator=(const PowerLaw3D& m);

  // PRIVATE MEMBERS ////////////////////////////////////////////////////
private :

  /// Viscous strain array.
  double_array _visStrain;

  /// Method to use for _calcElasticConsts().
  calcElasticConsts_fn_type _calcElasticConstsFn;

  /// Method to use for _calcStress().
  calcStress_fn_type _calcStressFn;

  /// Method to use for _updateProperties().
  updateProperties_fn_type _updatePropertiesFn;

}; // class PowerLaw3D

#include "PowerLaw3D.icc" // inline methods

#endif // pylith_materials_powerlaw3d_hh


// End of file 
