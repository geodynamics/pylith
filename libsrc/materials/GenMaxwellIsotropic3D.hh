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

/** @file libsrc/materials/GenMaxwellIsotropic3D.h
 *
 * @brief C++ GenMaxwellIsotropic3D object
 *
 * 3-D, isotropic, generalized linear Maxwell viscoelastic material.
 * This consists of several Maxwell models in parallel. At present,
 * the number of models is fixed at 3, but this will be changed in the
 * future. The physical properties are specified using density,
 * shear-wave speed, and compressional-wave speed. A viscosity and a
 * shear ratio are also given for each Maxwell model. The shear ratio
 * specifies how much of the total shear modulus is associated with
 * that model. The shear ratios must sum to a value less than one. If
 * the value is less than one, the remainder of the total shear
 * modulus is associated with a spring in parallel with the Maxwell
 * models.  The physical properties are stored internally using
 * density, lambdaTot, muTot, which are directly related to the
 * elasticity constants used in the finite-element integration. The
 * viscosity for each model is stored using Maxwell Time
 * (viscosity/mu), and the shear ratio is also stored for each Maxwell
 * model.
 */

// :TODO: Rewrite as template over the number of Maxwell models?
// We could instatiate for 2 and 3 models and provide example for how to
// instantiate over other numbers of Maxwell models.

#if !defined(pylith_materials_genmaxwellisotropic3d_hh)
#define pylith_materials_genmaxwellisotropic3d_hh

#include "ElasticMaterial.hh"

/// Namespace for pylith package
namespace pylith {
  namespace materials {
    class GenMaxwellIsotropic3D;
    class TestGenMaxwellIsotropic3D; // unit testing
  } // materials
} // pylith

/// 3-D, isotropic, generalized linear Maxwell viscoelastic material.
class pylith::materials::GenMaxwellIsotropic3D : public ElasticMaterial
{ // class GenMaxwellIsotropic3D
  friend class TestGenMaxwellIsotropic3D; // unit testing

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Default constructor
  GenMaxwellIsotropic3D(void);

  /// Destructor
  ~GenMaxwellIsotropic3D(void);

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
   * @param computeStateVars Flag indicating to compute updated state vars.
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
  typedef void (pylith::materials::GenMaxwellIsotropic3D::*calcStress_fn_type)
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
  typedef void (pylith::materials::GenMaxwellIsotropic3D::*calcElasticConsts_fn_type)
    (double* const,
     const int,
     const double*,
     const int,
     const double*,
     const int,
     const double*,
     const int);

  /// Member prototype for _updateProperties()
  typedef void (pylith::materials::GenMaxwellIsotropic3D::*updateProperties_fn_type)
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
  GenMaxwellIsotropic3D(const GenMaxwellIsotropic3D& m);

  /// Not implemented
  const GenMaxwellIsotropic3D& operator=(const GenMaxwellIsotropic3D& m);

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

}; // class GenMaxwellIsotropic3D

#include "GenMaxwellIsotropic3D.icc" // inline methods

#endif // pylith_materials_genmaxwellisotropic3d_hh


// End of file 
