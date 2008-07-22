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

/** @file libsrc/materials/ElasticPlaneStrain.h
 *
 * @brief C++ ElasticPlaneStrain object
 *
 * 2-D, isotropic, linear elastic material for plane strain. The
 * physical properties are specified using density, shear-wave speed,
 * and compressional-wave speed. The physical properties are stored
 * internally using density, lambda, and mu, which are directly
 * related to the elasticity constants used in the finite-element
 * integration.
 */

#if !defined(pylith_materials_elasticplanestrain_hh)
#define pylith_materials_elasticplanestrain_hh

#include "ElasticMaterial.hh"

/// Namespace for pylith package
namespace pylith {
  namespace materials {
    class ElasticPlaneStrain;
    class TestElasticPlaneStrain; // unit testing
  } // materials
} // pylith

/// 3-D, isotropic, linear elastic material.
class pylith::materials::ElasticPlaneStrain : public ElasticMaterial
{ // class ElasticPlaneStrain
  friend class TestElasticPlaneStrain; // unit testing

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Default constructor
  ElasticPlaneStrain(void);

  /// Destructor
  ~ElasticPlaneStrain(void);

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

  // NOT IMPLEMENTED ////////////////////////////////////////////////////
private :

  /// Not implemented
  ElasticPlaneStrain(const ElasticPlaneStrain& m);

  /// Not implemented
  const ElasticPlaneStrain& operator=(const ElasticPlaneStrain& m);

}; // class ElasticPlaneStrain

#endif // pylith_materials_elasticplanestrain_hh


// End of file 
