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

  /** Get names of values expected to be in database of parameters for
   *  physical properties.
   *
   * @returns Names of values
   */
  const char** _dbValues(void) const;

  /** Get number of values expected to be in database of parameters for
   *  physical properties.
   *
   * @returns Number of values
   */
  int _numDBValues(void) const;

  /** Compute parameters from values in spatial database.
   *
   * Order of values in arrays matches order used in dbValues() and
   * parameterNames().
   *
   * @param paramVals Array of parameters
   * @param dbValues Array of database values
   */
  void _dbToParameters(double* const paramVals,
		       const int numParams,
		       const double_array& dbValues) const;

  /** Get number of entries in stress tensor.
   *
   * 1-D = 1
   * 2-D = 3
   * 3-D = 6
   *
   * @returns Number of entries in stress tensor.
   */
  int _tensorSize(void) const;

  /** Get number of elastic constants for material.
   *
   * 1-D = 1
   * 2-D = 6
   * 3-D = 21
   *
   * @returns Number of elastic constants
   */
  int _numElasticConsts(void) const;

  /** Compute density from parameters.
   *
   * @param density Array for density.
   * @param parameters Parameters at location.
   * @param numParams Number of parameters.
   */
  void _calcDensity(double* const density,
		    const double* parameters,
		    const int numParams);

  /** Compute stress tensor from parameters.
   *
   * @param stress Array for stress tensor.
   * @param stressSize Size of stress tensor.
   * @param parameters Parameters at location.
   * @param numParams Number of parameters.
   * @param totalStrain Total strain at location.
   * @param strainSize Size of strain tensor.
   */
  void _calcStress(double* const stress,
		   const int stressSize,
		   const double* parameters,
		   const int numParams,
		   const double* totalStrain,
		   const int strainSize);

  /** Compute derivatives of elasticity matrix from parameters.
   *
   * @param elasticConsts Array for elastic constants.
   * @param numElasticConsts Number of elastic constants.
   * @param parameters Parameters at location.
   * @param numParams Number of parameters.
   * @param totalStrain Total strain at location.
   * @param strainSize Size of strain tensor.
   */
  void _calcElasticConsts(double* const elasticConsts,
			  const int numElasticConsts,
			  const double* parameters,
			  const int numParams,
			  const double* totalStrain,
			  const int strainSize);

  // NOT IMPLEMENTED ////////////////////////////////////////////////////
private :

  /// Not implemented
  ElasticPlaneStrain(const ElasticPlaneStrain& m);

  /// Not implemented
  const ElasticPlaneStrain& operator=(const ElasticPlaneStrain& m);

}; // class ElasticPlaneStrain

#endif // pylith_materials_elasticplanestrain_hh


// End of file 
