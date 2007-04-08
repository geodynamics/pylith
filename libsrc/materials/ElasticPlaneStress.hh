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

/** @file libsrc/materials/ElasticPlaneStress.h
 *
 * @brief C++ ElasticPlaneStress object
 *
 * 2-D, isotropic, linear elastic material for plane stress. The
 * physical properties are specified using density, shear-wave speed,
 * and compressional-wave speed. The physical properties are stored
 * internally using density, lambda, and mu, which are directly
 * related to the elasticity constants used in the finite-element
 * integration.
 */

#if !defined(pylith_materials_elasticplanestress_hh)
#define pylith_materials_elasticplanestress_hh

#include "ElasticMaterial.hh"

/// Namespace for pylith package
namespace pylith {
  namespace materials {
    class ElasticPlaneStress;
    class TestElasticPlaneStress; // unit testing
  } // materials
} // pylith

/// 2-D, isotropic, linear elastic material for plane stress.
class pylith::materials::ElasticPlaneStress : public ElasticMaterial
{ // class ElasticPlaneStress
  friend class TestElasticPlaneStress; // unit testing

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Default constructor
  ElasticPlaneStress(void);

  /// Destructor
  ~ElasticPlaneStress(void);

  /** Create a pointer to a copy of this.
   *
   * @returns Pointer to copy
   */
  ElasticMaterial* clone(void) const;

  /** Get number of entries in stress tensor.
   *
   * 1-D = 1
   * 2-D = 3
   * 3-D = 6
   *
   * @returns Number of entries in stress tensor.
   */
  int stressSize(void) const;

  /** Get number of elastic constants for material.
   *
   * 1-D = 1
   * 2-D = 6
   * 3-D = 21
   *
   * @returns Number of elastic constants
   */
  int numElasticConsts(void) const;

  // PROTECTED METHODS //////////////////////////////////////////////////
protected :

  /** Copy constructor
   *
   * @param m Material to copy
   */
  ElasticPlaneStress(const ElasticPlaneStress& m);

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

  /** Get names of parameters for physical properties.
   *
   * @returns Names of parameters
   */
  const char** _parameterNames(void) const;

  /** Get number of parameters for physical properties.
   *
   * @returns Number of parameters
   */
  int _numParameters(void) const;

  /** Compute parameters from values in spatial database.
   *
   * Order of values in arrays matches order used in dbValues() and
   * parameterNames().
   *
   * @param paramVals Array of parameters
   * @param numParams Number of parameters
   * @param dbValues Array of database values
   * @param numValues Number of database values
   */
  void _dbToParameters(double* paramVals,
		       const int numParams,
		       const double* dbValues,
		       const int numValues) const;

  /** Compute density from parameters.
   *
   * @param density Array for density
   * @param size Size of array for density
   * @param parameters Parameters at location
   * @param numParameters Number of parameters
   */
  void _calcDensity(double* const density,
		    const int size,
		    const double* parameters,
		    const int numParameters);

  /** Compute stress tensor from parameters.
   *
   * @param stress Array for stress tensor
   * @param size Size of array for stress tensor
   * @param parameters Parameters at locations.
   * @param numParameters Number of parameters.
   * @param totalStrain Total strain at locations.
   * @param spaceDim Spatial dimension for locations.
   */
  void _calcStress(double* const stress,
		   const int size,
		   const double* parameters,
		   const int numParameters,
		   const double* totalStrain,
		   const int spaceDim);

  /** Compute derivatives of elasticity matrix from parameters.
   *
   * @param elasticConsts Array for elastic constants
   * @param size Size of array
   * @param parameters Parameters at locations.
   * @param numParameters Number of parameters.
   * @param totalStrain Total strain at locations.
   * @param spaceDim Spatial dimension for locations.
   */
  void _calcElasticConsts(double* const elasticConsts,
			  const int size,
			  const double* parameters,
			  const int numParameters,
			  const double* totalStrain,
			  const int spaceDim);

  // NOT IMPLEMENTED ////////////////////////////////////////////////////
private :

  /// Not implemented
  const ElasticPlaneStress& operator=(const ElasticPlaneStress& m);

}; // class ElasticPlaneStress

#include "ElasticPlaneStress.icc" // inline methods

#endif // pylith_materials_elasticplanestress_hh


// End of file 
