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

/** @file libsrc/materials/ElasticIsotropic3D.h
 *
 * @brief C++ ElasticIsotropic3D object
 *
 * 3-D, isotropic, linear elastic material. The physical properties
 * are specified using density, shear-wave speed, and
 * compressional-wave speed. The physical properties stored are
 * internally using density, lambda, and mu, which are directly
 * related to the elasticity constants use in the finite-element
 * integration.
 */

#if !defined(pylith_materials_elasticisotropic3d_hh)
#define pylith_materials_elasticisotropic3d_hh

#include "ElasticMaterial.hh"

/// Namespace for pylith package
namespace pylith {
  namespace materials {
    class ElasticIsotropic3D;
    class TestElasticIsotropic3D; // unit testing
  } // materials
} // pylith

/// 3-D, isotropic, linear elastic material.
class pylith::materials::ElasticIsotropic3D : public ElasticMaterial
{ // class ElasticIsotropic3D
  friend class TestElasticIsotropic3D; // unit testing

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Default constructor
  ElasticIsotropic3D(void);

  /// Destructor
  ~ElasticIsotropic3D(void);

  /** Create a pointer to a copy of this.
   *
   * @returns Pointer to copy
   */
  ElasticMaterial* clone(void) const;

  /** Get number of elastic constants for material.
   *
   * 1-D = 1
   * 2-D = 6
   * 3-D = 21
   *
   * @returns Number of elastic constants
   */
  const int numElasticConsts(void) const;

  // PROTECTED METHODS //////////////////////////////////////////////////
protected :

  /** Copy constructor
   *
   * @param m Material to copy
   */
  ElasticIsotropic3D(const ElasticIsotropic3D& m);

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
  void _dbToParams(double* paramVals,
		   const int numParams,
		   double* dbValues,
		   const int numValues) const;

  /** Compute density at locations from parameters.
   *
   * Results are stored in _density.
   *
   * Index into parameters = iLoc*numParameters + iParam
   *
   * @param parameters Parameters at location
   * @param numParameters Number of parameters
   * @param numLocs Number of locations
   */
  void _calcDensity(const double* parameters,
		    const int numParameters,
		    const int numLocs);

  /** Compute density at locations from parameters.
   *
   * Results are stored in _elasticConsts.
   *
   * Index into parameters = iLoc*numParameters + iParam
   *
   * @param parameters Parameters at location
   * @param numParameters Number of parameters
   * @param numLocs Number of locations
   */
  void _calcElasticConsts(const double* parameters,
			  const int numParameters,
			  const int numLocs);

  // NOT IMPLEMENTED ////////////////////////////////////////////////////
private :

  /// Not implemented
  const ElasticIsotropic3D& operator=(const ElasticIsotropic3D& m);

}; // class ElasticIsotropic3D

#include "ElasticIsotropic3D.icc" // inline methods

#endif // pylith_materials_elasticisotropic3d_hh


// End of file 
