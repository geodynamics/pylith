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

#include "ElasticMaterial3D.hh"

/// Namespace for pylith package
namespace pylith {
  namespace materials {
    class ElasticIsotropic3D;
    class TestElasticIsotropic3D; // unit testing
  } // materials
} // pylith

/// 3-D, isotropic, linear elastic material.
class pylith::materials::ElasticIsotropic3D : public ElasticMaterial3D
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
  ElasticMaterial3D* clone(void) const;

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
  const char** dbValues(void) const;

  /** Get number of values expected to be in database of parameters for
   *  physical properties.
   *
   * @returns Number of values
   */
  int numDBValues(void) const;

  /** Get names of parameters for physical properties.
   *
   * @returns Names of parameters
   */
  const char** parameterNames(void) const;

  /** Get number of parameters for physical properties.
   *
   * @returns Number of parameters
   */
  int numParameters(void) const;

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
  void dbToParams(double* paramVals,
		  const int numParams,
		  double* dbValues,
		  const int numValues) const;

  /** Compute density at location from parameters.
   *
   * @param density Pointer to density at location
   * @param parameters Parameters at location
   * @param numParameters Number of parameters
   */
  void calcDensity(double* density,
		   const double* parameters,
		   const int numParameters);

  /** Compute density at location from parameters.
   *
   * @param elasticConsts Pointer to elastic constants at location
   * @param numConsts Number of elastic constants
   * @param parameters Parameters at location
   * @param numParameters Number of parameters
   */
  void calcElasticConsts(double* elasticConsts,
			 const int numConsts,
			 const double* parameters,
			 const int numParameters);

  // NOT IMPLEMENTED ////////////////////////////////////////////////////
private :

  /// Not implemented
  const ElasticIsotropic3D& operator=(const ElasticIsotropic3D& m);

  // PRIVATE MEMBERS ////////////////////////////////////////////////////
private :

  /// Number of values expected in database of physical property parameters
  static const int _NUMDBVALUES;

  /// Names of values expected in database of physical property parameters
  static const char* _NAMESDBVALUES[];

  /// Number of physical property parameters
  static const int _NUMPARAMETERS;

  /// Names of physical property parameters
  static const char* _NAMESPARAMETERS[];  

}; // class ElasticIsotropic3D

#include "ElasticIsotropic3D.icc" // inline methods

#endif // pylith_materials_elasticisotropic3d_hh


// End of file 
