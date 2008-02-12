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
 * compressional-wave speed. The physical properties are stored
 * internally using density, lambda, and mu, which are directly
 * related to the elasticity constants used in the finite-element
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

  /** Compute stress tensor from properties.
   *
   * @param stress Array for stress tensor.
   * @param stressSize Size of stress tensor.
   * @param properties Properties at location.
   * @param numProperties Number of properties.
   * @param totalStrain Total strain at location.
   * @param strainSize Size of strain tensor.
   */
  void _calcStress(double* const stress,
		   const int stressSize,
		   const double* properties,
		   const int numProperties,
		   const double* totalStrain,
		   const int strainSize);

  /** Compute derivatives of elasticity matrix from properties.
   *
   * @param elasticConsts Array for elastic constants.
   * @param numElasticConsts Number of elastic constants.
   * @param properties Properties at location.
   * @param numProperties Number of properties.
   * @param totalStrain Total strain at location.
   * @param strainSize Size of strain tensor.
   */
  void _calcElasticConsts(double* const elasticConsts,
			  const int numElasticConsts,
			  const double* properties,
			  const int numProperties,
			  const double* totalStrain,
			  const int strainSize);

  // NOT IMPLEMENTED ////////////////////////////////////////////////////
private :

  /// Not implemented
  ElasticIsotropic3D(const ElasticIsotropic3D&);

  /// Not implemented
  const ElasticIsotropic3D& operator=(const ElasticIsotropic3D&);

}; // class ElasticIsotropic3D

#endif // pylith_materials_elasticisotropic3d_hh


// End of file 
