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

/** @file libsrc/materials/MaxwellIsotropic3D.h
 *
 * @brief C++ MaxwellIsotropic3D object
 *
 * 3-D, isotropic, linear Maxwell viscoelastic material. The
 * physical properties are specified using density, shear-wave speed,
 * viscosity, and compressional-wave speed. The physical properties are
 * stored internally using density, lambda, mu, which are directly
 * related to the elasticity constants used in the finite-element
 * integration. The viscosity is stored using Maxwell Time (viscosity/mu).
 */

#if !defined(pylith_materials_maxwellisotropic3d_hh)
#define pylith_materials_maxwellisotropic3d_hh

#include "ElasticMaterial.hh"

/// Namespace for pylith package
namespace pylith {
  namespace materials {
    class MaxwellIsotropic3D;
    class TestMaxwellIsotropic3D; // unit testing
  } // materials
} // pylith

/// 3-D, isotropic, linear Maxwell viscoelastic material.
class pylith::materials::MaxwellIsotropic3D : public ElasticMaterial
{ // class MaxwellIsotropic3D
  friend class TestMaxwellIsotropic3D; // unit testing

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Default constructor
  MaxwellIsotropic3D(void);

  /// Destructor
  ~MaxwellIsotropic3D(void);

  /** Set current time step.
   *
   * @param dt Current time step.
   */
  void timestep(const double dt);

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

  /** Get names of parameters for physical properties.
   *
   * @returns Names of parameters
   */
  const char** _parameterNames(void) const;

  /** Get number of parameters for physical properties.
   *
   * @returns Number of parameters
   */
  void _numParamValues(int_array* numValues) const;

  /** Compute parameters from values in spatial database.
   *
   * Order of values in arrays matches order used in dbValues() and
   * parameterNames().
   *
   * @param paramVals Array of parameters
   * @param dbValues Array of database values
   */
  void _dbToParameters(std::vector<double_array>* paramVals,
		       const double_array& dbValues) const;

  /** Get number of entries in stress/strain tensors.
   *
   * 1-D = 1
   * 2-D = 3
   * 3-D = 6
   *
   * @returns Number of entries in stress/strain tensors.
   */
  int _tensorSize(void) const;

  /** Get number of entries in derivative of elasticity matrix.
   *
   * 1-D = 1
   * 2-D = 6
   * 3-D = 21
   *
   * @returns Number of entries in derivative of elasticity matrix.
   */
  int _numElasticConsts(void) const;

  /** Compute density from parameters.
   *
   * @param density Array for density
   * @param parameters Parameters at location
   */
  void _calcDensity(double_array* const density,
		    const std::vector<double_array>& parameters);

  /** Compute stress tensor from parameters.
   *
   * @param stress Array for stress tensor
   * @param parameters Parameters at locations.
   * @param totalStrain Total strain at locations.
   */
  void _calcStress(double_array* const stress,
		   const std::vector<double_array>& parameters,
		   const double_array& totalStrain);

  /** Compute derivatives of elasticity matrix from parameters.
   *
   * @param elasticConsts Array for elastic constants
   * @param parameters Parameters at locations.
   * @param totalStrain Total strain at locations.
   */
  void _calcElasticConsts(double_array* const elasticConsts,
			  const std::vector<double_array>& parameters,
			  const double_array& totalStrain);

  /** Update state variables after solve.
   *
   * @param parameters Parameters at locations.
   * @param totalStrain Total strain at locations.
   */
  void _updateState(std::vector<double_array>& parameters,
		    const double_array& totalStrain);

  // NOT IMPLEMENTED ////////////////////////////////////////////////////
private :

  /// Not implemented
  MaxwellIsotropic3D(const MaxwellIsotropic3D& m);

  /// Not implemented
  const MaxwellIsotropic3D& operator=(const MaxwellIsotropic3D& m);


}; // class MaxwellIsotropic3D

#include "MaxwellIsotropic3D.icc" // inline methods

#endif // pylith_materials_maxwellisotropic3d_hh


// End of file 
