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
  void timeStep(const double dt);

  /** Set whether elastic or inelastic constitutive relations are used.
   *
   * @param flag True to use elastic, false to use inelastic.
   */
  void useElasticBehavior(const bool flag);

  /** Get flag indicating whether material implements an empty
   * _updateState() method.
   *
   * @returns False if _updateState() is empty, true otherwise.
   */
  bool usesUpdateState(void) const;

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

  /** Update parameters (for next time step).
   *
   * @param parameters Parameters at location.
   * @param numParams Number of parameters.
   * @param totalStrain Total strain at location.
   * @param strainSize Size of strain tensor.
   */
  void _updateState(double* const parameters,
		    const int numParams,
		    const double* totalStrain,
		    const int strainSize);

  // PRIVATE TYPEDEFS ///////////////////////////////////////////////////
private :

  /// Member prototype for _calcStress()
  typedef void (pylith::materials::MaxwellIsotropic3D::*calcStress_fn_type)
    (double* const,
     const int,
     const double*,
     const int,
     const double*,
     const int);

  /// Member prototype for _calcElasticConsts()
  typedef void (pylith::materials::MaxwellIsotropic3D::*calcElasticConsts_fn_type)
    (double* const,
     const int,
     const double*,
     const int,
     const double*,
     const int);

  /// Member prototype for _updateState()
  typedef void (pylith::materials::MaxwellIsotropic3D::*updateState_fn_type)
    (double* const,
     const int,
     const double*,
     const int);

  // PRIVATE METHODS ////////////////////////////////////////////////////
private :

  /** Compute stress tensor from parameters as an elastic material.
   *
   * @param stress Array for stress tensor.
   * @param stressSize Size of stress tensor.
   * @param parameters Parameters at locations.
   * @param numParams Number of parameters.
   * @param totalStrain Total strain at locations.
   * @param strainSize Size of strain tensor.
   */
  void _calcStressElastic(double* const stress,
			  const int stressSize,
			  const double* parameters,
			  const int numParams,
			  const double* totalStrain,
			  const int strainSize);

  /** Compute stress tensor from parameters as an viscoelastic material.
   *
   * @param stress Array for stress tensor.
   * @param stressSize Size of stress tensor.
   * @param parameters Parameters at locations.
   * @param numParams Number of parameters.
   * @param totalStrain Total strain at locations.
   * @param strainSize Size of strain tensor.
   */
  void _calcStressViscoelastic(double* const stress,
			  const int stressSize,
			  const double* parameters,
			  const int numParams,
			  const double* totalStrain,
			  const int strainSize);

  /** Compute derivatives of elasticity matrix from parameters as an
   * elastic material.
   *
   * @param elasticConsts Array for elastic constants.
   * @param numElasticConsts Number of elastic constants.
   * @param parameters Parameters at location.
   * @param numParams Number of parameters.
   * @param totalStrain Total strain at location.
   * @param strainSize Size of strain tensor.
   */
  void _calcElasticConstsElastic(double* const elasticConsts,
				 const int numElasticConsts,
				 const double* parameters,
				 const int numParams,
				 const double* totalStrain,
				 const int strainSize);

  /** Compute derivatives of elasticity matrix from parameters as a
   * viscoelastic material.
   *
   * @param elasticConsts Array for elastic constants.
   * @param numElasticConsts Number of elastic constants.
   * @param parameters Parameters at location.
   * @param numParams Number of parameters.
   * @param totalStrain Total strain at location.
   * @param strainSize Size of strain tensor.
   */
  void _calcElasticConstsViscoelastic(double* const elasticConsts,
				      const int numElasticConsts,
				      const double* parameters,
				      const int numParams,
				      const double* totalStrain,
				      const int strainSize);

  /** Update state variables after solve as an elastic material.
   *
   * @param parameters Parameters at location.
   * @param numParams Number of parameters.
   * @param totalStrain Total strain at location.
   * @param strainSize Size of strain tensor.
   */
  void _updateStateElastic(double* const parameters,
			   const int numParams,
			   const double* totalStrain,
			   const int strainSize);

  /** Update state variables after solve as a viscoelastic material.
   *
   * @param parameters Parameters at location.
   * @param numParams Number of parameters.
   * @param totalStrain Total strain at location.
   * @param strainSize Size of strain tensor.
   */
  void _updateStateViscoelastic(double* const parameters,
				const int numParams,
				const double* totalStrain,
				const int strainSize);

  // NOT IMPLEMENTED ////////////////////////////////////////////////////
private :

  /// Not implemented
  MaxwellIsotropic3D(const MaxwellIsotropic3D& m);

  /// Not implemented
  const MaxwellIsotropic3D& operator=(const MaxwellIsotropic3D& m);

  // PRIVATE MEMBERS ////////////////////////////////////////////////////
private :

  /// Method to use for _calcElasticConsts().
  calcElasticConsts_fn_type _calcElasticConstsFn;

  /// Method to use for _calcStress().
  calcStress_fn_type _calcStressFn;

  /// Method to use for _updateState().
  updateState_fn_type _updateStateFn;

}; // class MaxwellIsotropic3D

#include "MaxwellIsotropic3D.icc" // inline methods

#endif // pylith_materials_maxwellisotropic3d_hh


// End of file 
