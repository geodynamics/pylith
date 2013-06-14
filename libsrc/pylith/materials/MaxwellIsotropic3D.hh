// -*- C++ -*-
//
// ----------------------------------------------------------------------
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University of Chicago
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2013 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

/** @file libsrc/materials/MaxwellIsotropic3D.hh
 *
 * @brief 3-D, isotropic, linear Maxwell viscoelastic material.
 */

#if !defined(pylith_materials_maxwellisotropic3d_hh)
#define pylith_materials_maxwellisotropic3d_hh

// Include directives ---------------------------------------------------
#include "ElasticMaterial.hh" // ISA ElasticMaterial

// MaxwellIsotropic3D ---------------------------------------------------
/** @brief 3-D, isotropic, linear Maxwell viscoelastic material.
 *
 * The physical properties are specified using density, shear-wave
 * speed, viscosity, and compressional-wave speed. The physical
 * properties are stored internally using density, lambda, mu, which
 * are directly related to the elasticity constants used in the
 * finite-element integration. The viscosity is stored using Maxwell
 * Time (viscosity/mu).
 */
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
  void timeStep(const PylithScalar dt);

  /** Set whether elastic or inelastic constitutive relations are used.
   *
   * @param flag True to use elastic, false to use inelastic.
   */
  void useElasticBehavior(const bool flag);

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
  void _dbToProperties(PylithScalar* const propValues,
		       const scalar_array& dbValues);

  /** Nondimensionalize properties.
   *
   * @param values Array of property values.
   * @param nvalues Number of values.
   */
  void _nondimProperties(PylithScalar* const values,
			 const int nvalues) const;

  /** Dimensionalize properties.
   *
   * @param values Array of property values.
   * @param nvalues Number of values.
   */
  void _dimProperties(PylithScalar* const values,
		      const int nvalues) const;

  /** Compute initial state variables from values in spatial database.
   *
   * @param stateValues Array of state variable values.
   * @param dbValues Array of database values.
   */
  void _dbToStateVars(PylithScalar* const stateValues,
		      const scalar_array& dbValues);

  // Note: We do not need to dimensionalize or nondimensionalize state
  // variables because there are strains, which are dimensionless.


  /** Compute density from properties.
   *
   * @param density Array for density.
   * @param properties Properties at location.
   * @param numProperties Number of properties.
   * @param stateVars Number of state variables.
   * @param numStateVars Number of state variables.
   */
  void _calcDensity(PylithScalar* const density,
		    const PylithScalar* properties,
		    const int numProperties,
		    const PylithScalar* stateVars,
		    const int numStateVars);

  /** Compute stress tensor from properties and state variables. If
   * the state variables are from the previous time step, then the
   * computeStateVars flag should be set to true so that the state
   * variables are updated (but not stored) when computing the
   * stresses.
   *
   * @param stress Array for stress tensor.
   * @param stressSize Size of stress tensor.
   * @param properties Properties at location.
   * @param numProperties Number of properties.
   * @param stateVars State variables at location.
   * @param numStateVars Number of state variables.
   * @param totalStrain Total strain at location.
   * @param strainSize Size of strain tensor.
   * @param initialStress Initial stress tensor at location.
   * @param initialStressSize Size of initial stress array.
   * @param initialStrain Initial strain tensor at location.
   * @param initialStrainSize Size of initial strain array.
   * @param computeStateVars Flag indicating to compute updated state variables.
   */
  void _calcStress(PylithScalar* const stress,
		   const int stressSize,
		   const PylithScalar* properties,
		   const int numProperties,
		   const PylithScalar* stateVars,
		   const int numStateVars,
		   const PylithScalar* totalStrain,
		   const int strainSize,
		   const PylithScalar* initialStress,
		   const int initialStressSize,
		   const PylithScalar* initialStrain,
		   const int initialStrainSize,
		   const bool computeStateVars);

  /** Compute derivatives of elasticity matrix from properties.
   *
   * @param elasticConsts Array for elastic constants.
   * @param numElasticConsts Number of elastic constants.
   * @param properties Properties at location.
   * @param numProperties Number of properties.
   * @param stateVars State variables at location.
   * @param numStateVars Number of state variables.
   * @param totalStrain Total strain at location.
   * @param strainSize Size of strain tensor.
   * @param initialStress Initial stress tensor at location.
   * @param initialStressSize Size of initial stress array.
   * @param initialStrain Initial strain tensor at location.
   * @param initialStrainSize Size of initial strain array.
   */
  void _calcElasticConsts(PylithScalar* const elasticConsts,
			  const int numElasticConsts,
			  const PylithScalar* properties,
			  const int numProperties,
			  const PylithScalar* stateVars,
			  const int numStateVars,
			  const PylithScalar* totalStrain,
			  const int strainSize,
			  const PylithScalar* initialStress,
			  const int initialStressSize,
			  const PylithScalar* initialStrain,
			  const int initialStrainSize);

  /** Update state variables (for next time step).
   *
   * @param stateVars State variables at location.
   * @param numStateVars Number of state variables.
   * @param properties Properties at location.
   * @param numProperties Number of properties.
   * @param totalStrain Total strain at location.
   * @param strainSize Size of strain tensor.
   * @param initialStress Initial stress tensor at location.
   * @param initialStressSize Size of initial stress array.
   * @param initialStrain Initial strain tensor at location.
   * @param initialStrainSize Size of initial strain array.
   */
  void _updateStateVars(PylithScalar* const stateVars,
			const int numStateVars,
			const PylithScalar* properties,
			const int numProperties,
			const PylithScalar* totalStrain,
			const int strainSize,
			const PylithScalar* initialStress,
			const int initialStressSize,
			const PylithScalar* initialStrain,
			const int initialStrainSize);

  /** Get stable time step for implicit time integration.
   *
   * @param properties Properties at location.
   * @param numProperties Number of properties.
   * @param stateVars State variables at location.
   * @param numStateVars Number of state variables.
   *
   * @returns Time step
   */
  PylithScalar _stableTimeStepImplicit(const PylithScalar* properties,
				 const int numProperties,
				 const PylithScalar* stateVars,
				 const int numStateVars) const;

  /** Get stable time step for explicit time integration.
   *
   * @param properties Properties at location.
   * @param numProperties Number of properties.
   * @param stateVars State variables at location.
   * @param numStateVars Number of state variables.
   * @param minCellWidth Minimum width across cell.
   *
   * @returns Time step
   */
  PylithScalar _stableTimeStepExplicit(const PylithScalar* properties,
				       const int numProperties,
				       const PylithScalar* stateVars,
				       const int numStateVars,
				       const double minCellWidth) const;
  
  // PRIVATE TYPEDEFS ///////////////////////////////////////////////////
private :

  /// Member prototype for _calcStress()
  typedef void (pylith::materials::MaxwellIsotropic3D::*calcStress_fn_type)
    (PylithScalar* const,
     const int,
     const PylithScalar*,
     const int,
     const PylithScalar*,
     const int,
     const PylithScalar*,
     const int,
     const PylithScalar*,
     const int,
     const PylithScalar*,
     const int,
     const bool);

  /// Member prototype for _calcElasticConsts()
  typedef void (pylith::materials::MaxwellIsotropic3D::*calcElasticConsts_fn_type)
    (PylithScalar* const,
     const int,
     const PylithScalar*,
     const int,
     const PylithScalar*,
     const int,
     const PylithScalar*,
     const int,
     const PylithScalar*,
     const int,
     const PylithScalar*,
     const int);

  /// Member prototype for _updateStateVars()
  typedef void (pylith::materials::MaxwellIsotropic3D::*updateStateVars_fn_type)
    (PylithScalar* const,
     const int,
     const PylithScalar*,
     const int,
     const PylithScalar*,
     const int,
     const PylithScalar*,
     const int,
     const PylithScalar*,
     const int);

  // PRIVATE METHODS ////////////////////////////////////////////////////
private :

  /** Compute stress tensor from properties as an elastic material.
   *
   * @param stress Array for stress tensor.
   * @param stressSize Size of stress tensor.
   * @param properties Properties at location.
   * @param numProperties Number of properties.
   * @param stateVars State variables at location.
   * @param numStateVars Number of state variables.
   * @param totalStrain Total strain at location.
   * @param strainSize Size of strain tensor.
   * @param initialStress Initial stress tensor at location.
   * @param initialStressSize Size of initial stress array.
   * @param initialStrain Initial strain tensor at location.
   * @param initialStrainSize Size of initial strain array.
   * @param computeStateVars Flag indicating to compute updated state variables.
   */
  void _calcStressElastic(PylithScalar* const stress,
			  const int stressSize,
			  const PylithScalar* properties,
			  const int numProperties,
			  const PylithScalar* stateVars,
			  const int numStateVars,
			  const PylithScalar* totalStrain,
			  const int strainSize,
			  const PylithScalar* initialStress,
			  const int initialStressSize,
			  const PylithScalar* initialStrain,
			  const int initialStrainSize,
			  const bool computeStateVars);

  /** Compute stress tensor from properties as a viscoelastic material.
   *
   * @param stress Array for stress tensor.
   * @param stressSize Size of stress tensor.
   * @param properties Properties at location.
   * @param numProperties Number of properties.
   * @param stateVars State variables at location.
   * @param numStateVars Number of state variables.
   * @param totalStrain Total strain at location.
   * @param strainSize Size of strain tensor.
   * @param initialStress Initial stress tensor at location.
   * @param initialStressSize Size of initial stress array.
   * @param initialStrain Initial strain tensor at location.
   * @param initialStrainSize Size of initial strain array.
   * @param computeStateVars Flag indicating to compute updated state variables.
   */
  void _calcStressViscoelastic(PylithScalar* const stress,
			       const int stressSize,
			       const PylithScalar* properties,
			       const int numProperties,
			       const PylithScalar* stateVars,
			       const int numStateVars,
			       const PylithScalar* totalStrain,
			       const int strainSize,
			       const PylithScalar* initialStress,
			       const int initialStressSize,
			       const PylithScalar* initialStrain,
			       const int initialStrainSize,
			       const bool computeStateVars);

  /** Compute derivatives of elasticity matrix from properties as an
   * elastic material.
   *
   * @param elasticConsts Array for elastic constants.
   * @param numElasticConsts Number of elastic constants.
   * @param properties Properties at location.
   * @param numProperties Number of properties.
   * @param stateVars State variables at location.
   * @param numStateVars Number of state variables.
   * @param totalStrain Total strain at location.
   * @param strainSize Size of strain tensor.
   * @param initialStress Initial stress tensor at location.
   * @param initialStressSize Size of initial stress array.
   * @param initialStrain Initial strain tensor at location.
   * @param initialStrainSize Size of initial strain array.
   */
  void _calcElasticConstsElastic(PylithScalar* const elasticConsts,
				 const int numElasticConsts,
				 const PylithScalar* properties,
				 const int numProperties,
				 const PylithScalar* stateVars,
				 const int numStateVars,
				 const PylithScalar* totalStrain,
				 const int strainSize,
				 const PylithScalar* initialStress,
				 const int initialStressSize,
				 const PylithScalar* initialStrain,
				 const int initialStrainSize);

  /** Compute derivatives of elasticity matrix from properties as a
   * viscoelastic material.
   *
   * @param elasticConsts Array for elastic constants.
   * @param numElasticConsts Number of elastic constants.
   * @param properties Properties at location.
   * @param numProperties Number of properties.
   * @param stateVars State variables at location.
   * @param numStateVars Number of state variables.
   * @param totalStrain Total strain at location.
   * @param strainSize Size of strain tensor.
   * @param initialStress Initial stress tensor at location.
   * @param initialStressSize Size of initial stress array.
   * @param initialStrain Initial strain tensor at location.
   * @param initialStrainSize Size of initial strain array.
   */
  void _calcElasticConstsViscoelastic(PylithScalar* const elasticConsts,
				      const int numElasticConsts,
				      const PylithScalar* properties,
				      const int numProperties,
				      const PylithScalar* stateVars,
				      const int numStateVars,
				      const PylithScalar* totalStrain,
				      const int strainSize,
				      const PylithScalar* initialStress,
				      const int initialStressSize,
				      const PylithScalar* initialStrain,
				      const int initialStrainSize);

  /** Update state variables after solve as an elastic material.
   *
   * @param stateVars State variables at location.
   * @param numStateVars Number of state variables.
   * @param properties Properties at location.
   * @param numProperties Number of properties.
   * @param totalStrain Total strain at location.
   * @param strainSize Size of strain tensor.
   * @param initialStress Initial stress tensor at location.
   * @param initialStressSize Size of initial stress array.
   * @param initialStrain Initial strain tensor at location.
   * @param initialStrainSize Size of initial strain array.
   */
  void _updateStateVarsElastic(PylithScalar* const stateVars,
			       const int numStateVars,
			       const PylithScalar* properties,
			       const int numProperties,
			       const PylithScalar* totalStrain,
			       const int strainSize,
			       const PylithScalar* initialStress,
			       const int initialStressSize,
			       const PylithScalar* initialStrain,
			       const int initialStrainSize);

  /** Update state variables after solve as a viscoelastic material.
   *
   * @param stateVars State variables at location.
   * @param numStateVars Number of state variables.
   * @param properties Properties at location.
   * @param numProperties Number of properties.
   * @param totalStrain Total strain at location.
   * @param strainSize Size of strain tensor.
   * @param initialStress Initial stress tensor at location.
   * @param initialStressSize Size of initial stress array.
   * @param initialStrain Initial strain tensor at location.
   * @param initialStrainSize Size of initial strain array.
   */
  void _updateStateVarsViscoelastic(PylithScalar* const stateVars,
				    const int numStateVars,
				    const PylithScalar* properties,
				    const int numProperties,
				    const PylithScalar* totalStrain,
				    const int strainSize,
				    const PylithScalar* initialStress,
				    const int initialStressSize,
				    const PylithScalar* initialStrain,
				    const int initialStrainSize);

  /** Compute viscous strains (state variables) for the current time
   * step.
   *
   * @param stateVars State variables at location.
   * @param numStateVars Number of state variables.
   * @param properties Properties at location.
   * @param numProperties Number of properties.
   * @param totalStrain Total strain at location.
   * @param strainSize Size of strain tensor.
   * @param initialStress Initial stress tensor at location.
   * @param initialStressSize Size of initial stress array.
   * @param initialStrain Initial strain tensor at location.
   * @param initialStrainSize Size of initial strain array.
   */
  void _computeStateVars(const PylithScalar* stateVars,
			 const int numStateVars,
			 const PylithScalar* properties,
			 const int numProperties,
			 const PylithScalar* totalStrain,
			 const int strainSize,
			 const PylithScalar* initialStress,
			 const int initialStressSize,
			 const PylithScalar* initialStrain,
			 const int initialStrainSize);

  // PRIVATE MEMBERS ////////////////////////////////////////////////////
private :

  scalar_array _viscousStrain; ///< Array for viscous strain tensor

  /// Method to use for _calcElasticConsts().
  calcElasticConsts_fn_type _calcElasticConstsFn;

  /// Method to use for _calcStress().
  calcStress_fn_type _calcStressFn;

  /// Method to use for _updateStateVars().
  updateStateVars_fn_type _updateStateVarsFn;

  // PRIVATE MEMBERS ////////////////////////////////////////////////////
private :

  static const int p_density;
  static const int p_mu;
  static const int p_lambda;
  static const int p_maxwellTime;
  static const int db_density;
  static const int db_vs;
  static const int db_vp;
  static const int db_viscosity;

  static const int s_totalStrain;
  static const int s_viscousStrain;
  static const int db_totalStrain;
  static const int db_viscousStrain;

  // NOT IMPLEMENTED ////////////////////////////////////////////////////
private :

  /// Not implemented
  MaxwellIsotropic3D(const MaxwellIsotropic3D&);

  /// Not implemented
  const MaxwellIsotropic3D& operator=(const MaxwellIsotropic3D&);

}; // class MaxwellIsotropic3D

#include "MaxwellIsotropic3D.icc" // inline methods

#endif // pylith_materials_maxwellisotropic3d_hh


// End of file 
