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
// Copyright (c) 2010-2015 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

/** @file libsrc/materials/DruckerPragerPlaneStrain.hh
 *
 * @brief 2-D, plane strain, Drucker-Prager elastic/perfectly plastic material. 
 */

#if !defined(pylith_materials_druckerpragerplanestrain_hh)
#define pylith_materials_druckerpragerplanestrain_hh

// Include directives ---------------------------------------------------
#include "ElasticMaterial.hh" // ISA ElasticMaterial

// DruckerPragerPlaneStrain ---------------------------------------------
/** @brief 2-D, plane strain, Drucker-Prager elastic/perfectly plastic material.
 *
 * The physical properties are specified using density, shear-wave
 * speed, friction angle, cohesion, dilatation angle, and
 * compressional-wave speed.  The physical properties are stored
 * internally using density, lambda, mu, which are directly related to
 * the elasticity constants used in the finite-element
 * integration. The plasticity information is retained as alpha_yield,
 * beta, and alpha_flow.
 */

class pylith::materials::DruckerPragerPlaneStrain : public ElasticMaterial
{ // class DruckerPragerPlaneStrain
  friend class TestDruckerPragerPlaneStrain; // unit testing

  // PUBLIC ENUMS ///////////////////////////////////////////////////////
public :

  enum FitMohrCoulombEnum {
    MOHR_COULOMB_CIRCUMSCRIBED=0, 
    MOHR_COULOMB_MIDDLE=1,
    MOHR_COULOMB_INSCRIBED=2,
  }; // FitMohrCoulombType

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Default constructor
  DruckerPragerPlaneStrain(void);

  /// Destructor
  ~DruckerPragerPlaneStrain(void);

  /** Set fit to Mohr-Coulomb surface.
   *
   * @param value Mohr-Coulomb surface match type.
   */
  void fitMohrCoulomb(FitMohrCoulombEnum value);

  /** Set flag for whether to allow tensile yield.
   *
   * @param flag True if tensile yield is allowed.
   */
  void allowTensileYield(const bool flag);

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

  /** Nondimensionalize state variables..
   *
   * @param values Array of state variables.
   * @param nvalues Number of values.
   */
  void _nondimStateVars(PylithScalar* const values,
			const int nvalues) const;

  /** Dimensionalize state variables.
   *
   * @param values Array of state variables.
   * @param nvalues Number of values.
   */
  void _dimStateVars(PylithScalar* const values,
		     const int nvalues) const;

  /** Compute density from properties.
   *
   * @param density Array for density.
   * @param properties Properties at location.
   * @param numProperties Number of properties.
   * @param stateVars State variables at location.
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
   * variables are updated (but not stored) when computing the stresses.
   *
   * @param stress Array for stress tensor.
   * @param stressSize Size of stress tensor.
   * @param properties Properties at location.
   * @param numProperties Number of properties.
   * @param stateVars State variables at location.
   * @param numStateVars Number of state variables.
   * @param totalStrain Total strain at location.
   * @param strainSize Size of strain tensor.
   * @param initialStress Initial stress values.
   * @param initialStressSize Size of initial stress array.
   * @param initialStrain Initial strain values.
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
   * @param initialStress Initial stress values.
   * @param initialStressSize Size of initial stress array.
   * @param initialStrain Initial strain values.
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
  
  /** Update state variables (for next time step).
   *
   * @param stateVars State variables at location.
   * @param numStateVars Number of state variables.
   * @param properties Properties at location.
   * @param numProperties Number of properties.
   * @param totalStrain Total strain at location.
   * @param strainSize Size of strain tensor.
   * @param initialStress Initial stress values.
   * @param initialStressSize Size of initial stress array.
   * @param initialStrain Initial strain values.
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

  // PRIVATE TYPEDEFS ///////////////////////////////////////////////////
private :

  /// Member prototype for _calcStress()
  typedef void (pylith::materials::DruckerPragerPlaneStrain::*calcStress_fn_type)
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
  typedef void (pylith::materials::DruckerPragerPlaneStrain::*calcElasticConsts_fn_type)
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
  typedef void (pylith::materials::DruckerPragerPlaneStrain::*updateStateVars_fn_type)
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
   * @param properties Properties at locations.
   * @param numProperties Number of properties.
   * @param stateVars State variables at locations.
   * @param numStateVars Number of state variables.
   * @param totalStrain Total strain at locations.
   * @param strainSize Size of strain tensor.
   * @param initialStress Initial stress values.
   * @param initialStressSize Size of initial stress array.
   * @param initialStrain Initial strain values.
   * @param initialStrainSize Size of initial strain array.
   * @param computeStateVars Flag indicating to compute updated state vars.
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

  /** Compute stress tensor from properties as an elastoplastic material.
   *
   * @param stress Array for stress tensor.
   * @param stressSize Size of stress tensor.
   * @param properties Properties at locations.
   * @param numProperties Number of properties.
   * @param stateVars State variables at locations.
   * @param numStateVars Number of state variables.
   * @param totalStrain Total strain at locations.
   * @param strainSize Size of strain tensor.
   * @param initialStress Initial stress values.
   * @param initialStressSize Size of initial stress array.
   * @param initialStrain Initial strain values.
   * @param initialStrainSize Size of initial strain array.
   * @param computeStateVars Flag indicating to compute updated state vars.
   */
  void _calcStressElastoplastic(PylithScalar* const stress,
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
   * @param stateVars State variables at locations.
   * @param numStateVars Number of state variables.
   * @param totalStrain Total strain at location.
   * @param strainSize Size of strain tensor.
   * @param initialStress Initial stress values.
   * @param initialStressSize Size of initial stress array.
   * @param initialStrain Initial strain values.
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

  /** Compute derivatives of elasticity matrix from properties as an
   * elastoplastic material.
   *
   * @param elasticConsts Array for elastic constants.
   * @param numElasticConsts Number of elastic constants.
   * @param properties Properties at location.
   * @param numProperties Number of properties.
   * @param stateVars State variables at locations.
   * @param numStateVars Number of state variables.
   * @param totalStrain Total strain at location.
   * @param strainSize Size of strain tensor.
   * @param initialStress Initial stress values.
   * @param initialStressSize Size of initial stress array.
   * @param initialStrain Initial strain values.
   * @param initialStrainSize Size of initial strain array.
   */
  void _calcElasticConstsElastoplastic(PylithScalar* const elasticConsts,
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
   * @param stateVars State variables at locations.
   * @param numStateVars Number of state variables.
   * @param properties Properties at location.
   * @param numProperties Number of properties.
   * @param totalStrain Total strain at location.
   * @param strainSize Size of strain tensor.
   * @param initialStress Initial stress values.
   * @param initialStressSize Size of initial stress array.
   * @param initialStrain Initial strain values.
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

  /** Update state variables after solve as an elastoplastic material.
   *
   * @param properties Properties at location.
   * @param numProperties Number of properties.
   * @param totalStrain Total strain at location.
   * @param strainSize Size of strain tensor.
   * @param initialState Initial state values.
   * @param initialStateSize Size of initial state array.
   */
  void _updateStateVarsElastoplastic(PylithScalar* const stateVars,
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

  /// Method to use for _calcElasticConsts().
  calcElasticConsts_fn_type _calcElasticConstsFn;

  /// Method to use for _calcStress().
  calcStress_fn_type _calcStressFn;

  /// Method to use for _updateStateVars().
  updateStateVars_fn_type _updateStateVarsFn;

  /// Fit to Mohr Coulomb surface
  FitMohrCoulombEnum _fitMohrCoulomb;

  /// Whether to allow tensile yield
  bool _allowTensileYield;

  static const int p_density;
  static const int p_mu;
  static const int p_lambda;
  static const int p_alphaYield;
  static const int p_beta;
  static const int p_alphaFlow;
  static const int db_density;
  static const int db_vs;
  static const int db_vp;
  static const int db_frictionAngle;
  static const int db_cohesion;
  static const int db_dilatationAngle;

  static const int s_stressZZInitial;
  static const int s_plasticStrain;
  static const int db_stressZZInitial;
  static const int db_plasticStrain;

  // NOT IMPLEMENTED ////////////////////////////////////////////////////
private :

  DruckerPragerPlaneStrain(const DruckerPragerPlaneStrain&); ///< Not implemented
  const DruckerPragerPlaneStrain& operator=(const DruckerPragerPlaneStrain&); ///< Not implemented

}; // class DruckerPragerPlaneStrain

#include "DruckerPragerPlaneStrain.icc" // inline methods

#endif // pylith_materials_druckerpragerplanestrain_hh


// End of file 
