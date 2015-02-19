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

/** @file libsrc/materials/PowerLawPlaneStrain.hh
 *
 * @brief 2-D, plane strain, power-law viscoelastic material. 
 */

#if !defined(pylith_materials_powerlawplanestrain_hh)
#define pylith_materials_powerlawplanestrain_hh

// Include directives ---------------------------------------------------
#include "ElasticMaterial.hh" // ISA ElasticMaterial

// PowerlawPlaneStrain----------------------------------------------------------
/** @brief 2-D, plane strain, power-law viscoelastic material. 
 *
 * The physical properties are specified using density, shear-wave
 * speed, viscosity coefficient, power-law exponent, and
 * compressional-wave speed.  The physical properties are stored
 * internally using density, lambda, mu, which are directly related to
 * the elasticity constants used in the finite-element
 * integration. The viscosity information is retained as specified.
 */

class pylith::materials::PowerLawPlaneStrain : public ElasticMaterial
{ // class PowerLawPlaneStrain
  friend class TestPowerLawPlaneStrain; // unit testing

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Default constructor
  PowerLawPlaneStrain(void);

  /// Destructor
  ~PowerLawPlaneStrain(void);

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

  /** Compute effective stress function.
   *
   * @param effStressTpdt Effective stress value.
   *
   * @returns Effective stress function value.
   */
  PylithScalar effStressFunc(const PylithScalar effStressTpdt);

  /** Compute effective stress function derivative.
   *
   * @param effStressTpdt Effective stress value.
   *
   * @returns Effective stress function derivative value.
   */
  PylithScalar effStressDerivFunc(const PylithScalar effStressTpdt);

  /** Compute effective stress function and derivative.
   *
   * @param func Returned effective stress function value.
   * @param dfunc Returned effective stress function derivative value.
   * @param effStressTpdt Effective stress value.
   *
   */
  void effStressFuncDerivFunc(PylithScalar* func,
			      PylithScalar* dfunc,
			      const PylithScalar effStressTpdt);

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
  typedef void (pylith::materials::PowerLawPlaneStrain::*calcStress_fn_type)
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
  typedef void (pylith::materials::PowerLawPlaneStrain::*calcElasticConsts_fn_type)
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
  typedef void (pylith::materials::PowerLawPlaneStrain::*updateStateVars_fn_type)
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

  /** Compute stress tensor from properties as an viscoelastic material.
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

  /** Compute derivatives of elasticity matrix from properties as a
   * viscoelastic material.
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

  /** Update state variables after solve as a viscoelastic material.
   *
   * @param properties Properties at location.
   * @param numProperties Number of properties.
   * @param totalStrain Total strain at location.
   * @param strainSize Size of strain tensor.
   * @param initialState Initial state values.
   * @param initialStateSize Size of initial state array.
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


  // PRIVATE STRUCTS ////////////////////////////////////////////////////
private :

  struct EffStressStruct {
    PylithScalar ae;
    PylithScalar b;
    PylithScalar c;
    PylithScalar d;
    PylithScalar alpha;
    PylithScalar dt;
    PylithScalar effStressT;
    PylithScalar powerLawExp;
    PylithScalar referenceStrainRate;
    PylithScalar referenceStress;
  };

  // PRIVATE MEMBERS ////////////////////////////////////////////////////
private :

  /// Structure to hold parameters for effective stress computation.
  EffStressStruct _effStressParams;

  /// Method to use for _calcElasticConsts().
  calcElasticConsts_fn_type _calcElasticConstsFn;

  /// Method to use for _calcStress().
  calcStress_fn_type _calcStressFn;

  /// Method to use for _updateStateVars().
  updateStateVars_fn_type _updateStateVarsFn;

  static const int p_density;
  static const int p_mu;
  static const int p_lambda;
  static const int p_referenceStrainRate;
  static const int p_referenceStress;
  static const int p_powerLawExponent;
  static const int db_density;
  static const int db_vs;
  static const int db_vp;
  static const int db_referenceStrainRate;
  static const int db_referenceStress;
  static const int db_powerLawExponent;

  static const int s_stressZZInitial;
  static const int s_viscousStrain;
  static const int s_stress4;
  static const int db_stressZZInitial;
  static const int db_viscousStrain;
  static const int db_stress4;

  // NOT IMPLEMENTED ////////////////////////////////////////////////////
private :

  PowerLawPlaneStrain(const PowerLawPlaneStrain&); ///< Not implemented
  const PowerLawPlaneStrain& operator=(const PowerLawPlaneStrain&); ///< Not implemented

}; // class PowerLawPlaneStrain

#include "PowerLawPlaneStrain.icc" // inline methods

#endif // pylith_materials_powerlawplanestrain_hh


// End of file 
