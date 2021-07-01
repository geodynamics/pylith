// -*- C++ -*-
//
// ----------------------------------------------------------------------
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University at Buffalo
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2021 University of California, Davis
//
// See LICENSE.md for license information.
//
// ----------------------------------------------------------------------
//

/* @brief C++ PlaneStrainState object that stores stress and strain as
 * state variables.
 *
 * This objects demonstrates how to extend PyLith by adding bulk
 * constitutive models. This PlaneStrainState object is identical to
 * the ElasticPlaneStrain object in PyLith except that it stores the
 * current stress and strain as state variables.
 *
 * 2-D, isotropic, linear elastic material for plane strain. The
 * physical properties are specified using density, shear-wave speed,
 * and compressional-wave speed. The physical properties are stored
 * internally using density, lambda, and mu, which are directly
 * related to the elasticity constants used in the finite-element
 * integration.
 *
 * $\sigma - \sigma_0 = C (\epsilon - \epsilon_0)
 *
 * This implies that when $\epsilon = \epsilon_0$, $\sigma =
 * \sigma_0$.
 */

#if !defined(pylith_materials_planestrainstate_hh)
#define pylith_materials_planestrainstate_hh

#include "pylith/materials/ElasticMaterial.hh" // ISA ElasticMaterial

namespace contrib {
  namespace materials {
    class PlaneStrainState;
  } // materials
} // contrib

// PlaneStrainState -----------------------------------------------------
class contrib::materials::PlaneStrainState : public pylith::materials::ElasticMaterial
{ // class PlaneStrainState
  friend class TestPlaneStrainState; // unit testing

  // --------------------------------------------------------------------
  // All of these functions are required to satisfy the PyLith
  // interface for a bulk constitutive model.
  // --------------------------------------------------------------------

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Default constructor
  PlaneStrainState(void);

  /// Destructor
  ~PlaneStrainState(void);

  /** Get stable time step for implicit time integration.
   *
   * Default is MAXDOUBLE (or 1.0e+30 if MAXFLOAT is not defined in math.h).
   *
   * This function is optional but provides an optimized
   * implementation of the more general
   * ElasticMaterial::stableTimeStepImplicit().
   *
   * @param mesh Finite-element mesh.
   * @returns Time step
   */
  PylithScalar stableTimeStepImplicit(const pylith::topology::Mesh& mesh);

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
		       const pylith::scalar_array& dbValues);

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

  /** Compute density from properties.
   *
   * @param density Array for density.
   * @param properties Properties at location.
   * @param numProperties Number of properties.
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
  
  // --------------------------------------------------------------------
  // Optional function in the PyLith interface for a bulk constitutive
  // model. Even though this function is optional, for it to be used
  // it the interface must exactly matched the one specified in
  // ElasticMaterial.
  // --------------------------------------------------------------------

  // PROTECTED METHODS //////////////////////////////////////////////////
protected :

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

  // PRIVATE MEMBERS ////////////////////////////////////////////////////
private :

  // --------------------------------------------------------------------
  // We use these constants for consistent access into the arrays of
  // physical properties and state variables.
  // --------------------------------------------------------------------

  static const int p_density;
  static const int p_mu;
  static const int p_lambda;
  static const int db_density;
  static const int db_vs;
  static const int db_vp;

  static const int s_totalStrain;
  static const int s_stress;

  // NOT IMPLEMENTED ////////////////////////////////////////////////////
private :

  /// Not implemented
  PlaneStrainState(const PlaneStrainState& m);

  /// Not implemented
  const PlaneStrainState& operator=(const PlaneStrainState& m);

}; // class PlaneStrainState

#endif // pylith_materials_planestrainstate_hh


// End of file 
