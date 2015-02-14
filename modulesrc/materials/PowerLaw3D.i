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

/** @file modulesrc/materials/PowerLaw3D.i
 *
 * @brief Python interface to C++ PowerLaw3D object.
 */

namespace pylith {
  namespace materials {

    class pylith::materials::PowerLaw3D : public ElasticMaterial
    { // class PowerLaw3D

      // PUBLIC METHODS /////////////////////////////////////////////////
    public :

      /// Default constructor
      PowerLaw3D(void);
      
      /// Destructor
      ~PowerLaw3D(void);
      
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
      
      // PROTECTED METHODS //////////////////////////////////////////////
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
      
      /** Compute initial state variables from values in spatial database.
       *
       * @param stateValues Array of state variable values.
       * @param dbValues Array of database values.
       */
      void _dbToStateVars(PylithScalar* const stateValues,
			  const pylith::scalar_array& dbValues);
      
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
       * @param computeStateVars Flag indicating to compute updated
       * state variables.
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

    }; // class PowerLaw3D

  } // materials
} // pylith

// End of file 
