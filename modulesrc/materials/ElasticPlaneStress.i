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
// Copyright (c) 2010-2011 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

/** @file modulesrc/materials/ElasticPlaneStress.i
 *
 * Python interface to C++ ElasticPlaneStress object.
 */

namespace pylith {
  namespace materials {

    class ElasticPlaneStress : public ElasticMaterial
    { // class ElasticPlaneStress

      // PUBLIC METHODS /////////////////////////////////////////////////
    public :
      
      /// Default constructor
      ElasticPlaneStress(void);
      
      /// Destructor
      ~ElasticPlaneStress(void);
      
      /** Get stable time step for implicit time integration.
       *
       * Default is MAXDOUBLE (or 1.0e+30 if MAXFLOAT is not defined
       * in math.h).
       *
       * @param mesh Finite-element mesh.
       * @returns Time step
       */
      double stableTimeStepImplicit(const pylith::topology::Mesh& mesh);

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
      void _dbToProperties(double* const propValues,
			   const pylith::double_array& dbValues);

      /** Nondimensionalize properties.
       *
       * @param values Array of property values.
       * @param nvalues Number of values.
       */
      void _nondimProperties(double* const values,
			     const int nvalues) const;

      /** Dimensionalize properties.
       *
       * @param values Array of property values.
       * @param nvalues Number of values.
       */
      void _dimProperties(double* const values,
			  const int nvalues) const;

      /** Compute density from properties.
       *
       * @param density Array for density.
       * @param properties Properties at location.
       * @param numProperties Number of properties.
       */
      void _calcDensity(double* const density,
			const double* properties,
			const int numProperties,
			const double* stateVars,
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
       * @param computeStateVars Flag indicating to compute updated
       * state variables.
       */
      void _calcStress(double* const stress,
		       const int stressSize,
		       const double* properties,
		       const int numProperties,
		       const double* stateVars,
		       const int numStateVars,
		       const double* totalStrain,
		       const int strainSize,
		       const double* initialStress,
		       const int initialStressSize,
		       const double* initialStrain,
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
      void _calcElasticConsts(double* const elasticConsts,
			      const int numElasticConsts,
			      const double* properties,
			      const int numProperties,
			      const double* stateVars,
			      const int numStateVars,
			      const double* totalStrain,
			      const int strainSize,
			      const double* initialStress,
			      const int initialStressSize,
			      const double* initialStrain,
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
      double _stableTimeStepImplicit(const double* properties,
				     const int numProperties,
				     const double* stateVars,
				     const int numStateVars) const;

    }; // class ElasticPlaneStress

  } // materials
} // pylith


// End of file 
