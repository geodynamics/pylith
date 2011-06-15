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

/** @file modulesrc/materials/ElasticMaterial.i
 *
 * Python interface to C++ abstract base ElasticMaterial.
 */

namespace pylith {
  namespace materials {

    class ElasticMaterial : public Material
    { // class ElasticMaterial

      // PUBLIC METHODS /////////////////////////////////////////////////
    public :

      /** Default constructor.
       *
       * @param dimension Spatial dimension associated with material.
       * @param tensorSize Array of names of database values for material.
       * @param numElasticConsts Number of elastic constants.
       * @param metadata Metadata for physical properties and state variables.
       */
      ElasticMaterial(const int dimension,
		      const int tensorSize,
		      const int numElasticConsts,
		      const Metadata& metadata);

      /// Destructor.
      virtual
      ~ElasticMaterial(void);

      /// Deallocate PETSc and local data structures.
      virtual
      void deallocate(void);
  
      /** Set database for initial stress state.
       *
       * @param db Spatial database.
       */
      void dbInitialStress(spatialdata::spatialdb::SpatialDB* db);
      
      /** Set database for initial strain state.
       *
       * @param db Spatial database.
       */
      void dbInitialStrain(spatialdata::spatialdb::SpatialDB* db);

      /** Get flag indicating whether material implements an empty
       * _updateProperties() method.
       *
       * @returns False if _updateProperties() is empty, true otherwise.
       */
      bool hasStateVars(void) const;

      /** Get stable time step for implicit time integration.
       *
       * @pre Must call retrievePropsAndVars for cell before calling
       * stableTimeStep().
       *
       * Default is MAXFLOAT (or 1.0e+30 if MAXFLOAT is not defined in math.h).
       *
       * @returns Time step
       */
      virtual
      double stableTimeStepImplicit(const pylith::topology::Mesh& mesh);

      /** Set whether elastic or inelastic constitutive relations are used.
       *
       * @param flag True to use elastic, false to use inelastic.
       */
      virtual
      void useElasticBehavior(const bool flag);

      /** Get initial stress/strain fields.
       *
       * @returns Initial stress field.
       */
      const pylith::topology::Fields<pylith::topology::Field<pylith::topology::Mesh> >* initialFields(void) const;

      // PROTECTED METHODS //////////////////////////////////////////////
    protected :

      /** Compute density from properties.
       *
       * @param density Array for density.
       * @param properties Properties at location.
       * @param numProperties Number of properties.
       */
      virtual
      void _calcDensity(double* const density,
			const double* properties,
			const int numProperties,
			const double* stateVars,
			const int numStateVars) = 0;
      
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
      virtual
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
		       const bool computeStateVars) = 0;
      
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
      virtual
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
			      const int initialStrainSize) = 0;
      
      /** Get stable time step for implicit time integration.
       *
       * @param properties Properties at location.
       * @param numProperties Number of properties.
       * @param stateVars State variables at location.
       * @param numStateVars Number of state variables.
       *
       * @returns Time step
       */
      virtual
      double _stableTimeStepImplicit(const double* properties,
				     const int numProperties,
				     const double* stateVars,
				     const int numStateVars) const = 0;
      
    }; // class ElasticMaterial

  } // materials
} // pylith


// End of file
