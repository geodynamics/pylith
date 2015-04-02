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
       * Default is MAXFLOAT (or 1.0e+30 if MAXFLOAT is not defined in math.h).
       *
       * @returns Time step
       */
      virtual
      PylithScalar stableTimeStepImplicit(const pylith::topology::Mesh& mesh);

      /** Get stable time step for explicit time integration.
       *
       * Default is MAXFLOAT (or 1.0e+30 if MAXFLOAT is not defined in math.h).
       *
       * @param mesh Finite-element mesh.
       * @param quadrature Quadrature for finite-element integration
       * @returns Time step
       */
      virtual
      PylithScalar stableTimeStepExplicit(const pylith::topology::Mesh& mesh,
					  pylith::feassemble::Quadrature* quadrature);
      
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
      const pylith::topology::Fields* initialFields(void) const;

      // PROTECTED METHODS //////////////////////////////////////////////
    protected :

      /** Compute density from properties.
       *
       * @param density Array for density.
       * @param properties Properties at location.
       * @param numProperties Number of properties.
       */
      virtual
      void _calcDensity(PylithScalar* const density,
			const PylithScalar* properties,
			const int numProperties,
			const PylithScalar* stateVars,
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
      PylithScalar _stableTimeStepImplicit(const PylithScalar* properties,
				     const int numProperties,
				     const PylithScalar* stateVars,
				     const int numStateVars) const = 0;
      
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
      virtual
      PylithScalar _stableTimeStepExplicit(const PylithScalar* properties,
					   const int numProperties,
					   const PylithScalar* stateVars,
					   const int numStateVars,
					   const double minCellWidth) const = 0;
  
    }; // class ElasticMaterial

  } // materials
} // pylith


// End of file
