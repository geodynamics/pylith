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

/** @file modulesrc/friction/FrictionModel.i
 *
 * Python interface to C++ abstract base FrictionModel.
 */

namespace pylith {
  namespace friction {

    class FrictionModel
    { // class FrictionModel

      // PUBLIC METHODS /////////////////////////////////////////////////
    public :
      
      /** Default constructor.
       *
       * @param metadata Metadata for physical properties and state variables.
       */
      FrictionModel(const pylith::materials::Metadata& metadata);

      /// Destructor.
      virtual
      ~FrictionModel(void);

      /// Deallocate PETSc and local data structures.
      virtual
      void deallocate(void);
  
      /** Set label of friction model.
       *
       * @param value Label of friction model.
       */
      void label(const char* value);

      /** Get label of friction model.
       *
       * @returns Label of friction model.
       */
      const char* label(void) const;

      /** Set current time step.
       *
       * @param dt Current time step.
       */
      virtual
      void timeStep(const float dt);

      /** Get current time step.
       *
       * @returns Current time step.
       */
      PylithScalar timeStep(void) const;

      /** Set database for physical property parameters.
       *
       * @param value Pointer to database.
       */
      void dbProperties(spatialdata::spatialdb::SpatialDB* value);

      /** Set database for initial state variables.
       *
       * @param value Pointer to database.
       */
      void dbInitialState(spatialdata::spatialdb::SpatialDB* value);

      /** Set scales used to nondimensionalize physical properties.
       *
       * @param dim Nondimensionalizer
       */
      void normalizer(const spatialdata::units::Nondimensional& dim);

      /** Initialize friction model by getting physical property
       * parameters from database.
       *
       * @pre Must call Quadrature::computeGeometry() before calling
       * initialize().
       *
       * @param mesh Finite-element mesh of subdomain.
       * @param quadrature Quadrature for finite-element integration
       */
      virtual
      void initialize(const pylith::topology::Mesh& mesh,
		      pylith::feassemble::Quadrature* quadrature);
  
      /** Check whether friction model has a field as a property or
       * state variable.
       *
       * @param name Name of field.
       *
       * @returns True if friction model has field as the given
       * property or state variable, false otherwise.
       */
      bool hasPropStateVar(const char* name);

      /** Get physical property or state variable field. Data is returned
       * via the argument.
       *
       * @param name Name of field to retrieve.
       * @returns Field over fault interface cells.
       */
      const pylith::topology::Field& getField(const char* name);

      /** Get the field with all properties and state variables.
       *
       * @returns Properties field.
       */
      const pylith::topology::Fields& fieldsPropsStateVars() const;

      /** Retrieve parameters for physical properties and state variables
       * for vertex.
       *
       * @param vertex Finite-element vertex on friction interface.
       */
      void retrievePropsStateVars(const int vertex);

      /** Compute friction at vertex.
       *
       * @pre Must call retrievePropsAndVars for cell before calling
       * calcFriction().
       *
       * @param t Current time.
       * @param slip Current slip at location.
       * @param slipRate Current slip rate at location.
       * @param normalTraction Normal traction at location.
       *
       * @returns Friction (magnitude of shear traction) at vertex.
       */
      PylithScalar calcFriction(const PylithScalar t,
				const PylithScalar slip,
				const PylithScalar slipRate,
				const PylithScalar normalTraction);
  

      /** Compute derivative of friction with slip at vertex.
       *
       * @pre Must call retrievePropsAndVars for cell before calling
       * calcFriction().
       *
       * @param t Time in simulation.
       * @param slip Current slip at location.
       * @param slipRate Current slip rate at location.
       * @param normalTraction Normal traction at location.
       *
       * @returns Derivative of friction (magnitude of shear traction).
       */
      PylithScalar calcFrictionDeriv(const PylithScalar t,
				     const PylithScalar slip,
				     const PylithScalar slipRate,
				     const PylithScalar normalTraction);
  
      /** Compute friction at vertex.
       *
       * @pre Must call retrievePropsAndVars for cell before calling
       * calcFriction().
       *
       * @param t Current time.
       * @param slip Current slip at location.
       * @param slipRate Current slip rate at location.
       * @param normalTraction Normal traction at location.
       * @param vertex Finite-element vertex on friction interface.
       */
      void updateStateVars(const PylithScalar t,
			   const PylithScalar slip,
			   const PylithScalar slipRate,
			   const PylithScalar normalTraction,
			   const int vertex);
  
      // PROTECTED METHODS //////////////////////////////////////////////
    protected :
      
      /// These methods should be implemented by every constitutive model.

      /** Compute properties from values in spatial database.
       *
       * @param propValues Array of property values.
       * @param dbValues Array of database values.
       */
      virtual
      void _dbToProperties(PylithScalar* const propValues,
			   const scalar_array& dbValues) const = 0;

      /** Nondimensionalize properties.
       *
       * @param values Array of property values.
       * @param nvalues Number of values.
       */
      virtual
      void _nondimProperties(PylithScalar* const values,
			     const int nvalues) const = 0;

      /** Dimensionalize properties.
       *
       * @param values Array of property values.
       * @param nvalues Number of values.
       */
      virtual
      void _dimProperties(PylithScalar* const values,
			  const int nvalues) const = 0;

      /** Compute initial state variables from values in spatial database.
       *
       * @param stateValues Array of state variable values.
       * @param dbValues Array of database values.
       */
      virtual
      void _dbToStateVars(PylithScalar* const stateValues,
			  const scalar_array& dbValues) const;

      /** Nondimensionalize state variables.
       *
       * @param values Array of initial state values.
       * @param nvalues Number of values.
       */
      virtual
      void _nondimStateVars(PylithScalar* const values,
			    const int nvalues) const;
  
      /** Dimensionalize state variables.
       *
       * @param values Array of initial state values.
       * @param nvalues Number of values.
       */
      virtual
      void _dimStateVars(PylithScalar* const values,
			 const int nvalues) const;

      /** Compute friction from properties and state variables.
       *
       * @param t Current time.
       * @param slip Current slip at location.
       * @param slipRate Current slip rate at location.
       * @param normalTraction Normal traction at location.
       * @param properties Properties at location.
       * @param numProperties Number of properties.
       * @param stateVars State variables at location.
       * @param numStateVars Number of state variables.
       */
      virtual
      PylithScalar _calcFriction(const PylithScalar t,
				 const PylithScalar slip,
				 const PylithScalar slipRate,
				 const PylithScalar normalTraction,
				 const PylithScalar* properties,
				 const int numProperties,
				 const PylithScalar* stateVars,
				 const int numStateVars) = 0;

      /** Compute derivative friction with slip from properties and state
       * variables.
       *
       * @param t Time in simulation.
       * @param slip Current slip at location.
       * @param slipRate Current slip rate at location.
       * @param normalTraction Normal traction at location.
       * @param properties Properties at location.
       * @param numProperties Number of properties.
       * @param stateVars State variables at location.
       * @param numStateVars Number of state variables.
       *
       * @returns Derivative of friction (magnitude of shear traction) at vertex.
       */
      virtual
      PylithScalar _calcFrictionDeriv(const PylithScalar t,
				      const PylithScalar slip,
				      const PylithScalar slipRate,
				      const PylithScalar normalTraction,
				      const PylithScalar* properties,
				      const int numProperties,
				      const PylithScalar* stateVars,
				      const int numStateVars) = 0;
  
      /** Update state variables (for next time step).
       *
       * @param t Current time.
       * @param stateVars State variables at location.
       * @param numStateVars Number of state variables.
       * @param properties Properties at location.
       * @param numProperties Number of properties.
       */
      virtual
      void _updateStateVars(const PylithScalar t,
			    const PylithScalar slip,
			    const PylithScalar slipRate,
			    PylithScalar* const stateVars,
			    const int numStateVars,
			    const PylithScalar* properties,
			    const int numProperties);

    }; // class FrictionModel

  } // friction
} // pylith


// End of file 
