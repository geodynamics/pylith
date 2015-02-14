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

/** @file libsrc/friction/FrictionModel.hh
 *
 * @brief C++ abstract base class for FrictionModel object.
 */

#if !defined(pylith_friction_frictionmodel_hh)
#define pylith_friction_frictionmodel_hh

// Include directives ---------------------------------------------------
#include "frictionfwd.hh" // forward declarations

#include "pylith/topology/topologyfwd.hh" // forward declarations
#include "pylith/feassemble/feassemblefwd.hh" // forward declarations
#include "spatialdata/spatialdb/spatialdbfwd.hh" // forward declarations
#include "spatialdata/units/unitsfwd.hh" // forward declarations

#include "pylith/materials/Metadata.hh" // HASA Metadata

#include <string> // HASA std::string

// FrictionModel --------------------------------------------------------
/** @brief C++ abstract base class for FrictionModel object.
 *
 * Interface definition for a friction model. The physical properties
 * for the friction model are associated with the constants in the
 * constitutive model.
 */

class pylith::friction::FrictionModel
{ // class FrictionModel
  friend class TestFrictionModel; // unit testing

  // PUBLIC METHODS /////////////////////////////////////////////////////
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
  void timeStep(const PylithScalar dt);

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
  void initialize(const topology::Mesh& mesh,
		  feassemble::Quadrature* quadrature);
  
  /** Check whether friction model has a field as a property or state
   * variable.
   *
   * @param name Name of field.
   *
   * @returns True if friction model has field as the given property
   * or state variable, false otherwise.
   */
  bool hasPropStateVar(const char* name);
  
  /** Return the property and state variable metadata.
   *
   * @returns Metadata for properties and state variables.
   */
  const pylith::materials::Metadata& getMetadata();

  /** Get physical property or state variable field. Data is returned
   * via the argument.
   *
   * @param name Name of field to retrieve.
   * @returns Field over fault interface cells.
   */
  const topology::Field& getField(const char* name);

  /** Get the field with all properties and state variables.
   *
   * @returns Properties field.
   */
  const topology::Fields& fieldsPropsStateVars() const;

  /** Retrieve properties and state variables for a point.
   *
   * @param point Finite-element point.
   */
  void retrievePropsStateVars(const int point);

  /** Compute friction at vertex.
   *
   * @pre Must call retrievePropsAndVars for cell before calling
   * calcFriction().
   *
   * @param t Time in simulation.
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
  
  /** Compute update to state variables at vertex.
   *
   * @pre Must call retrievePropsAndVars for cell before calling
   * calcFriction().
   *
   * @param t Time in simulation.
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
  
  // PROTECTED METHODS //////////////////////////////////////////////////
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
   * @param t Time in simulation.
   * @param slip Current slip at location.
   * @param slipRate Current slip rate at location.
   * @param normalTraction Normal traction at location.
   * @param properties Properties at location.
   * @param numProperties Number of properties.
   * @param stateVars State variables at location.
   * @param numStateVars Number of state variables.
   *
   * @returns Friction (magnitude of shear traction) at vertex.
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
   * @param t Time in simulation.
   * @param slip Current slip at location.
   * @param slipRate Current slip rate at location.
   * @param normalTraction Normal traction at location.
   * @param stateVars State variables at location.
   * @param numStateVars Number of state variables.
   * @param properties Properties at location.
   * @param numProperties Number of properties.
   */
  virtual
  void _updateStateVars(const PylithScalar t,
			const PylithScalar slip,
			const PylithScalar slipRate,
			const PylithScalar normalTraction,
			PylithScalar* const stateVars,
			const int numStateVars,
			const PylithScalar* properties,
			const int numProperties);

  // PRIVATE METHODS ////////////////////////////////////////////////////
private :

  /// Setup fields for physical properties and state variables.
  void _setupPropsStateVars(void);

  // PROTECTED MEMBERS //////////////////////////////////////////////////
protected :

  PylithScalar _dt; ///< Current time step

  spatialdata::units::Nondimensional* _normalizer; ///< Nondimensionalizer
  
  /// Property and state variable metadata.
  const pylith::materials::Metadata _metadata;

  // PRIVATE MEMBERS ////////////////////////////////////////////////////
private :

  std::string _label; ///< Label of friction model.

  /// Database of parameters for physical properties of friction model.
  spatialdata::spatialdb::SpatialDB* _dbProperties;

  /// Database of initial state variables for the friction model.
  spatialdata::spatialdb::SpatialDB* _dbInitialState;

  /// Field containing physical properties and state variables of
  /// friction model.
  topology::Fields* _fieldsPropsStateVars;

  /// Buffer for properties and state variables at vertex.
  scalar_array _propsStateVarsVertex;

  int _propsFiberDim; ///< Number of properties per point.
  int _varsFiberDim; ///< Number of state variables per point.

  // NOT IMPLEMENTED ////////////////////////////////////////////////////
private :

  FrictionModel(const FrictionModel&); ///< Not implemented.
  const FrictionModel& operator=(const FrictionModel&); ///< Not implemented

}; // class FrictionModel

#include "FrictionModel.icc" // inline methods

#endif // pylith_friction_frictionmodel_hh


// End of file 
