// -*- C++ -*-
//
// ----------------------------------------------------------------------
//
//                           Brad T. Aagaard
//                        U.S. Geological Survey
//
// {LicenseText}
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
  void timeStep(const double dt);

  /** Get current time step.
   *
   * @returns Current time step.
   */
  double timeStep(void) const;

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
   * @param area Area at vertices of subdomain.
   */
  virtual
  void initialize(const topology::SubMesh& mesh,
		  feassemble::Quadrature<topology::SubMesh>* quadrature,
		  const topology::Field<topology::SubMesh>& area);
  
  /** Check whether friction model has a field as a property.
   *
   * @param name Name of field.
   *
   * @returns True if friction model has field as a property, false
   * otherwise.
   */
  bool hasProperty(const char* name);

  /** Check whether friction model has a field as a state variable.
   *
   * @param name Name of field.
   *
   * @returns True if friction model has field as a state variable,
   * false otherwise.
   */
  bool hasStateVar(const char* name);

  /** Get physical property or state variable field. Data is returned
   * via the argument.
   *
   * @param field Field over fault interface cells.
   * @param name Name of field to retrieve.
   */
  void getField(topology::Field<topology::SubMesh> *field,
		const char* name) const;

  /** Get the field with all properties.
   *
   * @returns Properties field.
   */
  const topology::Field<topology::SubMesh>* propertiesField() const;

  /** Get the field with all of the state variables.
   *
   * @returns State variables field.
   */
  const topology::Field<topology::SubMesh>* stateVarsField() const;

  /** Retrieve parameters for physical properties and state variables
   * for vertex.
   *
   * @param vertex Finite-element vertex on friction interface.
   */
  void retrievePropsAndVars(const int vertex);

  /** Compute friction at vertex.
   *
   * @pre Must call retrievePropsAndVars for cell before calling
   * calcFriction().
   *
   * @param slip Current slip at location.
   * @param slipRate Current slip rate at location.
   * @param normalTraction Normal traction at location.
   *
   * @returns Friction (magnitude of shear traction) at vertex.
   */
  double calcFriction(const double slip,
                      const double slipRate,
                      const double normalTraction);
  
  // PROTECTED METHODS //////////////////////////////////////////////////
protected :

  /// These methods should be implemented by every constitutive model.

  /** Compute properties from values in spatial database.
   *
   * @param propValues Array of property values.
   * @param dbValues Array of database values.
   */
  virtual
  void _dbToProperties(double* const propValues,
		       const double_array& dbValues) const = 0;

  /** Nondimensionalize properties.
   *
   * @param values Array of property values.
   * @param nvalues Number of values.
   */
  virtual
  void _nondimProperties(double* const values,
			 const int nvalues) const = 0;

  /** Dimensionalize properties.
   *
   * @param values Array of property values.
   * @param nvalues Number of values.
   */
  virtual
  void _dimProperties(double* const values,
		      const int nvalues) const = 0;

  /** Compute initial state variables from values in spatial database.
   *
   * @param stateValues Array of state variable values.
   * @param dbValues Array of database values.
   */
  virtual
  void _dbToStateVars(double* const stateValues,
		      const double_array& dbValues) const;

  /** Nondimensionalize state variables.
   *
   * @param values Array of initial state values.
   * @param nvalues Number of values.
   */
  virtual
  void _nondimStateVars(double* const values,
			   const int nvalues) const;
  
  /** Dimensionalize state variables.
   *
   * @param values Array of initial state values.
   * @param nvalues Number of values.
   */
  virtual
  void _dimStateVars(double* const values,
			const int nvalues) const;

  /** Compute friction from properties and state variables.
   *
   * @param slip Current slip at location.
   * @param slipRate Current slip rate at location.
   * @param normalTraction Normal traction at location.
   * @param properties Properties at location.
   * @param numProperties Number of properties.
   * @param stateVars State variables at location.
   * @param numStateVars Number of state variables.
   */
  virtual
  double _calcFriction(const double slip,
		       const double slipRate,
		       const double normalTraction,
		       const double* properties,
		       const int numProperties,
		       const double* stateVars,
		       const int numStateVars) = 0;

  /** Update state variables (for next time step).
   *
   * @param stateVars State variables at location.
   * @param numStateVars Number of state variables.
   * @param properties Properties at location.
   * @param numProperties Number of properties.
   */
  virtual
  void _updateStateVars(const double slip,
			const double slipRate,
			double* const stateVars,
			const int numStateVars,
			const double* properties,
			const int numProperties);

  // PROTECTED MEMBERS //////////////////////////////////////////////////
protected :

  double _dt; ///< Current time step

  /// Field containing physical properties of friction model.
  topology::Field<topology::SubMesh> *_properties;

  /// Field containing the state variables for the friction model.
  topology::Field<topology::SubMesh> *_stateVars;

  spatialdata::units::Nondimensional* _normalizer; ///< Nondimensionalizer
  
  int _numPropsVertex; ///< Number of properties per vertex.
  int _numVarsVertex; ///< Number of state variables per vertex.

  // PRIVATE METHODS ////////////////////////////////////////////////////
private :
  
  /** Get indices for physical property or state variable field. Index
   * of physical property or state variable is set, unknown values are
   * -1.
   *
   * @param propertyIndex Index of field in properties array.
   * @param stateVarIndex Index of field in state variables array.
   * @param name Name of field.
   */
  void _findField(int* propertyIndex,
		  int* stateVarIndex,
		  const char* name) const;

  // PRIVATE MEMBERS ////////////////////////////////////////////////////
private :

  /// Database of parameters for physical properties of friction model.
  spatialdata::spatialdb::SpatialDB* _dbProperties;

  /// Database of initial state variables for the friction model.
  spatialdata::spatialdb::SpatialDB* _dbInitialState;

  std::string _label; ///< Label of friction model.

  /// Property and state variable metadata.
  const pylith::materials::Metadata _metadata;

  /** Properties for current vertex.
   *
   * size = numProps
   * index = iProp
   */
  double_array _propertiesVertex;

  /** State variables for current vertex.
   *
   * size = numVars
   * index = iStateVar
   */
  double_array _stateVarsVertex;

  // NOT IMPLEMENTED ////////////////////////////////////////////////////
private :

  FrictionModel(const FrictionModel&); ///< Not implemented.
  const FrictionModel& operator=(const FrictionModel&); ///< Not implemented

}; // class FrictionModel

#include "FrictionModel.icc" // inline methods

#endif // pylith_friction_frictionmodel_hh


// End of file 
