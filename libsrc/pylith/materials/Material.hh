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

/** @file libsrc/materials/Material.hh
 *
 * @brief C++ abstract base class for Material object.
 */

#if !defined(pylith_materials_material_hh)
#define pylith_materials_material_hh

// Include directives ---------------------------------------------------
#include "materialsfwd.hh" // forward declarations

#include "pylith/topology/topologyfwd.hh" // forward declarations
#include "pylith/feassemble/feassemblefwd.hh" // forward declarations
#include "spatialdata/spatialdb/spatialdbfwd.hh" // forward declarations
#include "spatialdata/units/unitsfwd.hh" // forward declarations

#include "Metadata.hh" // HASA Metadata

#include <string> // HASA std::string

// Material -------------------------------------------------------------
/** @brief C++ abstract base class for Material object.
 *
 * Interface definition for a material. The physical properties for
 * the material are associated with the constants in the constitutive
 * model.
 */

class pylith::materials::Material
{ // class Material
  friend class TestMaterial; // unit testing

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /** Default constructor.
   *
   * @param dimension Spatial dimension associated with material.
   * @param tensorSize Array of names of database values for material.
   * @param metadata Metadata for physical properties and state variables.
   */
  Material(const int dimension,
	   const int tensorSize,
	   const Metadata& metadata);

  /// Destructor.
  virtual
  ~Material(void);

  /// Deallocate PETSc and local data structures.
  virtual
  void deallocate(void);
  
  /** Get spatial dimension of material.
   *
   * @returns Spatial dimension.
   */
  int dimension(void) const;

  /** Set identifier of material.
   *
   * @param value Material identifier
   */
  void id(const int value);

  /** Get identifier of material.
   *
   * @returns Material identifier
   */
  int id(void) const;

  /** Set label of material.
   *
   * @param value Label of material
   */
  void label(const char* value);

  /** Get label of material.
   *
   * @returns Label of material
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

  /** Initialize material by getting physical property parameters from
   * database.
   *
   * @pre Must call Quadrature::computeGeometry() before calling
   * initialize().
   *
   * @param mesh Finite-element mesh.
   * @param quadrature Quadrature for finite-element integration
   */
  virtual
  void initialize(const topology::Mesh& mesh,
		  feassemble::Quadrature* quadrature);
  
  /** Get size of stress/strain tensor associated with material.
   *
   * @returns Size of array holding stress/strain tensor.
   */
  int tensorSize(void) const;

  /** Get flag indicating whether Jacobian matrix must be reformed for
   * current state.
   *
   * @returns True if Jacobian matrix must be reformed, false otherwise.
   */
  bool needNewJacobian(void) const;

   /** Check whether material generates a symmetric Jacobian.
    *
    * @returns True if material generates symmetric Jacobian.
    */
   bool isJacobianSymmetric(void) const;

  /// Reset flag indicating whether Jacobian matrix must be reformed for
  /// current state.
  void resetNeedNewJacobian(void);

  /** Set whether elastic or inelastic constitutive relations are used.
   *
   * @param flag True to use elastic, false to use inelastic.
   */
  virtual
  void useElasticBehavior(const bool flag);

  /** Check whether material has a field as a property.
   *
   * @param name Name of field.
   *
   * @returns True if material has field as a property, false otherwise.
   */
  bool hasProperty(const char* name);

  /** Check whether material has a field as a state variable.
   *
   * @param name Name of field.
   *
   * @returns True if material has field as a state variable, false otherwise.
   */
  bool hasStateVar(const char* name);

  /** Get physical property or state variable field. Data is returned
   * via the argument.
   *
   * @param field Field over material cells.
   * @param name Name of field to retrieve.
   */
  void getField(topology::Field *field,
		const char* name) const;

  /** Get the field with all properties.
   *
   * @returns Properties field.
   */
  const topology::Field* propertiesField() const;

  /** Get the field with all of the state variables.
   *
   * @returns State variables field.
   */
  const topology::Field* stateVarsField() const;

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
		       const scalar_array& dbValues) = 0;

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
		      const scalar_array& dbValues);

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

  // PROTECTED MEMBERS //////////////////////////////////////////////////
protected :

  PylithScalar _dt; ///< Current time step

  /// Field containing physical properties of material.
  topology::Field *_properties;

  /// Field containing the state variables for the material.
  topology::Field *_stateVars;

  spatialdata::units::Nondimensional* _normalizer; ///< Nondimensionalizer
  
  topology::StratumIS* _materialIS; ///< Index set for material cells.

  int _numPropsQuadPt; ///< Number of properties per quad point.
  int _numVarsQuadPt; ///< Number of state variables per quad point.
  const int _dimension; ///< Spatial dimension associated with material.
  const int _tensorSize; ///< Tensor size for material.
  bool _needNewJacobian; ///< True if need to reform Jacobian, false otherwise.
  bool _isJacobianSymmetric; ///< True if Jacobian is symmetric;

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

  /// Database of parameters for physical properties of material.
  spatialdata::spatialdb::SpatialDB* _dbProperties;

  /// Database of initial state variables for the material.
  spatialdata::spatialdb::SpatialDB* _dbInitialState;

  int _id; ///< Material identifier.
  std::string _label; ///< Label of material.

  const Metadata _metadata; ///< Property and state variable metadata.

  // NOT IMPLEMENTED ////////////////////////////////////////////////////
private :

  Material(const Material&); ///< Not implemented.
  const Material& operator=(const Material&); ///< Not implemented

}; // class Material

#include "Material.icc" // inline methods

#endif // pylith_materials_material_hh


// End of file 
