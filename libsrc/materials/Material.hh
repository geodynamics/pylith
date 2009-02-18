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

/** @file libsrc/materials/Material.hh
 *
 * @brief C++ abstract base class for Material object.
 *
 * Interface definition for a material. The physical properties for
 * the material are associated with the constants in the constitutive
 * model.
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

  /** Get spatial dimension of material.
   *
   * @returns Spatial dimension.
   */
  int dimension(void) const;

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
  const std::string& label(void) const;

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

  /** Set scales used to nondimensionalize physical properties.
   *
   * @param dim Nondimensionalizer
   */
  void normalizer(const spatialdata::units::Nondimensional& dim);

  /** Initialize material by getting physical property parameters from
   * database.
   *
   * @param mesh Finite-element mesh.
   * @param quadrature Quadrature for finite-element integration
   */
  virtual
  void initialize(const topology::Mesh& mesh,
		  feassemble::Quadrature<topology::Mesh>* quadrature);
  
  /** Get flag indicating whether Jacobian matrix must be reformed for
   * current state.
   *
   * @returns True if Jacobian matrix must be reformed, false otherwise.
   */
  bool needNewJacobian(void) const;

  /// Reset flag indicating whether Jacobian matrix must be reformed for
  /// current state.
  void resetNeedNewJacobian(void);

#if 0
  /** Get physical property field. Data is returned via the
   * argument.
   *
   * @param field Proeprty field.
   * @param name Name of physical property.
   * @param mesh Finite-element mesh.
   * @param numQuadPoints Number of quadrature points.
   */
  void propertyField(ALE::Obj<real_section_type>* field,
		     const char* name,
		     const topology::Mesh& mesh,
		     const int numQuadPts) const;
#endif

  // PROTECTED METHODS //////////////////////////////////////////////////
protected :

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
		      const double_array& dbValues) const = 0;

  /** Nondimensionalize state variables.
   *
   * @param values Array of initial state values.
   * @param nvalues Number of values.
   */
  virtual
  void _nondimStateVars(double* const values,
			   const int nvalues) const = 0;
  
  /** Dimensionalize state variables.
   *
   * @param values Array of initial state values.
   * @param nvalues Number of values.
   */
  virtual
  void _dimStateVars(double* const values,
			const int nvalues) const = 0;

  // PROTECTED MEMBERS //////////////////////////////////////////////////
protected :

  double _dt; ///< Current time step

  /// Field containing physical properties of material.
  topology::Field<topology::Mesh>* _properties;

  /// Field containing the state variables for the material.
  topology::Field<topology::Mesh>* _stateVars;

  spatialdata::units::Nondimensional* _normalizer; ///< Nondimensionalizer
  
  int _numPropsQuadPt; ///< Number of properties per quad point.
  int _numVarsQuadPt; ///< Number of state variables per quad point.
  const int _dimension; ///< Spatial dimension associated with material.
  const int _tensorSize; ///< Tensor size for material.
  bool _needNewJacobian; ///< True if need to reform Jacobian, false otherwise.

  // PRIVATE MEMBERS ////////////////////////////////////////////////////
private :

  /// Database of parameters for physical properties of material.
  spatialdata::spatialdb::SpatialDB* _dbProperties;

  /// Database of initial state variables for the material.
  spatialdata::spatialdb::SpatialDB* _dbInitialState;

  int _id; ///< Material identifier.
  std::string _label; ///< Label of material.

  const Metadata& _metadata; ///< Property and state variable metadata.

  // NOT IMPLEMENTED ////////////////////////////////////////////////////
private :

  Material(const Material&); ///< Not implemented.
  const Material& operator=(const Material&); ///< Not implemented

}; // class Material

#include "Material.icc" // inline methods

#endif // pylith_materials_material_hh


// End of file 
