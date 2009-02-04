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
 * the material include BOTH parameters for the physical properties
 * AND state variables associated with the material constitutive
 * model.
 */

#if !defined(pylith_materials_material_hh)
#define pylith_materials_material_hh

#include "pylith/utils/array.hh" // USES double_array
#include <string> // HASA std::string
#include "pylith/utils/sievetypes.hh" // USES real_section_type

/// Namespace for pylith package
namespace pylith {
  namespace materials {
    class Material;
    class TestMaterial; // unit testing
  } // materials

  namespace feassemble {
    class Quadrature; // USES Quadrature
  } // feassemble
} // pylith

/// Namespace for spatialdata package
namespace spatialdata {
  namespace spatialdb {
    class SpatialDB; // forward declaration
  } // spatialdb
  namespace geocoords {
    class CoordSys; // forward declaration
  } // geocoords
  namespace units {
    class Nondimensional; // forward declaration
  } // units
} // spatialdata

/// C++ abstract base class for Material object.
class pylith::materials::Material
{ // class Material
  friend class TestMaterial; // unit testing

  // PUBLIC STRUCTURES //////////////////////////////////////////////////
public :

  struct PropMetaData {
    const char* name;
    int fiberDim;
    VectorFieldEnum fieldType;
  }; // PropMetaData

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /** Default constructor.
   *
   * @param tensorSize Array of names of database values for material.
   * @param initialStateDBValues Names of initial state database values for material.
   * @param dbValues Array of names of database values for material.
   * @param numDBValues Number of database values.
   * @param properties Array of physical property meta data.
   * @param numProperties Number of physical properties for material.
   */
  Material(const int tensorSize,
	   const char** dbValues,
	   const char** initialStateDBValues,
	   const int numDBValues,
	   const PropMetaData* properties,
	   const int numProperties);

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
  void db(spatialdata::spatialdb::SpatialDB* value);

  /** Set database for initial state variables.
   *
   * @param value Pointer to database.
   */
  void initialStateDB(spatialdata::spatialdb::SpatialDB* value);

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
  void initialize(const topology::Mesh& mesh,
		  pylith::feassemble::Quadrature* quadrature);
  
  /** Get flag indicating whether Jacobian matrix must be reformed for
   * current state.
   *
   * @returns True if Jacobian matrix must be reformed, false otherwise.
   */
  bool needNewJacobian(void) const;

  /// Reset flag indicating whether Jacobian matrix must be reformed for
  /// current state.
  void resetNeedNewJacobian(void);

  /** Get type of field associated with physical property.
   *
   * @param name Name of physical property.
   *
   * @returns Type of vector field associated with property.
   */
  VectorFieldEnum propertyFieldType(const char* name) const;

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

  // PROTECTED METHODS //////////////////////////////////////////////////
protected :

  /** Compute properties from values in spatial database.
   *
   * @param propVals Array of property values.
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

  /** Nondimensionalize initial state.
   *
   * @param values Array of initial state values.
   * @param nvalues Number of values.
   */
  virtual
  void _nondimInitState(double* const values,
			const int nvalues) const = 0;

  /** Dimensionalize initial state.
   *
   * @param values Array of initial state values.
   * @param nvalues Number of values.
   */
  virtual
  void _dimInitState(double* const values,
		     const int nvalues) const = 0;

  // NOT IMPLEMENTED ////////////////////////////////////////////////////
private :

  /// Not implemented
  Material(const Material& m);

  /// Not implemented
  const Material& operator=(const Material& m);

  // PROTECTED MEMBERS //////////////////////////////////////////////////
protected :

  double _dt; ///< Current time step

  /// Section containing physical properties of material.
  ALE::Obj<real_section_type> _properties;

  /// Section containing the initial state variables for the material.
  ALE::Obj<real_section_type> _initialState;

  spatialdata::units::Nondimensional* _normalizer; ///< Nondimensionalizer
  
  int _totalPropsQuadPt; ///< Total # of property values per quad point.
  int _dimension; ///< Spatial dimension associated with material.
  int _tensorSize; ///< Tensor size for material.
  int _initialStateSize; ///< Initial state size for material.
  bool _needNewJacobian; ///< True if need to reform Jacobian, false otherwise.

  // PRIVATE MEMBERS ////////////////////////////////////////////////////
private :

  /// Database of parameters for physical properties of material
  spatialdata::spatialdb::SpatialDB* _db;

  /// Database of initial state values for the material
  spatialdata::spatialdb::SpatialDB* _initialStateDB;

  int _id; ///< Material identifier
  std::string _label; ///< Label of material

  const PropMetaData* _propMetaData; ///< Property meta data.
  const int _numProperties; ///< Number of properties

  const char** _dbValues; ///< Names of database values
  const int _numDBValues; ///< Number of database values

  const char** _initialStateDBValues; ///< Names of initial state database values

}; // class Material

#include "Material.icc" // inline methods

#endif // pylith_materials_material_hh


// End of file 
