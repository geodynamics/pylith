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
#include "pylith/utils/vectorfields.hh" // USES VectorFieldEnum

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
   * @param dbValues Array of names of database values for material.
   * @param numDBValues Number of database values.
   * @param properties Array of physical property meta data.
   * @param numProperties Number of physical properties for material.
   */
  Material(const char** dbValues,
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

  /** Initialize material by getting physical property parameters from
   * database.
   *
   * @param mesh PETSc mesh
   * @param cs Coordinate system associated with mesh
   * @param quadrature Quadrature for finite-element integration
   */
  void initialize(const ALE::Obj<Mesh>& mesh,
		  const spatialdata::geocoords::CoordSys* cs,
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

  /** Get metadata for physical property. Values are returned through
   * the arguments.
   *
   * @param space Subspace index.
   * @param fiberDim Fiber dimension.
   * @param fieldType Vector field type.
   * @param name Name of physical property.
   */
  void propertyInfo(int* space,
		    int* fiberDim,
		    VectorFieldEnum* fieldType,
		    const char* name) const;

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
  
  int _totalPropsQuadPt; ///< Total # of property values per quad point.
  int _dimension; ///< Spatial dimension associated with material.
  bool _needNewJacobian; ///< True if need to reform Jacobian, false otherwise.

  // PRIVATE MEMBERS ////////////////////////////////////////////////////
private :

  /// Database of parameters for physical properties of material
  spatialdata::spatialdb::SpatialDB* _db;

  int _id; ///< Material identifier
  std::string _label; ///< Label of material

  const PropMetaData* _propMetaData; ///< Property meta data.
  const int _numProperties; ///< Number of properties

  const char** _dbValues; ///< Names of database values
  const int _numDBValues; ///< Number of database values

}; // class Material

#include "Material.icc" // inline methods

#endif // pylith_materials_material_hh


// End of file 
