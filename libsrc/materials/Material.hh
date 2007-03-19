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
 * Interface definition for material.
 */

#if !defined(pylith_materials_material_hh)
#define pylith_materials_material_hh

#include <string> // HASA std::string

/// Namespace for pylith package
namespace pylith {
  namespace materials {
    class Material;
    class TestMaterial; // unit testing
  } // materials

  namespace feassemble {
    class ParameterManager; // HOLDSA ParameterManager
    class Quadrature; // USES Quadrature
  } // feassemble
} // pylith

/// Namespace for spatialdata package
namespace spatialdata {
  namespace spatialdb {
    class SpatialDB; // forward declaration
  } // spatialdb
  namespace geocoords {
    class CoordSys;
  } // geocoords
} // spatialdata

/// Namespace for spatialdata package
namespace ALE {
  class Mesh;
  template<class T> class Obj;
} // ALE

/// C++ abstract base class for Material object.
class pylith::materials::Material
{ // class Material
  friend class TestMaterial; // unit testing

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Default constructor.
  Material(void);

  /// Destructor.
  virtual
  ~Material(void);

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

  /** Initialize material by getting physical property parameters from
   * database.
   *
   * @param mesh PETSc mesh
   * @param cs Coordinate system associated with mesh
   * @param quadrature Quadrature for finite-element integration
   */
  void initialize(const ALE::Obj<ALE::Mesh>& mesh,
		  const spatialdata::geocoords::CoordSys* cs,
		  pylith::feassemble::Quadrature* quadrature);
  
  // PROTECTED METHODS //////////////////////////////////////////////////
protected :

  /** Copy constructor.
   *
   * @param m Material to copy
   */
  Material(const Material& m);

  /** Get names of values expected to be in database of parameters for
   *  physical properties.
   *
   * @returns Names of values
   */
  virtual
  const char** _dbValues(void) const = 0;

  /** Get number of values expected to be in database of parameters for
   *  physical properties.
   *
   * @returns Number of values
   */
  virtual
  int _numDBValues(void) const = 0;

  /** Get names of parameters for physical properties.
   *
   * @returns Names of parameters
   */
  virtual
  const char** _parameterNames(void) const = 0;

  /** Get number of parameters for physical properties.
   *
   * @returns Number of parameters
   */
  virtual
  int _numParameters(void) const = 0;

  /** Compute parameters from values in spatial database.
   *
   * @param paramVals Array of parameters
   * @param numParams Number of parameters
   * @param dbValues Array of database values
   * @param numValues Number of database values
   */
  virtual
  void _dbToParameters(double* paramVals,
		       const int numParams,
		       const double* dbValues,
		       const int numValues) const = 0;

  // NOT IMPLEMENTED ////////////////////////////////////////////////////
private :

  /// Not implemented
  const Material& operator=(const Material& m);

  // PROTECTED MEMBERS //////////////////////////////////////////////////
protected :

  ///< Manager of parameters for physical properties of material
  pylith::feassemble::ParameterManager* _parameters;

  // PRIVATE MEMBERS ////////////////////////////////////////////////////
private :

  /// Database of parameters for physical properties of material
  spatialdata::spatialdb::SpatialDB* _db;

  int _id; ///< Material identifier
  std::string _label; ///< Label of material

}; // class Material

#include "Material.icc" // inline methods

#endif // pylith_materials_material_hh


// End of file 
