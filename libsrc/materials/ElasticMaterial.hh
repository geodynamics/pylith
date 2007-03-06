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

/** @file libsrc/materials/ElasticMaterial.hh
 *
 * @brief C++ ElasticMaterial object
 *
 * Interface definition for material constitutive model.
 */

#if !defined(pylith_materials_elasticmaterial_hh)
#define pylith_materials_elasticmaterial_hh

/// Namespace for pylith package
namespace pylith {
  namespace materials {
    class ElasticMaterial;
  } // materials
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

/// C++ object for material constitutive model.
class pylith::materials::ElasticMaterial
{ // class ElasticMaterial

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Default constructor.
  ElasticMaterial(void);

  /// Destructor.
  virtual
  ~ElasticMaterial(void);

  /** Copy constructor.
   *
   * @param m Material to copy
   */
  ElasticMaterial(const ElasticMaterial& m);

  /** Create a pointer to a copy of this.
   *
   * @returns Pointer to copy
   */
  virtual
  ElasticMaterial* clone(void) const = 0;

  /** Set database for physical property parameters.
   *
   * @param pDB Pointer to database.
   */
  void parametersDB(spatialdata::spatialdb::SpatialDB* pDB);

  /** Get physical property parameters from database.
   *
   * @param mesh PETSc mesh
   * @param label Label identifying material
   */
  void createParameters(const ALE::Obj<ALE::Mesh>& mesh);
  
  /** Get section with inertia values.
   *
   */
  void inertia(void);

  /** Get section with elasticity constants.
   *
   * Index into array of elasticity constants:
   * index = iPoint*NUMELASTCONSTS + iConstant
   *
   * Order of elasticity constants:
   *  0: C1111,  1: C1122,  2: C1133,  3: C1112,  4: C1123,  5: C1113,
   *             6: C2222,  7: C2233,  8: C2212,  9: C2223, 10: C2213,
   *                       11: C3333, 12: C3312, 13: C3323, 14: C3313,
   *                                  15: C1212, 16: C1223, 17: C1213,
   *                                             18: C2323, 19: C2313,
   *                                                        20: C1313
   *
   */
  void elasticityConstants(void);

  // PUBLIC MEMBERS /////////////////////////////////////////////////////
public :

  static const int NUMELASTCONSTS; ///< Number of elastic constants

  // NOT IMPLEMENTED ////////////////////////////////////////////////////
private :

  /// Not implemented
  const ElasticMaterial& operator=(const ElasticMaterial& m);

  // PRIVATE MEMBERS ////////////////////////////////////////////////////
private :

  /// Database of parameters for material
  spatialdata::spatialdb::SpatialDB* _pDB;

}; // class ElasticMaterial

#endif // pylith_materials_elasticmaterial_hh

// End of file 
