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

  /** Set database for material parameters.
   *
   * @param pDB Pointer to database.
   */
  void parametersDB(const spatialdata::spatialdb::SpatialDB* pDB);

  /// Initialize material.
  void initialize(void);

  /// Open database.
  void openDB(void);

  /// Close database.
  void closeDB(void);

  /** Get inertia at points.
   *
   * The values are returned through the parameters.
   *
   * Index into array of densities:
   * index = ipt
   *
   * @param ppInertia Array of inertia values.
   * @param pSize Size of mass densities array.
   * @param pPts Array of coordinates for points [npts x 3].
   * @param npts Number of points.
   * @param pCS Coordinate system associated with points.
   * @param pState Pointer to system state at points
   */
  void inertia(double** pDensity,
	       int* pSize,
	       const double* pPts,
	       const int npts,
	       const spatialdata::geocoords::CoordSys* pCS,
	       const void* pState) const;
  
  /** Get elasticity constants at points.
   *
   * The values are returned through the parameters and are grouped by
   * point.
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
   * @param ppElasticityC Array of elasticity constants
   * @param pSize Size of elastiticy constants array
   * @param pPts Array of coordinates for points [npts x 3]
   * @param npts Number of points
   * @param pCS Coordinate system associated with points
   * @param pState Pointer to system state at points
   */
  void elasticityConsts(double** ppElasticityC,
			int* pSize,
			const double* pPts,
			const int npts,
			const spatialdata::geocoords::CoordSys* pCS,
			const void* pState) const;

  // PUBLIC MEMBERS /////////////////////////////////////////////////////
public :

  static const int NUMELASTCONSTS; ///< Number of elastic constants

  // PROTECTED METHODS //////////////////////////////////////////////////
protected :

  /** Get names of parameters for material.
   *
   * @returns Names of parameters.
   */
  virtual
  const char** namesParams(void) const = 0;

  /** Get number of parameters for material.
   *
   * @returns Number of parameters.
   */
  virtual
  int numParams(void) const = 0;

  /** Compute inertia at points using material parameters.
   *
   * The values are returned through the parameters.
   *
   * Index into array of inertia values:
   * index = iPoint
   *
   * @param ppIntertia Array of mass densities
   * @param pSize Size of mass densities array
   * @param pParams Array of material parameters [npts x numParams]
   * @param npts Number of points
   * @param pState Pointer to system state at points
   */
  virtual
  void calcInertia(double** ppInertia,
		   int* pSize,
		   const double* pParams,
		   const int npts,
		   const void* pState) const = 0;
  
  /** Compute elasticity constants at points using material parameters.
   *
   * The values are returned through the parameters and are grouped by
   * point.
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
   * @param pElasticityC Array of elasticity constants
   * @param pSize Size of elastiticy constants array
   * @param pParams Array of material parameters [npts x numParams]
   * @param npts Number of points
   * @param pState Pointer to system state at points
   */
  virtual
  void calcElasticityConsts(double** pElasticityC,
			    int* pSize,
			    const double* pParams,
			    const int npts,
			    const void* pState) const = 0;

  // PRIVATE METHODS ////////////////////////////////////////////////////
private :

  /** Get material parameters at a set of points.
   *
   * @param ppParams Pointer to array of paramters.
   * @param pSize Size of array of parameters
   * @param pPts Array of coordinates for points [npts x 3]
   * @param npts Number of points
   * @param pCS Coordinate system associated with points.
   */
  void
  _getParameters(double** ppParams,
		 int* pSize,
		 const double* pPts,
		 const int npts,
		 const spatialdata::geocoords::CoordSys* pCS) const;

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

// version
// $Id$

// End of file 
