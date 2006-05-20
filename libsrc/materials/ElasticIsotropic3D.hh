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

/** @file libsrc/materials/ElasticIsotropic3D.h
 *
 * @brief C++ ElasticIsotropic3D object
 *
 * 3-D, isotropic, linear elastic material.
 */

#if !defined(pylith_materials_elasticisotropic3d_hh)
#define pylith_materials_elasticisotropic3d_hh

/// Namespace for pylith package
namespace pylith {
  namespace materials {
    class ElasticMaterial;
    class ElasticIsotropic3D;
    class TestElasticIsotropic3D; // unit testing
  } // materials
} // pylith

/// 3-D, isotropic, linear elastic material.
class pylith::materials::ElasticIsotropic3D : 
  public pylith::materials::ElasticMaterial
{ // class ElasticIsotropic3D
  friend class TestElasticIsotropic3D; // unit testing

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Default constructor
  ElasticIsotropic3D(void);

  /// Destructor
  ~ElasticIsotropic3D(void);

  /** Copy constructor
   *
   * @param m Material to copy
   */
  ElasticIsotropic3D(const ElasticIsotropic3D& m);

  /** Create a pointer to a copy of this.
   *
   * @returns Pointer to copy
   */
  ElasticMaterial* clone(void) const;

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
  void calcInertia(double** ppInertia,
		   int* pSize,
		   const double* pParams,
		   const int npts,
		   const void* pState) const;
  
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
  void calcElasticityConsts(double** pElasticityC,
			    int* pSize,
			    const double* pParams,
			    const int npts,
			    const void* pState) const;

  // PROTECTED METHODS //////////////////////////////////////////////////
protected :

  /** Get names of parameters for material.
   *
   * @returns Array of names of parameters
   */
  const char** namesParams(void) const;
  
  /** Get number of parameters for material.
   *
   * @returns Number of parameters
   */
  int numParams(void) const;
  
  // PRIVATE METHODS ////////////////////////////////////////////////////
private :

  /// Not implemented
  const ElasticIsotropic3D& operator=(const ElasticIsotropic3D& m);

}; // class ElasticIsotropic3D

#include "ElasticIsotropic3D.icc" // inline methods

#endif // pylith_materials_elasticisotropic3d_hh

// version
// $Id$

// End of file 
