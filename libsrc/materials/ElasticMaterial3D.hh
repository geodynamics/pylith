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

/** @file libsrc/materials/ElasticMaterial3D.hh
 *
 * @brief C++ ElasticMaterial object
 *
 * Interface definition for 3-D linear and nonlinear elastic materials.
 */

#if !defined(pylith_materials_elasticmaterial3d_hh)
#define pylith_materials_elasticmaterial3d_hh

#include "Material.hh" // ISA Material

#include <petscmesh.h> // USES Mesh


/// Namespace for pylith package
namespace pylith {
  namespace materials {
    class ElasticMaterial3D;
  } // materials
} // pylith

/// C++ object for material constitutive model.
class pylith::materials::ElasticMaterial3D : public Material
{ // class ElasticMaterial3D

  typedef ALE::Mesh::real_section_type real_section_type;


  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Default constructor.
  ElasticMaterial3D(void);

  /// Destructor.
  virtual
  ~ElasticMaterial3D(void);

  /** Create a pointer to a copy of this.
   *
   * @returns Pointer to copy
   */
  virtual
  ElasticMaterial3D* clone(void) const = 0;

  /** Get section with density values.
   *
   * @returns Section with density values.
   */
  const ALE::Obj<real_section_type>& density(void);

  /** Get section with elasticity constants.
   *
   * Index into array of elasticity constants:
   * index = iQuadPt*NUMELASTCONSTS + iConstant
   *
   * Order of elasticity constants:
   *  0: C1111,  1: C1122,  2: C1133,  3: C1112,  4: C1123,  5: C1113,
   *             6: C2222,  7: C2233,  8: C2212,  9: C2223, 10: C2213,
   *                       11: C3333, 12: C3312, 13: C3323, 14: C3313,
   *                                  15: C1212, 16: C1223, 17: C1213,
   *                                             18: C2323, 19: C2313,
   *                                                        20: C1313
   *
   * @returns Section with elasticity constants.
   */
  const ALE::Obj<real_section_type>& elasticityConstants(void);

  // PUBLIC MEMBERS /////////////////////////////////////////////////////
public :

  static const int NUMELASTCONSTS; ///< Number of elastic constants

  // PROTECTED METHODS //////////////////////////////////////////////////
protected :

  /** Copy constructor.
   *
   * @param m Material to copy
   */
  ElasticMaterial3D(const ElasticMaterial3D& m);


  /** Compute density at location from parameters.
   *
   * @param density Pointer to density at location
   * @param parameters Parameters at location
   * @param numParameters Number of parameters
   */
  virtual
  void calcDensity(double* density,
		   const double* parameters,
		   const int numParameters) = 0;

  /** Compute density at location from parameters.
   *
   * @param elasticConsts Pointer to elastic constants at location
   * @param numConsts Number of elastic constants
   * @param parameters Parameters at location
   * @param numParameters Number of parameters
   */
  virtual
  void calcElasticConsts(double* elasticConsts,
			 const int numConstants,
			 const double* parameters,
			 const int numParameters) = 0;

  // NOT IMPLEMENTED ////////////////////////////////////////////////////
private :

  /// Not implemented
  const ElasticMaterial3D& operator=(const ElasticMaterial3D& m);

}; // class ElasticMaterial3D

#endif // pylith_materials_elasticmaterial3d_hh


// End of file 
