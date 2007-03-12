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
 * Interface definition for 3-D linear and nonlinear elastic materials.
 */

#if !defined(pylith_materials_elasticmaterial_hh)
#define pylith_materials_elasticmaterial_hh

#include "Material.hh" // ISA Material

#include <petscmesh.h> // USES Mesh

/// Namespace for pylith package
namespace pylith {
  namespace materials {
    class ElasticMaterial;
  } // materials
} // pylith

/// C++ object for material constitutive model.
class pylith::materials::ElasticMaterial : public Material
{ // class ElasticMaterial

  // PUBLIC TYPEDEFS ////////////////////////////////////////////////////
public :

  typedef ALE::Mesh Mesh;
  typedef Mesh::topology_type topology_type;
  typedef Mesh::real_section_type real_section_type;

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Default constructor.
  ElasticMaterial(void);

  /// Destructor.
  virtual
  ~ElasticMaterial(void);

  /** Create a pointer to a copy of this.
   *
   * @returns Pointer to copy
   */
  virtual
  ElasticMaterial* clone(void) const = 0;

  /** Get density values for cell at quadrature points.
   *
   * @returns Array of density values at cell's quadrature points.
   */
  const double* density(void) const;

  /** Get elasticity constants for cell at quadrature points.
   *
   * Index into array of elasticity constants:
   * index = iQuadPt*NUMELASTCONSTS + iConstant
   *
   * Order of elasticity constants for 3-D:
   *  0: C1111,  1: C1122,  2: C1133,  3: C1112,  4: C1123,  5: C1113,
   *             6: C2222,  7: C2233,  8: C2212,  9: C2223, 10: C2213,
   *                       11: C3333, 12: C3312, 13: C3323, 14: C3313,
   *                                  15: C1212, 16: C1223, 17: C1213,
   *                                             18: C2323, 19: C2313,
   *                                                        20: C1313
   *
   * Order of elasticity constants for 2-D:
   *  0: C1111,  1: C1122,  2: C1112,
   *             3: C2222,  4: C2212,
   *                        5: C1212
   *
   * Order of elasticity constants for 1-D:
   *  0: C1111
   *
   * @returns Array of elasticity constants at cell's quadrature points.
   */
  const double* elasticConsts(void) const;

  /** Get number of elastic constants for material.
   *
   * 1-D = 1
   * 2-D = 6
   * 3-D = 21
   *
   * @returns Number of elastic constants
   */
  virtual
  const int numElasticConsts(void) const = 0;

  /** Compute physical properties of cell at quadrature points.
   *
   * @param cell Finite-element cell
   * @param patch Finite-element patch
   * @param numQuadPts Number of quadrature points (consistency check)
   */
  void calcProperties(const topology_type::point_type& cell,
		      const topology_type::patch_type& patch,
		      const int numQuadPts);

  // PROTECTED METHODS //////////////////////////////////////////////////
protected :

  /** Copy constructor.
   *
   * @param m Material to copy
   */
  ElasticMaterial(const ElasticMaterial& m);

  /** Initialize arrays holding cell data.
   *
   * @param numQuadPts Number of quadrature points
   */
  void _initCellData(const int numQuadPts);

  /** Compute density at locations from parameters.
   *
   * Results are stored in _density.
   *
   * @param parameters Parameters at location
   * @param numParameters Number of parameters
   * @param numLocs Number of locations
   */
  virtual
  void _calcDensity(const double* parameters,
		    const int numParameters,
		    const int numLocs) = 0;

  /** Compute density at locations from parameters.
   *
   * Results are stored in _elasticConsts.
   *
   * @param parameters Parameters at location
   * @param numParameters Number of parameters
   * @param numLocs Number of locations
   */
  virtual
  void _calcElasticConsts(const double* parameters,
			  const int numParameters,
			  const int numLocs) = 0;

  // NOT IMPLEMENTED ////////////////////////////////////////////////////
private :

  /// Not implemented
  const ElasticMaterial& operator=(const ElasticMaterial& m);

  // PROTECTED MEMBERS //////////////////////////////////////////////////
protected :

  /** Density value at quadrature points for current cell.
   *
   * size = numQuadPts
   * index = iQuadPt
   */
  double* _density;

  /** Array of elasticity constants at quadrature points for current cell.
   *
   * size = numElasticConsts*numQuadPts
   * index = iQuadPt*numElasticConsts+iElasticConst
   */
  double* _elasticConsts;

}; // class ElasticMaterial

#include "ElasticMaterial.icc" // inline methods

#endif // pylith_materials_elasticmaterial_hh


// End of file 
