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
 * Interface definition for linear and nonlinear elastic materials.
 */

#if !defined(pylith_materials_elasticmaterial_hh)
#define pylith_materials_elasticmaterial_hh

#include "Material.hh" // ISA Material

#include <petscmesh.h> // USES Mesh
#include "pylith/utils/arrayfwd.hh" // USES double_array
#include <vector> // USES std::vector

/// Namespace for pylith package
namespace pylith {
  namespace materials {
    class ElasticMaterial;

    class TestElasticMaterial; // unit testing
  } // materials
} // pylith

/// C++ object for material constitutive model.
class pylith::materials::ElasticMaterial : public Material
{ // class ElasticMaterial
  friend class TestElasticMaterial; ///< unit testing

  // PUBLIC TYPEDEFS ////////////////////////////////////////////////////
public :

  typedef ALE::Mesh Mesh;
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

  /** Initialize arrays holding cell data.
   *
   * @param cell Finite element cell
   * @param numQuadPts Number of quadrature points
   */
  void initCellData(const Mesh::point_type& cell,
		    const int numQuadPts);

  /** Compute density for cell at quadrature points.
   *
   * @returns Array of density values at cell's quadrature points.
   */
  const std::vector<double_array>& calcDensity(void);
  
  /** Get stress tensor at quadrature points.
   *
   * Size of array of stress tensors = [numQuadPts][tensorSize].
   *
   * Order of stresses for 3-D:
   *  0: S11,  1: S22,  2: S33,  3: S12,  4: S23,  5: S13,
   *
   * Order of stresses for 2-D:
   *  0: S11,  1: S22,  2: S12,
   *
   * Order of elasticity constants for 1-D:
   *  0: S11
   *
   * @param totalStrain Total strain tensor at quadrature points
   *    [numQuadPts][tensorSize]
   *
   * @returns Array of stresses at cell's quadrature points.
   */
  const std::vector<double_array>&
  calcStress(const std::vector<double_array>& totalStrain);

  /** Compute derivative of elasticity matrix for cell at quadrature points.
   *
   * Size of array of elasticity constants = [numQuadPts][numElasticConsts]
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
   * @param totalStrain Total strain tensor at quadrature points
   *    [numQuadPts][tensorSize]
   */
  const std::vector<double_array>&
  calcDerivElastic(const std::vector<double_array>& totalStrain);

  // PROTECTED METHODS //////////////////////////////////////////////////
protected :

  /** Copy constructor.
   *
   * @param m Material to copy
   */
  ElasticMaterial(const ElasticMaterial& m);

  /** Get number of entries in stress/strain tensors.
   *
   * @returns Size of stress/strain tensors.
   */
  virtual
  int _tensorSize(void) const = 0;

  /** Get number of entries in elastic constants array.
   *
   * @returns Number of entries in array of elastic constants.
   */
  virtual
  int _numElasticConsts(void) const = 0;

  /** Compute density from parameters.
   *
   * @param density Array for density
   * @param parameters Parameters at location
   */
  virtual
  void _calcDensity(double_array* const density,
		    const double_array& parameters) = 0;

  /** Compute stress tensor from parameters.
   *
   * @param stress Array for stress tensor
   * @param parameters Parameters at locations.
   * @param totalStrain Total strain at locations.
   */
  virtual
  void _calcStress(double_array* const stress,
		   const double_array& parameters,
		   const double_array& totalStrain) = 0;

  /** Compute derivatives of elasticity matrix from parameters.
   *
   * @param elasticConsts Array for elastic constants
   * @param parameters Parameters at locations.
   * @param totalStrain Total strain at locations.
   */
  virtual
  void _calcElasticConsts(double_array* const elasticConsts,
			  const double_array& parameters,
			  const double_array& totalStrain) = 0;

  // NOT IMPLEMENTED ////////////////////////////////////////////////////
private :

  /// Not implemented
  const ElasticMaterial& operator=(const ElasticMaterial& m);

  // PRIVATE METHODS ////////////////////////////////////////////////////
private :

  /** Get parameters for cell.
   *
   * @param cell Finite-element cell
   */
  void _getParameters(const Mesh::point_type& cell);

  // PRIVATE MEMBERS ////////////////////////////////////////////////////
private :

  int _numQuadPts; ///< Number of quadrature points

  /** Parameters at quadrature points for current cell.
   *
   * size = [numQuadPts][numParams]
   */
  std::vector<double_array> _paramsCell;

  /** Density value at quadrature points for current cell.
   *
   * size = [numQuadPts][1]
   * index = [iQuadPt][0]
   */
  std::vector<double_array> _density;

  /** Stress tensor at quadrature points for current cell.
   *
   * size = [numQuadPts][tensorSize]
   * index = [iQuadPt][iStress]
   */
  std::vector<double_array> _stress;

  /** Elasticity matrix at quadrature points for current cell.
   *
   * size = [numQuadPts][numElasticConsts]
   * index = [iQuadPt][iConstant]
   */
  std::vector<double_array> _elasticConsts;

}; // class ElasticMaterial

#endif // pylith_materials_elasticmaterial_hh


// End of file 
