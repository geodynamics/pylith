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

  /** Compute density for cell at quadrature points.
   *
   * @param cell Finite-element cell
   * @param patch Finite-element patch
   * @param numQuadPts Number of quadrature points (consistency check)
   *
   * @returns Array of density values at cell's quadrature points.
   */
  const double* calcDensity(const Mesh::point_type& cell,
			    const int numQuadPts);
  
  /** Get stress tensor at quadrature points.
   *
   * Index into array of elasticity constants:
   * index = iQuadPt*STRESSSIZE + iStress
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
   * @param cell Finite-element cell
   * @param patch Finite-element patch
   * @param totalStrain Total strain tensor at quadrature points
   * @param numQuadPts Number of quadrature points (consistency check)
   * @param spaceDim Spatial dimension (consistency check)
   *
   * @returns Array of stresses at cell's quadrature points.
   */
  const double* calcStress(const Mesh::point_type& cell,
			   const double* totalStrain,
			   const int numQuadPts,
			   const int spaceDim);

  /** Get number of entries in stress tensor for material.
   *
   * 1-D = 1
   * 2-D = 3
   * 3-D = 6
   *
   * @returns Number of entries in stress tensor
   */
  virtual
  int stressSize(void) const = 0;

  /** Compute derivative of elasticity matrix for cell at quadrature points.
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
   * @param cell Finite-element cell
   * @param patch Finite-element patch
   * @param totalStrain Total strain tensor at quadrature points
   * @param numQuadPts Number of quadrature points (consistency check)
   * @param spaceDim Spatial dimension (consistency check)
   */
  const double* calcDerivElastic(const Mesh::point_type& cell,
				 const double* totalStrain,
				 const int numQuadPts,
				 const int spaceDim);

  /** Get number of entries in derivatives of elasticity matrix for material.
   *
   * 1-D = 1
   * 2-D = 6
   * 3-D = 21
   *
   * @returns Number of entries in derivative of elasticity matrix
   */
  virtual
  int numElasticConsts(void) const = 0;

  /** Initialize arrays holding cell data.
   *
   * @param numQuadPts Number of quadrature points
   */
  void initCellData(const int numQuadPts);

  // PROTECTED METHODS //////////////////////////////////////////////////
protected :

  /** Copy constructor.
   *
   * @param m Material to copy
   */
  ElasticMaterial(const ElasticMaterial& m);

  /** Compute density from parameters.
   *
   * @param density Array for density
   * @param size Size of array for density
   * @param parameters Parameters at location
   * @param numParameters Number of parameters
   */
  virtual
  void _calcDensity(double* const density,
		    const int size,
		    const double* parameters,
		    const int numParameters) = 0;

  /** Compute stress tensor from parameters.
   *
   * @param stress Array for stress tensor
   * @param size Size of array for stress tensor
   * @param parameters Parameters at locations.
   * @param numParameters Number of parameters.
   * @param totalStrain Total strain at locations.
   * @param spaceDim Spatial dimension for locations.
   */
  virtual
  void _calcStress(double* const stress,
		   const int size,
		   const double* parameters,
		   const int numParameters,
		   const double* totalStrain,
		   const int spaceDim) = 0;

  /** Compute derivatives of elasticity matrix from parameters.
   *
   * @param elasticConsts Array for elastic constants
   * @param size Size of array
   * @param parameters Parameters at locations.
   * @param numParameters Number of parameters.
   * @param totalStrain Total strain at locations.
   * @param spaceDim Spatial dimension for locations.
   */
  virtual
  void _calcElasticConsts(double* const elasticConsts,
			  const int size,
			  const double* parameters,
			  const int numParameters,
			  const double* totalStrain,
			  const int spaceDim) = 0;

  // NOT IMPLEMENTED ////////////////////////////////////////////////////
private :

  /// Not implemented
  const ElasticMaterial& operator=(const ElasticMaterial& m);

  // PRIVATE METHODS ////////////////////////////////////////////////////
private :

  /** Get parameters for cell.
   *
   * Parameters are returned in paramsCells.
   *
   * size = numQuadPts * numParams
   * index = iQuad*numParams + iParam
   *
   * @param paramsCell Array of parameters for cell
   * @param cell Finite-element cell
   * @param numQuadPts Number of quadrature points (consistency check)
   */
  void _getParameters(double** paramsCells,
		      const Mesh::point_type& cell,
		      const int numQuadPts);

  // PRIVATE MEMBERS ////////////////////////////////////////////////////
private :

  /** Density value at quadrature points for current cell.
   *
   * size = numQuadPts
   * index = iQuadPt
   */
  double* _density;

  /** Stress tensor at quadrature points for current cell.
   *
   * size = stressSize*numQuadPts
   * index = iQuadPt*stressSize+iStress
   */
  double* _stress;

  /** Elasticity matrix at quadrature points for current cell.
   *
   * size = numElasticConsts*numQuadPts
   * index = iQuadPt*numElasticConsts+iConstant
   */
  double* _elasticConsts;

}; // class ElasticMaterial

#endif // pylith_materials_elasticmaterial_hh


// End of file 
