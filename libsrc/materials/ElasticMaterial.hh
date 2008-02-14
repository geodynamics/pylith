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

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /** Default constructor.
   *
   * @param tensorSize Number of entries in stress tensor.
   * @param numElasticConsts Number of elastic constants.
   * @param dbValues Array of names of database values for material.
   * @param numDBValues Number of database values.
   * @param properties Array of physical property meta data.
   * @param numProperties Number of physical properties for material.
   */
  ElasticMaterial(const int tensorSize,
		  const int numElasticConsts,
		  const char** dbValues,
		  const int numDBValues,
		  const PropMetaData* properties,
		  const int numProperties);

  /// Destructor.
  virtual
  ~ElasticMaterial(void);

  /** Get cell's property information from material's section.
   *
   * @param cell Finite element cell
   * @param numQuadPts Number of quadrature points
   */
  void getPropertiesCell(const Mesh::point_type& cell,
			 const int numQuadPts);

  /** Compute density for cell at quadrature points.
   *
   * @returns Array of density values at cell's quadrature points.
   */
  const double_array& calcDensity(void);
  
  /** Get stress tensor at quadrature points. If the state variables
   * are from the previous time step, then the computeStateVars flag
   * should be set to true so that the state variables are updated
   * (but not stored) when computing the stresses.
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
   * @param computeStateVars Flag indicating to compute updated state vars.
   *
   * @returns Array of stresses at cell's quadrature points.
   */
  const double_array&
  calcStress(const double_array& totalStrain,
	     const bool computeStateVars =false);

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
  const double_array&
  calcDerivElastic(const double_array& totalStrain);

  /** Update properties (for next time step).
   *
   * @param totalStrain Total strain tensor at quadrature points
   *    [numQuadPts][tensorSize]
   * @param cell Finite element cell
   */
  void updateProperties(const double_array& totalStrain,
			const Mesh::point_type& cell);

  /** Set whether elastic or inelastic constitutive relations are used.
   *
   * @param flag True to use elastic, false to use inelastic.
   */
  virtual
  void useElasticBehavior(const bool flag);

  /** Get flag indicating whether material implements an empty
   * _updateProperties() method.
   *
   * @returns False if _updateProperties() is empty, true otherwise.
   */
  virtual
  bool usesUpdateProperties(void) const;

  // PROTECTED METHODS //////////////////////////////////////////////////
protected :

  /** Compute density from properties.
   *
   * @param density Array for density.
   * @param properties Properties at location.
   * @param numProperties Number of properties.
   */
  virtual
  void _calcDensity(double* const density,
		    const double* properties,
		    const int numProperties) = 0;

  /** Compute stress tensor from properties. If the state variables
   * are from the previous time step, then the computeStateVars flag
   * should be set to true so that the state variables are updated
   * (but not stored) when computing the stresses.
   *
   * @param stress Array for stress tensor.
   * @param stressSize Size of stress tensor.
   * @param properties Properties at location.
   * @param numProperties Number of properties.
   * @param totalStrain Total strain at location.
   * @param strainSize Size of strain tensor.
   * @param computeStateVars Flag indicating to compute updated state vars.
   */
  virtual
  void _calcStress(double* const stress,
		   const int stressSize,
		   const double* properties,
		   const int numProperties,
		   const double* totalStrain,
		   const int strainSize,
		   const bool computeStateVars) = 0;

  /** Compute derivatives of elasticity matrix from properties.
   *
   * @param elasticConsts Array for elastic constants.
   * @param numElasticConsts Number of elastic constants.
   * @param properties Properties at location.
   * @param numProperties Number of properties.
   * @param totalStrain Total strain at location.
   * @param strainSize Size of strain tensor.
   */
  virtual
  void _calcElasticConsts(double* const elasticConsts,
			  const int numElasticConsts,
			  const double* properties,
			  const int numProperties,
			  const double* totalStrain,
			  const int strainSize) = 0;

  /** Update properties (for next time step).
   *
   * @param properties Properties at location.
   * @param numProperties Number of properties.
   * @param totalStrain Total strain at location.
   * @param strainSize Size of strain tensor.
   */
  virtual
  void _updateProperties(double* const properties,
			 const int numProperties,
			 const double* totalStrain,
			 const int strainSize);

  // NOT IMPLEMENTED ////////////////////////////////////////////////////
private :

  /// Not implemented
  ElasticMaterial(const ElasticMaterial& m);

  /// Not implemented
  const ElasticMaterial& operator=(const ElasticMaterial& m);

  // PRIVATE METHODS ////////////////////////////////////////////////////
private :

  /** Get properties for cell.
   *
   * @param cell Finite-element cell
   */
  void _getProperties(const Mesh::point_type& cell);

  // PRIVATE MEMBERS ////////////////////////////////////////////////////
private :

  int _numQuadPts; ///< Number of quadrature points
  int _tensorSize; ///< Number of entries in stress tensor.
  int _numElasticConsts; ///< Number of elastic constants.

  /** Properties at quadrature points for current cell.
   *
   * size = numQuadPts*numPropsQuadPt
   * index = iQuadPt*iParam*iValue
   */
  double_array _propertiesCell;

  /** Density value at quadrature points for current cell.
   *
   * size = numQuadPts
   * index = iQuadPt
   */
  double_array _density;

  /** Stress tensor at quadrature points for current cell.
   *
   * size = numQuadPts*tensorSize
   * index = *iQuadPt*tensorSize + iStress
   */
  double_array _stress;

  /** Elasticity matrix at quadrature points for current cell.
   *
   * size = numQuadPts*numElasticConsts
   * index = iQuadPt*numElasticConsts+iConstant
   */
  double_array _elasticConsts;

}; // class ElasticMaterial

#include "ElasticMaterial.icc" // inline methods

#endif // pylith_materials_elasticmaterial_hh


// End of file 
