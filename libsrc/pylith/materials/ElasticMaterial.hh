// -*- C++ -*-
//
// ----------------------------------------------------------------------
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University of Chicago
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2015 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

/** @file libsrc/materials/ElasticMaterial.hh
 *
 * @brief Interface definition for linear and nonlinear elastic
 * materials.
 */

#if !defined(pylith_materials_elasticmaterial_hh)
#define pylith_materials_elasticmaterial_hh

// Include directives ---------------------------------------------------
#include "Material.hh" // ISA Material

// ElasticMaterial ------------------------------------------------------
/** @brief Interface definition for linear and nonlinear elastic
 *  materials.
 */
class pylith::materials::ElasticMaterial : public Material
{ // class ElasticMaterial
  friend class TestElasticMaterial; ///< unit testing

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /** Default constructor.
   *
   * @param dimension Spatial dimension associated with material.
   * @param tensorSize Array of names of database values for material.
   * @param numElasticConsts Number of elastic constants.
   * @param metadata Metadata for physical properties and state variables.
   */
  ElasticMaterial(const int dimension,
		  const int tensorSize,
		  const int numElasticConsts,
		  const Metadata& metadata);

  /// Destructor.
  virtual
  ~ElasticMaterial(void);

  /// Deallocate PETSc and local data structures.
  virtual
  void deallocate(void);
  
  /** Set database for initial stress state.
   *
   * @param db Spatial database.
   */
  void dbInitialStress(spatialdata::spatialdb::SpatialDB* db);

  /** Set database for initial strain state.
   *
   * @param db Spatial database.
   */
  void dbInitialStrain(spatialdata::spatialdb::SpatialDB* db);

  /** Initialize material by getting physical property parameters from
   * database.
   *
   * @param mesh Finite-element mesh.
   * @param quadrature Quadrature for finite-element integration
   */
  void initialize(const topology::Mesh& mesh,
		  feassemble::Quadrature* quadrature);
  
  /// Create visitors for properties and state variables.
  void createPropsAndVarsVisitors(void);

  /// Destroy visitors for properties and state variables.
  void destroyPropsAndVarsVisitors(void);

  /** Retrieve parameters for physical properties and state variables
   * for cell.
   *
   * @param cell Finite-element cell
   */
  void retrievePropsAndVars(const int cell);

  /** Compute density for cell at quadrature points.
   *
   * @pre Must call retrievePropsAndVars for cell before calling
   * calcDensity().
   *
   * @returns Array of density values at cell's quadrature points.
   */
  const scalar_array& calcDensity(void);
  
  /** Get stress tensor at quadrature points. If the state variables
   * are from the previous time step, then the computeStateVars flag
   * should be set to true so that the state variables are updated
   * (but not stored) when computing the stresses.
   *
   * @pre Must call retrievePropsAndVars for cell before calling
   * calcStress().
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
  const scalar_array&
  calcStress(const scalar_array& totalStrain,
	     const bool computeStateVars =false);

  /** Compute derivative of elasticity matrix for cell at quadrature points.
   *
   * @pre Must call retrievePropsAndVars for cell before calling
   * calcDerivElastic().
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
  const scalar_array&
  calcDerivElastic(const scalar_array& totalStrain);

  /** Update state variables (for next time step).
   *
   * @param totalStrain Total strain tensor at quadrature points
   *    [numQuadPts][tensorSize]
   * @param cell Finite element cell
   */
  void updateStateVars(const scalar_array& totalStrain,
		       const int cell);

  /** Get flag indicating whether material implements an empty
   * _updateProperties() method.
   *
   * @returns False if _updateProperties() is empty, true otherwise.
   */
  bool hasStateVars(void) const;

  /** Get stable time step for implicit time integration.
   *
   * Default is MAXFLOAT (or 1.0e+30 if MAXFLOAT is not defined in math.h).
   *
   * @param mesh Finite-element mesh.
   * @param field Field for storing min stable time step for each cell.
   *
   * @returns Time step
   */
  virtual
  PylithScalar stableTimeStepImplicit(const topology::Mesh& mesh,
				      topology::Field* field =0);

  /** Get stable time step for explicit time integration.
   *
   * Default is MAXFLOAT (or 1.0e+30 if MAXFLOAT is not defined in math.h).
   *
   * @param mesh Finite-element mesh.
   * @param quadrature Quadrature for finite-element integration
   * @param field Field for storing min stable time step for each cell.
   *
   * @returns Time step
   */
  virtual
  PylithScalar stableTimeStepExplicit(const topology::Mesh& mesh,
				      feassemble::Quadrature* quadrature,
				      topology::Field* field =0);

  /** Set whether elastic or inelastic constitutive relations are used.
   *
   * @param flag True to use elastic, false to use inelastic.
   */
  virtual
  void useElasticBehavior(const bool flag);

  /** Get initial stress/strain fields.
   *
   * @returns Initial stress field.
   */
  const topology::Fields* initialFields(void) const;

  // PROTECTED METHODS //////////////////////////////////////////////////
protected :

  /// These methods must be implemented by every elasticity
  /// constitutive model.

  /** Compute density from properties.
   *
   * @param density Array for density.
   * @param properties Properties at location.
   * @param numProperties Number of properties.
   */
  virtual
  void _calcDensity(PylithScalar* const density,
		    const PylithScalar* properties,
		    const int numProperties,
		    const PylithScalar* stateVars,
		    const int numStateVars) = 0;

  /** Compute stress tensor from properties and state variables. If
   * the state variables are from the previous time step, then the
   * computeStateVars flag should be set to true so that the state
   * variables are updated (but not stored) when computing the
   * stresses.
   *
   * @param stress Array for stress tensor.
   * @param stressSize Size of stress tensor.
   * @param properties Properties at location.
   * @param numProperties Number of properties.
   * @param stateVars State variables at location.
   * @param numStateVars Number of state variables.
   * @param totalStrain Total strain at location.
   * @param strainSize Size of strain tensor.
   * @param initialStress Initial stress tensor at location.
   * @param initialStressSize Size of initial stress array.
   * @param initialStrain Initial strain tensor at location.
   * @param initialStrainSize Size of initial strain array.
   * @param computeStateVars Flag indicating to compute updated state variables.
   */
  virtual
  void _calcStress(PylithScalar* const stress,
		   const int stressSize,
		   const PylithScalar* properties,
		   const int numProperties,
		   const PylithScalar* stateVars,
		   const int numStateVars,
		   const PylithScalar* totalStrain,
		   const int strainSize,
		   const PylithScalar* initialStress,
		   const int initialStressSize,
		   const PylithScalar* initialStrain,
		   const int initialStrainSize,
		   const bool computeStateVars) = 0;

  /** Compute derivatives of elasticity matrix from properties.
   *
   * @param elasticConsts Array for elastic constants.
   * @param numElasticConsts Number of elastic constants.
   * @param properties Properties at location.
   * @param numProperties Number of properties.
   * @param stateVars State variables at location.
   * @param numStateVars Number of state variables.
   * @param totalStrain Total strain at location.
   * @param strainSize Size of strain tensor.
   * @param initialStress Initial stress tensor at location.
   * @param initialStressSize Size of initial stress array.
   * @param initialStrain Initial strain tensor at location.
   * @param initialStrainSize Size of initial strain array.
   */
  virtual
  void _calcElasticConsts(PylithScalar* const elasticConsts,
			  const int numElasticConsts,
			  const PylithScalar* properties,
			  const int numProperties,
			  const PylithScalar* stateVars,
			  const int numStateVars,
			  const PylithScalar* totalStrain,
			  const int strainSize,
			  const PylithScalar* initialStress,
			  const int initialStressSize,
			  const PylithScalar* initialStrain,
			  const int initialStrainSize) = 0;

  /** Update state variables (for next time step).
   *
   * @param stateVars State variables at location.
   * @param numStateVars Number of state variables.
   * @param properties Properties at location.
   * @param numProperties Number of properties.
   * @param totalStrain Total strain at location.
   * @param strainSize Size of strain tensor.
   * @param initialStress Initial stress tensor at location.
   * @param initialStressSize Size of initial stress array.
   * @param initialStrain Initial strain tensor at location.
   * @param initialStrainSize Size of initial strain array.
   */
  virtual
  void _updateStateVars(PylithScalar* const stateVars,
			const int numStateVars,
			const PylithScalar* properties,
			const int numProperties,
			const PylithScalar* totalStrain,
			const int strainSize,
			const PylithScalar* initialStress,
			const int initialStressSize,
			const PylithScalar* initialStrain,
			const int initialStrainSize);

  /** Get stable time step for implicit time integration.
   *
   * @param properties Properties at location.
   * @param numProperties Number of properties.
   * @param stateVars State variables at location.
   * @param numStateVars Number of state variables.
   *
   * @returns Time step
   */
  virtual
  PylithScalar _stableTimeStepImplicit(const PylithScalar* properties,
				       const int numProperties,
				       const PylithScalar* stateVars,
				       const int numStateVars) const = 0;

  /** Get stable time step for explicit time integration.
   *
   * @param properties Properties at location.
   * @param numProperties Number of properties.
   * @param stateVars State variables at location.
   * @param numStateVars Number of state variables.
   * @param minCellWidth Minimum width across cell.
   *
   * @returns Time step
   */
  virtual
  PylithScalar _stableTimeStepExplicit(const PylithScalar* properties,
				       const int numProperties,
				       const PylithScalar* stateVars,
				       const int numStateVars,
				       const double minCellWidth) const = 0;
  
  // PROTECTED METHODS //////////////////////////////////////////////////
protected :

  /** Get stable time step for implicit time integration for a
   * material where the stable time step is infinite.
   *
   * @param mesh Finite-element mesh.
   * @param field Field for storing min stable time step for each cell.
   *
   * @returns PYLITH::MAX_SCALAR;
   */
  PylithScalar _stableTimeStepImplicitMax(const topology::Mesh& mesh,
					  topology::Field* field =0);

  /** Compute 2D deviatoric stress/strain from vector and mean value.
   *
   * @param deviatoric Array for deviatoric tensor.
   * @param vec Input tensor (as vector).
   * @param vecMean Tensor trace divided by spatial_dimension.
   */
  static
  void calcDeviatoric2D(PylithScalar* const deviatoric,
			const PylithScalar* vec,
			const PylithScalar vecMean);

  /** Compute 2D deviatoric stress/strain from vector (length 4) and mean value.
   *
   * @param deviatoric Array for deviatoric tensor.
   * @param vec Input tensor (as vector).
   * @param vecMean Tensor trace divided by spatial_dimension.
   */
  static
  void calcDeviatoric2DPS(PylithScalar* const deviatoric,
			  const PylithScalar* vec,
			  const PylithScalar vecMean);
  
  /** Compute 3D deviatoric stress/strain from vector and mean value.
   *
   * @param deviatoric Array for deviatoric tensor.
   * @param vec Input tensor (as vector).
   * @param vecMean Tensor trace divided by spatial_dimension.
   */
  static
  void calcDeviatoric3D(PylithScalar* const deviatoric,
			const PylithScalar* vec,
			const PylithScalar vecMean);
  
  /** Compute 2D scalar product of two tensors represented as vectors.
   *
   * @param tensor1 First tensor.
   * @param tensor2 Second tensor.
   */
  static
  PylithScalar scalarProduct2D(const PylithScalar* tensor1,
			       const PylithScalar* tensor2);
  
  /** Compute 2D scalar product of two tensors represented as vectors of
   * length 4.
   *
   * @param tensor1 First tensor.
   * @param tensor2 Second tensor.
   */
  static
  PylithScalar scalarProduct2DPS(const PylithScalar* tensor1,
				 const PylithScalar* tensor2);
  
  /** Compute 3D scalar product of two tensors represented as vectors.
   *
   * @param tensor1 First tensor.
   * @param tensor2 Second tensor.
   */
  static
  PylithScalar scalarProduct3D(const PylithScalar* tensor1,
			       const PylithScalar* tensor2);
  
  // PRIVATE METHODS ////////////////////////////////////////////////////
private :

  /** Allocate cell arrays.
   *
   * @param numQuadPts Number of quadrature points.
   */
  void _allocateCellArrays(void);

  /** Initialize initial stress field.
   *
   * @param mesh Finite-element mesh.
   * @param quadrature Quadrature for finite-element integration
   */
  void _initializeInitialStress(const topology::Mesh& mesh,
				feassemble::Quadrature* quadrature);

  /** Initialize initial strain field.
   *
   * @param mesh Finite-element mesh.
   * @param quadrature Quadrature for finite-element integration
   */
  void _initializeInitialStrain(const topology::Mesh& mesh,
				feassemble::Quadrature* quadrature);

  // PRIVATE MEMBERS ////////////////////////////////////////////////////
private :

  /// Database for initial stress tensor;
  spatialdata::spatialdb::SpatialDB* _dbInitialStress;

  /// Database for initial strain tensor;
  spatialdata::spatialdb::SpatialDB* _dbInitialStrain;

  /// Initial stress/strain fields.
  topology::Fields* _initialFields;
  
  /** Properties at quadrature points for current cell.
   *
   * size = numQuadPts * numPropsQuadPt
   * index = iQuadPt * numPropsQuadPt + iPropQuadPt
   */
  scalar_array _propertiesCell;

  /** State variables at quadrature points for current cell.
   *
   * size = numQuadPts * numVarsQuadPt
   * index = iQuadPt * numVarsQuadPt + iStateVar
   */
  scalar_array _stateVarsCell;

  /** Initial stress state for current cell.
   *
   * size = numQuadPts * tensorSize
   * index = iQuadPt * tensorSize + iComponent
   */
  scalar_array _initialStressCell;

  /** Initial strain state for current cell.
   *
   * size = numQuadPts * tensorSize
   * index = iQuadPt * tensorSize + iComponent
   */
  scalar_array _initialStrainCell;

  /** Density value at quadrature points for current cell.
   *
   * size = numQuadPts
   * index = iQuadPt
   */
  scalar_array _densityCell;

  /** Stress tensor at quadrature points for current cell.
   *
   * size = numQuadPts * tensorSize
   * index = iQuadPt * tensorSize + iStress
   */
  scalar_array _stressCell;

  /** Elasticity matrix at quadrature points for current cell.
   *
   * size = numQuadPts * numElasticConsts
   * index = iQuadPt * numElasticConsts + iConstant
   */
  scalar_array _elasticConstsCell;

  int _numQuadPts; ///< Number of quadrature points
  const int _numElasticConsts; ///< Number of elastic constants.

  pylith::topology::VecVisitorMesh* _propertiesVisitor; ///< Visitor for properties field.
  pylith::topology::VecVisitorMesh* _stateVarsVisitor; ///< Visitor for stateVars field.
  pylith::topology::VecVisitorMesh* _stressVisitor; ///< Visitor for initial stress field.
  pylith::topology::VecVisitorMesh* _strainVisitor; ///< Visitor for initial strain field.

  // NOT IMPLEMENTED ////////////////////////////////////////////////////
private :

  ElasticMaterial(const ElasticMaterial&); ///< Not implemented.
  const ElasticMaterial& operator=(const ElasticMaterial&); ///< Not implemented

}; // class ElasticMaterial

#include "ElasticMaterial.icc" // inline methods

#endif // pylith_materials_elasticmaterial_hh


// End of file 
