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

/** @file libsrc/bc/DirichletBC.hh
 *
 * @brief C++ implementation of Dirichlet (prescribed values at
 * degrees of freedom) boundary conditions with a set of points.
 */

#if !defined(pylith_bc_dirichletbc_hh)
#define pylith_bc_dirichletbc_hh

// Include directives ---------------------------------------------------
#include "BoundaryCondition.hh" // ISA BoundaryCondition
#include "pylith/feassemble/Constraint.hh" // ISA Constraint

#include "pylith/utils/array.hh" // HASA std::vector, double_array, int_array
#define NEWPYLITHMESH 1
#include "pylith/utils/sievetypes.hh" // HASA SieveMesh::point_type

// DirichletBC ------------------------------------------------------
class pylith::bc::DirichletBC : public BoundaryCondition, 
				public feassemble::Constraint
{ // class DirichletBC
  friend class TestDirichletBC; // unit testing

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Default constructor.
  DirichletBC(void);

  /// Destructor.
  ~DirichletBC(void);

  /** Set database for rate of change of values.
   *
   * @param db Spatial database
   */
  void dbRate(spatialdata::spatialdb::SpatialDB* const db);

  /** Set indices of fixed degrees of freedom. 
   *
   * Note: all points associated with boundary condition has same
   * degrees of freedom fixed.
   *
   * Example: [0, 1] to fix x and y degrees of freedom in Cartesian system.
   *
   * @param flags Array of indices of fixed degrees of freedom.
   * @param size Size of array
   */
  void fixedDOF(const int* flags,
		const int size);

  /** Set time at which rate of change begins.
   *
   * @param t Reference time.
   */
  void referenceTime(const double t);

  /** Initialize boundary condition.
   *
   * @param mesh PETSc mesh
   * @param upDir Vertical direction (somtimes used in 3-D problems).
   */
  void initialize(const topology::Mesh& mesh,
		  const double upDir[3]);

  /** Set number of degrees of freedom that are constrained at points in field.
   *
   * @param field Solution field
   */
  void setConstraintSizes(const topology::Field& field);

  /** Set which degrees of freedom are constrained at points in field.
   *
   * @param field Solution field
   */
  void setConstraints(const topology::Field& field);

  /** Set values in field.
   *
   * @param t Current time
   * @param field Solution field
   */
  void setField(const double t,
		const topology::Field& field);

  // PROTECTED METHODS //////////////////////////////////////////////////
protected :

  /** Get mesh labels for points associated with Dirichlet BC.
   *
   * @param mesh Finite-element mesh.
   */
  void _getPoints(const topology::Mesh& mesh);

  /// Setup initial and rate of change databases for querying.
  void _setupQueryDatabases(void);

  /** Query initial and rate of change databases for values.
   *
   * @param mesh Finite-element mesh.
   */
  void _queryDatabases(const topology::Mesh& mesh);

  // NOT IMPLEMENTED ////////////////////////////////////////////////////
private :

  /// Not implemented
  DirichletBC(const DirichletBC& m);

  /// Not implemented
  const DirichletBC& operator=(const DirichletBC& m);

  // PROTECTED MEMBERS //////////////////////////////////////////////////
protected :

  double _tRef; /// Time when rate of change for values begins.
  double_array _valuesInitial; ///< Initial values at points.
  double_array _valuesRate; ///< Rate of change of values at points.

  std::vector<SieveMesh::point_type> _points; ///< Points for BC
  int_array _fixedDOF; ///< Indices of fixed degrees of freedom

  /// Offset in list of fixed DOF at point to get to fixed DOF
  /// associated with this DirichletBC boundary condition.
  int_array _offsetLocal;

  /// Spatial database with parameters for rate of change values.
  spatialdata::spatialdb::SpatialDB* _dbRate;

}; // class DirichletBC

#include "DirichletBC.icc" // inline methods

#endif // pylith_bc_dirichletbc_hh


// End of file 
