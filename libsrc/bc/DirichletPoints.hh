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

/** @file libsrc/bc/DirichletPoints.hh
 *
 * @brief C++ implementation of Dirichlet (prescribed values at
 * degrees of freedom) boundary conditions with a set of points.
 */

#if !defined(pylith_bc_dirichletpoints_hh)
#define pylith_bc_dirichletpoints_hh

#include "BoundaryCondition.hh" // ISA BoundaryCondition
#include "pylith/feassemble/Constraint.hh" // ISA Constraint

#include "pylith/utils/array.hh" // USES std::vector, double_array, int_array

/// Namespace for pylith package
namespace pylith {
  namespace bc {
    class DirichletPoints;
    class TestDirichletPoints; // unit testing
  } // bc
} // pylith


/// C++ implementation of DirichletPoints boundary conditions.
class pylith::bc::DirichletPoints : public BoundaryCondition, 
				    public feassemble::Constraint
{ // class DirichletPoints
  friend class TestDirichletPoints; // unit testing

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Default constructor.
  DirichletPoints(void);

  /// Destructor.
  ~DirichletPoints(void);

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
   * @param flags Indices of fixed degrees of freedom.
   */
  void fixedDOF(const int_array& flags);

  /** Set time at which rate of change begins.
   *
   * @param t Reference time.
   */
  void referenceTime(const double t);

  /** Initialize boundary condition.
   *
   * @param mesh PETSc mesh
   * @param cs Coordinate system for mesh
   */
  void initialize(const ALE::Obj<ALE::Mesh>& mesh,
		  const spatialdata::geocoords::CoordSys* cs,
		  const double_array& upDir);

  /** Set number of degrees of freedom that are constrained at points in field.
   *
   * @param field Solution field
   * @param mesh PETSc mesh
   */
  void setConstraintSizes(const ALE::Obj<real_section_type>& field,
			  const ALE::Obj<ALE::Mesh>& mesh);

  /** Set which degrees of freedom are constrained at points in field.
   *
   * @param field Solution field
   * @param mesh PETSc mesh
   */
  void setConstraints(const ALE::Obj<real_section_type>& field,
		      const ALE::Obj<ALE::Mesh>& mesh);

  /** Set values in field.
   *
   * @param t Current time
   * @param field Solution field
   * @param mesh PETSc mesh
   */
  void setField(const double t,
		const ALE::Obj<real_section_type>& field,
		const ALE::Obj<ALE::Mesh>& mesh);

  // NOT IMPLEMENTED ////////////////////////////////////////////////////
private :

  /// Not implemented
  DirichletPoints(const DirichletPoints& m);

  /// Not implemented
  const DirichletPoints& operator=(const DirichletPoints& m);

  // PRIVATE MEMBERS ////////////////////////////////////////////////////
private :

  double _tRef; /// Time when rate of change for values begins
  double_array _valuesInitial; ///< Initial values at degrees of freedom
  double_array _valuesRate; ///< Rate of change of Values at degrees of freedom

  std::vector<Mesh::point_type> _points; ///< Locations of boundary condition
  int_array _fixedDOF; ///< Indices of fixed degrees of freedom

  /// Offset in list of fixed DOF at point to get to fixed DOF
  /// associated with this DirichletPoints boundary condition.
  int_array _offsetLocal;

  /// Spatial database with parameters for rate of change values.
  spatialdata::spatialdb::SpatialDB* _dbRate;

}; // class DirichletPoints

#include "DirichletPoints.icc" // inline methods

#endif // pylith_bc_dirichletpoints_hh


// End of file 
