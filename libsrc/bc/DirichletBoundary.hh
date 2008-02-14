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

/** @file libsrc/bc/DirichletBoundary.hh
 *
 * @brief C++ implementation of Dirichlet (prescribed values at
 * degrees of freedom) boundary conditions with points on a boundary.
 */

#if !defined(pylith_bc_dirichletboundary_hh)
#define pylith_bc_dirichletboundary_hh

#include "BoundaryCondition.hh" // ISA BoundaryCondition
#include "pylith/feassemble/Constraint.hh" // ISA Constraint

#include "pylith/utils/array.hh" // USES std::vector, double_array, int_array

/// Namespace for pylith package
namespace pylith {
  namespace bc {
    class DirichletBoundary;
    class TestDirichletBoundary; // unit testing
  } // bc
} // pylith


/// C++ implementation of DirichletBoundary boundary conditions.
class pylith::bc::DirichletBoundary : public BoundaryCondition, 
				    public feassemble::Constraint
{ // class DirichletBoundary
  friend class TestDirichletBoundary; // unit testing

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Default constructor.
  DirichletBoundary(void);

  /// Destructor.
  ~DirichletBoundary(void);

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

  /** Get data mesh.
   *
   * @return Boundary mesh.
   */
  const ALE::Obj<Mesh>& dataMesh(void) const;

  /** Get vertex field of BC initial or rate of change of values.
   *
   * @param name Name of field {'initial', 'rate-of-change'}
   */
  const ALE::Obj<real_section_type>&
  getVertexField(const char* name);

  // NOT IMPLEMENTED ////////////////////////////////////////////////////
private :

  /// Not implemented
  DirichletBoundary(const DirichletBoundary& m);

  /// Not implemented
  const DirichletBoundary& operator=(const DirichletBoundary& m);

  // PRIVATE MEMBERS ////////////////////////////////////////////////////
private :

  double _tRef; /// Time when rate of change for values begins

  /// Initial values and rate of change of values at DOF.
  ALE::Obj<real_section_type> _values;
  ALE::Obj<real_section_type> _buffer; ///< Buffer for values.

  ALE::Obj<Mesh> _boundaryMesh; ///< Boundary mesh.
  int_array _fixedDOF; ///< Indices of fixed degrees of freedom

  /// Offset in list of fixed DOF at point to get to fixed DOF
  /// associated with this DirichletBoundary boundary condition.
  ALE::Obj<int_section_type> _offsetLocal;

  /// Spatial database with parameters for rate of change values.
  spatialdata::spatialdb::SpatialDB* _dbRate;


}; // class DirichletBoundary

#include "DirichletBoundary.icc" // inline methods

#endif // pylith_bc_dirichletboundary_hh


// End of file 
