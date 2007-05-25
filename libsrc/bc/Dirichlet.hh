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

/** @file libsrc/bc/Dirichlet.hh
 *
 * @brief C++ implementation of Dirichlet (prescribed values at
 * degrees of freedom) boundary conditions.
 */

#if !defined(pylith_bc_dirichlet_hh)
#define pylith_bc_dirichlet_hh

#include "BoundaryCondition.hh" // ISA BoundaryCondition
#include "pylith/feassemble/Constraint.hh" // ISA Constraint

#include "pylith/utils/array.hh" // USES std::vector, double_array, int_array

/// Namespace for pylith package
namespace pylith {
  namespace bc {
    class Dirichlet;
    class TestDirichlet; // unit testing
  } // bc
} // pylith


/// C++ implementation of Dirichlet boundary conditions.
class pylith::bc::Dirichlet : public BoundaryCondition, 
			      public feassemble::Constraint
{ // class Dirichlet
  friend class TestDirichlet; // unit testing

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Default constructor.
  Dirichlet(void);

  /// Destructor.
  ~Dirichlet(void);

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

  /** Initialize boundary condition.
   *
   * @param mesh PETSc mesh
   * @param cs Coordinate system for mesh
   */
  void initialize(const ALE::Obj<ALE::Mesh>& mesh,
		  const spatialdata::geocoords::CoordSys* cs);

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
  Dirichlet(const Dirichlet& m);

  /// Not implemented
  const Dirichlet& operator=(const Dirichlet& m);

  // PRIVATE MEMBERS ////////////////////////////////////////////////////
private :

  std::vector<Mesh::point_type> _points; ///< Locations of boundary condition
  double_array _values; ///< Values at degrees of freedom
  int_array _fixedDOF; ///< Indices of fixed degrees of freedom

}; // class Dirichlet

#include "Dirichlet.icc" // inline methods

#endif // pylith_bc_dirichlet_hh


// End of file 
