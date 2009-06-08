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
#include "TimeDependentPoints.hh" // ISA TimeDependentPoints
#include "pylith/feassemble/Constraint.hh" // ISA Constraint

#include "pylith/utils/array.hh" // HASA int_array

// DirichletBC ------------------------------------------------------
class pylith::bc::DirichletBC : public TimeDependentPoints, 
				public feassemble::Constraint
{ // class DirichletBC
  friend class TestDirichletBC; // unit testing

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Default constructor.
  DirichletBC(void);

  /// Destructor.
  ~DirichletBC(void);

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
  void setConstraintSizes(const topology::Field<topology::Mesh>& field);

  /** Set which degrees of freedom are constrained at points in field.
   *
   * @param field Solution field
   */
  void setConstraints(const topology::Field<topology::Mesh>& field);

  /** Set values in field.
   *
   * @param t Current time
   * @param field Solution field
   */
  void setField(const double t,
		const topology::Field<topology::Mesh>& field);

  /** Set increment in values from t0 to t1 in field.
   *
   * @param t0 Time t.
   * @param t1 Time t+dt.
   * @param field Solution field
   */
  void setFieldIncr(const double t0,
		    const double t1,
		    const topology::Field<topology::Mesh>& field);

  /** Verify configuration is acceptable.
   *
   * @param mesh Finite-element mesh
   */
  void verifyConfiguration(const topology::Mesh& mesh) const;

  // PROTECTED METHODS //////////////////////////////////////////////////
protected :

  /** Get manager of scales used to nondimensionalize problem.
   *
   * @returns Nondimensionalizer.
   */
  const spatialdata::units::Nondimensional& _getNormalizer(void) const;

  // PROTECTED MEMBERS //////////////////////////////////////////////////
protected :

  /// Offset in list of fixed DOF at point to get to fixed DOF
  /// associated with this DirichletBC boundary condition.
  int_array _offsetLocal;

  // NOT IMPLEMENTED ////////////////////////////////////////////////////
private :

  DirichletBC(const DirichletBC&); ///< Not implemented
  const DirichletBC& operator=(const DirichletBC&); ///< Not implemented

}; // class DirichletBC

#include "DirichletBC.icc" // inline methods

#endif // pylith_bc_dirichletbc_hh


// End of file 
