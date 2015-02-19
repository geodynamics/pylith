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

/** @file libsrc/feassemble/Constraint.hh
 *
 * @brief C++ abstract base class defining interface for constraints
 * applied to finite-elements.
 */

#if !defined(pylith_feassemble_constraint_hh)
#define pylith_feassemble_constraint_hh

// Include directives ---------------------------------------------------
#include "feassemblefwd.hh"

#include "pylith/topology/topologyfwd.hh" // USES Mesh, Field
#include "spatialdata/units/unitsfwd.hh" // USES Nondimensional

// Constraint -----------------------------------------------------------
/** @brief Abstract base class for defining constraints applied to
 *  vertices of finite-elements.
 */
class pylith::feassemble::Constraint
{ // class Constraint
  friend class TestConstraint; // unit testing

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Default constructor.
  Constraint(void);

  /// Destructor.
  virtual
  ~Constraint(void);

  /// Deallocate PETSc and local data structures.
  virtual
  void deallocate(void);
  
  /** Set manager of scales used to nondimensionalize problem.
   *
   * @param dim Nondimensionalizer.
   */
  void normalizer(const spatialdata::units::Nondimensional& dim);

  /** Get number of constraints per location.
   *
   * @returns Number of constraints per location.
   */
  virtual
  int numDimConstrained(void) const = 0;

  /** Set number of degrees of freedom that are constrained at points in field.
   *
   * @param field Solution field
   */
  virtual
  void setConstraintSizes(const topology::Field& field) = 0;

  /** Set which degrees of freedom are constrained at points in field.
   *
   * @param field Solution field
   */
  virtual
  void setConstraints(const topology::Field& field) = 0;

  /** Set values in field.
   *
   * @param t Current time
   * @param field Solution field
   */
  virtual
  void setField(const PylithScalar t,
		const topology::Field& field) = 0;

  /** Set increment in values from t0 to t1 in field.
   *
   * @param t0 Time t.
   * @param t1 Time t+dt.
   * @param field Solution field
   */
  virtual
  void setFieldIncr(const PylithScalar t0,
		    const PylithScalar t1,
		    const topology::Field& field) = 0;

  // PROTECTED MEMBERS //////////////////////////////////////////////////
protected :

  spatialdata::units::Nondimensional* _normalizer; ///< Nondimensionalizer.

  // NOT IMPLEMENTED ////////////////////////////////////////////////////
private :

  /// Not implemented
  Constraint(const Constraint& m);

  /// Not implemented
  const Constraint& operator=(const Constraint& m);

}; // class Constraint

#include "Constraint.icc" // inline methods

#endif // pylith_feassemble_constraint_hh


// End of file 
