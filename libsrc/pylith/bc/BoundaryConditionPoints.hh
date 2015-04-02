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

/** @file libsrc/bc/BoundaryConditionPoints.hh
 *
 * @brief C++ abstract base class for BoundaryCondition object with
 * boundary condition applied at a set of points.
 *
 * Interface definition for boundary conditions applied at a set of points.
 */

#if !defined(pylith_bc_boundaryconditionpoints_hh)
#define pylith_bc_boundaryconditionpoints_hh

// Include directives ---------------------------------------------------
#include "BoundaryCondition.hh" // ISA BoundaryCondition

#include "pylith/topology/topologyfwd.hh" // HOLDSA Fields
#include "pylith/utils/array.hh" // HASA int_array

// BoundaryConditionPoints ----------------------------------------------
/** @brief C++ abstract base class for BoundaryCondition object with
 * boundary condition applied at a set of points.
 */
class pylith::bc::BoundaryConditionPoints : public BoundaryCondition
{ // class BoundaryCondition
  friend class TestBoundaryConditionPoints; // unit testing

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Default constructor.
  BoundaryConditionPoints(void);

  /// Destructor.
  virtual
  ~BoundaryConditionPoints(void);
  
  /// Deallocate PETSc and local data structures.
  virtual
  void deallocate(void);
  
  /** Get parameter fields.
   *
   * @returns Parameter fields.
   */
  const topology::Fields* parameterFields(void) const;

  // PROTECTED METHODS //////////////////////////////////////////////////
protected :

  /** Get mesh labels for points associated with boundary condition.
   *
   * @param mesh Finite-element mesh.
   */
  void _getPoints(const topology::Mesh& mesh);

  // PROTECTED MEMBERS //////////////////////////////////////////////////
protected :

  /// Parameters for boundary condition.
  topology::Fields* _parameters;

  int_array _points; ///< Points for boundary condition.

  // NOT IMPLEMENTED ////////////////////////////////////////////////////
private :

  /// Not implemented
  BoundaryConditionPoints(const BoundaryConditionPoints&);

  /// Not implemented
  const BoundaryConditionPoints& operator=(const BoundaryConditionPoints&);

}; // class BoundaryConditionPoints

#endif // pylith_bc_boundaryconditionpoints_hh


// End of file 
