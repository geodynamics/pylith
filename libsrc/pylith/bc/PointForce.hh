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

/** @file libsrc/bc/PointForce.hh
 *
 * @brief C++ implementation of point force on vertices.
 */

#if !defined(pylith_bc_pointforce_hh)
#define pylith_bc_pointforce_hh

// Include directives ---------------------------------------------------
#include "TimeDependentPoints.hh" // ISA TimeDependentPoints
#include "pylith/feassemble/Integrator.hh" // ISA Integrator

// PointForce ------------------------------------------------------
/// Point forces applied at vertices.
class pylith::bc::PointForce : public TimeDependentPoints,
			       public feassemble::Integrator
{ // class PointForce
  friend class TestPointForce; // unit testing

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Default constructor.
  PointForce(void);

  /// Destructor.
  ~PointForce(void);

  /// Deallocate PETSc and local data structures.
  void deallocate(void);
  
  /** Initialize boundary condition.
   *
   * @param mesh PETSc mesh
   * @param upDir Vertical direction (somtimes used in 3-D problems).
   */
  void initialize(const topology::Mesh& mesh,
		  const PylithScalar upDir[3]);

  /** Integrate contributions to residual term (r) for operator.
   *
   * @param residual Field containing values for residual
   * @param t Current time
   * @param fields Solution fields
   */
  void integrateResidual(const topology::Field& residual,
			 const PylithScalar t,
			 topology::SolutionFields* const fields);

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

  // NOT IMPLEMENTED ////////////////////////////////////////////////////
private :

  PointForce(const PointForce&); ///< Not implemented.
  const PointForce& operator=(const PointForce&); ///< Not implemented.

}; // class PointForce

#include "PointForce.icc" // inline methods

#endif // pylith_bc_pointforce_hh


// End of file 
