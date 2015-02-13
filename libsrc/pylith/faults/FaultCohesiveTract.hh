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

/** @file libsrc/faults/FaultCohesiveTract.hh
 *
 * @brief C++ implementation for a fault surface with tractions
 * applied to the fault surface using cohesive cells.
 */

#if !defined(pylith_faults_faultcohesivetract_hh)
#define pylith_faults_faultcohesivetract_hh

// Include directives ---------------------------------------------------
#include "FaultCohesive.hh" // ISA FaultCohesive

// FaultCohesiveTract -----------------------------------------------------
/** 
 * @brief C++ implementation for a fault surface with tractions
 * applied to the fault surface using cohesive cells.
 */
class pylith::faults::FaultCohesiveTract : public FaultCohesive
{ // class FaultCohesiveTract
  friend class TestFaultCohesiveTract; // unit testing

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Default constructor.
  FaultCohesiveTract(void);

  /// Destructor.
  virtual
  ~FaultCohesiveTract(void);

  /// Deallocate PETSc and local data structures.
  virtual
  void deallocate(void);

  /** Initialize fault. Determine orientation and setup boundary
   * condition parameters.
   *
   * @param mesh Finite-element mesh.
   * @param upDir Direction perpendicular to along-strike direction that is 
   *   not collinear with fault normal (usually "up" direction but could 
   *   be up-dip direction; applies to fault surfaces in 2-D and 3-D).
   */
  void initialize(const topology::Mesh& mesh,
		  const PylithScalar upDir[3]);

  /** Integrate contribution of cohesive cells to residual term.
   *
   * @param residual Field containing values for residual
   * @param t Current time
   * @param fields Solution fields
   */
  void integrateResidual(const topology::Field& residual,
			 const PylithScalar t,
			 topology::SolutionFields* const fields);

  /** Integrate contributions to Jacobian matrix (A) associated with
   * operator.
   *
   * @param jacobian Sparse matrix for Jacobian of system.
   * @param t Current time
   * @param fields Solution fields
   */
  void integrateJacobian(topology::Jacobian* jacobian,
			 const PylithScalar t,
			 topology::SolutionFields* const fields);
  
  /** Verify configuration is acceptable.
   *
   * @param mesh Finite-element mesh
   */
  void verifyConfiguration(const topology::Mesh& mesh) const;

  /** Get vertex field associated with integrator.
   *
   * @param name Name of vertex field.
   * @param fields Solution fields.
   *
   * @returns Vertex field.
   */
  const topology::Field& vertexField(const char* name,
				     const topology::SolutionFields* fields =0);
  
  /** Get cell field associated with integrator.
   *
   * @param name Name of cell field.
   * @param fields Solution fields.
   *
   * @returns Cell field.
   */
  const topology::Field& cellField(const char* name,
				   const topology::SolutionFields* fields =0);

  // NOT IMPLEMENTED ////////////////////////////////////////////////////
private :

  /// Not implemented
  FaultCohesiveTract(const FaultCohesiveTract&);

  /// Not implemented
  const FaultCohesiveTract& operator=(const FaultCohesiveTract&);

}; // class FaultCohesiveTract

#endif // pylith_faults_faultcohesivetract_hh


// End of file 
