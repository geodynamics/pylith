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
// Copyright (c) 2010 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

/** @file libsrc/faults/FaultCohesiveKin.hh
 *
 * @brief C++ implementation for a fault surface with kinematic
 * (prescribed) slip implemented with cohesive elements.
 */

#if !defined(pylith_faults_faultcohesivekin_hh)
#define pylith_faults_faultcohesivekin_hh

// Include directives ---------------------------------------------------
#include "FaultCohesiveLagrange.hh" // ISA FaultCohesive

#include <string> // HASA std::string
#include <map> // HASA std::map

// FaultCohesiveKin -----------------------------------------------------
/**
 * @brief C++ implementation for a fault surface with kinematic
 * (prescribed) slip implemented with cohesive elements.
 */
class pylith::faults::FaultCohesiveKin : public FaultCohesiveLagrange
{ // class FaultCohesiveKin
  friend class TestFaultCohesiveKin; // unit testing

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Default constructor.
  FaultCohesiveKin(void);

  /// Destructor.
  ~FaultCohesiveKin(void);

  /// Deallocate PETSc and local data structures.
  void deallocate(void);
  
  /** Set kinematic earthquake sources.
   *
   * @param names Array of kinematic earthquake source names.
   * @param numNames Number of earthquake sources.
   * @param sources Array of kinematic earthquake sources.
   * @param numSources Number of earthquake sources.
   */
  void eqsrcs(const char* const* names,
	      const int numNames,
	      EqKinSrc** sources,
	      const int numSources);

  /** Initialize fault. Determine orientation and setup boundary
   * condition parameters.
   *
   * @param mesh Finite-element mesh.
   * @param upDir Direction perpendicular to along-strike direction that is 
   *   not collinear with fault normal (usually "up" direction but could 
   *   be up-dip direction; applies to fault surfaces in 2-D and 3-D).
   */
  void initialize(const topology::Mesh& mesh,
		  const double upDir[3]);

  /** Integrate contributions to residual term (r) for operator that
   * do not require assembly across cells, vertices, or processors.
   *
   * @param residual Field containing values for residual
   * @param t Current time
   * @param fields Solution fields
   */
  void integrateResidual(const topology::Field<topology::Mesh>& residual,
			 const double t,
			 topology::SolutionFields* const fields);

  /** Get vertex field associated with integrator.
   *
   * @param name Name of cell field.
   * @param fields Solution fields.
   * @returns Vertex field.
   */
  const topology::Field<topology::SubMesh>&
  vertexField(const char* name,
	      const topology::SolutionFields* fields =0);

  /** Get cell field associated with integrator.
   *
   * @param name Name of cell field.
   * @param fields Solution fields.
   * @returns Cell field.
   */
  const topology::Field<topology::SubMesh>&
  cellField(const char* name,
	    const topology::SolutionFields* fields =0);

  // PRIVATE TYPEDEFS ///////////////////////////////////////////////////
private :

  typedef std::map<std::string, EqKinSrc*> srcs_type;

  // PRIVATE MEMBERS ////////////////////////////////////////////////////
private :

  srcs_type _eqSrcs; ///< Array of kinematic earthquake sources.

  // NOT IMPLEMENTED ////////////////////////////////////////////////////
private :

  /// Not implemented
  FaultCohesiveKin(const FaultCohesiveKin&);

  /// Not implemented
  const FaultCohesiveKin& operator=(const FaultCohesiveKin&);

}; // class FaultCohesiveKin

#endif // pylith_faults_faultcohesivekin_hh


// End of file 
