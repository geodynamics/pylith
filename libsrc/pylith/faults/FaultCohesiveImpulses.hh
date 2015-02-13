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

/** @file libsrc/faults/FaultCohesiveImpulses.hh
 *
 * @brief C++ implementation for Green's functions impulses
 * implemented with cohesive elements.
 */

#if !defined(pylith_faults_faultcohesiveimpulses_hh)
#define pylith_faults_faultcohesiveimpulses_hh

// Include directives ---------------------------------------------------
#include "FaultCohesiveLagrange.hh" // ISA FaultCohesive

#include <map> // HASA std::map

// FaultCohesiveImpulses -----------------------------------------------------
/**
 * @brief C++ implementation for Green's functions impulses
 * implemented with cohesive elements.
 */
class pylith::faults::FaultCohesiveImpulses : public FaultCohesiveLagrange
{ // class FaultCohesiveImpulses
  friend class TestFaultCohesiveImpulses; // unit testing

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Default constructor.
  FaultCohesiveImpulses(void);

  /// Destructor.
  ~FaultCohesiveImpulses(void);

  /// Deallocate PETSc and local data structures.
  void deallocate(void);
  
  /** Sets the spatial database for amplitudes of the impulses.
   *
   * @param db spatial database for amplitudes of impulses.
   */
  void dbImpulseAmp(spatialdata::spatialdb::SpatialDB* db);
  
  /** Set indices of fault degrees of freedom associated with
   * impulses.
   *
   * @param flags Array of indices for degrees of freedom.
   * @param size Size of array
   */
  void impulseDOF(const int* flags,
		  const int size);  

  /** Set threshold for nonzero impulse amplitude.
   *
   * @param value Threshold for detecting nonzero amplitude.
   */
  void threshold(const PylithScalar value);

  /** Get number of impulses.
   *
   * Multiply by number of components to get total number of impulses.
   *
   * @returns Number of points with impulses.
   */
  int numImpulses(void) const;

  /** Get number of components for impulses at each point.
   *
   * Multiply by number of components to get total number of impulses.
   *
   * @returns Number of points with impulses.
   */
  int numComponents(void) const;

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

  /** Integrate contributions to residual term (r) for operator that
   * do not require assembly across cells, vertices, or processors.
   *
   * @param residual Field containing values for residual
   * @param t Current time
   * @param fields Solution fields
   */
  void integrateResidual(const topology::Field& residual,
			 const PylithScalar t,
			 topology::SolutionFields* const fields);

  /** Get vertex field associated with integrator.
   *
   * @param name Name of cell field.
   * @param fields Solution fields.
   * @returns Vertex field.
   */
  const topology::Field& vertexField(const char* name,
				     const topology::SolutionFields* fields =0);

  // PRIVATE METHODS ////////////////////////////////////////////////////
private :

  /// Setup impulses.
  void _setupImpulses(void);

  /** Setup impulse order.
   *
   * @param pointOrder Map from point to impulse (local) index.
   */
  void _setupImpulseOrder(const std::map<int, int>& pointOrder);

  /** Set relative displacemet associated with impulse.
   *
   * @param dispRel Relative displacement field.
   * @parm impulse Index of impulse.
   */
  void _setRelativeDisp(const topology::Field& dispRel,
			const int impulse);

  // PRIVATE TYPEDEFS ///////////////////////////////////////////////////
private :

  struct ImpulseInfoStruct {
    int indexCohesive;
    int indexDOF;
  }; // ImpulseInfoStruct
  typedef std::map<int, ImpulseInfoStruct> srcs_type;

  // PRIVATE MEMBERS ////////////////////////////////////////////////////
private :

  /// Threshold for nonzero impulse amplitude.
  PylithScalar _threshold;

  /// Database for amplitudes of impulses.
  spatialdata::spatialdb::SpatialDB* _dbImpulseAmp;

  /// Map from impulse index to corresponding point
  srcs_type _impulsePoints;

  int_array _impulseDOF; ///< Degrees of freedom associated with impulses.

  // NOT IMPLEMENTED ////////////////////////////////////////////////////
private :

  /// Not implemented
  FaultCohesiveImpulses(const FaultCohesiveImpulses&);

  /// Not implemented
  const FaultCohesiveImpulses& operator=(const FaultCohesiveImpulses&);

}; // class FaultCohesiveImpulses

#endif // pylith_faults_faultcohesiveimpulses_hh


// End of file 
