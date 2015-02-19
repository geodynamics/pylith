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

/** @file libsrc/faults/FaultCohesive.hh
 *
 * @brief C++ abstract base class for a fault surface implemented with
 * cohesive elements.
 */

#if !defined(pylith_faults_faultcohesive_hh)
#define pylith_faults_faultcohesive_hh

// Include directives ---------------------------------------------------
#include "Fault.hh" // ISA Fault

#include "pylith/feassemble/Integrator.hh" // ISA Integrator

#include <map>

// FaultCohesive --------------------------------------------------------
/// Absract base class for fault surface implemented with cohesive cells.
class pylith::faults::FaultCohesive : public Fault,
				      public feassemble::Integrator
{ // class FaultCohesive
  friend class TestFaultCohesive; // unit testing

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Default constructor.
  FaultCohesive(void);

  /// Destructor.
  virtual
  ~FaultCohesive(void);

  /// Deallocate PETSc and local data structures.
  virtual
  void deallocate(void);
  
  /** Set flag for using fault mesh or group of vertices to define
   * fault surface.
   *
   * This method is part of a KLUDGE to allow creation of cohesive
   * cells in cases where domain cells have more than one face (edge
   * for 2-D problems) on the fault.
   *
   * @param flag True if using fault mesh, false if using vertices.
   */
  void useFaultMesh(const bool flag);

  /** Get the number of vertices associated with the fault (before
   * fault mesh exists).
   *
   * @param mesh PETSc mesh
   * @return Number of vertices on the fault.
   */
  int numVerticesNoMesh(const topology::Mesh& mesh) const;

  /** Adjust mesh topology for fault implementation.
   *
   * If firstFaultVertex == 0, then firstFaultVertex is set to the first point
   * not currently used in the mesh, and firstLagrangeVertex/firstFaultCell are
   * incremented with this point. These values are updated as new fault vertices
   * and cells are added.
   *
   * @param mesh PETSc mesh.
   * @param firstFaultVertex The first point eligible to become a new fault vertex
   * @param firstLagrangeVertex The first point eligible to become a new Lagrange vertex
   * @param firstFaultCell The first point eligible to become a new fault cell
   */
  void adjustTopology(topology::Mesh* const mesh,
                      int *firstFaultVertex,
                      int *firstLagrangeVertex,
                      int *firstFaultCell);

  /** Cohesive cells use Lagrange multiplier constraints?
   *
   * @returns True if implementation using Lagrange multiplier
   * constraints, false otherwise.
   */
  bool useLagrangeConstraints(void) const;

  /** Get fields associated with fault.
   *
   * @returns Fields associated with fault.
   */
  const topology::Fields* fields(void) const;

  // PROTECTED MEMBERS //////////////////////////////////////////////////
protected :

  /// Fields for fault information.
  topology::Fields* _fields;

  bool _useLagrangeConstraints; ///< True if uses Lagrange multipliers.

  /// Map label of cohesive cell to label of fault cell.
  std::map<PetscInt, PetscInt> _cohesiveToFault;

// PRIVATE MEMBERS ////////////////////////////////////////////////////
private :

  /// If true, use fault mesh to define fault; otherwise, use group of
  /// vertices to define fault.
  bool _useFaultMesh;

  // NOT IMPLEMENTED ////////////////////////////////////////////////////
private :

  FaultCohesive(const FaultCohesive&); ///< Not implemented
  const FaultCohesive& operator=(const FaultCohesive&); ///< Not implemented

}; // class FaultCohesive

#include "FaultCohesive.icc" // inline methods

#endif // pylith_faults_faultcohesive_hh


// End of file 
