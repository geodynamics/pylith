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

/** @file libsrc/faults/Fault.hh
 *
 */

#if !defined(pylith_faults_fault_hh)
#define pylith_faults_fault_hh

// Include directives ---------------------------------------------------
#include "faultsfwd.hh" // forward declarations

#include "pylith/topology/topologyfwd.hh" // USES Field
#include "pylith/utils/arrayfwd.hh" // USES scalar_array

#include "spatialdata/spatialdb/spatialdbfwd.hh" // USES SpatialDB

#include <string> // HASA std::string

// Fault ----------------------------------------------------------------
/**
 * @brief C++ abstract base class for Fault object.
 *
 * Interface definition for fault.
 *
 * The fault id is associated with the material-id for the fault and
 * the label is associated with the group of vertices that define the
 * fault surface.
 */
class pylith::faults::Fault
{ // class Fault
  friend class TestFault; // unit testing

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Default constructor.
  Fault(void);

  /// Destructor.
  virtual
  ~Fault(void);

  /// Deallocate PETSc and local data structures.
  virtual
  void deallocate(void);
  
  /** Set material identifier of fault.
   *
   * @param value Fault identifier
   */
  void id(const int value);

  /** Get material identifier of fault.
   *
   * @returns Fault identifier
   */
  int id(void) const;

  /** Set label of group of vertices associated with fault.
   *
   * @param value Label of fault
   */
  void label(const char* value);

  /** Get label of group of vertices associated with fault.
   *
   * @returns Label of fault
   */
  const char* label(void) const;

  /** Set label of group of vertices defining buried edge of fault.
   *
   * @param value Label of fault
   */
  void edge(const char* value);

  /** Get label of group of vertices defining buried edge of fault.
   *
   * @returns Label of fault
   */
  const char* edge(void) const;

  /** Get dimension of mesh.
   *
   * @returns Dimension of mesh.
   */
  int dimension(void) const;

  /** Get number of vertices per cell for mesh.
   *
   * @returns Number of vertices per cell for mesh.
   */
  int numCorners(void) const;
  
  /** Get number of vertices in mesh.
   *
   * @returns Number of vertices in mesh.
   */
  int numVertices(void) const;
  
  /** Get number of cells in mesh.
   *
   * @returns Number of cells in mesh.
   */
  int numCells(void) const;

  /** Get the number of vertices associated with the fault (before
   * fault mesh exists).
   *
   * @param mesh PETSc mesh
   * @return Number of vertices on the fault.
   */
  virtual
  int numVerticesNoMesh(const topology::Mesh& mesh) const = 0;

  /** Adjust mesh topology for fault implementation.
   *
   * @param mesh PETSc mesh
   * @param firstFaultVertex The first point eligible to become a new fault vertex
   * @param firstLagrangeVertex The first point eligible to become a new Lagrange vertex
   * @param firstFaultCell The first point eligible to become a new fault cell
   */
  virtual
  void adjustTopology(topology::Mesh* const mesh,
                      int *firstFaultVertex,
                      int *firstLagrangeVertex,
                      int *firstFaultCell) = 0;

  /** Initialize fault. Determine orientation and setup boundary
   * condition parameters.
   *
   * @param mesh PETSc mesh
   * @param cs Coordinate system for mesh
   * @param upDir Direction perpendicular to along-strike direction that is 
   *   not collinear with fault normal (usually "up" direction but could 
   *   be up-dip direction; only applies to fault surfaces 2-D and 3-D).
   */
  virtual
  void initialize(const topology::Mesh& mesh,
		  const PylithScalar upDir[3]) = 0;

  /** Get mesh associated with fault fields.
   *
   * @returns PETSc mesh object
   */
  const topology::Mesh& faultMesh(void) const;

  /** Get vertex field associated with integrator.
   *
   * @param name Name of vertex field.
   * @param fields Solution fields.
   * @returns Vertex field.
   */
  virtual
  const topology::Field& vertexField(const char* name,
				     const topology::SolutionFields* fields =0) = 0;

  /** Get cell field associated with integrator.
   *
   * @param name Name of cell field.
   * @param fields Solution fields.
   * @returns Cell field.
   */
  virtual
  const topology::Field& cellField(const char* name,
				   const topology::SolutionFields* fields =0) = 0;

  // NOT IMPLEMENTED ////////////////////////////////////////////////////
private :

  Fault(const Fault&); ///< Not implemented
  const Fault& operator=(const Fault&); ///< Not implemented

// PROTECTED MEMBERS ////////////////////////////////////////////////////
protected :

  topology::Mesh* _faultMesh; ///< Mesh over fault surface

  // PRIVATE MEMBERS ////////////////////////////////////////////////////
private :

  int _id; ///< Fault identifier
  std::string _label; ///< Label for points associated with fault.
  std::string _edge; ///< Label for points defining edge of fault.

}; // class Fault

#include "Fault.icc" // inline methods

#endif // pylith_faults_fault_hh


// End of file 
