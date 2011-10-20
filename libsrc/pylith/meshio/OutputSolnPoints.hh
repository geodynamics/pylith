// -*- C++ -*-
//
// ======================================================================
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University of Chicago
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2011 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

/**
 * @file libsrc/meshio/OutputSolnPoints.hh
 *
 * @brief C++ object for managing output of finite-element data over a
 * subdomain.
 */

#if !defined(pylith_meshio_outputsolnpoints_hh)
#define pylith_meshio_outputsolnpoints_hh

// Include directives ---------------------------------------------------
#include "meshiofwd.hh" // forward declarations

#include "pylith/topology/SubMesh.hh" // ISA OutputManager<Mesh>
#include "pylith/topology/Field.hh" // ISA OutputManager<Field<Mesh>>
#include "OutputManager.hh" // ISA OutputManager

#include <string> // HASA std::string

// OutputSolnPoints -----------------------------------------------------
/** @brief C++ object for managing output of finite-element data over
 * a subdomain.
 */
class pylith::meshio::OutputSolnPoints : 
  public OutputManager<topology::Mesh, topology::Field<topology::Mesh> >
{ // OutputSolnPoints
  friend class TestOutputSolnPoints; // unit testing

// PUBLIC METHODS ///////////////////////////////////////////////////////
public :

  /// Constructor
  OutputSolnPoints(void);

  /// Destructor
  ~OutputSolnPoints(void);

  /// Deallocate PETSc and local data structures.
  void deallocate(void);
  
  /** Get mesh associated with points.
   *
   * @returns Mesh associated with points.
   */
  const topology::Mesh& createPointsMesh(const PylithScalar* points,
					 const int numPoints,
					 const int spaceDim);
  
// NOT IMPLEMENTED //////////////////////////////////////////////////////
private :

  OutputSolnPoints(const OutputSolnPoints&); ///< Not implemented.
  const OutputSolnPoints& operator=(const OutputSolnPoints&); ///< Not implemented

// PRIVATE TYPEDEFS /////////////////////////////////////////////////////
private :

  typedef DMMeshInterpolationInfo PetscDMMeshInterpolationInfo;

// PRIVATE MEMBERS //////////////////////////////////////////////////////
private :

  topology::Mesh* _mesh; ///< Domain mesh.
  topology::Mesh* _pointsMesh; ///< Mesh for points (no cells).
  PetscDMMeshInterpolationInfo* _interpolator;

}; // OutputSolnPoints

#endif // pylith_meshio_outputsolnpoints_hh

// End of file 
