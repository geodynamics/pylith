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

/** @file libsrc/bc/BoundaryCondition.hh
 *
 * @brief C++ abstract base class for BoundaryCondition object.
 *
 * Interface definition for boundary conditions.
 */

#if !defined(pylith_bc_boundarycondition_hh)
#define pylith_bc_boundarycondition_hh

// Include directives ---------------------------------------------------
#include "bcfwd.hh" // forward declarations

#include "pylith/topology/topologyfwd.hh" // USES Mesh
#include "pylith/utils/arrayfwd.hh" // USES scalar_array

#include "spatialdata/geocoords/geocoordsfwd.hh" // USES CoordSys
#include "spatialdata/spatialdb/spatialdbfwd.hh" // USES SpatialDB

#include <string> // HASA std::string

// BoundaryCondition ----------------------------------------------------
/// Abstract base class for boundary conditions.
class pylith::bc::BoundaryCondition
{ // class BoundaryCondition
  friend class TestBoundaryCondition; // unit testing

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Default constructor.
  BoundaryCondition(void);

  /// Destructor.
  virtual
  ~BoundaryCondition(void);

  /// Deallocate PETSc and local data structures.
  virtual
  void deallocate(void);
  
  /** Set label of boundary condition surface.
   *
   * @param value Label of surface (from mesh generator).
   */
  void label(const char* value);

  /** Get label of boundary condition surface.
   *
   * @returns Label of surface (from mesh generator).
   */
  const char* label(void) const;

  /** Verify configuration.
   *
   * @param mesh Finite-element mesh.
   */
  virtual
  void verifyConfiguration(const topology::Mesh& mesh) const;

  /** Initialize boundary condition.
   *
   * @param mesh Finite-element mesh.
   * @param upDir Vertical direction (somtimes used in 3-D problems).
   */
  virtual
  void initialize(const topology::Mesh& mesh,
		  const PylithScalar upDir[3]) = 0;

  // PROTECTED MEMBERS //////////////////////////////////////////////////
protected :

  std::string _label; ///< Label of boundary condition

  // NOT IMPLEMENTED ////////////////////////////////////////////////////
private :

  /// Not implemented
  BoundaryCondition(const BoundaryCondition&);

  /// Not implemented
  const BoundaryCondition& operator=(const BoundaryCondition&);

}; // class BoundaryCondition

#include "BoundaryCondition.icc" // inline methods

#endif // pylith_bc_boundarycondition_hh


// End of file 
