// -*- C++ -*-
//
// ----------------------------------------------------------------------
//
//                           Brad T. Aagaard
//                        U.S. Geological Survey
//
// {LicenseText}
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
#include "pylith/utils/arrayfwd.hh" // USES double_array

#include "spatialdata/geocoords/geocoordsfwd.hh" // USES CoordSys
#include "spatialdata/spatialdb/spatialdbfwd.hh" // USES SpatialDB

#include <string> // HASA std::string

// BoundaryCondition ----------------------------------------------------
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

  /** Set database for boundary condition parameters.
   *
   * @param db Spatial database
   */
  void db(spatialdata::spatialdb::SpatialDB* const db);

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
		  const double upDir[3]) = 0;

  // NOT IMPLEMENTED ////////////////////////////////////////////////////
private :

  /// Not implemented
  BoundaryCondition(const BoundaryCondition& m);

  /// Not implemented
  const BoundaryCondition& operator=(const BoundaryCondition& m);

  // PROTECTED MEMBERS //////////////////////////////////////////////////
protected :

  std::string _label; ///< Label of boundary condition
  spatialdata::spatialdb::SpatialDB* _db; ///< Spatial database w/parameters

}; // class BoundaryCondition

#include "BoundaryCondition.icc" // inline methods

#endif // pylith_bc_boundarycondition_hh


// End of file 
