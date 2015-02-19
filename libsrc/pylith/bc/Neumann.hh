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

/** @file libsrc/bc/Neumann.hh
 *
 * @brief C++ implementation of time dependent Neumann (traction)
 * boundary conditions applied to a simply-connected boundary.
 */

#if !defined(pylith_bc_neumann_hh)
#define pylith_bc_neumann_hh

// Include directives ---------------------------------------------------
#include "BCIntegratorSubMesh.hh" // ISA BCIntegratorSubMesh
#include "TimeDependent.hh" // ISA TimeDependent

// Neumann ------------------------------------------------------
/// @brief Time dependent Neumann (traction) boundary conditions
/// applied to a simply-connected boundary.
class pylith::bc::Neumann : public BCIntegratorSubMesh, 
			    public TimeDependent
{ // class Neumann
  friend class TestNeumann; // unit testing

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Default constructor.
  Neumann(void);

  /// Destructor.
  ~Neumann(void);

  /// Deallocate PETSc and local data structures.
  void deallocate(void);
  
  /** Initialize boundary condition.
   *
   * @param mesh Finite-element mesh.
   * @param upDir Direction perpendicular to horizontal surface tangent 
   *   direction that is not collinear with surface normal.
   */
  void initialize(const topology::Mesh& mesh,
		  const PylithScalar upDir[3]);

  /** Integrate contributions to residual term (r) for operator.
   *
   * @param residual Field containing values for residual.
   * @param t Current time.
   * @param fields Solution fields.
   */
  void integrateResidual(const topology::Field& residual,
			 const PylithScalar t,
			 topology::SolutionFields* const fields);

  /** Verify configuration is acceptable.
   *
   * @param mesh Finite-element mesh
   */
  void verifyConfiguration(const topology::Mesh& mesh) const;

  /** Get cell field with BC information.
   *
   * @param name Name of field.
   * @param fields Solution fields.
   *
   * @returns Traction vector at integration points.
   */
  const topology::Field& cellField(const char* name,
				   topology::SolutionFields* const fields =0);

  // PROTECTED METHODS //////////////////////////////////////////////////
protected :

  /** Get label of boundary condition surface.
   *
   * @returns Label of surface (from mesh generator).
   */
  const char* _getLabel(void) const;

  /** Get manager of scales used to nondimensionalize problem.
   *
   * @returns Nondimensionalizer.
   */
  const spatialdata::units::Nondimensional& _getNormalizer(void) const;

  /// Query databases for time dependent parameters.
  void _queryDatabases(void);

  /** Query database for values.
   *
   * @param name Name of field associated with database.
   * @param db Spatial database with values.
   * @param querySize Number of values at each location.
   * @param scale Dimension scale associated with values.
   */
  void _queryDB(const char* name,
		spatialdata::spatialdb::SpatialDB* const db,
		const int querySize,
		const PylithScalar scale);

  /** Convert parameters in local coordinates to global coordinates.
   *
   * @param upDir Direction perpendicular to horizontal surface tangent 
   *   direction that is not collinear with surface normal.
   */
  void _paramsLocalToGlobal(const PylithScalar upDir[3]);

  /** Calculate spatial and temporal variation of value over the list
   *  of submesh.
   *
   * @param t Current time.
   */
  void _calculateValue(const PylithScalar t);

  // NOT IMPLEMENTED ////////////////////////////////////////////////////
private :

  Neumann(const Neumann&); ///< Not implemented.
  const Neumann& operator=(const Neumann&); ///< Not implemented.

}; // class Neumann

#include "Neumann.icc"

#endif // pylith_bc_neumann_hh


// End of file 
