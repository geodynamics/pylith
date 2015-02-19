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

/** @file libsrc/faults/TractPerturbation.hh
 *
 * @brief C++ implementation of a spatial and temporal perturbation in
 * tractions.
 */

#if !defined(pylith_faults_tractperturbation_hh)
#define pylith_faults_tractperturbation_hh

// Include directives ---------------------------------------------------
#include "faultsfwd.hh" // forward declarations
#include "pylith/bc/TimeDependent.hh" // ISA TimeDependent

#include "spatialdata/units/unitsfwd.hh" // USES Nondimensional

// TractPerturbation -----------------------------------------------------------
/**
 * @brief C++ implementation of a spatial and temporal perturbation in
 * tractions.
 */
class pylith::faults::TractPerturbation : public pylith::bc::TimeDependent
{ // class TractPerturbation
  friend class TestTractPerturbation; // unit testing

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Default constructor.
  TractPerturbation(void);

  /// Destructor.
  virtual
  ~TractPerturbation(void);

  /// Deallocate PETSc and local data structures.
  virtual
  void deallocate(void);
  
  /** Set label for traction perturbation.
   *
   * @param value Label.
   */
  void label(const char* value);

  /** Get parameter fields.
   *
   * @returns Parameter fields.
   */
  const topology::Fields* parameterFields(void) const;
  
  /** Initialize slip time function.
   *
   * @param faultMesh Finite-element mesh of fault.
   * @param faultOrientation Orientation of fault.
   * @param normalizer Nondimensionalization of scales.
   */
  void initialize(const topology::Mesh& faultMesh,
		  const topology::Field& faultOrientation,
		  const spatialdata::units::Nondimensional& normalizer);

  /** Calculate spatial and temporal variation of value.
   *
   * @param t Current time.
   */
  void calculate(const PylithScalar t);

  /** Determine if perturbation has a given parameter.
   *
   * @param name Name of parameter field.
   * @returns True if perturbation has parameter field, false otherwise.
   */
  bool hasParameter(const char* name) const;

  /** Get vertex field with traction perturbation information.
   *
   * @param name Name of field.
   * @param fields Solution fields.
   *
   * @returns Traction vector field.
   */
  const topology::Field& vertexField(const char* name,
				     const topology::SolutionFields* const fields =0);
  
  // PROTECTED METHODS //////////////////////////////////////////////////
protected :

  /** Get label of boundary condition surface.
   *
   * @returns Label of surface (from mesh generator).
   */
  const char* _getLabel(void) const;

  /** Query database for values.
   *
   * @param name Name of field associated with database.
   * @param db Spatial database with values.
   * @param querySize Number of values at each location.
   * @param scale Dimension scale associated with values.
   * @param normalizer Nondimensionalization of scales.
   */
  void _queryDB(const char* name,
		spatialdata::spatialdb::SpatialDB* const db,
		const int querySize,
		const PylithScalar scale,
		const spatialdata::units::Nondimensional& normalizer);

// PRIVATE MEMBERS //////////////////////////////////////////////////////
private :
  
  
  topology::Fields* _parameters; ///< Parameters for perturbations.

  /// Time scale for current time.
  PylithScalar _timeScale;

  /// Label for traction perturbation.
  std::string _label;

  // NOT IMPLEMENTED ////////////////////////////////////////////////////
private :

  TractPerturbation(const TractPerturbation&); ///< Not implemented
  const TractPerturbation& operator=(const TractPerturbation&); ///< Not implemented

}; // class TractPerturbation

#endif // pylith_faults_tractperturbation_hh


// End of file 
