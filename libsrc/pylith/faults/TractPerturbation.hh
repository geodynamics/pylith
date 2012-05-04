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
// Copyright (c) 2010-2012 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

/** @file libsrc/faults/TractPerturbation.hh
 *
 * @brief C++ abstract base class for traction perturbation function.
 */

#if !defined(pylith_faults_tractperturbation_hh)
#define pylith_faults_tractperturbation_hh

// Include directives ---------------------------------------------------
#include "faultsfwd.hh" // forward declarations

#include "pylith/topology/topologyfwd.hh" // USES Fields<SubMesh>

#include "spatialdata/units/unitsfwd.hh" // USES Nondimensional

// TractPerturbation -----------------------------------------------------------
/**
 * @brief Abstract base class for traction perturbation function.
 *
 * Interface definition for traction perturbation function.
 */
class pylith::faults::TractPerturbation
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
  
  /** Initialize slip time function.
   *
   * @param faultMesh Finite-element mesh of fault.
   * @param cs Coordinate system for mesh
   * @param normalizer Nondimensionalization of scales.
   * @param originTime Origin time for earthquake source.
   */
  virtual
  void initialize(const topology::SubMesh& faultMesh,
		  const spatialdata::units::Nondimensional& normalizer) = 0;

  /** Get traction perturbation on fault surface at time t.
   *
   * @param tractionField Traction field over fault surface.
   * @param t Time t.
   */
  virtual
  void traction(topology::Field<topology::SubMesh>* const tractionField,
		const PylithScalar t) = 0;
  
  /** Get amplitude of traction perturbation.
   *
   * @returns Traction field.
   */
  virtual
  const topology::Field<topology::SubMesh>& amplitude(void) = 0;

  /** Get parameter fields.
   *
   * @returns Parameter fields.
   */
  const topology::Fields<topology::Field<topology::SubMesh> >*
  parameterFields(void) const;

// PROTECTED MEMBERS ////////////////////////////////////////////////////
protected :

  /// Parameters for slip time function.
  topology::Fields<topology::Field<topology::SubMesh> >* _parameters;

  // NOT IMPLEMENTED ////////////////////////////////////////////////////
private :

  TractPerturbation(const TractPerturbation&); ///< Not implemented
  const TractPerturbation& operator=(const TractPerturbation&); ///< Not implemented

}; // class TractPerturbation

#endif // pylith_faults_tractperturbation_hh


// End of file 
