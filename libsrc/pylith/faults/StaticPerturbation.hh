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

/** @file libsrc/faults/StaticPerturbation.hh
 *
 * @brief C++ implementation of a static perturbation in tractions.
 */

#if !defined(pylith_faults_staticperturbation_hh)
#define pylith_faults_staticperturbation_hh

// Include directives ---------------------------------------------------
#include "TractPerturbation.hh"

#include "pylith/utils/array.hh" // HASA scalar_array

#include "spatialdata/spatialdb/spatialdbfwd.hh"

// StaticPerturbation -----------------------------------------------------------
/**
 * @brief C++ implementation of a static perturbation in tractions.
 *
 * T = F(x)
*/
class pylith::faults::StaticPerturbation : public TractPerturbation
{ // class StaticPerturbation
  friend class TestStaticPerturbation; // unit testing

// PUBLIC METHODS ///////////////////////////////////////////////////////
public :

  /// Default constructor.
  StaticPerturbation(void);

  /// Destructor.
  ~StaticPerturbation(void);

  /// Deallocate PETSc and local data structures.
  virtual
  void deallocate(void);
  
  /** Set spatial database for traction amplitude.
   *
   * @param db Spatial database
   */
  void dbAmplitude(spatialdata::spatialdb::SpatialDB* const db);
  
  /** Initialize static perturbation.
   *
   * @param faultMesh Finite-element mesh of fault.
   * @param normalizer Nondimensionalization of scales.
   */
  void initialize(const topology::SubMesh& faultMesh,
		  const spatialdata::units::Nondimensional& normalizer);
  
  /** Get traction on fault surface at time t.
   *
   * @param tractionField Traction field over fault surface.
   * @param t Time t.
   */
  void traction(topology::Field<topology::SubMesh>* const tractionField,
		const PylithScalar t);
  
  /** Get amplitude of traction perturbation.
   *
   * @returns Final slip.
   */
  const topology::Field<topology::SubMesh>& amplitude(void);

// NOT IMPLEMENTED //////////////////////////////////////////////////////
private :

  StaticPerturbation(const StaticPerturbation&); ///< Not implemented.
  const StaticPerturbation& operator=(const StaticPerturbation&); ///< Not implemented

// PRIVATE MEMBERS //////////////////////////////////////////////////////
private :

  /// Spatial database for traction perturbation.
  spatialdata::spatialdb::SpatialDB* _dbAmplitude;

}; // class StaticPerturbation

#endif // pylith_faults_staticperturbation_hh


// End of file 
