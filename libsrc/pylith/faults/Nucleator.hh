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

/** @file libsrc/faults/Nucleator.hh
 *
 * @brief C++ object for managing traction perturbations on faults
 * with spontaneous rupture.
 */

#if !defined(pylith_faults_nucleator_hh)
#define pylith_faults_nucleator_hh

// Include directives ---------------------------------------------------
#include "faultsfwd.hh" // forward declarations

#include "pylith/topology/topologyfwd.hh" // USES Field<SubMesh>

#include "spatialdata/units/unitsfwd.hh" // USES Nondimensional

// Nucleator -------------------------------------------------------------
/** @brief Traction perturbations on faults with spontaneous rupture.
 *
 * Nucleator is responsible for providing perturbations in the
 * traction field over a fault surface at time t.
 */
class pylith::faults::Nucleator
{ // class Nucleator
  friend class TestNucleator; // unit testing

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Default constructor.
  Nucleator(void);

  /// Destructor.
  ~Nucleator(void);

  /// Deallocate PETSc and local data structures.
  virtual
  void deallocate(void);
  
  /** Set origin time for perturbation.
   *
   * @param value Origin time for perturbation.
   */
  void originTime(const PylithScalar value);

  /** Get origin time for perturbation.
   *
   * @returns Origin time for perturbation.
   */
  PylithScalar originTime(void) const;

  /** Set traction perturbation.
   *
   * @param tractfn Traction perturbation function.
   */
  void perturbationFn(TractPerturbation* tractfn);

  /** Initialize traction perturbation function.
   *
   * @param faultMesh Finite-element mesh of fault.
   * @param normalizer Nondimensionalization of scales.
   */
  void initialize(const topology::SubMesh& faultMesh,
		  const spatialdata::units::Nondimensional& normalizer);

  /** Get traction perturbation on fault surface at time t.
   *
   * @param tractionField Traction field over fault mesh.
   * @param t Time t.
   */
  void traction(topology::Field<topology::SubMesh>* const tractionField,
		const PylithScalar t);

  /** Get amplitude of spatial variation of traction.
   *
   * @returns Spatial variation of traction amplitude.
   */
  const topology::Field<topology::SubMesh>& tractionAmp(void) const;

  // NOT IMPLEMENTED ////////////////////////////////////////////////////
private :

  Nucleator(const Nucleator&); ///< Not implemented
  const Nucleator& operator=(const Nucleator&); ///< Not implemented

  // PRIVATE MEMBERS ////////////////////////////////////////////////////
private :

  PylithScalar _originTime; ///< Origin time for traction perturbation.
  TractPerturbation* _tractfn; ///< Traction perturbation function.

}; // class Nucleator

#endif // pylith_faults_nucleator_hh


// End of file 
