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

/** @file libsrc/faults/SlipTimeFn.hh
 *
 * @brief C++ abstract base class for kinematic slip time function.
 */

#if !defined(pylith_faults_sliptimefn_hh)
#define pylith_faults_sliptimefn_hh

// Include directives ---------------------------------------------------
#include "faultsfwd.hh" // forward declarations

#include "pylith/topology/topologyfwd.hh" // USES Fields<Mesh>

#include "spatialdata/units/unitsfwd.hh" // USES Nondimensional

// SlipTimeFn -----------------------------------------------------------
/**
 * @brief Abstract base class for kinematic slip time function.
 *
 * Interface definition for slip time function.
 */
class pylith::faults::SlipTimeFn
{ // class SlipTimeFn
  friend class TestSlipTimeFn; // unit testing

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Default constructor.
  SlipTimeFn(void);

  /// Destructor.
  virtual
  ~SlipTimeFn(void);

  /// Deallocate PETSc and local data structures.
  virtual
  void deallocate(void);
  
  /** Initialize slip time function.
   *
   * @param faultMesh Finite-element mesh of fault.
   * @param normalizer Nondimensionalization of scales.
   * @param originTime Origin time for earthquake source.
   */
  virtual
  void initialize(const topology::Mesh& faultMesh,
		  const spatialdata::units::Nondimensional& normalizer,
		  const PylithScalar originTime =0.0) = 0;

  /** Get slip on fault surface at time t.
   *
   * @param slipField Slip field over fault surface.
   * @param t Time t.
   */
  virtual
  void slip(topology::Field* const slipField,
	    const PylithScalar t) = 0;
  
  /** Get final slip.
   *
   * @returns Final slip.
   */
  virtual
  const topology::Field& finalSlip(void) = 0;

  /** Get time when slip begins at each point.
   *
   * @returns Time when slip begins.
   */
  virtual
  const topology::Field& slipTime(void) = 0;

  /** Get parameter fields.
   *
   * @returns Parameter fields.
   */
  const topology::Fields* parameterFields(void) const;

// PROTECTED MEMBERS ////////////////////////////////////////////////////
protected :

  topology::Fields* _parameters; ///< Parameters for slip time function.

  // NOT IMPLEMENTED ////////////////////////////////////////////////////
private :

  SlipTimeFn(const SlipTimeFn&); ///< Not implemented
  const SlipTimeFn& operator=(const SlipTimeFn&); ///< Not implemented

}; // class SlipTimeFn

#endif // pylith_faults_sliptimefn_hh


// End of file 
