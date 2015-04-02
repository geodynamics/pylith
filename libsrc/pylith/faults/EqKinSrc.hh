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

/** @file libsrc/faults/EqKinSrc.hh
 *
 * @brief C++ object for managing parameters for a kinematic
 * earthquake source.
 */

#if !defined(pylith_faults_eqkinsrc_hh)
#define pylith_faults_eqkinsrc_hh

// Include directives ---------------------------------------------------
#include "faultsfwd.hh" // forward declarations

#include "pylith/topology/topologyfwd.hh" // USES Field<Mesh>

#include "spatialdata/units/unitsfwd.hh" // USES Nondimensional

// EqKinSrc -------------------------------------------------------------
/** @brief Kinematic earthquake source.
 *
 * EqKinSrc is responsible for providing the value of slip at time t
 * over a fault surface.
 */
class pylith::faults::EqKinSrc
{ // class EqKinSrc
  friend class TestEqKinSrc; // unit testing

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Default constructor.
  EqKinSrc(void);

  /// Destructor.
  ~EqKinSrc(void);

  /// Deallocate PETSc and local data structures.
  virtual
  void deallocate(void);
  
  /** Set origin time for earthquake source.
   *
   * @param value Origin time for earthquake source.
   */
  void originTime(const PylithScalar value);

  /** Get origin time for earthquake source.
   *
   * @returns Origin time for earthquake source.
   */
  PylithScalar originTime(void) const;

  /** Set slip time function.
   *
   * @param slipfn Slip time function.
   */
  void slipfn(SlipTimeFn* slipfn);

  /** Initialize slip time function.
   *
   * @param faultMesh Finite-element mesh of fault.
   * @param normalizer Nondimensionalization of scales.
   */
  void initialize(const topology::Mesh& faultMesh,
		  const spatialdata::units::Nondimensional& normalizer);

  /** Get slip on fault surface at time t.
   *
   * @param slipField Slip field over fault mesh.
   * @param t Time t.
   */
  void slip(topology::Field* const slipField,
	    const PylithScalar t);

  /** Get final slip.
   *
   * @returns Final slip.
   */
  const topology::Field& finalSlip(void) const;

  /** Get time when slip begins at each point.
   *
   * @returns Time when slip begins.
   */
  const topology::Field& slipTime(void) const;

  // NOT IMPLEMENTED ////////////////////////////////////////////////////
private :

  EqKinSrc(const EqKinSrc&); ///< Not implemented
  const EqKinSrc& operator=(const EqKinSrc&); ///< Not implemented

  // PRIVATE MEMBERS ////////////////////////////////////////////////////
private :

  PylithScalar _originTime; ///< Origin time for earthquake source
  SlipTimeFn* _slipfn; ///< Slip time function

}; // class EqKinSrc

#endif // pylith_faults_eqkinsrc_hh


// End of file 
