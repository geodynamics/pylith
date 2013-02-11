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

/** @file libsrc/faults/DKSelector.hh
 *
 * @brief C++ implementation of the dynamic kinematic selector
 */

#if !defined(pylith_faults_dkselector_hh)
#define pylith_faults_dkselector_hh

// Include directives ---------------------------------------------------
#include "faultsfwd.hh" // forward declaration


#include "pylith/topology/topologyfwd.hh" // USES Fields<Field<SubMesh> >
#include "spatialdata/spatialdb/spatialdbfwd.hh" // USES SpatialDB
#include "spatialdata/units/unitsfwd.hh" // USES Nondimensional

#include "pylith/topology/SubMesh.hh" // USES Mesh

#include "pylith/utils/array.hh" // HASA scalar_array

// DKSelector -----------------------------------------------------------
/** @brief Dynamic-Kinematic Selector
 *
 * If value is over .5 the fault vertex has its slip value controlled by EQsrc
 *
 */

// ----------------------------------------------------------------------
typedef pylith::topology::SubMesh::RealUniformSection RealUniformSection;

// ----------------------------------------------------------------------
class pylith::faults::DKSelector 
{ // class DKSelector

// PUBLIC METHODS ///////////////////////////////////////////////////////
public :

  /// Default constructor.
  DKSelector(void);

  /// Destructor.
  ~DKSelector(void);

  /// Deallocate PETSc and local data structures.
  virtual
  void deallocate(void);
  
  /** Set spatial database for dkselector
   *
   * @param db dksel
   */
  void dbdksel(spatialdata::spatialdb::SpatialDB* const db);

  /** Initialize dkselector
   *
   * @param faultMesh Finite-element mesh of fault.
   * @param normalizer Nondimensionalization of scales.
   */
  void initialize(const topology::SubMesh& faultMesh,
		  const spatialdata::units::Nondimensional& normalizer);

  /** Get dk on fault surface (time will be implemented through this guy)
   *
   * @params faultMesh Finite-element mesh of fault.
   * 
   * @returns a section for the given time (future)
   */
  void dk(const topology::Field<topology::SubMesh>* dk);
  
// NOT IMPLEMENTED //////////////////////////////////////////////////////
private :

  DKSelector(const DKSelector&); ///< Not implemented
  const DKSelector& operator=(const DKSelector&); ///< Not implemented

// PRIVATE MEMBERS //////////////////////////////////////////////////////
private :

  scalar_array _dkselVertex; ///< dk sel at vertex.
  PylithScalar _dkselvv; ///< used to check 

  /// Spatial database for dk selector.
  spatialdata::spatialdb::SpatialDB* _dbdksel;

  /// Parameters for perturbations.
  topology::Fields<topology::Field<topology::SubMesh> >* _parameters;

}; // class DKSelector

#endif // pylith_faults_dkselector_hh


// End of file 
