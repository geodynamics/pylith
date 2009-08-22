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

/** @file libsrc/bc/TimeDependent.hh
 *
 * @brief C++ Abstract base class for time-dependent boundary
 * conditions.
 */

#if !defined(pylith_bc_timedependent_hh)
#define pylith_bc_timedependent_hh

// Include directives ---------------------------------------------------
#include "bcfwd.hh" // forward declarations

#include "pylith/utils/array.hh" // HASA int_array
#include "pylith/topology/topologyfwd.hh" // USES Mesh
#include "spatialdata/spatialdb/spatialdbfwd.hh" // USES SpatialDB
#include "spatialdata/units/unitsfwd.hh" // USES Nondimensional

// TimeDependent ------------------------------------------------------
/// Abstract base class for time-dependent boundary conditions.
class pylith::bc::TimeDependent
{ // class TimeDependent
  friend class TestTimeDependent; // unit testing

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Default constructor.
  TimeDependent(void);

  /// Destructor.
  ~TimeDependent(void);

  /// Deallocate PETSc and local data structures.
  virtual
  void deallocate(void);
  
  /** Set database for initial values.
   *
   * @param db Spatial database
   */
  void dbInitial(spatialdata::spatialdb::SpatialDB* const db);

  /** Set database for rate of change of values.
   *
   * @param db Spatial database
   */
  void dbRate(spatialdata::spatialdb::SpatialDB* const db);

  /** Set database for change in values.
   *
   * @param db Spatial database
   */
  void dbChange(spatialdata::spatialdb::SpatialDB* const db);

  /** Set database for temporal evolution of change in value.
   *
   * @param db Time history database.
   */
  void dbTimeHistory(spatialdata::spatialdb::TimeHistory* const db);

  /** Verify configuration is acceptable.
   *
   * @param mesh Finite-element mesh
   */
  virtual
  void verifyConfiguration(const topology::Mesh& mesh) const;

  // PROTECTED METHODS //////////////////////////////////////////////////
protected :

  /** Get label of boundary condition surface.
   *
   * @returns Label of surface (from mesh generator).
   */
  virtual
  const char* _getLabel(void) const = 0;

  /** Get manager of scales used to nondimensionalize problem.
   *
   * @returns Nondimensionalizer.
   */
  virtual
  const spatialdata::units::Nondimensional& _getNormalizer(void) const = 0;

  // PROTECTED MEMBERS //////////////////////////////////////////////////
protected :

  /// Spatial database for initial values.
  spatialdata::spatialdb::SpatialDB* _dbInitial;

  /// Spatial database for rate of change of values.
  spatialdata::spatialdb::SpatialDB* _dbRate;

  /// Spatial database for change in value.
  spatialdata::spatialdb::SpatialDB* _dbChange;

  /// Temporal evolution of amplitude for change in value;
  spatialdata::spatialdb::TimeHistory* _dbTimeHistory;
  
  // NOT IMPLEMENTED ////////////////////////////////////////////////////
private :

  TimeDependent(const TimeDependent&); ///< Not implemented.
  const TimeDependent& operator=(const TimeDependent&); ///< Not implemented.

}; // class TimeDependent

#include "TimeDependent.icc" // inline methods

#endif // pylith_bc_timedependent_hh


// End of file 
