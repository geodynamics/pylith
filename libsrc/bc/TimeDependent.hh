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
 * @brief C++ implementation of point force on vertices.
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
class pylith::bc::TimeDependent
{ // class TimeDependent
  friend class TestTimeDependent; // unit testing

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Default constructor.
  TimeDependent(void);

  /// Destructor.
  ~TimeDependent(void);

  /** Set indices of degrees of freedom associated with BC.
   *
   * Note: Forces at all points are applied to the same degrees of freedom.
   *
   * Example: [0, 1] to apply forces to x and y degrees of freedom in
   * Cartesian system.
   *
   * @param flags Array of indices for degrees of freedom for forces.
   * @param size Size of array
   */
  void bcDOF(const int* flags,
	     const int size);  

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

  /** Set database for start time of rate change.
   *
   * @param db Spatial database
   */
  void dbRateTime(spatialdata::spatialdb::SpatialDB* const db);

  /** Set database for change in values.
   *
   * @param db Spatial database
   */
  void dbChange(spatialdata::spatialdb::SpatialDB* const db);

  /** Set database for start time of change in values.
   *
   * @param db Spatial database
   */
  void dbChangeTime(spatialdata::spatialdb::SpatialDB* const db);

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
  const char* getLabel(void) const = 0;

  /** Get manager of scales used to nondimensionalize problem.
   *
   * @returns Nondimensionalizer.
   */
  virtual
  const spatialdata::units::Nondimensional& getNormalizer(void) const = 0;

  // PROTECTED MEMBERS //////////////////////////////////////////////////
protected :

  /// Spatial database for initial values.
  spatialdata::spatialdb::SpatialDB* _dbInitial;

  /// Spatial database for rate of change of values.
  spatialdata::spatialdb::SpatialDB* _dbRate;

  /// Spatial database for start time of rate change.
  spatialdata::spatialdb::SpatialDB* _dbRateTime;

  /// Spatial database for change in value.
  spatialdata::spatialdb::SpatialDB* _dbChange;

  /// Spatial database for start time of change in value.
  spatialdata::spatialdb::SpatialDB* _dbChangeTime;

  /// Temporal evolution of amplitude for change in value;
  spatialdata::spatialdb::TimeHistory* _dbTimeHistory;
  
  int_array _bcDOF; ///< Degrees of freedom associated with BC.

  // NOT IMPLEMENTED ////////////////////////////////////////////////////
private :

  TimeDependent(const TimeDependent&); ///< Not implemented.
  const TimeDependent& operator=(const TimeDependent&); ///< Not implemented.

}; // class TimeDependent

#include "TimeDependent.icc" // inline methods

#endif // pylith_bc_timedependent_hh


// End of file 
