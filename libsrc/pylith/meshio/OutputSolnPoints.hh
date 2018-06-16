// -*- C++ -*-
//
// ======================================================================
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University of Chicago
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2017 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

/**
 * @file libsrc/meshio/OutputSolnPoints.hh
 *
 * @brief C++ object for managing output of finite-element data over a
 * an arbitrary set of points.
 */

#if !defined(pylith_meshio_outputsolnpoints_hh)
#define pylith_meshio_outputsolnpoints_hh

// Include directives ---------------------------------------------------
#include "meshiofwd.hh" // forward declarations
#include "pylith/meshio/OutputSoln.hh" // ISA OutputSoln

#include "pylith/utils/petscfwd.h"

#include "spatialdata/geocoords/geocoordsfwd.hh" // USES CoordSys

// OutputSolnPoints -----------------------------------------------------
/** @brief C++ object for managing output of finite-element data over
 * an arbitrary set of points.
 */
class pylith::meshio::OutputSolnPoints : public pylith::meshio::OutputSoln {
    friend class TestOutputSolnPoints;   // unit testing

    // PUBLIC METHODS ///////////////////////////////////////////////////////
public:

    /** Constructor
     *
     * @param[in] problem Problem to observe.
     */
    OutputSolnPoints(pylith::problems::Problem* const problem);

    /// Destructor
    ~OutputSolnPoints(void);

    /// Deallocate PETSc and local data structures.
    void deallocate(void);

    /** Set station coordinates and names.
     *
     * @param[in] points Array of station coordinates [numStations * spaceDim].
     * @param[in] numStations Number of stations.
     * @param[in] spaceDim Spatial dimension for coordinates.
     * @param[in] stationNames Array with station names.
     * @param[in] numStationNames Number of station banes.
     */
    void stations(const PylithReal* stationCoords,
                  const PylithInt numStations,
                  const PylithInt spaceDim,
                  const char* const* stationNames,
                  const PylithInt numStationNames);

    // PROTECTED MEMBERS ////////////////////////////////////////////////////
protected:

    /** Write solution at time step.
     *
     * @param[in] t Current time.
     * @param[in] tindex Current time step.
     * @param[in] solution Solution at time t.
     */
    void _writeDataStep(const PylithReal t,
                        const PylithInt tindex,
                        const pylith::topology::Field& solution);

    /** Setup interpolatior.
     *
     * @param[in] mesh Finite-element mesh of domain.
     */
    void _setupInterpolator(const pylith::topology::Mesh& mesh);

    // Interpolate field.
    void _interpolateField(void);

    /// Write dataset with names of points to file.
    void _writePointNames(void);

    // PRIVATE MEMBERS //////////////////////////////////////////////////////
private:

    pylith::scalar_array _stationCoords; ///< Array of station coordinates.
    pylith::string_vector _stationNames; ///< Array of station names.
    pylith::topology::Mesh* _stationsMesh;   ///< Mesh for stations (no cells).
    DMInterpolationInfo _interpolator;   ///< Field interpolator.

    static const char* _pyreComponent; ///< Name of Pyre component.

    // NOT IMPLEMENTED //////////////////////////////////////////////////////
private:

    OutputSolnPoints(const OutputSolnPoints&);   ///< Not implemented.
    const OutputSolnPoints& operator=(const OutputSolnPoints&);   ///< Not implemented

}; // OutputSolnPoints

#endif // pylith_meshio_outputsolnpoints_hh

// End of file
