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
 * subdomain.
 */

#if !defined(pylith_meshio_outputsolnpoints_hh)
#define pylith_meshio_outputsolnpoints_hh

// Include directives ---------------------------------------------------
#include "meshiofwd.hh" // forward declarations
#include "pylith/utils/petscfwd.h"

#include "spatialdata/units/unitsfwd.hh" // USES Nondimensional

#include "pylith/topology/Mesh.hh" // ISA OutputManager<Mesh>
#include "pylith/topology/Field.hh" // ISA OutputManager<Field<Mesh>>
#include "OutputManager.hh" // ISA OutputManager

// OutputSolnPoints -----------------------------------------------------
/** @brief C++ object for managing output of finite-element data over
 * a subdomain.
 */
class pylith::meshio::OutputSolnPoints : public OutputManager {
    friend class TestOutputSolnPoints;   // unit testing

    // PUBLIC METHODS ///////////////////////////////////////////////////////
public:

    /// Constructor
    OutputSolnPoints(void);

    /// Destructor
    ~OutputSolnPoints(void);

    /// Deallocate PETSc and local data structures.
    void deallocate(void);

    /** Setup interpolator.
     *
     * @param mesh Domain mesh.
     * @param points Array of dimensioned coordinates for points [numPoints*spaceDim].
     * @param numPoints Number of points.
     * @param spaceDim Spatial dimension for coordinates.
     * @param names Array with name for each point, e.g., station name.
     * @param nunNames Number of names in array.
     * @param normalizer Nondimensionalizer.
     */
    void setupInterpolator(pylith::topology::Mesh* mesh,
                           const PylithScalar* points,
                           const int numPoints,
                           const int spaceDim,
                           const char* const* names,
                           const int numNames,
                           const spatialdata::units::Nondimensional& normalizer);

    // PROTECTED MEMBERS ////////////////////////////////////////////////////
protected:

    /** Prepare for output.
     *
     * @param[in] mesh Finite-element mesh object.
     * @param[in] isInfo True if only writing info values.
     */
    void _open(const pylith::topology::Mesh& mesh,
               const bool isInfo);

    /** Append finite-element vertex field to file.
     *
     * @param t Time associated with field.
     * @param field Vertex field.
     * @param mesh Mesh for output.
     */
    void _appendField(const PylithReal t,
                      pylith::topology::Field& field,
                      const topology::Mesh& mesh);

    /// Write dataset with names of points to file.
    void _writePointNames(void);

    // NOT IMPLEMENTED //////////////////////////////////////////////////////
private:

    OutputSolnPoints(const OutputSolnPoints&);   ///< Not implemented.
    const OutputSolnPoints& operator=(const OutputSolnPoints&);   ///< Not implemented

    // PRIVATE MEMBERS //////////////////////////////////////////////////////
private:

    pylith::topology::Mesh* _mesh;   ///< Domain mesh.
    pylith::topology::Mesh* _pointsMesh;   ///< Mesh for points (no cells).
    DMInterpolationInfo _interpolator;   ///< Field interpolator.
    pylith::string_vector _stations; ///< Array of station names.

}; // OutputSolnPoints

#endif // pylith_meshio_outputsolnpoints_hh

// End of file
