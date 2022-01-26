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

/** @file libsrc/sources/Source.hh
 *
 * @brief C++ abstract base class for sources.
 */

#if !defined(pylith_sources_source_hh)
#define pylith_sources_source_hh

#include "sourcesfwd.hh" // forward declarations
#include "pylith/topology/topologyfwd.hh"
#include "pylith/utils/petscfwd.h"
#include "spatialdata/geocoords/geocoordsfwd.hh" // USES CoordSys
#include "pylith/problems/Physics.hh" // ISA Physics

#include <string> // HASA std::string

// Source -------------------------------------------------------------
/** @brief C++ abstract base class for sources.
 *
 */

class pylith::sources::Source : public pylith::problems::Physics {
    friend class TestSource; // unit testing

    // PUBLIC METHODS //////////////////////////////////////////////////////////////////////////////////////////////////
public:

    /// Default constructor.
    Source(void);

    /// Destructor.
    virtual ~Source(void);

    /// Deallocate PETSc and local data structures.
    virtual
    void deallocate(void);

    /** Set value of label source-id used to identify source cells.
     *
     * @param value Source identifier
     */
    void setSourceId(const int value);

    /** Get value of label source-id used to identify source cells.
     *
     * @returns Source identifier
     */
    int getSourceId(void) const;

    /** Set descriptive label for source.
     *
     * @param value Label of source.
     */
    void setDescriptiveLabel(const char* value);

    /** Get descruptive label of source.
     *
     * @returns Label of source
     */
    const char* getDescriptiveLabel(void) const;

    /** Create constraint and set kernels.
     *
     * @param[in] solution Solution field.
     * @returns Constraint if applicable, otherwise NULL.
     */
    virtual
    std::vector<pylith::feassemble::Constraint*> createConstraints(const pylith::topology::Field& solution);

    /** Set coordinates and names of points.
     *
     * @param[in] points Array of coordinates [numPoints * spaceDim].
     * @param[in] numPoints Number of points.
     * @param[in] spaceDim Spatial dimension for coordinates.
     * @param[in] pointNames Array with point names.
     * @param[in] numPointNames Number of point banes.
     */
    void setPoints(const PylithReal* pointCoords,
                   const int numPoints,
                   const int spaceDim,
                   const char* const* pointNames,
                   const int numPointNames);

    // PROTECTED MEMBERS ///////////////////////////////////////////////////////////////////////////////////////////////
protected:

    pylith::scalar_array _pointCoords; ///< Array of point coordinates.
    pylith::string_vector _pointNames; ///< Array of point names.
    int _sourceId; ///< Value of source-id label in mesh.
    std::string _descriptiveLabel; ///< Descriptive label for source.

    // PRIVATE MEMBERS /////////////////////////////////////////////////////////////////////////////////////////////////
private:

    pylith::topology::Mesh* _pointMesh; ///< Mesh for points (no cells).
    pylith::topology::Field* _pointSoln; ///< Solution field at points.

    // NOT IMPLEMENTED /////////////////////////////////////////////////////////////////////////////////////////////////
private:

    Source(const Source&); ///< Not implemented.
    const Source& operator=(const Source&); ///< Not implemented

}; // Source

#endif // pylith_sources_source_hh

// End of file
