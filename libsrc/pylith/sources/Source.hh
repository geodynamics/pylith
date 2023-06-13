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

    /** Set descriptive label for source.
     *
     * @param value Label of source.
     */
    void setDescription(const char* value);

    /** Get descruptive label of source.
     *
     * @returns Label of source
     */
    const char* getDescription(void) const;

    /** Set name of label marking boundary associated with boundary condition surface.
     *
     * @param[in] value Name of label for surface (from mesh generator).
     */
    void setLabelName(const char* value);

    /** Get name of label marking boundary associated with source location.
     *
     * @returns Name of label for surface (from mesh generator).
     */
    const char* getLabelName(void) const;

    /** Set value of label marking boundary associated with source location.
     *
     * @param[in] value Value of label for surface (from mesh generator).
     */
    void setLabelValue(const int value);

    /** Get value of label marking boundary associated with source location.
     *
     * @returns Value of label for surface (from mesh generator).
     */
    int getLabelValue(void) const;

    /** Set first choice for reference direction to discriminate among tangential directions in 3-D.
     *
     * @param vec Reference direction unit vector.
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
    std::string _description; ///< Descriptive label for source.
    std::string _labelName; ///< Name of label to identify source points in mesh.
    int _labelValue; ///< Value of label to identify source points in mesh.
    
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
