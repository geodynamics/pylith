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
// Copyright (c) 2010-2015 University of California, Davis
//
// See COPYING for license information.
//
// ======================================================================
//

/**
 * @file libsrc/problems/InitialConditionsPatch.hh
 *
 * @brief C++ object for specifying initial conditions over a portion of the domain (patch).
 */
#if !defined(pylith_problems_initialconditionspatch_hh)
#define pylith_problems_initialconditionspatch_hh

#include "InitialConditions.hh" // ISA InitialConditions

#include "spatialdata/spatialdb/spatialdbfwd.hh" // HOLDSA SpatialDB

class pylith::problems::InitialConditionsPatch : public pylith::problems::InitialConditions {
    friend class TestInitialConditionsPatch; // unit testing

    // PUBLIC METHODS //////////////////////////////////////////////////////////////////////////////////////////////////
public:

    /// Constructor
    InitialConditionsPatch(void);

    /// Destructor
    virtual ~InitialConditionsPatch(void);

    /// Deallocate PETSc and local data structures.
    void deallocate(void);

    /** Set label for marker associated with patch.
     *
     * @param[in] value Label for marker associated with patch.
     */
    void setMarkerLabel(const char* value);

    /** Get label for marker associated with patch.
     *
     * @returns Label for marker associated with patch.
     */
    const char* getMarkerLabel(void) const;

    /** Set spatial database holding initial conditions.
     *
     * @param[in] db Spatial database holding initial conditions.
     */
    void setDB(spatialdata::spatialdb::SpatialDB* db);

    /** Set solver type.
     *
     * @param[out] solution Solution field.
     * @param[in] normalizer Nondimensionalization.
     */
    void setValues(pylith::topology::Field* solution,
                   const spatialdata::units::Nondimensional& normalizer);

    // PRIVATE MEMEBRS
    // //////////////////////////////////////////////////////////////////////////////////////////////////
private:

    std::string _patchLabel; ///< Marker label associated with patch.
    spatialdata::spatialdb::SpatialDB* _db; ///< Spatial database with values for initial condition.

    // NOT IMPLEMENTED /////////////////////////////////////////////////////////////////////////////////////////////////
private:

    InitialConditionsPatch(const InitialConditionsPatch&); ///< Not implemented
    const InitialConditionsPatch& operator=(const InitialConditionsPatch&); ///< Not implemented

}; // InitialConditionsPatch

#endif // pylith_problems_initialconditionspatch_hh

// End of file
