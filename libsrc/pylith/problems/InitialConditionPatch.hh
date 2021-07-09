// -*- C++ -*-
//
// ======================================================================
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University at Buffalo
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2021 University of California, Davis
//
// See LICENSE.md for license information.
//
// ======================================================================
//

/**
 * @file libsrc/problems/InitialConditionPatch.hh
 *
 * @brief C++ object for specifying initial conditions over a portion of the domain (patch).
 */
#if !defined(pylith_problems_initialconditionspatch_hh)
#define pylith_problems_initialconditionspatch_hh

#include "InitialCondition.hh" // ISA InitialCondition

#include "spatialdata/spatialdb/spatialdbfwd.hh" // HOLDSA SpatialDB

class pylith::problems::InitialConditionPatch : public pylith::problems::InitialCondition {
    friend class TestInitialConditionPatch; // unit testing

    // PUBLIC METHODS //////////////////////////////////////////////////////////////////////////////////////////////////
public:

    /// Constructor
    InitialConditionPatch(void);

    /// Destructor
    virtual ~InitialConditionPatch(void);

    /// Deallocate PETSc and local data structures.
    void deallocate(void);

    /** Set material id associated with patch.
     *
     * @param[in] value Material id associated with patch.
     */
    void setMaterialId(const int value);

    /** Get material id associated with patch.
     *
     * @returns Material id associated with patch.
     */
    int getMaterialId(void) const;

    /** Verify configuration is acceptable.
     *
     * @param[in] solution Solution field.
     */
    void verifyConfiguration(const pylith::topology::Field& solution) const;

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

    PylithInt _patchId; ///< Material id associated with patch.
    spatialdata::spatialdb::SpatialDB* _db; ///< Spatial database with values for initial condition.

    // NOT IMPLEMENTED /////////////////////////////////////////////////////////////////////////////////////////////////
private:

    InitialConditionPatch(const InitialConditionPatch&); ///< Not implemented
    const InitialConditionPatch& operator=(const InitialConditionPatch&); ///< Not implemented

}; // InitialConditionPatch

#endif // pylith_problems_initialconditionspatch_hh

// End of file
