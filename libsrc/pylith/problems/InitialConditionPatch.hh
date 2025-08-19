// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2025, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================
#pragma once

#include "pylith/problems/InitialCondition.hh" // ISA InitialCondition

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

    /** Set name of label marking material.
     *
     * @param[in] value Name of label for material (from mesh generator).
     */
    void setLabelName(const char* value);

    /** Get name of label marking material.
     *
     * @returns Name of label for material (from mesh generator).
     */
    const char* getLabelName(void) const;

    /** Set value of label marking material.
     *
     * @param[in] value Value of label for material (from mesh generator).
     */
    void setLabelValue(const int value);

    /** Get value of label marking material.
     *
     * @returns Value of label for material (from mesh generator).
     */
    int getLabelValue(void) const;

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
     * @param[in] scales Nondimensionalization.
     */
    void setValues(pylith::topology::Field* solution,
                   const spatialdata::units::Scales& scales);

    // PRIVATE MEMEBRS
    // //////////////////////////////////////////////////////////////////////////////////////////////////
private:

    std::string _labelName; ///< Name of label associated with patch.
    int _labelValue; ///< Value of label associated with patch.
    spatialdata::spatialdb::SpatialDB* _db; ///< Spatial database with values for initial condition.

    // NOT IMPLEMENTED /////////////////////////////////////////////////////////////////////////////////////////////////
private:

    InitialConditionPatch(const InitialConditionPatch&); ///< Not implemented
    const InitialConditionPatch& operator=(const InitialConditionPatch&); ///< Not implemented

}; // InitialConditionPatch

// End of file
