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
// Copyright (c) 2010-2017 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//
#pragma once

/** @file libsrc/sources/AuxiliaryFactorySquarePulseSource.hh
 *
 * @brief C++ helper class for setting up auxiliary subfields for the squarepulse source equation.
 */

#if !defined(pylith_sources_auxiliaryfactorysquarepulsesource_hh)
#define pylith_sources_auxiliaryfactorysquarepulsesource_hh

#include "pylith/sources/sourcesfwd.hh"                         // forward declarations
#include "pylith/feassemble/AuxiliaryFactory.hh" // ISA AuxiliaryFactory

#include "spatialdata/spatialdb/spatialdbfwd.hh" // USES GravityField

class pylith::sources::AuxiliaryFactorySquarePulseSource : public pylith::feassemble::AuxiliaryFactory {
    friend class TestAuxiliaryFactorySquarePulseSource; // unit testing

    // PUBLIC ENUMS ////////////////////////////////////////////////////////////////////////////////////////////////////
public:

    enum ReferenceEnum {
        XYZ=0, ///< Coordinate directions (x, y, z).
        TANGENTIAL_NORMAL=1, ///< Directions tangential and normal to the boundary (tangential_1, tangential_2, normal).
    };

    // PUBLIC METHODS //////////////////////////////////////////////////////////////////////////////////////////////////
public:

    /// Default constructor.
    AuxiliaryFactorySquarePulseSource(const ReferenceEnum reference=XYZ);

    /// Destructor.
    virtual ~AuxiliaryFactorySquarePulseSource(void);

    /// Add fluid density subfield to auxiliary subfields.
    void addVolumeFlowRate(void);

    /// Add time delay subfield to auxiliary subfields.
    void addTimeDelay(void);

    /// Add initial amplitude field to auxiliary fields.
    void addInitialAmplitude(void);

    /// Add rate amplitude field to auxiliary fields.
    void addRateAmplitude(void);

    /// Add rate start time amplitude field to auxiliary fields.
    void addRateStartTime(void);

    /// Add time history amplitude field to auxiliary fields.
    void addTimeHistoryAmplitude(void);

    /// Add time history start time field to auxiliary fields.
    void addTimeHistoryStartTime(void);

    /// Add time history value field to auxiliary fields.
    void addTimeHistoryValue(void);

    /** Update auxiliary field for current time.
     *
     * @param[inout] auxiliaryField Auxiliary field to update.
     * @param[in] t Current time.
     * @param[in] timeScale Time scale for nondimensionalization.
     * @param[in] dbTimeHistory Time history database.
     */
    static
    void updateAuxiliaryField(pylith::topology::Field* auxiliaryField,
                              const PylithReal t,
                              const PylithReal timeScale,
                              spatialdata::spatialdb::TimeHistory* const dbTimeHistory);

    // PRIVATE METHODS /////////////////////////////////////////////////////////////////////////////////////////////////
private:

    /** Set names of vector field components in auxiliary subfield.
     *
     * @param[in] description Subfield description.
     */
    void _setVectorFieldComponentNames(pylith::topology::FieldBase::Description* description);

    // PRIVATE MEMBERS /////////////////////////////////////////////////////////////////////////////////////////////////
private:

    ReferenceEnum _auxComponents; ///< Coordinate system reference for field components.

    // NOT IMPLEMENTED /////////////////////////////////////////////////////////////////////////////////////////////////
private:

    AuxiliaryFactorySquarePulseSource(const AuxiliaryFactorySquarePulseSource &);                  ///< Not implemented.
    const AuxiliaryFactorySquarePulseSource &operator=(const AuxiliaryFactorySquarePulseSource &); ///< Not implemented

}; // class AuxiliaryFactorySquarePulseSource

#endif // pylith_sources_auxiliaryfactorysquarepulsesource_hh

// End of file
