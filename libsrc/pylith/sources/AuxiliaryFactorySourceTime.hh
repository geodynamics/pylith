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

/** @file libsrc/sources/AuxiliaryFactorySourceTime.hh
 *
 * @brief C++ helper class for setting up auxiliary subfields for the wellbore source equation.
 */

#if !defined(pylith_sources_auxiliaryfactorysourcetime_hh)
#define pylith_sources_auxiliaryfactorysourcetime_hh

#include "sourcesfwd.hh" // forward declarations
#include "pylith/sources/AuxiliaryFactoryMomentTensorForce.hh" // ISA AuxiliaryFactoryMomentTensorForce

class pylith::sources::AuxiliaryFactorySourceTime : public pylith::sources::AuxiliaryFactoryMomentTensorForce {
    friend class TestAuxiliaryFactorySourceTime; // unit testing

    // PUBLIC METHODS //////////////////////////////////////////////////////////////////////////////////////////////////
public:

    /// Default constructor.
    AuxiliaryFactorySourceTime(void);

    /// Destructor.
    virtual ~AuxiliaryFactorySourceTime(void);

    /// Add center frequency subfield to auxiliary subfields.
    void addCenterFrequency(void);

    // /// Add time history amplitude field to auxiliary fields.
    // void addTimeHistoryAmplitude(void);

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

    // NOT IMPLEMENTED /////////////////////////////////////////////////////////////////////////////////////////////////
private:

    AuxiliaryFactorySourceTime(const AuxiliaryFactorySourceTime &); ///< Not implemented.
    const AuxiliaryFactorySourceTime& operator=(const AuxiliaryFactorySourceTime&); ///< Not implemented

}; // class AuxiliaryFactorySourceTime

#endif // pylith_sources_auxiliaryfactorysourcetime_hh

// End of file
