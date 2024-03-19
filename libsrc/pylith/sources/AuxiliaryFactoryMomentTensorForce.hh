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

/** @file libsrc/sources/AuxiliaryFactoryMomentTensorForce.hh
 *
 * @brief C++ helper class for setting up auxiliary subfields for the moment tensor force equation.
 */

#if !defined(pylith_sources_auxiliaryfactorymomenttensorforce_hh)
#define pylith_sources_auxiliaryfactorymomenttensorforce_hh

#include "sourcesfwd.hh" // forward declarations
#include "pylith/feassemble/AuxiliaryFactory.hh" // ISA AuxiliaryFactory

#include "spatialdata/spatialdb/spatialdbfwd.hh" // USES GravityField

class pylith::sources::AuxiliaryFactoryMomentTensorForce : public pylith::feassemble::AuxiliaryFactory {
    friend class TestAuxiliaryFactoryMomentTensorForce; // unit testing

    // PUBLIC METHODS //////////////////////////////////////////////////////////////////////////////////////////////////
public:

    /// Default constructor.
    AuxiliaryFactoryMomentTensorForce(void);

    /// Destructor.
    virtual ~AuxiliaryFactoryMomentTensorForce(void);

    /// Add moment tensor subfield to auxiliary subfields.
    void addMomentTensor(void);

    /// Add time delay subfield to auxiliary subfields.
    void addTimeDelay(void);

    /// Add center frequency subfield to auxiliary subfields.
    void addCenterFrequency(void);

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

    // NOT IMPLEMENTED /////////////////////////////////////////////////////////////////////////////////////////////////
private:

    AuxiliaryFactoryMomentTensorForce(const AuxiliaryFactoryMomentTensorForce &); ///< Not implemented.
    const AuxiliaryFactoryMomentTensorForce& operator=(const AuxiliaryFactoryMomentTensorForce&); ///< Not implemented

}; // class AuxiliaryFactoryMomentTensorForce

#endif // pylith_sources_auxiliaryfactorymomenttensorforce_hh

// End of file
