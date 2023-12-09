// =================================================================================================
// This code is part of PyLith, developed through the Computational Infrastructure
// for Geodynamics (https://github.com/geodynamics/pylith).
//
// Copyright (c) 2010-2023, University of California, Davis and the PyLith Development Team.
// All rights reserved.
//
// See https://mit-license.org/ and LICENSE.md and for license information.
// =================================================================================================

/** @file libsrc/faults/KinSrcAuxiliaryFactory.hh
 *
 * @brief C++ helper class for setting up auxiliary fields for slip time functions.
 */

#if !defined(pylith_faults_kinsrcauxiliaryfactory_hh)
#define pylith_faults_kinsrcauxiliaryfactory_hh

#include "faultsfwd.hh" // forward declarations
#include "pylith/feassemble/AuxiliaryFactory.hh" // ISA AuxiliaryFactory

class pylith::faults::KinSrcAuxiliaryFactory : public pylith::feassemble::AuxiliaryFactory {
    friend class TestSlipFnAuxiliaryFactory; // unit testing

    // PUBLIC METHODS //////////////////////////////////////////////////////////////////////////////////////////////////
public:

    /// Default constructor.
    KinSrcAuxiliaryFactory(void);

    /// Destructor.
    virtual ~KinSrcAuxiliaryFactory(void);

    /// Add slip initiation time (relative to origin time) subfield to auxiliary field.
    void addInitiationTime(void);

    /// Add rise time subfield to auxiliary field.
    void addRiseTime(void);

    /// Add impulse duration subfield to auxiliary field.
    void addImpulseDuration(void);

    /// Add final_slip subfield to auxiliary field.
    void addFinalSlip(void);

    /// Add slip_rate subfield to auxiliary field.
    void addSlipRate(void);

    /// Add time_history_value subfield to auxiliary field.
    void addTimeHistoryValue(void);

    /** Update time history value subfield for current time.
     *
     * @param[inout] auxiliaryField Auxiliary field to update.
     * @param[in] t Current time.
     * @param[in] timeScale Time scale for nondimensionalization.
     * @param[in] dbTimeHistory Time history database.
     */
    static
    void updateTimeHistoryValue(pylith::topology::Field* auxiliaryField,
                                const PylithReal t,
                                const PylithReal timeScale,
                                spatialdata::spatialdb::TimeHistory* const dbTimeHistory);

    // NOT IMPLEMENTED /////////////////////////////////////////////////////////////////////////////////////////////////
private:

    KinSrcAuxiliaryFactory(const KinSrcAuxiliaryFactory&); ///< Not implemented.
    const KinSrcAuxiliaryFactory& operator=(const KinSrcAuxiliaryFactory&); ///< Not implemented

}; // class KinSrcAuxiliaryFactory

#endif // pylith_faults_kinsrcauxiliaryfactory_hh

// End of file
