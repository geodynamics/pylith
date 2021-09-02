// -*- C++ -*-
//
// ----------------------------------------------------------------------
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
// See LICENSE.md.md for license information.
//
// ----------------------------------------------------------------------
//

/** @file libsrc/bc/TimeDependentAuxiliaryFactory.hh
 *
 * @brief C++ helper class for setting up auxiliary fields for time-dependent boundary conditions.
 */

#if !defined(pylith_bc_timedependentauxiliaryfactory_hh)
#define pylith_bc_timedependentauxiliaryfactory_hh

#include "bcfwd.hh" // forward declarations
#include "pylith/feassemble/AuxiliaryFactory.hh" // ISA AuxiliaryFactory

#include "spatialdata/spatialdb/spatialdbfwd.hh" // USES TimeHistory

class pylith::bc::TimeDependentAuxiliaryFactory : public pylith::feassemble::AuxiliaryFactory {
    friend class TestDirichletAuxiliaryFactory; // unit testing

    // PUBLIC ENUMS ////////////////////////////////////////////////////////////////////////////////////////////////////
public:

    enum ReferenceEnum {
        XYZ=0, ///< Coordinate directions (x, y, z).
        TANGENTIAL_NORMAL=1, ///< Directions tangential and normal to the boundary (tangential_1, tangential_2, normal).
    };

    // PUBLIC METHODS //////////////////////////////////////////////////////////////////////////////////////////////////
public:

    /** Default constructor.
     *
     * @param[in] reference Reference for coordinate directions in auxiliary subfield.s
     */
    TimeDependentAuxiliaryFactory(const ReferenceEnum reference=XYZ);

    /// Destructor.
    ~TimeDependentAuxiliaryFactory(void);

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
     * @param[in] dt Current time step.
     * @param[in] timeScale Time scale for nondimensionalization.
     * @param[in] dbTimeHistory Time history database.
     */
    static
    void updateAuxiliaryField(pylith::topology::Field* auxiliaryField,
                              const PylithReal t,
                              const PylithReal dt,
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

    TimeDependentAuxiliaryFactory(const TimeDependentAuxiliaryFactory &); ///< Not implemented.
    const TimeDependentAuxiliaryFactory& operator=(const TimeDependentAuxiliaryFactory&); ///< Not implemented

}; // class TimeDependentAuxiliaryFactory

#endif // pylith_bc_timedependentauxiliaryfactory_hh

// End of file
