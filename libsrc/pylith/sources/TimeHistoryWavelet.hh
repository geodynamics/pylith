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
// See LICENSE.md for license information.
//
// ----------------------------------------------------------------------
//

/** @file libsrc/sources/TimeHistoryWavelet.hh
 *
 * @brief C++ class for Ricker source time function.
 */

#if !defined(pylith_materials_timehistorywavelet_hh)
#define pylith_materials_timehistorywavelet_hh

#include "sourcesfwd.hh" // forward declarations

#include "pylith/sources/SourceTimeFunctionMomentTensorForce.hh" // ISA SourceTimeFunctionMomentTensorForce

class pylith::sources::TimeHistoryWavelet : public pylith::sources::SourceTimeFunctionMomentTensorForce {
    friend class TestTimeHistoryWavelet; // unit testing

    // PUBLIC METHODS //////////////////////////////////////////////////////////////////////////////////////////////////
public:

    /// Default constructor.
    TimeHistoryWavelet(void);

    /// Destructor.
    ~TimeHistoryWavelet(void);

    /// Deallocate PETSc and local data structures.
    void deallocate(void);

    /** Set time history database.
     *
     * @param[in] db Time history database.
     */
    void setTimeHistoryDB(spatialdata::spatialdb::TimeHistory* th);

    /** Get time history database.
     *
     * @preturns Time history database.
     */
    const spatialdata::spatialdb::TimeHistory* getTimeHistoryDB(void);

    /** Use time history term in time history expression.
     *
     * @param[in] value True if using time history term in expression.
     */
    void useTimeHistory(const bool value);

    /** Get flag associated with using time history term in time history expression.
     *
     * @returns True if using time history term in expression, false otherwise.
     */
    bool useTimeHistory(void) const;

    /** Get auxiliary factory associated with physics.
     *
     * @return Auxiliary factory for physics object.
     */
    pylith::sources::AuxiliaryFactoryMomentTensorForce* getAuxiliaryFactory(void);

    /** Add source time subfields to auxiliary field.
     *
     * @param[inout] auxiliaryField Auxiliary field.
     */
    void addAuxiliarySubfields(void);

    /** Get g1v kernel for residual, G(t,s).
     *
     * @param[in] coordsys Coordinate system.
     *
     * @return residual kernel for g1v.
     */
    PetscPointFunc getKernelg1v_explicit(const spatialdata::geocoords::CoordSys* coordsys) const;

    /** Update auxiliary field for current time.
     *
     * @param[inout] auxiliaryField Auxiliary field to update.
     * @param[in] t Current time.
     * @param[in] timeScale Time scale for nondimensionalization.
     */
    PetscPointFunc updateAuxiliaryField(pylith::topology::Field* auxiliaryField,
                              const PylithReal t,
                              const PylithReal timeScale);
    
    // PRIVATE MEMBERS /////////////////////////////////////////////////////////////////////////////////////////////////
private:

    pylith::sources::AuxiliaryFactorySourceTime* _auxiliaryFactory; ///< Factory for creating auxiliary subfields.
    spatialdata::spatialdb::TimeHistory* _dbTimeHistory; ///< Time history database.
    bool _useTimeHistory; ///< Use time history term.

    // NOT IMPLEMENTED /////////////////////////////////////////////////////////////////////////////////////////////////
private:

    TimeHistoryWavelet(const TimeHistoryWavelet&); ///< Not implemented.
    const TimeHistoryWavelet& operator=(const TimeHistoryWavelet&); /// Not implemented.

}; // class TimeHistoryWavelet

#endif // pylith_materials_timehistorywavelet_hh

// End of file
