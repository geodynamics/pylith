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

    // PRIVATE MEMBERS /////////////////////////////////////////////////////////////////////////////////////////////////
private:

    spatialdata::spatialdb::TimeHistory *_dbTimeHistory;                   ///< Time history database.
    pylith::sources::AuxiliaryFactorySourceTime* _auxiliaryFactory; ///< Factory for creating auxiliary subfields.

    // NOT IMPLEMENTED /////////////////////////////////////////////////////////////////////////////////////////////////
private:

    TimeHistoryWavelet(const TimeHistoryWavelet&); ///< Not implemented.
    const TimeHistoryWavelet& operator=(const TimeHistoryWavelet&); /// Not implemented.

}; // class TimeHistoryWavelet

#endif // pylith_materials_timehistorywavelet_hh

// End of file
