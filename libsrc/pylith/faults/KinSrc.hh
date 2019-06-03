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

/** @file libsrc/faults/KinSrc.hh
 *
 * @brief C++ object for managing parameters for a kinematic
 * earthquake source.
 */

#if !defined(pylith_faults_kinsrc_hh)
#define pylith_faults_kinsrc_hh

// Include directives ---------------------------------------------------
#include "faultsfwd.hh" // forward declarations

#include "pylith/utils/PyreComponent.hh" // ISA PyreComponent

#include "pylith/topology/FieldBase.hh" // USES FieldBase::SpaceDim
#include "pylith/feassemble/feassemblefwd.hh" // USES AuxiliaryFactory
#include "pylith/utils/petscfwd.h" // HASA PetscVec

#include "spatialdata/geocoords/geocoordsfwd.hh" // USES CoordSys
#include "spatialdata/spatialdb/spatialdbfwd.hh" // USES SpatialDB
#include "spatialdata/units/unitsfwd.hh" // USES Nondimensional

// KinSrc -------------------------------------------------------------
/** @brief Kinematic earthquake source.
 *
 * KinSrc is responsible for providing the value of slip at time t
 * over a fault surface.
 *
 * The fault integrator's auxiliary field has the 'slip' subfield. The
 * auxiliary subfields in this object use the discretization from the
 * 'slip' subfield.
 */
class pylith::faults::KinSrc : public pylith::utils::PyreComponent {
    friend class TestKinSrc; // unit testing

    // PUBLIC METHODS /////////////////////////////////////////////////////
public:

    /// Default constructor.
    KinSrc(void);

    /// Destructor.
    virtual ~KinSrc(void);

    /// Deallocate PETSc and local data structures.
    virtual
    void deallocate(void);

    /** Set origin time for earthquake source.
     *
     * @param value Origin time for earthquake source.
     */
    void originTime(const PylithReal value);

    /** Get origin time for earthquake source.
     *
     * @returns Origin time for earthquake source.
     */
    PylithReal originTime(void) const;

    /** Get auxiliary field associated with the kinematic source.
     *
     * @return field Auxiliary field for the kinematic source.
     */
    const pylith::topology::Field& auxField(void) const;

    /** Set the spatial database for filling auxiliary subfields.
     *
     * @param[in] value Pointer to database.
     */
    void auxFieldDB(spatialdata::spatialdb::SpatialDB* value);

    /** Initialize kinematic (prescribed slip) earthquake source.
     *
     * @param[in] auxField Auxiliary field associated with fault finite-element integration.
     * @param[in] normalizer Normalizer for nondimensionalizing values.
     * @param[in] cs Coordinate system for problem.
     */
    void initialize(const pylith::topology::Field& auxField,
                    const spatialdata::units::Nondimensional& normalizer,
                    const spatialdata::geocoords::CoordSys* cs);

    /** Set slip subfield in fault integrator's auxiliary field at time t.
     *
     * @param[out] auxField Fault integrator's auxiliary field.
     * @param[in] t Time t.
     */
    void slip(pylith::topology::Field* const auxField,
              const PylithScalar t);

    // PROTECTED METHODS //////////////////////////////////////////////////
protected:

    /** Setup auxiliary subfields (discretization and query fns).
     *
     * Create subfields in auxiliary fields (includes name of the field,
     * vector field type, discretization, and scale for
     * nondimensionalization) and set query functions for filling them
     * from a spatial database.
     *
     * @attention The order of the calls to subfieldAdd() must match the
     * order of the auxiliary fields in the FE kernels.
     *
     * @param[in] normalizer Normalizer for nondimensionalizing values.
     * @param[in] cs Coordinate system for problem.
     */
    virtual
    void _auxFieldSetup(const spatialdata::units::Nondimensional& normalizer,
                        const spatialdata::geocoords::CoordSys* cs) = 0;

    /** Set constants used in finite-element integrations.
     *
     * @param[in] auxField Auxiliary field associated with fault finite-element integration.
     */
    void _setFEConstants(const pylith::topology::Field& auxField) const;

    // PROTECTED MEMBERS //////////////////////////////////////////////////
protected:

    pylith::faults::KinSrcAuxiliaryFactory* _auxFactory; ///< Factory for auxiliary subfields.
    PetscPointFunc _slipFnKernel; ///< Kernel for slip time function.
    pylith::topology::Field* _auxField; ///< Auxiliary field for this integrator.
    PetscVec _slipLocalVec; ///< PETSc local Vec to hold slip for this rupture.

    // PRIVATE MEMBERS ////////////////////////////////////////////////////
private:

    PylithReal _originTime; ///< Origin time for earthquake source

    // NOT IMPLEMENTED ////////////////////////////////////////////////////
private:

    KinSrc(const KinSrc&); ///< Not implemented
    const KinSrc& operator=(const KinSrc&); ///< Not implemented

}; // class KinSrc

#endif // pylith_faults_kinsrc_hh

// End of file
