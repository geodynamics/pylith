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
// Copyright (c) 2010-2022 University of California, Davis
//
// See LICENSE.md.md for license information.
//
// ----------------------------------------------------------------------
//

/** @file libsrc/faults/TractionPerturbation.hh
 *
 * @brief C++ implementation of traction perturbation with time-dependent expression for value.
 *
 * f(x,t) = f_0(x) + \dot{f}_1(x)*(t-t_1(x)) + f_2(x)*a(t-t_2(x)).
 *
 * Auxiliary fields:
 *     if _useInitial
 *         initial amplitude (vector) f_0(x)
 *    if _useRate
 *        rate amplitude (vector) \dot{f}_1(x)
 *        rate start (scalar) t_1(x)
 *    if _useTimeHistory
 *        time history amplitude (vector) f_2(x)
 *        time history start (scalar) t_2(x)
 *        time history value (scalar) a(t-t_2(x))
 */

#if !defined(pylith_faults_tractionperturbation_hh)
#define pylith_faults_tractionperturbation_hh

#include "pylith/utils/PyreComponent.hh" // ISA PyreComponent
#include "pylith/faults/faultsfwd.hh" // forward declaration

#include "pylith/topology/topologyfwd.hh" // USES Field
#include "pylith/bc/bcfwd.hh" // USES TimeDependentAuxiliaryFactory
#include "pylith/utils/petscfwd.h" // USES PetscVec

#include "spatialdata/units/Nondimensional.hh" // USGES Nondimensional
#include "spatialdata/geocoords/geocoordsfwd.hh" // USES CoordSys
#include "spatialdata/spatialdb/spatialdbfwd.hh" // HOLDSA TimeHistory

class pylith::faults::TractionPerturbation : pylith::utils::PyreComponent {
    friend class TestTractionPerturbation; // unit testing

    // PUBLIC METHODS /////////////////////////////////////////////////////////////////////////////
public:

    /// Default constructor.
    TractionPerturbation(void);

    /// Destructor.
    ~TractionPerturbation(void);

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

    /** Use initial value term in time history expression.
     *
     * @param[in] value True if using initial value term in expression.
     */
    void useInitial(const bool value);

    /** Get flag associated with using initial value term in time history expression.
     *
     * @returns True if using initial value term in expression, false otherwise.
     */
    bool useInitial(void) const;

    /** Use rate value term in time history expression.
     *
     * @param[in] value True if using rate value term in expression.
     */
    void useRate(const bool value);

    /** Get flag associated with using rate value term in time history expression.
     *
     * @returns True if using rate value term in expression, false otherwise.
     */
    bool useRate(void) const;

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

    /** Set origin time for earthquake source.
     *
     * @param value Origin time for earthquake source.
     */
    void setOriginTime(const PylithReal value);

    /** Get origin time for earthquake source.
     *
     * @returns Origin time for earthquake source.
     */
    PylithReal getOriginTime(void) const;

    /** Get auxiliary field associated with the kinematic source.
     *
     * @return field Auxiliary field for the kinematic source.
     */
    const pylith::topology::Field& getAuxiliaryField(void) const;

    /** Set the spatial database for filling auxiliary subfields.
     *
     * @param[in] value Pointer to database.
     */
    void setAuxiliaryFieldDB(spatialdata::spatialdb::SpatialDB* value);

    /** Initialize kinematic (prescribed slip) earthquake source.
     *
     * @param[in] auxField Auxiliary field associated with fault finite-element integration.
     * @param[in] normalizer Normalizer for nondimensionalizing values.
     * @param[in] cs Coordinate system for problem.
     */
    void initialize(const pylith::topology::Field& auxField,
                    const spatialdata::units::Nondimensional& normalizer,
                    const spatialdata::geocoords::CoordSys* cs);

    /** Get requested slip subfields at time t.
     *
     * @param[inout] tractionLocalVec Local PETSc vector for traction values.
     * @param[in] faultAuxiliaryField Auxiliary field for fault.
     * @param[in] t Time t.
     * @param[in] timeScale Time scale for nondimensionalization.
     */
    virtual
    void getTraction(PetscVec tractionLocalVec,
                     pylith::topology::Field* faultAuxiliaryField,
                     const PylithScalar t,
                     const PylithScalar timeScale);

    // PROTECTED METHODS //////////////////////////////////////////////////////////////////////////
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
    void _auxiliaryFieldSetup(const spatialdata::units::Nondimensional& normalizer,
                              const spatialdata::geocoords::CoordSys* cs) = 0;

    /// Set kernels for computing traction perturbation.
    void _setTractionKernel(void);

    /** Set constants used in finite-element integrations.
     *
     * @param[in] auxField Auxiliary field associated with fault finite-element integration.
     */
    void _setFEConstants(const pylith::topology::Field& auxField) const;

    // PRIVATE MEMBERS ////////////////////////////////////////////////////////////////////////////
private:

    spatialdata::spatialdb::TimeHistory* _dbTimeHistory; ///< Time history database.
    pylith::bc::TimeDependentAuxiliaryFactory* _auxiliaryFactory; ///< Factory for auxiliary subfields.
    pylith::topology::Field* _auxiliaryField; ///< Auxiliary field for this integrator.
    PetscPointFunc _tractionKernel; ///< Kernel for traction perturbation.

    bool _useInitial; ///< Use initial value term.
    bool _useRate; ///< Use rate term.
    bool _useTimeHistory; ///< Use time history term.
    PylithReal _originTime; ///< Origin time for earthquake source

    // NOT IMPLEMENTED ////////////////////////////////////////////////////////////////////////////
private:

    TractionPerturbation(const TractionPerturbation&); ///< Not implemented.
    const TractionPerturbation& operator=(const TractionPerturbation&); ///< Not implemented.

}; // class TractionPerturbation

#endif // pylith_faults_tractionperturbation_hh

// End of file
