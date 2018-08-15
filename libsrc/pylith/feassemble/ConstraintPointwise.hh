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
// Copyright (c) 2010-2016 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

/** @file libsrc/feassemble/Constraint.hh
 *
 * @brief C++ abstract base class defining interface for constraints
 * applied to finite-elements.
 */

#if !defined(pylith_feassemble_constraintpointwise_hh)
#define pylith_feassemble_constraintpointwise_hh

// Include directives ---------------------------------------------------
#include "pylith/feassemble/feassemblefwd.hh"

#include "pylith/feassemble/Observers.hh" // ISA Observers

#include "pylith/topology/FieldBase.hh" // USES FieldBase::discretizations_map
#include "pylith/utils/array.hh" // HASA int_array
#include "pylith/utils/utilsfwd.hh" // HOLDSA Logger

#include "spatialdata/units/unitsfwd.hh" // USES Nondimensional
#include "spatialdata/spatialdb/spatialdbfwd.hh" // USES SpatialDB

// ConstraintPointwise---------------------------------------------------
/** @brief Abstract base class for defining constraints on boundaries.
 */
class pylith::feassemble::ConstraintPointwise :
    public pylith::feassemble::Observers {
    friend class TestConstraintPointwise;   // unit testing

    // PUBLIC METHODS /////////////////////////////////////////////////////
public:

    /// Default constructor.
    ConstraintPointwise(void);

    /// Destructor.
    virtual ~ConstraintPointwise(void);

    /// Deallocate PETSc and local data structures.
    virtual
    void deallocate(void);

    /** Set indices of constrained degrees of freedom at each location.
     *
     * Example: [0, 1] to apply forces to x and y degrees of freedom in
     * a Cartesian coordinate system.
     *
     * @param[in] dof Array of indices for constrained degrees of freedom.
     * @param[in] size Size of array
     */
    void constrainedDOF(const int* flags,
                        const int size);

    /** Get indices of constrained degrees of freedom.
     *
     * @returns Array of indices for constrained degrees of freedom.
     */
    const pylith::int_array& constrainedDOF(void) const;

    /** Get mesh associated with integrator domain.
     *
     * @returns Mesh associated with integrator domain.
     */
    virtual
    const pylith::topology::Mesh& domainMesh(void) const = 0;

    /** Get auxiliary field.
     *
     * @returns field Field over boundary.
     */
    const pylith::topology::Field* auxField(void) const;

    /** Set spatial database for filling auxiliary subfields.
     *
     * @param[in] value Pointer to database.
     */
    void auxFieldDB(spatialdata::spatialdb::SpatialDB* value);

    /** Set discretization information for auxiliary subfield.
     *
     * @param[in] name Name of auxiliary subfield.
     * @param[in] basisOrder Polynomial order for basis.
     * @param[in] quadOrder Order of quadrature rule.
     * @param[in] isBasisContinuous True if basis is continuous.
     * @param[in] feSpace Finite-element space.
     */
    void auxSubfieldDiscretization(const char* name,
                                   const int basisOrder,
                                   const int quadOrder,
                                   const bool isBasisContinuous,
                                   const pylith::topology::FieldBase::SpaceEnum feSpace);

    /** Set manager of scales used to nondimensionalize problem.
     *
     * @param dim Nondimensionalizer.
     */
    void normalizer(const spatialdata::units::Nondimensional& dim);

    /** Verify configuration is acceptable.
     *
     * @param[in] solution Solution field.
     */
    virtual
    void verifyConfiguration(const pylith::topology::Field& solution) const = 0;

    /** Initialize integrator.
     *
     * @param[in] solution Solution field (layout).
     */
    virtual
    void initialize(const pylith::topology::Field& solution) = 0;

    /** Update auxiliary fields at beginning of time step.
     *
     * @param[in] t Current time.
     * @param[in] dt Current time step.
     */
    virtual
    void prestep(const double t,
                 const double dt);

    /** Update at end of time step.
     *
     * @param[in] t Current time.
     * @param[in] tindex Current time step.
     * @param[in] dt Current time step.
     * @param[in] solution Solution at time t.
     */
    virtual
    void poststep(const PylithReal t,
                  const PylithInt tindex,
                  const PylithReal dt,
                  const pylith::topology::Field& solution);

    /** Set constrained values in solution field.
     *
     * @param[out] solution Solution field.
     * @param[in] t Current time.
     */
    virtual
    void setSolution(pylith::topology::Field* solution,
                     const double t) = 0;

    // PROTECTED METHODS //////////////////////////////////////////////////
protected:

    /** Get factory for setting up auxliary fields.
     *
     * @returns Factor for auxiliary fields.
     */
    virtual
    pylith::feassemble::AuxiliaryFactory* _auxFactory(void) = 0;

    // PROTECTED MEMBERS //////////////////////////////////////////////////
protected:

    spatialdata::units::Nondimensional* _normalizer;   ///< Nondimensionalizer.
    int_array _constrainedDOF; ///< List of constrained degrees of freedom at each location.

    pylith::topology::Field *_auxField; ///< Auxiliary field for this constraint.

    pylith::utils::EventLogger* _logger;   ///< Event logger.

    // NOT IMPLEMENTED ////////////////////////////////////////////////////
private:

    ConstraintPointwise(const ConstraintPointwise &m); ///< Not implemented
    const ConstraintPointwise& operator=(const ConstraintPointwise& m); ///< Not implemented

}; // class Constraint

#endif // pylith_feassemble_constraintpointwise_hh


// End of file
