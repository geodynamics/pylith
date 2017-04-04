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
#include "feassemblefwd.hh"

#include "pylith/utils/PyreComponent.hh" // ISA PyreComponent

#include "pylith/topology/topologyfwd.hh" // HOLDSA Field
#include "pylith/utils/utilsfwd.hh" // HOLDSA Logger
#include "pylith/utils/array.hh" // HASA int_array
#include "spatialdata/units/unitsfwd.hh" // USES Nondimensional
#include "spatialdata/spatialdb/spatialdbfwd.hh" // USES SpatialDB

#include "pylith/topology/FieldBase.hh" // USES DiscretizeInfo


#include <map> // HOLDSA std::map

// ConstraintPointwise---------------------------------------------------
/** @brief Abstract base class for defining constraints on boundaries.
 */
class pylith::feassemble::ConstraintPointwise : public pylith::utils::PyreComponent
{ // class ConstraintPointwise
    friend class TestConstraintPointwise;   // unit testing

    // PUBLIC METHODS /////////////////////////////////////////////////////
public:

    /// Default constructor.
    ConstraintPointwise(void);

    /// Destructor.
    virtual
    ~ConstraintPointwise(void);

    /// Deallocate PETSc and local data structures.
    virtual
    void deallocate(void);

    /** Set name of field in solution to constrain.
     *
     * @param[in] value Name of field in solution to constrain.
     */
    void field(const char* value);

    /** Get name of field in solution to constrain.
     *
     * @returns Name of field in solution to constrain.
     */
    const char* field(void) const;

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

    /** Get auxiliary fields.
     *
     * @returns field Field over boundary.
     */
    const pylith::topology::Field& auxFields(void) const;

    /** Check whether constraint has a given auxiliary field.
     *
     * @param[in] name Name of field.
     *
     * @returns True if constraint has auxiliary field, false otherwise.
     */
    bool hasAuxField(const char* name);

    /** Get auxiliary field.
     *
     * @param[out] field Pointer to field.
     * @param[in] name Name of field to retrieve.
     */
    void getAuxField(pylith::topology::Field *field,
                     const char* name) const;

    /** Set spatial database for auxiliary fields.
     *
     * @param[in] value Pointer to database.
     */
    void auxFieldsDB(spatialdata::spatialdb::SpatialDB* value);

    /** Set discretization information for auxiliary subfield.
     *
     * @param[in] name Name of auxiliary subfield.
     * @param[in] basisOrder Polynomial order for basis.
     * @param[in] quadOrder Order of quadrature rule.
     * @param[in] isBasisContinuous True if basis is continuous.
     */
    void auxFieldDiscretization(const char* name,
                                const int basisOrder,
                                const int quadOrder,
                                const bool isBasisContinuous);

    /** Get discretization information for auxiliary subfield.
     *
     * @param[in] name Name of subfield.
     * @return Discretization information for auxiliary subfield. If
     * discretization information was not set, then use "default".
     */
    const pylith::topology::FieldBase::DiscretizeInfo& auxFieldDiscretization(const char* name) const;

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
    void verifyConfiguration(const pylith::topology::Field& solution) const;

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

    /** Set constrained values in solution field.
     *
     * @param[out] solution Solution field.
     * @param[in] t Current time.
     */
    virtual
    void setSolution(pylith::topology::Field* solution,
                     const double t) = 0;

    // PROTECTED TYPEDEFS /////////////////////////////////////////////////
protected:

    typedef std::map<std::string, pylith::topology::FieldBase::DiscretizeInfo> discretizations_type;

    // PROTECTED MEMBERS //////////////////////////////////////////////////
protected:

    spatialdata::units::Nondimensional* _normalizer;   ///< Nondimensionalizer.
    std::string _field; ///< Name of solution field to constrain.
    int_array _constrainedDOF; ///< List of constrained degrees of freedom at each location.

    /// Auxiliary fields for this problem
    pylith::topology::Field *_auxFields;

    /// Database of values for auxiliary fields.
    spatialdata::spatialdb::SpatialDB* _auxFieldsDB;

    /// Set auxiliary fields via query.
    pylith::topology::FieldQuery* _auxFieldsQuery;

    /// Map from auxiliary field to discretization.
    discretizations_type _auxFieldsFEInfo;

    pylith::utils::EventLogger* _logger;   ///< Event logger.

    // NOT IMPLEMENTED ////////////////////////////////////////////////////
private:

    /// Not implemented
    ConstraintPointwise(const ConstraintPointwise &m);

    /// Not implemented
    const ConstraintPointwise& operator=(const ConstraintPointwise& m);

}; // class Constraint

#endif // pylith_feassemble_constraintpointwise_hh


// End of file
