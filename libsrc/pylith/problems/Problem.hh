// -*- C++ -*-
//
// ======================================================================
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
// ======================================================================
//

/**
 * @file libsrc/problems/Problem.hh
 *
 * @brief C++ object that manages the solution of a problem.formulating the equations.
 *
 * We cast the problem in terms of F(t,s,\dot{s}) = G(t,s), s(t0) = s0.
 *
 * In PETSc time stepping (TS) notation, G is the RHS, and F is the I
 * function (which we call the LHS).
 */
#if !defined(pylith_problems_problem_hh)
#define pylith_problems_problem_hh

#include "problemsfwd.hh" // forward declarations

#include "pylith/utils/PyreComponent.hh" // ISA PyreComponent

#include "pylith/feassemble/feassemblefwd.hh" // HOLDSA Integrator, Constraint, Observers
#include "pylith/materials/materialsfwd.hh" // HOLDSA Material
#include "pylith/bc/bcfwd.hh" // HOLDSA BoundaryCondition
#include "pylith/sources/sourcesfwd.hh" // HOLDSA Source
#include "pylith/faults/faultsfwd.hh" // HOLDSA FaultCohesive
#include "spatialdata/spatialdb/spatialdbfwd.hh" // HASA GravityField

#include "pylith/topology/topologyfwd.hh" // USES Mesh, Field
#include "spatialdata/units/unitsfwd.hh" // HASA Nondimensional

#include "pylith/utils/petscfwd.h" // USES PetscVec, PetscMat

#include "pylith/problems/Physics.hh" // USES Problem::Formulation

#include "pylith/utils/array.hh" // HASA std::vector

class pylith::problems::Problem : public pylith::utils::PyreComponent {
    friend class TestProblem; // unit testing

    // PUBLIC ENUM /////////////////////////////////////////////////////////////////////////////////////////////////////
public:

    enum SolverTypeEnum {
        LINEAR, // Linear solver.
        NONLINEAR, // Nonlinear solver.
    }; // SolverType

    // PUBLIC METHODS //////////////////////////////////////////////////////////////////////////////////////////////////
public:

    /// Constructor
    Problem(void);

    /// Destructor
    virtual ~Problem(void);

    /// Deallocate PETSc and local data structures.
    void deallocate(void);

    /** Set formulation for equations.
     *
     * @param[in] value Formulation type.
     */
    void setFormulation(const pylith::problems::Physics::FormulationEnum value);

    /** Get formulation for equations.
     *
     * @returns Formulation type.
     */
    pylith::problems::Physics::FormulationEnum getFormulation(void) const;

    /** Set solver type.
     *
     * @param[in] value Solver type.
     */
    void setSolverType(const SolverTypeEnum value);

    /** Get solver type.
     *
     * @returns Solver type.
     */
    SolverTypeEnum getSolverType(void) const;

    /** Set manager of scales used to nondimensionalize problem.
     *
     * @param[in] dim Nondimensionalizer.
     */
    void setNormalizer(const spatialdata::units::Nondimensional& dim);

    /** Set gravity field.
     *
     * @param[in] g Gravity field.
     */
    void setGravityField(spatialdata::spatialdb::GravityField* const g);

    /** Register observer to receive notifications.
     *
     * Observers are used for output.
     *
     * @param[in] observer Observer to receive notifications.
     */
    void registerObserver(pylith::problems::ObserverSoln* observer);

    /** Remove observer from receiving notifications.
     *
     * @param[in] observer Observer to remove.
     */
    void removeObserver(pylith::problems::ObserverSoln* observer);

    /** Set solution field.
     *
     * @param[in] field Solution field.
     */
    void setSolution(pylith::topology::Field* field);

    /** Set materials.
     *
     * @param[in] materials Array of materials.
     * @param[in] numMaterials Number of materials.
     */
    void setMaterials(pylith::materials::Material* materials[],
                      const int numMaterials);

    /** Set sources.
     *
     * @param[in] sources Array of sources.
     * @param[in] numSource Number of sources.
     */
    void setSources(pylith::sources::Source* sources[],
                    const int numSources);

    /** Set boundary conditions.
     *
     * @param[in] bc Array of boundary conditions.
     * @param[in] numBC Number of boundary conditions.
     */
    void setBoundaryConditions(pylith::bc::BoundaryCondition* bc[],
                               const int numBC);

    /** Set interior interface conditions.
     *
     * @param[in] interfaces Array of interior interfaces.
     * @param[in] numInterfaces Number of interior interfaces.
     */
    void setInterfaces(pylith::faults::FaultCohesive* faults[],
                       const int numFaults);

    /** Do minimal initialization.
     *
     * @param mesh Finite-element mesh.
     */
    virtual
    void preinitialize(const pylith::topology::Mesh& mesh);

    /// Verify configuration.
    virtual
    void verifyConfiguration(void) const;

    /// Initialize problem.
    virtual
    void initialize(void);

    // PROTECTED MEMBERS ///////////////////////////////////////////////////////////////////////////////////////////////
protected:

    pylith::topology::Field* _solution; ///< Solution field.

    spatialdata::units::Nondimensional* _normalizer; ///< Nondimensionalization of scales.
    spatialdata::spatialdb::GravityField* _gravityField; ///< Gravity field.

    std::vector<pylith::materials::Material*> _materials; ///< Array of materials.
    std::vector<pylith::sources::Source*> _sources; ///< Array of sources.
    std::vector<pylith::bc::BoundaryCondition*> _bc; ///< Array of boundary conditions.
    std::vector<pylith::faults::FaultCohesive*> _interfaces; ///< Array of interior interfaces.

    std::vector<pylith::feassemble::Integrator*> _integrators; ///< Array of integrators.
    std::vector<pylith::feassemble::Constraint*> _constraints; ///< Array of constraints.
    pylith::problems::ObserversSoln* _observers; ///< Subscribers of solution updates.

    pylith::problems::Physics::FormulationEnum _formulation; ///< Formulation for equations.
    SolverTypeEnum _solverType; ///< Problem (solver) type.

    // PRIVATE METHODS /////////////////////////////////////////////////////////////////////////////////////////////////
private:

    /// Check material and interface ids.
    void _checkMaterialIds(void) const;

    /// Create array of integrators from materials, interfaces, and boundary conditions.
    void _createIntegrators(void);

    /// Create array of constraints from materials, interfaces, and boundary conditions.
    void _createConstraints(void);

    /// Setup solution subfields and discretization.
    void _setupSolution(void);

    // Setup field so Lagrange multiplier subfield is limited to degrees of freedom associated with the cohesive cells.
    void _setupLagrangeMultiplier(void);

    // NOT IMPLEMENTED /////////////////////////////////////////////////////////////////////////////////////////////////
private:

    Problem(const Problem&); ///< Not implemented
    const Problem& operator=(const Problem&); ///< Not implemented

}; // Problem

#endif // pylith_problems_problem_hh

// End of file
