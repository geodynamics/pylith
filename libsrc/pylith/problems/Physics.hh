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

/*** @file libsrc/problems/Physics.hh
 *
 * @brief Abstract base class for objects defining physics, such as behavior
 * of a bulk makterial, boundary condition, interface, or constraint.
 */

#if !defined(pylith_problems_physics_hh)
#define pylith_problems_physics_hh

#include "problemsfwd.hh" // forward declarations

#include "pylith/utils/PyreComponent.hh" // ISA PyreComponent

#include "pylith/topology/Field.hh" // USES FieldBase
#include "pylith/feassemble/feassemblefwd.hh" // USES Integrator, Constraint

#include "spatialdata/spatialdb/spatialdbfwd.hh" // USES SpatialDB
#include "spatialdata/units/unitsfwd.hh" // HASA Nondimensional

class pylith::problems::Physics : public pylith::utils::PyreComponent {
    friend class TestPhysics; // unit testing

    // PUBLIC ENUM /////////////////////////////////////////////////////////////////////////////////////////////////////
public:

    enum FormulationEnum {
        QUASISTATIC, // Without inertia; implicit time stepping.
        DYNAMIC, // With inertia; explicit time stepping).
        DYNAMIC_IMEX, // With inertia; implicit+explicit time stepping).
    }; // FormulationEnum

    // PUBLIC MEMBERS //////////////////////////////////////////////////////////////////////////////////////////////////
public:

    /// Constructor
    Physics(void);

    /// Destructor
    virtual ~Physics(void);

    /// Deallocate PETSc and local data structures.
    virtual
    void deallocate(void);

    /** Set manager of scales used to nondimensionalize problem.
     *
     * @param dim Nondimensionalizer.
     */
    void setNormalizer(const spatialdata::units::Nondimensional& dim);

    /** Get manager of scales used to nondimensionalize problem.
     *
     * @param dim Nondimensionalizer.
     */
    const spatialdata::units::Nondimensional& getNormalizer(void) const;

    /** Set formulation for equations.
     *
     * @param[in] value Formulation for equations.
     */
    void setFormulation(const FormulationEnum value);

    /** Set spatial database for populating auxiliary field.
     *
     * @param[in] value Spatial database with iniital values for auxiliary field.
     */
    void setAuxiliaryFieldDB(spatialdata::spatialdb::SpatialDB* const value);

    /** Set discretization information for auxiliary subfield.
     *
     * @param[in] subfieldName Name of auxiliary subfield.
     * @param[in] basisOrder Polynomial order for basis.
     * @param[in] quadOrder Order of quadrature rule.
     * @param[in] dimension Dimension of points for discretization.
     * @param[in] cellBasis Type of basis functions to use (e.g., simplex, tensor, or default).
     * @param[in] feSpace Finite-element space.
     * @param[in] isBasisContinuous True if basis is continuous.
     */
    void setAuxiliarySubfieldDiscretization(const char* subfieldName,
                                            const int basisOrder,
                                            const int quadOrder,
                                            const int dimension,
                                            const pylith::topology::FieldBase::CellBasis cellBasis,
                                            const pylith::topology::FieldBase::SpaceEnum feSpace,
                                            const bool isBasisContinuous);

    /** Set discretization information for derived subfield.
     *
     * @param[in] subfieldName Name of auxiliary subfield.
     * @param[in] basisOrder Polynomial order for basis.
     * @param[in] quadOrder Order of quadrature rule.
     * @param[in] dimension Dimension of points for discretization.
     * @param[in] cellBasis Type of basis functions to use (e.g., simplex, tensor, or default).
     * @param[in] feSpace Finite-element space.
     * @param[in] isBasisContinuous True if basis is continuous.
     */
    void setDerivedSubfieldDiscretization(const char* subfieldName,
                                          const int basisOrder,
                                          const int quadOrder,
                                          const int dimension,
                                          const pylith::topology::FieldBase::CellBasis cellBasis,
                                          const pylith::topology::FieldBase::SpaceEnum feSpace,
                                          const bool isBasisContinuous);

    /** Register observer to receive notifications.
     *
     * Observers are used for output.
     *
     * @param[in] observer Observer to receive notifications.
     */
    void registerObserver(pylith::problems::ObserverPhysics* observer);

    /** Remove observer from receiving notifications.
     *
     * @param[in] observer Observer to remove.
     */
    void removeObserver(pylith::problems::ObserverPhysics* observer);

    /** Get observers receiving notifications of physics updates.
     *
     * @returns Observers receiving notifications.
     */
    pylith::problems::ObserversPhysics* getObservers(void);

    /** Get constants used in kernels (point-wise functions).
     *
     * @param[in] solution Solution field.
     * @param[in] dt Current time step.
     *
     * @return Array of constants.
     */
    const pylith::real_array& getKernelConstants(const PylithReal dt);

    /** Verify configuration is acceptable.
     *
     * @param[in] solution Solution field.
     */
    virtual
    void verifyConfiguration(const pylith::topology::Field& solution) const = 0;

    /** Create integrator and set kernels.
     *
     * @param[in] solution Solution field.
     * @returns Integrator if applicable, otherwise NULL.
     */
    virtual
    pylith::feassemble::Integrator* createIntegrator(const pylith::topology::Field& solution) = 0;

    /** Create constraint and set kernels.
     *
     * @param[in] solution Solution field.
     * @returns Constraint if applicable, otherwise NULL.
     */
    virtual
    pylith::feassemble::Constraint* createConstraint(const pylith::topology::Field& solution) = 0;

    /** Create auxiliary field.
     *
     * @param[in] solution Solution field.
     * @param[in\ physicsMesh Finite-element mesh associated with physics.
     *
     * @returns Auxiliary field if applicable, otherwise NULL.
     */
    virtual
    pylith::topology::Field* createAuxiliaryField(const pylith::topology::Field& solution,
                                                  const pylith::topology::Mesh& physicsMesh) = 0;

    /** Create derived field.
     *
     * @param[in] solution Solution field.
     * @param[in] physicsMesh Finite-element mesh associated with physics.
     *
     * @returns Derived field if applicable, otherwise NULL.
     */
    virtual
    pylith::topology::Field* createDerivedField(const pylith::topology::Field& solution,
                                                const pylith::topology::Mesh& physicsMesh) = 0;

    /** Update time-dependent auxiliary field.
     *
     * @param[inout] auxiliaryField Auxiliary field.
     * @param[in] t Current time.
     */
    virtual
    void updateAuxiliaryField(pylith::topology::Field* auxiliaryField,
                              const double t);

    // PROTECTED METHODS ///////////////////////////////////////////////////////////////////////////////////////////////
protected:

    /** Get auxiliary factory associated with physics.
     *
     * @return Auxiliary factory for physics object.
     */
    virtual
    pylith::feassemble::AuxiliaryFactory* _getAuxiliaryFactory(void) = 0;

    /** Get derived factory associated with physics.
     *
     * @return Derived factory for physics object.
     */
    virtual
    pylith::topology::FieldFactory* _getDerivedFactory(void);

    /** Update kernel constants.
     *
     * @param[in] dt Current time step.
     */
    virtual
    void _updateKernelConstants(const PylithReal dt);

    // PROTECTED MEMBERS ///////////////////////////////////////////////////////////////////////////////////////////////
protected:

    spatialdata::units::Nondimensional* _normalizer; ///< Nondimensionalizer.
    FormulationEnum _formulation; ///< Formulation for equations.
    pylith::real_array _kernelConstants; ///< Constants used in finite-element kernels (point-wise functions).

    // PRIVATE MEMBERS //////////////////////////////////////////////////////
private:

    pylith::problems::ObserversPhysics* _observers; ///< Subscribers of updates.

    // NOT IMPLEMENTED /////////////////////////////////////////////////////////////////////////////////////////////////
private:

    Physics(const Physics&); ///< Not implemented.
    const Physics& operator=(const Physics&); ///< Not implemented.

}; // Physics

#endif // pylith_problems_physics_hh

// End of file
